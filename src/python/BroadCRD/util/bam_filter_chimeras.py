#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2013) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

import sys
import argparse
import re
from BroadCRD.util.Samtools import openbam
from BroadCRD.util.Samtools import SamLine
from BroadCRD.util.pushback_wrapper import pushback_wrapper


def copy_header_lines( input, output, readgroup = None ):
    for line in input:
        if line[:3] == '@RG' and readgroup != None:
            print >>sys.stderr, "warning: replacing the readgroup"
        elif line[0] == '@': output.write(line)
        else:
            input.pushback(line)
            if readgroup:
                output.write("@RG\t{}\n".format("\t".join(readgroup.split())))
            break

def grab_one_query( input ):
    last = None
    output_lines = []

    for line in input:
        sline = SamLine( line )
        if not last or sline.qname == last:
            output_lines.append( sline )
            last = sline.qname
        else:
            input.pushback(line)
            break

    return output_lines


def is_chimera( slines ):
    for sline in slines:
        if sline.flag & 2048:
            return True
    else:
        return False

def check_closest( one, two, checkme ):
    min_one=float("inf")
    min_two=float("inf")
    for c in checkme:
        if one.rname == c.rname:
            one_dist = abs(one.pos - c.pos )
            if one_dist < min_one and one_dist < 1000:
                min_one = one_dist
        if two.rname == c.rname:
            two_dist = abs(two.pos - c.pos )
            if two_dist < min_two and two_dist < 1000:
                min_two = two_dist

    return min_one <= min_two


def filter_chimerae( slines, checkme, debug = False ):
    if len(slines) == 1 or not is_chimera(slines):
            return slines
    elif len(slines) > 2:
        if debug: print >>sys.stderr, "WARNING: many reads, choose none:\n{}\n".format(slines)
        return []
    elif len(slines) < 1:
        raise Exception("BUG - no records to filter_chimeras")
    else:
        one = slines[0]
        two = slines[1]

        if debug: print "CHIMERIC READS:\n\t{}\n\t{}".format(one,two)

        # if RC and clipped at the end or FW and clipped at the front,
        # then take the other one
        if one.flag & 0x10 and re.match( ".*\d+[HS]$", one.cigar ) or \
            not one.flag & 0x10 and re.match( "^\d+[HS].*", one.cigar ):
            two.flag = two.flag & ~2048         # clear the supplemental flag
            if debug:
                print "RETURNING:\n{}".format(two)
                if check_closest( one, two, checkme ): print "GOOD\n"
                else:
                    print "BAD, compare against:\n{}\n".format(checkme)
            return [two]
        elif two.flag & 0x10 and re.match(".*\d+[HS]$", two.cigar ) or \
            not two.flag&0x10 and re.match("^\d+[HS].*", two.cigar ):
            one.flag = one.flag & ~2048         # clear the supplemental flag
            if debug:
                print "RETURNING:\n{}".format(one)
                if check_closest( two, one, checkme ): print "GOOD\n"
                else:
                    print "BAD, compare against:\n{}\n".format(checkme)
            return [one]
        else:
            raise Exception("neither of pair seems to be clipped as expected:\n{}\n{}\n".format(one,two) )


def find_closest( pone, ptwo ):
    """
    find elements in pone and ptwo that are closest, on the same
    chromosome and closer than 1000bp
    """
    min_dist=float("inf")
    min_winner=(None,None)
    for o in range(len(pone)):
        for t in range(len(ptwo)):
            if pone[o].rname == ptwo[t].rname:
                dist = abs(pone[o].pos - ptwo[t].pos)
                if dist < 1000 and dist < min_dist:
                    min_winner=( o, t )
                    min_dist = dist
    return min_winner


def filter_a_pair( pone, ptwo, lgoods, debug ):
    c_one = is_chimera(pone)
    c_two = is_chimera(ptwo)

    if c_one and c_two:
        raise Exception("caught a double chimera in 'filter-a-pair'")
    if not c_one and not c_two:
        return True

    # from the division of unnecessary parentheses
    if ( not c_one and len(pone) > 1 ) or len(pone) > 2:
        return False
    if ( not c_two and len(ptwo) > 1 ) or len(ptwo) > 2:
        return False

    closest=find_closest( pone, ptwo )

    if closest != (None, None):
        if c_one:
            if debug:
                print "REMOVING CHIMERA:"
                print "{}".format( pone[ closest[0] ] )
                print "CLOSE TO OTHER:"
                print "{}".format( ptwo )
            del pone[closest[0]]
            pone[0].flag &= ~2048        # clear supplemental flag
            if debug:
                print "KEEPING:"
                print "{}\n".format( pone )
        else:
            if debug:
                print "REMOVING CHIMERA:"
                print "{}".format( ptwo[ closest[1] ] )
                print "CLOSE TO OTHER:"
                print "{}".format( pone )
            del ptwo[closest[1]]
            ptwo[0].flag &= ~2048        # clear supplemental flag
            if debug:
                print "KEEPING:"
                print "{}\n".format( ptwo )
        lgoods[0] += 1
        return True
    elif debug:
        print "NO SOLUTION:"
        print "ONE:\n{}".format("\n".join(map(str,pone)))
        print "TWO:\n{}\n".format("\n".join(map(str,ptwo)))

    return False




def filter_paired_chimeras( f1, f2, o, readgroup, debug=False ):

    doubles = 0
    multis = 0
    bads = 0
    lgoods = [0]
    total = 0

    input1=pushback_wrapper( openbam( f1, 'r', header=True) )
    input2=pushback_wrapper( openbam( f2, 'r', header=True) )

    output=openbam( o, 'w' )

    readgroup_id = None
    if readgroup:
        # extract the ID: tag from the readgroup line
        readgroup_id=dict( [ tuple( x.split( ':' ) ) for x in readgroup.split() ]).get('ID')
        if not readgroup_id:
            raise Exception('readgroup line seems to be missing an ID: field {}'.format(readgroup ) )

    copy_header_lines( input1, output, readgroup )
    copy_header_lines( input2, open('/dev/null','w') )  # seek past the header

    while True:

        one = grab_one_query(input1)
        two = grab_one_query(input2)


        # double check
        if set( [ q.qname for q in one] ) != set( [ q.qname for q in two ] ):
            raise Exception("the input is not paired properly")

        if not one and not two: break

        total += 1

        if not one and two or not two and one:
            raise Exception("one file ran out of reads before the other")

        if is_chimera(one) and is_chimera(two):
            doubles += 1
            continue

#        one = filter_chimerae(one,two,debug)
#        two = filter_chimerae(two,one,debug)

#        if not one or not two:
#            bads += 1
#            continue

        if not filter_a_pair(one, two, lgoods, debug):
            bads += 1
            continue

        if len(one) > 1:
            print >>sys.stderr, "warning: multiple records in query one:"
            print >>sys.stderr, one
            multis += 1
            continue
        if len(two) > 1:
            print >>sys.stderr, "warning: multiple records in query two:"
            print >>sys.stderr, two
            multis += 1
            continue

        if debug and ( one[0].flag & 0x800 or two[0].flag & 0x800 ):
            print "secondary aligment flag was not cleared!"
            print "ONE: {}".format(one[0])
            print "TWO: {}".format(two[0])

        one[0].flag |= 0x41             # first in pair plus is_paired
        two[0].flag |= 0x81             # second in pair plus is_paired

        one[0].flag |= (two[0].flag & 0x10) << 1        # note if the OTHER is RC
        two[0].flag |= (one[0].flag & 0x10) << 1

        one[0].flag |= (two[0].flag & 0x4) << 1         # note if the OTHER is unmapped
        two[0].flag |= (one[0].flag & 0x4) << 1

        # set the "both are mapped" bit (0x2) if neither is unmapped (0x4)
        one[0].flag |= (( one[0].flag & 0x4 | two[0].flag & 0x4 ) >> 1) ^ 0x2
        two[0].flag |= (( one[0].flag & 0x4 | two[0].flag & 0x4 ) >> 1) ^ 0x2

        # set the contig name and position for the other
        one[0].rnext = two[0].rname
        one[0].pnext = two[0].pos

        two[0].rnext = one[0].rname
        two[0].pnext = one[0].pos

        if readgroup_id:
            one[0].rest = [ x for x in one[0].rest if x[:3] != 'ID:' ]
            one[0].rest.append( 'RG:Z:{}'.format( readgroup_id ) )
            two[0].rest = [ x for x in two[0].rest if x[:3] != 'ID:' ]
            two[0].rest.append( 'RG:Z:{}'.format( readgroup_id ) )

        output.write( str(one[0])+"\n" );
        output.write( str(two[0])+"\n" );

    print "total read pairs processed = {}".format(total)
    print "pairs with chimeras properly processed = {}%".format( lgoods[0]*100.0/total)
    print "pairs with chimeras on both reads = {}%".format(doubles*100.0/total)
    print "unfilterable or multiple placements = {}%".format(bads*100.0/total )
    print "ignored multiple alignments = {}%".format(multis*100.0/total)


def main( argv=[__name__] ):
    parser=argparse.ArgumentParser(description="filter chimeras in SAM/BAM files in various ways")
    parser.add_argument( 'f1', help='input for first of read pair' )
    parser.add_argument( 'f2', help='input for second of read pair' )
    parser.add_argument( 'o', help='output for paired SAM/BAM file' )
    parser.add_argument( '-r', '--rg', help='read group line for the reads minus the @RG (must have ID: and SM:' )
    parser.add_argument( '-d', '--debug', action='store_true' )
    args=parser.parse_args( argv[1:] )

    filter_paired_chimeras( args.f1, args.f2, args.o, args.rg, args.debug )


if __name__ == "__main__":
    sys.exit(main(sys.argv))
