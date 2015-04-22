#!/usr/bin/env python
###############################################################################
##                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     ##
##       This software and its documentation are copyright (2013) by the     ##
##   Broad Institute.  All rights are reserved.  This software is supplied   ##
##   without any warranty or guaranteed support whatsoever. The Broad        ##
##   Institute is not responsible for its use, misuse, or functionality.     ##
###############################################################################

# filter_vcfs - filter variant lines where phred score in either REFP or
# ALTP exceeds a threshold.
#
# The following rules apply:
#
# - already filtered lines are left alone (lines not marked as "PASS" or
# ".")
#
# - lines that are "PASS" or "." with a phred score for any element in
# ALTP or REFP which is >0 AND <threshold are marked as filtered (e.g.
# "PhredFilter10" for lines filtered with a threshold of 10.
#
# - lines that are "." and pass the phred threshold are marked as
# "PASS"
#

import sys
import argparse
import math
import os

def filter_vcf( input, output ):
    debug = False
    phred90 = int(-10.0*math.log10( 1.0-0.9 ))
    phred995 = int(-10.0*math.log10( 1.0-0.995 ))
    print "input={}, output={}, phred90={}, phred995={}".format(input, output, phred90, phred995)

    # new filters must be added below as symbol=("name","description")
    # this is just to enforce consistency between what's in the head and
    # in the body of the vcf
    filters = dict(
        Filter090=( "PhredFilter{}".format(phred90), "Some allele was below prob 0.90 = phred {}".format(phred90) ),
        NoAllelesPass0995=( \
            "NoAllelesPassPhred{}".format(phred995), "Not one allele was above prob 0.995 = phred {}".format(phred995) ),
        RefCallOnly=( "RefCallOnly", "Non-variant line; only the reference was called" ),
        TooManyCalls=( "TooManyCalls", "Multi-allelic site (greater than two alleles)" )
    )


    with open(output, 'w') as fd:
        for line in open(input, 'r'):
            if line[0] == '#':
                if line[:len("#CHROM")] == "#CHROM":
                    # dump out the ##FILTER lines before the #CHROM
                    for id, descr in filters.viewvalues():
                        fd.write("##FILTER=<ID={id},Description=\"{descr}\">\n".format( id=id, descr=descr ))
                fd.write(line)
            else:
                debug_line = False

                fields=line.split()
                (chrom,pos,id,ref,alt,qual,filter,info,format)=fields[:9]
                samples = fields[9:]

                # don't filter unless we're already passing or not yet filtered
                if filter == '.' or filter == 'PASS':
                    names = format.split(':')
                    if (not "REFP" in names) or (not "ALTP" in names):
                        raise Exception("missing REFP and ALTP tags in line {}".format(line.strip() ) )

                    new_samples=[]
                    new_filters=set()

                    for sample in samples:
                        vals = sample.split(':')
                        if len(vals) != len(names):
                            raise Exception("sample values {} doesn't match format {} for line {}".format(
                                sample, format, line.strip()))
                        sample_info = dict( zip( names, vals ) )

                        probs = [ float( sample_info['REFP'] ) ]
                        probs.extend( map(float, sample_info['ALTP'].split(',') ) )

                        if True in map( lambda x: x>0 and x<phred90, probs ):
                            new_samples.append( sample )
                            new_filters.add( filters['Filter090'][0])
                            continue

                        calls = map( lambda x: x > phred995, probs )

                        if sum(calls) == 0:
                            new_samples.append( sample )
                            new_filters.add( filters['NoAllelesPass0995'][0] )
                            continue

                        if sum(calls) == 1 and calls[0] == True:
                            new_samples.append( sample )
                            new_filters.add( filters['RefCallOnly'][0] )
                            continue

                        if sum(calls) > 2:
                            new_samples.append( sample )
                            new_filters.add( filters['TooManyCalls'][0] )
                            continue

                        allele1 = calls.index(True)
                        if sum(calls) == 1:
                            allele2 = allele1
                        else:
                            allele2 = calls.index(True, allele1+1)

                        gt_idx = names.index('GT')
                        old_gt = vals[gt_idx]
                        vals[gt_idx] = '{}/{}'.format(allele1, allele2)

                        if debug and old_gt != vals[gt_idx]:
                            debug_line = True

                        new_samples.append( ':'.join(vals) )
                        new_filters.add( 'PASS' )

                    # okay, now re-write the filter tag based on
                    # "new_filters"

                    # RefCallOnly doesn't matter if there are other objections
                    if len(new_filters) > 1: new_filters.discard( filters['RefCallOnly'][0] )

                    # PASS only counts if there are no other objections
                    # (besides RefCallOnly removed above)
                    if len(new_filters) > 1: new_filters.discard('PASS')

                    filter= ",".join( new_filters )

                else:           # re-using previous filter becase line wasn't . or PASS
                    new_samples = samples


                if debug_line:
                    print "re-wrote genotypes:"
                    print "\told line:"
                    print "\t",line
                    print "\tnew genotypes:"
                    print "\t","\t".join( new_samples )
                    print "\n"

                fd.write( "\t".join( [chrom, pos, id, ref, alt, qual, filter, info, format] ) )
                fd.write( "\t" )
                fd.write( "\t".join( new_samples ) )
                fd.write( "\n" )



def main( argv=[__name__] ):
    parser=argparse.ArgumentParser(description='filter Discovar-generated VCF based on REFP/ALTP')
    parser.add_argument( '-o', '--output', help='VCF output file', required=True)
    parser.add_argument( 'input', help='VCF input file' )
    args = parser.parse_args(argv[1:])

    if os.path.exists( args.output ):
        print >>sys.stderr, \
            "Output file {} already exists.  Move it out of the way first.".format( args.output )
        return 1

    try:
        return(filter_vcf( args.input, args.output ) )
    except:
        print >>sys.stderr,"removing partial output file {}".format( args.output )
        os.unlink(args.output)
        raise

if __name__ == "__main__":
    sys.exit(main(sys.argv))
