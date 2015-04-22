#!/usr/bin/env python

# annotate_variantsall - takes a variants.all(.all) file and adds
# transition (Ti) and transversion (Tv) markers to it.

import sys
import argparse

def load_dbsnp( file ):
    print "loading SNPs from {}".format(file)
    values={}
    nlines=0
    for line in open(file, 'r'):
        nlines = nlines+1
        if nlines % 10000 == 0: print "passed {} lines".format(nlines)
        if line[0] != '#':
            fields=line.split()
            (chr, coord, id, ref, alt, qual, filter, info) = fields
            if len(ref) == 1 and 1 in map(len,alt.split(',')):
                values[(chr,int(coord))] = (ref,alt)

    return values

def markup( input, output, discordance ):

    if discordance:
        dbsnp = load_dbsnp( discordance )

    print "reading variants.all {}".format(input)
    if output: outfd = open(output, 'w')
    else: outfd = sys.stdout

    for line in open( input, 'r' ):
        fields=line.split()
        (fosid,context,dunno,chrcoord,ref,alt) = fields[:6]
        (chr,coord)=chrcoord.split(':')
        marker=''
        if len(ref) == 1 and len(alt) == 1 and ref != alt:
            bases=sorted([ref,alt])
            # c.f. http://en.wikipedia.org/wiki/Transversion
            if bases==['A','G'] or bases==['C','T']:
                marker=marker+',Ti'
            else:
                marker=marker+',Tv'

            print chrcoord+":"
            if discordance:
                if dbsnp.has_key( (chr,int(coord)) ):
                    snp=dbsnp[ ( chr, int(coord) ) ]
                    print "\tdbsnp: {}".format(snp)
                    if snp[0] != ref:
                        raise Exception('dbsnp say ref is {}, we say {} at {}'.format( snp[0], ref, chrcoord ))
                    if alt in snp[1].split(','):
                        marker=marker+",dbsnp"
                else:
                    print "\tdnsnp: NONE"

        print >>outfd,"{}{}".format( line.strip(), marker)






def main( argv=[__name__] ):
    parser=argparse.ArgumentParser(description='annotate variants.all files e.g. with Ti/Tv')
    parser.add_argument( '-d', '--discordance', help='report SNP novelty with respect to this VCF', required=False)
    parser.add_argument( '-o', '--output', help='variants.all output file')
    parser.add_argument( 'input', help='variants.all input file' )
    args = parser.parse_args(argv[1:])

    return(markup( args.input, args.output, args.discordance ) )

if __name__ == "__main__":
    sys.exit(main(sys.argv))
