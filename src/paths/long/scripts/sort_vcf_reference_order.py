#!/usr/bin/env python

#
# sort_vcf_reference_order -- primarily made to patch up the cortex
# calls, but should work with other data.  Essentially reads the
# reference .dict file and then re-writes the VCF in the same contig
# order as the reference .dict file.  Does NOT rearrange the metadata
# "contig" lines -- not sure if that's a format violation, but doesn't
# matter for cortex.
#


import sys
import argparse
import math
import os
from VCF import *


def get_contig_list_from_fai(file):
    return [ line.split()[0] for line in open( file, 'r' ) ]

def get_contig_list_from_fasta(reference):
    contigs = []
    for line in open(reference, 'r'):
        if line[0] == '>':
            contig = line[1:].split()[0]
            if len(contigs) == 0 or contigs[-1] != contig:
                contigs.append(contig)
                print "found {}".format(contig)

    return contigs

def get_contig_list( reference ):
    d1 = reference + ".fai"
    d2 = os.path.splitext(reference)[0] + ".fai"

    if os.path.exists(d1):
        contigs = get_contig_list_from_fai(d1)
    elif os.path.exists(d2):
        contigs = get_contig_list_from_fai(d2)
    else:
        print "no {} or {} file (e.g. from samtools faidx)\n"\
                "must read reference the slow way".format( d1, d2 )
        contigs = get_contig_list_from_fasta(reference)

    return contigs



def sort_vcf( input, reference, output ):
    contig_list = get_contig_list( reference )

    print "read {} contigs".format(len(contig_list))

    v = VCF( input, True )      # request index of VCF upon open

    with open( output, 'w') as fd_out:
        fd_out.writelines( [ line+"\n" for line in v.metadata ] )

        for contig in contig_list:
            print "writing entries for contig {}".format(contig)
            # filter lines from vcf by contig
            count = 0
            if v.seek(contig) < 0:
                print "skipped {} because it's not in the VCF".format( contig )
            else:
                for line in v.lines( True, \
                        lambda raw_line: VCFLine(raw_line).chr == contig ):
                    fd_out.write( line.line+"\n" )
                    count += 1
                print "wrote {} entries for {}".format( count, contig )


def main( argv=[__name__] ):
    parser=argparse.ArgumentParser(description='rewrite VCF file in reference contig order')
    parser.add_argument( '-o', '--output', help='VCF output file', required=True)
    parser.add_argument( 'input_vcf', help='VCF input file' )
    parser.add_argument( 'reference_fasta', help='reference fasta file')
    args = parser.parse_args(argv[1:])

    return sort_vcf( args.input_vcf, args.reference_fasta, args.output )

if __name__ == "__main__":
    sys.exit(main(sys.argv))
