#!/usr/bin/env python

# nr_concordance -- calculate non-reference call concordance between
# VCFs based on a truth set.

HapMap3="hapmap3.3_r27_na12878_nonref_snps.b37.fwd.vcf"

import sys
import argparse
import math
from VCF import *

def check_sort_order( filename, last_chr_pos, chr, pos ):
    if chr == last_chr_pos[0] and pos < last_chr_pos[1]:
        raise Exception( """
        file {} doesn't seem to be sorted.
        Saw {}:{} and then {}:{}.
        """.format( filename, last_chr_pos[0], last_chr_pos[1], chr, pos) )

def nr_sensitivity( input, sample, truth, minqual=0, misses=None, debug=False ):
    print >>sys.stderr,"input={}\ntruth={}\nminqual={}".format(input,truth,minqual)

    truth_vcf = VCF( truth )
    eval_vcf = VCF( input )

    if misses:
        misses_fd = open( misses, 'w' )
        found = False
        for line in truth_vcf.metadata:
            if line[:len("#CHROM")]=="#CHROM":
                misses_fd.write( '##nr_concordance="comment={subset of '\
                    'missed sites created by DISCOVAR release bundle Python program '\
                    'nr_concordance.py}"\n')
                found = True
            misses_fd.write(line+"\n")
        if not found:
            raise Exception("program bug? made it this far without #CHROM in truth VCF?")

    else: misses_fd = None

    eval_sample_index = eval_vcf.sample_names.index(sample)
    truth_sample_index = truth_vcf.sample_names.index(sample)

    eval_chr=[]
    truth_chr=[]

    eval_gen = eval_vcf.lines()
    eval_line = eval_gen.next()
    eval_chr.append( eval_line.chr )

    eval_done = False

    n_truth_lines = 0
    n_site_hits = 0
    n_site_concords = 0

    last_truth=(None,None)
    last_eval=(None,None)

    for truth_line in truth_vcf.lines():

        check_sort_order( truth, last_truth, truth_line.chr, truth_line.pos )

        truth_genotype=truth_line.get_sample_dict(truth_sample_index)["GT"]
        if truth_genotype == "0/0" or truth_genotype=="0|0" or truth_genotype == ".": continue
        if debug: print >>sys.stderr,"seeking {}:{}".format( truth_line.chr, truth_line.pos )
        n_truth_lines += 1
        if n_truth_lines % 1000 == 0: print n_truth_lines

        # skip to the correct chromosome
        if len(truth_chr) == 0 or truth_chr[-1] != truth_line.chr:
            truth_chr.append( truth_line.chr )

        # if we've already passed this chr in the eval file, then spin
        if not eval_done and truth_line.chr != eval_line.chr and truth_line.chr in eval_chr:
            continue

        # if we've not already passed this chr, then find it
        while not eval_done and eval_line.chr != truth_line.chr:
            try:
                eval_line = eval_gen.next()
                check_sort_order( input, last_eval, eval_line.chr, eval_line.pos )
                if debug: print >>sys.stderr,"...next chr={}".format( eval_line.chr )
                if eval_chr[-1] != eval_line.chr:
                    eval_chr.append(eval_line.chr)
                    print eval_line.chr
            except StopIteration: eval_done=True

        # try to find the correct position
        while not eval_done and eval_line.pos < truth_line.pos \
                and eval_line.chr == truth_line.chr:
            try:
                eval_line = eval_gen.next()
                check_sort_order( input, last_eval, eval_line.chr, eval_line.pos )
                if debug: print >>sys.stderr,"...next chr:pos={}:{}".format( eval_line.chr, eval_line.pos )
                if eval_chr[-1] != eval_line.chr:
                    eval_chr.append(eval_line.chr)
                    print eval_line.chr
            except StopIteration: eval_done=True

        if minqual > 0 and eval_line.qual == '.':
            raise Exception("not sure what to do here, we're filtering on qual, but qual is '.'")

        if eval_done or eval_line.pos != truth_line.pos \
                or  eval_line.chr != truth_line.chr \
                or  ( eval_line.qual != '.' and float(eval_line.qual) < minqual ):
            if misses_fd: misses_fd.write(truth_line.line+"\n")
        else:
            if truth_line.ref != eval_line.ref:
                raise Exception("""
                    Your truth set does not seem to be called on the
                    same reference as your call set.  We're done here.
                    truth={}
                    truth_pos={}:{}
                    truth_ref={}

                    input={}
                    input_pos={}:{}
                    input_ref={}
                    """.format( truth, truth_line.chr, truth_line.pos,
                        truth_line.ref, input, eval_line.chr,
                        eval_line.pos, eval_line.ref ) )

            if debug:
                print >>sys.stderr,"""
                Evaluating:
                truth_pos={}:{}
                truth_ref={}
                truth_alt={}

                input_pos={}:{}
                input_ref={}
                input_alt={}

                """.format( truth_line.chr, truth_line.pos,
                        truth_line.ref, truth_line.alt, eval_line.chr,
                        eval_line.pos, eval_line.ref, eval_line.alt )

            # grab truth NR bases and eval NR bases
            eval_genotype=eval_line.get_sample_dict( eval_sample_index)["GT"]
            if eval_genotype != ".":
                eval_calls_idx = eval_genotype.split("/")
                if eval_calls_idx[0] == eval_genotype: eval_calls_idx = eval_genotype.split("|")
                if '0' in eval_calls_idx: eval_calls_idx.remove('0')
                eval_calls_idx = map(int, eval_calls_idx )
                truth_calls_idx = truth_genotype.split("/")
                if truth_calls_idx[0] == truth_genotype: truth_calls_idx = truth_genotype.split("|")
                if '0' in truth_calls_idx: truth_calls_idx.remove('0')
                truth_calls_idx = map(int, truth_calls_idx )

                if len(eval_calls_idx) > 0:
                    n_site_hits += 1
                    if debug: print >>sys.stderr, "accepting site hit at {}:{}".format(eval_line.chr, eval_line.pos)
                elif debug:
                    print >>sys.stderr, "no NR calls"

                for truth_nr in [ truth_line.alt_list[i-1] for i in truth_calls_idx ]:
                    # if any truth non-reference call is not found on
                    # the eval line, then we break without counting the
                    # concordance
                    if not truth_nr in eval_line.alt_list:
                        if debug: print >>sys.stderr,"""
                            non-concordant:
                            pos={}:{}
                            truth_ref={}
                            eval_ref={}
                            truth_alt={}
                            eval_alt={}
                            """.format( truth_line.chr, truth_line.pos,
                                    truth_line.ref, eval_line.ref,
                                    truth_line.alt_list,
                                    eval_line.alt_list )
                        break

                else:
                    # normal termination of the for loop, so count the
                    # concordance
                    if debug: print >>sys.stderr,"""
                            CONCORDANT:
                            pos={}:{}
                            truth_ref={}
                            eval_ref={}
                            truth_alt={}
                            eval_alt={}
                            """.format( truth_line.chr, truth_line.pos,
                                    truth_line.ref, eval_line.ref,
                                    truth_line.alt_list,
                                    eval_line.alt_list )
                    n_site_concords += 1

    print "n_truth_lines={}, n_site_hits={}, n_site_concords={}, site_hit_frac={}, site_concord_frac={}".format(
            n_truth_lines, n_site_hits, n_site_concords,
            n_site_hits/float(n_truth_lines),
            n_site_concords/float(n_truth_lines)
            )

    # check that for the eval chromosomes that are also in the truth
    # set, that they come in the same order

    # first form intersection set
    overlap_chr = set(truth_chr).intersection(set(eval_chr))

    truth_chr_rev = [ chr for chr in truth_chr if chr in overlap_chr ]
    truth_chr_rev.reverse()

    for chr in eval_chr:
        if chr not in overlap_chr: continue
        if chr != truth_chr_rev[-1]:
            raise Exception("""
            input chromosome ordering doesn't match truth chromosome ordering:
            input={}
            truth={}
            """.format( eval_chr, truth_chr ) )
        truth_chr_rev.pop()

    if misses_fd: misses_fd.close()




def main( argv=[__name__] ):
    parser=argparse.ArgumentParser(description='calculate non-reference call concordance between an evaluation VCF and a truth set VCF')
    parser.add_argument( 'input', help='VCF file to evaluate', nargs=1)
    parser.add_argument( 'sample_name', help='sample to evaluate -- must be in both VCFs', nargs=1)
    parser.add_argument( 'truth', help='VCF to use for comparison', nargs='?', default=HapMap3 )
    parser.add_argument( '-d','--debug', help='detail debugging output',action='store_true')
    parser.add_argument( '-m', '--misses', help='pseudo-vcf file (excerpt of truth) of misses')
    parser.add_argument( '-q', '--minqual', help='minimumal QUAL value in input (not truth) to accept', default=0, type=float)
    args = parser.parse_args(argv[1:])

    nr_sensitivity( args.input[0], args.sample_name[0], args.truth, args.minqual, args.misses, args.debug )

    return(0)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
