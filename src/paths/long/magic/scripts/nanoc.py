#!/usr/bin/env python

import BroadCRD.util.Nanopore as Nanopore
import sys

class ReadPrinter:
    CONSENSUS = 0
    TEMPLATE = 1
    COMPLEMENT = 2
    FASTQ = 0
    FASTA = 1

    def __init__(self, filename):
	self.r = Nanopore.Read( filename )

    def split( self, string, columns ):
	return [ string[i:i+columns] \
		for i in range(0,len(string),columns) ]

    def text_fastq(self, rtype, out_type ):
	if rtype == self.CONSENSUS: fq = self.r.fastq
	elif rtype == self.TEMPLATE: fq = self.r.tfastq
	else: fq = self.r.cfastq 

	cols=80
	out=str()
	if out_type == self.FASTA:
	    out = out + ">"
            out = out + fq[0] + "\n"
            for f in self.split( fq[1], cols ): out = out + f + "\n"
        elif out_type == self.FASTQ:
            out = "\n".join(fq)
        else:
            raise("Bug: unrecognized out_type in text_fastq");

	return out

    def events(self, rtype):
        if rtype == self.CONSENSUS: return self.r.events
        elif rtype == self.COMPLEMENT: return self.r.cevents
        elif rtype == self.TEMPLATE: return self.r.tevents

    def model(self, rtype):
        if rtype == self.CONSENSUS:
            raise Exception("no model for consensus")
        elif rtype == self.COMPLEMENT:
            return self.r.cmodel
        elif rtype == self.TEMPLATE:
            return self.r.tmodel


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Extract Oxford Nanopore reads', \
	epilog="""
actions to perform:
	events - extract events
        model - extract models
	fastq - print fastq file
	fasta - print fasta file
    """, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-H', action='store_true', help='include header lines')
    parser.add_argument( 'action', help='action to perform -- see below')
    parser.add_argument( 'component', help='consensus, template, or complement')
    parser.add_argument( 'filenames', nargs='+', help='the .fast5 file to read')
    args=parser.parse_args()

    for filename in args.filenames:
        rp=ReadPrinter(filename)

        comp=ReadPrinter.CONSENSUS
        if args.component == 'template':
            comp=ReadPrinter.TEMPLATE
        elif args.component == 'complement':
            comp=ReadPrinter.COMPLEMENT

        if args.action=='fastq':
            sys.stdout.write(ReadPrinter(filename).text_fastq(comp, ReadPrinter.FASTQ))
        elif args.action=='fasta':
            sys.stdout.write(ReadPrinter(filename).text_fastq(comp, ReadPrinter.FASTA))
        elif args.action=='events':
            print ReadPrinter(filename).events(comp).to_csv(path_or_buf=None, header=args.H, index=False)
        elif args.action=='model':
            print ReadPrinter(filename).model(comp).to_csv(path_or_buf=None, header=args.H, index=False)


