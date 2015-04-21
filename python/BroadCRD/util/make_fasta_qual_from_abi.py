#!/usr/bin/env python


import sys
import os.path
from Bio import SeqIO

class BadFileError(Exception): 
    def __init__(self,msg): 
        Exception.__init__(self,msg)

def main( argv=[__name__] ):
    if ( len(argv) < 3 ):
        print "usage: {0} OUTPUT_HEAD file.abi".format( argv[0] )
        print "\nAssumes one record per .abi file in its current form"
        return 1

    output_seq = argv[1] + ".fasta"
    output_qual = argv[1] + ".qual"

    if os.path.exists( output_seq ) or os.path.exists( output_qual ):
        raise BadFileError("output file " + output_seq + " (or .qual) already exists ")


    dup_ids=dict()

    with open( output_seq, 'w' ) as fdout_seq:
        with open( output_qual, 'w' ) as fdout_qual:
            for input in argv[2:]:
                print input
                # assumes one record per file
                record=SeqIO.read( open(input, 'rb') , "abi" )

                id = record.id
                if not dup_ids.has_key( id ): dup_ids[id]=1
                record.id = "{}.{}".format( id, dup_ids[id] )
                dup_ids[id]+=1

                SeqIO.write(record, fdout_seq, 'fasta' )
                SeqIO.write(record, fdout_qual, 'qual' )

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))


# vim:tw=1000:
