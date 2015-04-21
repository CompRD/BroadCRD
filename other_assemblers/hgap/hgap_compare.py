#!/usr/bin/env python
###############################################################################
##                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     ##
##       This software and its documentation are copyright (2013) by the     ##
##   Broad Institute.  All rights are reserved.  This software is supplied   ##
##   without any warranty or guaranteed support whatsoever. The Broad        ##
##   Institute is not responsible for its use, misuse, or functionality.     ##
###############################################################################


import re
import sys
import shutil
import os.path


def main( argv = [__name__]  ):


    fosmids = range(0,107)

    finished_path = '/wga/dev/references/Homo_sapiens/NA12878_Fosmid_Pool.regions.fin'
    pb_assembly = 'polished_assembly.fasta'

    format_str = '{},\t'*7
    print format_str.format('fosmid', 'f_len', 'start', 'end','fraction', 'mismatches', 'indels')


    for fosmid in fosmids:

        # copy finished fosmid
        finished_filename = 'fos.{}.fasta'.format(fosmid)
        finished_source = finished_path + '/' + finished_filename
        if (os.path.isfile(finished_source) ):
            shutil.copy(finished_source, finished_filename);
        else:
            continue;

        lookup_table = finished_filename + '.lookup'
        lookup_table_log  = lookup_table + '.log'
#        os.system('MakeLookupTable SOURCE={} > {}'.format(finished_filename, lookup_table_log));
            
        qlt_log = '{}.visual_aligns'.format(fosmid)
 #       os.system('QueryLookupTable K=12 MM=12 MC=0.15 SEQS={} L={} SMITH_WAT=True VISUAL=True > {}'.format(pb_assembly, lookup_table, qlt_log))

        output = os.popen('grep vs {}'.format(qlt_log)).readlines()
#        os.system('grep vs {}'.format(qlt_log))
        if (len(output) == 0): 
            print '{},\tNo alignments found'.format(fosmid)
        else:
            coverage = 0
            for line in output :

                m = re.match( r'.*, (\d+) mismatches/(\d+) indels \(of (\d+)\), from (\d+)-(\d+) to (\d+)-(\d+) \(of (\d+)\)', line, re.M|re.I)

                if m:
                    mismatches = m.group(1)
                    indels = m.group(2)
                    seq_length = m.group(3)
                    seq_start = m.group(4)
                    seq_end = m.group(5)
                    target_start = m.group(6)
                    target_end = m.group(7)
                    target_length = m.group(8)
                    format_str = '{},\t'*7
                    fraction = '{0:.2f}'.format((int(target_end) - int( target_start)) / float(target_length))
                    print format_str.format(fosmid, target_length, target_start, target_end, fraction, mismatches, indels)



        print '*'*75
               
if __name__ == "__main__":
    sys.exit(main(sys.argv))
