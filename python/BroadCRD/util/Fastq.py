###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2012) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# Simple functions for reading and writing FASTQ data.

# Original author: Michael G. Ross <mgross@broadinstitute.org>

class FastqFormatError(Exception):
    def __init__(self, state):
        self.state = state
        return
    def __str__(self):
        return repr(self.state)

def readFastq(filename):
    name_read_quality_triples = []
    with open(filename, 'r') as fastq_data:
        current_name = ''
        current_seq = ''
        current_plus = ''
        current_qual = ''

        for line_num, fastq_line in enumerate(fastq_data):
            if (not current_name and not current_seq and not current_plus and 
                not current_qual and fastq_line[0] == '@'):
                current_name = fastq_line[1:].rstrip()
            elif (current_name and not current_plus and not current_qual and 
                fastq_line[0] != '+'):
                current_seq += fastq_line.rstrip()
            elif (current_name and current_seq and not current_plus and
                not current_qual and fastq_line[0] == '+'):
                current_plus = fastq_line.rstrip()
            elif (current_name and current_seq and current_plus and
                len(current_qual) < len(current_seq)):
                current_qual += fastq_line.rstrip()
            else:
                raise FastqFormatError, '{0}:{1}'.format(filename, line_num)
            if (current_name and current_seq and current_plus and
                len(current_qual) == len(current_seq)):
                name_read_quality_triples.append((current_name, current_seq,
                    current_qual))
                current_name = ''
                current_seq = ''
                current_plus = ''
                current_qual = ''

    return name_read_quality_triples

def writeFastq(name_read_quality_triples, filename):
    with open(filename, 'w') as fastq_out:
        for nrq in name_read_quality_triples:
            print >>fastq_out, '@{0}\n{1}\n+\n{2}\n'.format(nrq[0], nrq[1],
                nrq[2]),
    return
    
def readFastqDict(filename):
    fastq_list = readFastq(filename)
    fastq_dict = {}
    for f in fastq_list:
        fastq_dict[f[0]] = (f[1], f[2])
    return fastq_dict

def writeFastqDict(fastq_dict, filename):
    contig_names = list(fastq_dict.keys())
    contig_names.sort()
    name_read_quality_triples = []
    for c in contig_names:
        (read, qual) = fastq_dict[c]
        name_read_pairs.append((c, read, qual))
    writeFastq(name_read_quality_triples, filename)
    return
