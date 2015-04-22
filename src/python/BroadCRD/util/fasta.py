###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# Python functions for reading and writing FASTA files.

from __future__ import print_function

import re

def deNewline(inString):
    if (' ' in inString):
        return inString.replace('\n',' ')
    else:
        return inString.replace('\n','')
    
def cleanUpNewlines(inTuple):
    return [deNewline(s) for s in inTuple]

fastaContigRE = r'^>(.+)\n((?:[^>]+)+)'
fastaContigPattern = re.compile(fastaContigRE, re.MULTILINE)

fastaContigNoCommentRE = r'^>(\S+).*\n((?:[^>]+)+)'
fastaContigNoCommentPattern = re.compile(fastaContigNoCommentRE, re.MULTILINE)


def parseFasta(data, comments=True):
    if comments:
        return [cleanUpNewlines(t) for t in fastaContigPattern.findall(data)]
    else:
        return [cleanUpNewlines(t) for t in 
            fastaContigNoCommentPattern.findall(data)]

def readFasta(filename, comments=False):
    '''Read the FASTA file in as a list of pairs, first element name
    second element sequence, line breaks stripped out.'''
    f = open(filename)
    data = f.read()
    f.close()
    return parseFasta(data, comments)

def readFastaDict(filename, comments=False):
    '''Read the FASTA file, returning it in dictionary format with the
    contig names as keys.'''
    fasta_list = readFasta(filename, comments)
    fasta_dict = {}
    for f in fasta_list:
        fasta_dict[f[0]] = f[1]
    return fasta_dict

def writeFasta(name_read_pairs, filename, indices=None):
    '''Write the FASTA file, line-breaking the reads every 80 characters.'''
    if not indices:
        indices = xrange(0, len(name_read_pairs))
    with open(filename, 'w') as output_file:
        for nri in indices:
            print('>' + name_read_pairs[nri][0], file=output_file)
            contig_length = len(name_read_pairs[nri][1])
            lines = contig_length / 80
            if contig_length % 80 != 0:
                lines += 1
            for li in xrange(lines):
                start_line = li * 80
                end_line = min((li + 1) * 80, contig_length)
                print(name_read_pairs[nri][1][start_line:end_line],
                    file=output_file)
    return

def writeFastaDict(name_read_dict, filename):
    '''Write the FASTA file using a dictionary, keys will be alphabetically
    ordered.'''
    contig_names = list(name_read_dict.keys())
    contig_names.sort()
    name_read_pairs = []
    for c in contig_names:
        name_read_pairs.append((c, name_read_dict[c]))
    writeFasta(name_read_pairs, filename)
    return
