###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

import subprocess
import shlex

class SamtoolsError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class SamLine:
    def __init__(self, line):
        fields = line.strip().split()
        try:
            ( qname, flag, rname, pos, mapq,
                cigar, rnext, pnext, tlen, seq,
                qual ) = fields[:11]
        except ValueError as e:
            print "Error parsing SAM format line -- too few fields?"
            print line
            raise

        self.qname = qname
        self.flag = int(flag)
        self.rname = rname
        self.pos = int(pos)
        self.mapq = int(mapq)
        self.cigar = cigar
        self.rnext = rnext
        self.pnext = int(pnext)
        self.tlen = int(tlen)
        self.seq = seq
        self.qual = qual
        self.rest = fields[11:]

    def __repr__(self):
        return "\t".join( [ self.qname, str(self.flag), self.rname, str(self.pos), \
                str(self.mapq), self.cigar, self.rnext, str(self.pnext), \
                str(self.tlen), self.seq, self.qual ] ) + "\t" +"\t".join(self.rest)


def openbam(bamfile,  mode='r', header=False, records=True):
    '''Pipe "samtools view"'s output to create a line-buffered file-like
    object of SAM header lines and records from a SAM or BAM file (if mode="r") 
    or to write BAM files from a stream of SAM records (if mode="w").'''

    if mode == 'r':
        samtools_cmd = 'samtools view '    

        if header and records:
            samtools_cmd += '-h '
        elif header and not records:
            samtools_cmd += '-H '
        elif not header and records:
            samtools_cmd += ''
        else:
            raise SamtoolsError('when reading BAM/SAMs, either the header or '
                'records argument must be true')

        if bamfile[-4:] == '.sam':
            samtools_cmd += '-S '
        samtools_cmd += bamfile
        return subprocess.Popen(shlex.split(samtools_cmd), shell=False, bufsize=-1,
            stdout=subprocess.PIPE).stdout
    elif mode == 'w':
        if (bamfile[-4:] == '.bam'):
            return subprocess.Popen(shlex.split('samtools view -S -b -'), shell=False, 
                bufsize=-1, stdin=subprocess.PIPE,
                stdout=open(bamfile, 'w')).stdin
        elif (bamfile[-4:] == '.sam'):
            return open(bamfile, 'w')
        else:
            raise SamtoolsError('when writing BAM/SAMS, the filename must '
                'end with either ".sam" or ".bam" to indicate the output ' 
                'format')
    else:
        raise SamtoolsError('mode must be "r" or "w"')
