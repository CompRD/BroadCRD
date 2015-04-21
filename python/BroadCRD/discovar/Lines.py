import re
import sys

class Lines(list):
    """Load an a.s.lines.map file as a list of 8-tuples

    constants for fields CHR, START, STOP, LINE, ESTART, ESTOP, LEN, COV"""

    CHR=0
    START=1
    STOP=2
    LINE=3
    ESTART=4
    ESTOP=5
    LEN=6
    COV=7

    def __init__(self, filename=None):
        super(Lines,self).__init__()
        if filename == None: return
        for line in open(filename,'r'):
            m=re.match('(\d+)\s+([\d,]+)\s+([\d,]+)\s+line\[(\d+)\]\s+(\d+)..(\d+)\s+len=(\d+)\s+(cov=)?([\d.]+)?', line)
            if m:
                g=m.groups()
                chr=int(g[Lines.CHR])
                start=int(g[Lines.START].translate(None,','))
                stop=int(g[Lines.STOP].translate(None,','))
                line=int(g[Lines.LINE])
                start_edge=int(g[Lines.ESTART])
                stop_edge=int(g[Lines.ESTOP])
                line_len=int(g[Lines.LEN])
                cov=g[8]        # ugly...skipped one
                if cov: cov=float(cov)
                self.append( (chr,start,stop,line,start_edge,stop_edge,line_len,cov))
