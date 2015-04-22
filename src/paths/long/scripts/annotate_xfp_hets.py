#!/usr/bin/env python

# "use Python-2.7" or newer


import locus_masks as lm
import subprocess
import re

class CounterDict(dict):
    def __init__(self):
        super(CounterDict,self).__init__(self)

    def incr(self,key):
        if self.has_key(key):
            self[key]+=1
        else:
            self[key]=1

print "loading dust mask..."
fDust="/home/unix/blau/wga/dustmasker/hg19/ofn"
dustmask=lm.dustmask(fDust)
print "loading segdup mask..."
fSegDup="/home/unix/blau/wga/segdup/GRCh37GenomicSuperDup.tab"
segdupmask=lm.segdupmask(fSegDup)

counts=CounterDict()
callers=set()
#for line in open('/home/unix/neilw/l/dev/BroadCRD/table.txt','r'):
print "running XFP..."
for line in subprocess.check_output(["XFP","x"]).split('\n'):
    m=re.match('^(X[XY]) ([^ ]+) (\d+)$', line.strip())
    if m:
        (mf,caller,locus)=m.groups()
        locus=int(locus)

        if dustmask.isMasked('X',locus,locus):
            counts.incr((mf,caller,'DUST'))
        if segdupmask.isMasked('X',locus,locus):
            counts.incr((mf,caller,'SEG_DUP'))

        counts.incr((mf,caller,'ALL'))

        callers.add(caller)


for key,val in counts.items():
    print key,val


print "sanity check:"
for caller in callers:
    print "x-het FP for {} is {:.2f}".format(caller, counts[('XY',caller,'ALL')] \
                    *100.0 / ( 1.8 * counts[('XX',caller,'ALL')] ) )


for genome in ['XY','XX']:
    print "{} table ===================".format(genome)
    titles=['CALLER','ALL','DUST','SEG_DUP']
    for title in titles:
        print "| {:12}".format(title),
    print "|"
    for title in titles:
        print "| {:12}".format(':---'),
    print "|"


    for caller in callers:
        print "| {:12}".format(caller),
        for measure in ['ALL','DUST','SEG_DUP']:
            print "| {:12}".format(counts[(genome,caller,measure)]),
        print "|"


