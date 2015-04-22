#!/usr/bin/env python
import h5py
import sys
import pandas as pd

if len(sys.argv) != 3:
    print >>sys.stderr, "usage: {} tag filename".format(sys.argv[0])
    sys.exit(2)

tag = sys.argv[1]
filename = sys.argv[2]

h=h5py.File(filename,'r')

try:
    if h[tag].shape:
        data=h[tag].value
        print pd.DataFrame(data).to_csv(index=False)
    else:
        print h[tag].value
except:
    print "{} does not seem to exist in {}".format(tag,filename)
    sys.exit(1)

