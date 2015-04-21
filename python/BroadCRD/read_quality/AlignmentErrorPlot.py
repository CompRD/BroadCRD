#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2010) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# A program for plotting the output of alignment_error_stats.py.

# Original author: Michael G. Ross <mgross@broadinstitute.org>

import optparse
import matplotlib.pyplot
import numpy
import os
import sys

oparser = optparse.OptionParser(usage='usage: %prog [options] '
    'alignment_errors1.txt [alignment_errors2.txt] ...')
oparser.add_option('--min_readlength', type='int', default=0)
oparser.add_option('--max_readlength', type='int', default=sys.maxint)
oparser.add_option('--output_file', type='string', default=None)
oparser.add_option('--nodisplay', action='store_true')
(options, args) = oparser.parse_args()

print ('plotting reads with lengths between ' + str(options.min_readlength) +
    ' and ' + str(options.max_readlength))

overallmax = -1
overallcountmax = -1
num_datasets = None
num_errorfiles = len(args)
for p in xrange(num_errorfiles):
    errorfile = args[p]

    # how many alignments are in this file?    
    errordat = open(errorfile, 'r')
    if not num_datasets:
        num_datasets = len(errordat.readline().rstrip().split('\t')) / 5
    elif num_datasets != len(errordat.readline().rstrip().split('\t')) / 5:
        print >> sys.stderr, 'mismatch in number of alignments per file'
        exit(1)
    errordat.close()

    # load alignment data
    data_cols = []
    for n in range(num_datasets):
        data_cols.extend(range(n * 5 + 2, n * 5 + 6))
    data = numpy.loadtxt(errorfile, usecols=data_cols, skiprows=1,
        dtype='int')
        
    # exclude too short or too long reads
    long_enough = numpy.all(data[:,range(0, num_datasets * 4, 4)] >= 
        options.min_readlength, axis=1)
    short_enough = numpy.all(data[:,range(0, num_datasets * 4, 4)] <= 
        options.max_readlength, axis=1)

    histdata = (data[:,range(3, num_datasets * 4, 4)] /
        numpy.array(data[:,range(0, num_datasets * 4, 4)], dtype='float'))
    histdata = histdata[long_enough & short_enough,:]

    # update the cross-file max
    maxdata = histdata[:].max()
    overallmax = max(overallmax, maxdata)
    
    # plot the data
    for n in range(num_datasets):
        matplotlib.pyplot.subplot(num_errorfiles, num_datasets,
            p * num_datasets + 1 + n)
        (hcounts, hbins, hpatches) = matplotlib.pyplot.hist(histdata[:,n],
            bins=numpy.arange(-0.01, 1.01, .01))
        overallcountmax = max(overallcountmax, max(hcounts))
    matplotlib.pyplot.subplot(num_errorfiles, num_datasets,
        p * num_datasets + 1)
    matplotlib.pyplot.title(errorfile)

# label axes in an aesthetically pleasing way
for p in xrange(num_errorfiles):
    matplotlib.pyplot.subplot(num_errorfiles, num_datasets,
        p * num_datasets + 1)
    matplotlib.pyplot.ylabel('# reads')
for p in xrange(num_errorfiles * num_datasets):
    matplotlib.pyplot.subplot(num_errorfiles, num_datasets, p + 1)
    matplotlib.pyplot.xlim(xmin=-0.01, xmax=1.01)
    matplotlib.pyplot.ylim(ymin=0, ymax=overallcountmax)
for p in xrange(num_datasets):
    matplotlib.pyplot.subplot(num_errorfiles, num_datasets,
        (num_errorfiles - 1) * num_datasets + 1 + p)
    matplotlib.pyplot.xlabel('errors per base')
if options.output_file:
    matplotlib.pyplot.savefig(options.output_file)
if not options.nodisplay:
    matplotlib.pyplot.show()

exit(0)
