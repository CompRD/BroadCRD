#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

import optparse
import os
import re
import sys

oparser = optparse.OptionParser(usage='%prog [options] first.badcoverage '
    '[second.badcoverage] ...')
oparser.add_option('--field_regex', help='fields to include', default=None)
oparser.add_option('--transpose', help='swap rows and columns in output',
    action="store_true", default=False)
(options, args) = oparser.parse_args()
    
bc_files = args[0:]
bc_files.sort()

bc_reports = {}
all_keys = set()
all_keys_list = []
motifs = {}
for bcf in bc_files:
    with open(bcf, 'r') as bc_data:
        bc_reports[bcf] = {}
        for bc_line in bc_data:
            data_match = re.match('^(.+)\s=\s(.+)', bc_line)
            if data_match:
                key = data_match.group(1)
                if not options.field_regex or re.match(options.field_regex, 
                    key):
                    value = data_match.group(2)
                    if key in bc_reports[bcf]:
                        print >>sys.stderr, ('Error - {0} appears twice in {1}'
                            .format(key, bcf))
                        exit(1)
                    bc_reports[bcf][key] = value
                    if key not in all_keys:
                        all_keys.add(key)
                        all_keys_list.append(key)
                    if re.match('^MOTIF\d+$', key):
                        if key not in motifs:
                            motifs[key] = value
                        elif motifs[key] != value:
                            print >>sys.stderr, ('Error - {0} is defined as {1} '
                                ' in {2}, but previously was defined as {3}'
                                .format(key, value, bcf, motifs[key]))
                            exit(1)
                        

if not options.transpose:
    # print header
    print 'BC FILE\t',
    for k in all_keys_list:
        print k + '\t',
    print
    
    # print data
    for bcf in bc_files:
        print bcf + '\t',
        for k in all_keys_list:
            if k in bc_reports[bcf]:
                print bc_reports[bcf][k],
            print '\t',
        print
else:
    # print header
    print 'FIELD\t',
    for bcf in bc_files:
        print bcf + '\t',
    print
    
    # print data
    for k in all_keys_list:
        print k + '\t',
        for bcf in bc_files:
            if k in bc_reports[bcf]:
                print bc_reports[bcf][k],
            print '\t',
        print
    
exit(0)