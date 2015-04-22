#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# Find SMRT cells with samples matching a specified regular expression.

# Original author: Michael G. Ross <mgross@broadinstitute.org>


import glob
import os
import re
import sys

def check_sample(metadata_xml, sample_name_regexp):
    with open(metadata_xml, 'r') as metadata_file:
        metadata = metadata_file.read()
    
    sample_regexp = re.compile('<Sample>.*<Name>(.*)</Name>.*</Sample>',
        re.MULTILINE)
    
    sample_name = sample_regexp.search(metadata).group(1)
    return re.match(sample_name_regexp, sample_name) != None

def check_smrtpack(cell_dir, sample_name_regexp):
    smrtcells = os.listdir(cell_dir)
    for s in smrtcells:
        full_s = os.path.join(cell_dir, s)
        metadata_xmls = glob.glob(os.path.join(full_s, '*metadata.xml'))
        for m in metadata_xmls:
            if check_sample(os.path.join(full_s, m), sample_name_regexp):
                return True
    return False

def is_smrtpack(cell_dir):
    if not os.path.isdir(cell_dir):
        return False
        
    subdirs = os.listdir(cell_dir)

    if not all(map(lambda x : re.match('[A-Z]\d+_\d+', x), subdirs)):
        return False
    else:
        return True
    
def main(sample_name_regexp):
    smrtcell_root = '/seq/pacbio_results'
    
    dirs = os.listdir(smrtcell_root)
    desired_packs = set()
    
    for d in dirs:
        full_d = os.path.join(smrtcell_root, d)
        if is_smrtpack(full_d) and check_smrtpack(full_d, sample_name_regexp):
            desired_packs.add(full_d)

    out_packs = list(desired_packs)
    out_packs.sort()
    
    for o in out_packs:
        print o
    
    return 0

if __name__ == '__main__':
    exit(main(sys.argv[1]))
