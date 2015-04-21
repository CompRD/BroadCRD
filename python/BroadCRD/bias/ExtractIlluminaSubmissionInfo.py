#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2012) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# This program extracts sequence submission info from BadCoverage output
# files and the Picard filesystem and outputs the relevant data in a
# spreadsheet.

import os
import re
import shlex
import subprocess
import sys

def extract_fields(params_file):
    output = {'PROJECT':None, 'WORK_REQUEST_ID':None, 'FLOWCELL_BARCODE':None,
        'LANE':None, 'LIBRARY_NAME':None, 'SAMPLE_BARCODE':None,
        'SAMPLE_ALIAS':None}
    with open(params_file, 'r') as params:
        for param_value in params:
            (param, value) = param_value.rstrip().split('=', 1)
            if param in output:
                output[param] = value
    return output

def main():
    badcov_file = sys.argv[1]

    with open(badcov_file, 'r') as badcov_data:
        for badcov_line in badcov_data:
            pu_match = re.match('^production units = (.*)$', badcov_line)
            if pu_match:
                production_units = pu_match.group(1).split(', ')

    field_names = None
    for punit in production_units:
        punit_elements = punit.split('.')
        flowcell = punit_elements[0]
        lane = punit_elements[1]
        barcode = None
        if len(punit_elements) == 3:
            barcode = punit_elements[2]
        elif len(punit_elements) > 3:
            print >>sys.stderr, 'failed to comprehend ' + punit
            exit(1)
        flowcell = flowcell[0:-6]
        picard_path = os.path.join('/seq/picard', flowcell)

        analysis_dir = [a for a in os.listdir(picard_path) if
            os.path.isdir(os.path.join(picard_path, a))]
        analysis_dir.sort()
        picard_path = os.path.join(picard_path, analysis_dir[-1], lane)

        library_dir = [l for l in os.listdir(picard_path) if
            os.path.isdir(os.path.join(picard_path, l))]

        if len(library_dir) > 1 and barcode:
            pfile = subprocess.check_output(shlex.split('find ' +
                picard_path + ' -name params.txt ' +
                '-exec grep -q \'MOLECULAR_BARCODE_SEQUENCE=' +
                barcode + '\' {} \\; -print'))
            library_dir = pfile.split()
            library_dir[0] = os.path.dirname(library_dir[0])

        if len(library_dir) != 1:
            print >>sys.stderr, ('failed to find a unique library dir for '
                + punit)
            exit(1)

        picard_path = os.path.join(picard_path, library_dir[0])

        seq_info = extract_fields(os.path.join(picard_path, 'params.txt'))

        if not field_names:
            field_names = seq_info.keys()
            field_names.sort()
            print '\t'.join(field_names)

        print '\t'.join([seq_info[f] for f in field_names])



    return 0

if __name__ == '__main__':
    exit(main())
