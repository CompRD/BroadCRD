#!/usr/bin/env python

###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2010) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# Execute a Picard program using the currently released version
# and automatically setting the TMP_DIR to /tmp/<username>/ (which can
# be over-riden by specifying a 'tmp_dir=/alt/tmp/dir' keyword-style argument).
# There is no need to append the path to the JAR file, for example:
# picard_exec('SortSam', 'I=unsorted.bam', 'O=sorted.bam', 'SO=coordinate').

# Original author: Michael G. Ross <mgross@broadinstitute.org>

import getpass
import os
import subprocess

def picard_exec(command, *args, **keywords):
    if 'tmp_dir' in keywords:
        my_temp_dir = keywords[tmp_dir]
    else:
        my_temp_dir = os.path.join('/tmp/', getpass.getuser())
    if not os.path.exists(my_temp_dir):
        os.mkdir(my_temp_dir)
    elif not os.path.isdir(my_temp_dir):
        raise IOError, (my_temp_dir + ' exists, but is not a directory.')

    subprocess.check_call(('java', '-jar',
        os.path.join('/seq/software/picard/current/bin',
                   command + '.jar'),
        'TMP_DIR=' + my_temp_dir)
        + args)
    return
