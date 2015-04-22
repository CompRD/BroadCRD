#!/usr/bin/env python
#
# R.I.P. W.R.S. - rest in peace, W. Richard Stevens
#

import os
import signal
import resource
import sys

def init_daemon(cwd="/", out_file = None, verbose=True):
    # fork
    try:
        if os.fork(): sys.exit(0)
    except OSError:
        print "failed forking initially"
        raise

    # set session id and make process process group leader
    os.setsid()

    # ignore SIGHUP
    signal.signal( signal.SIGHUP, signal.SIG_IGN )

    # fork again
    try:
        pid = os.fork()
        if pid: os._exit(0)
    except OSError:
        print "failed second fork"
        raise

    # chdir /
    os.chdir(cwd)

    # set up outfile
    if out_file:
        out_file = "{}_{}.log".format( out_file, os.getpid() )
    else:
        out_file = "/dev/null"
    if verbose:
        print "daemon has PID {}, output to {}".format( os.getpid(), out_file )

    # close all FDs and open stdin, stdout, and stderr
    try:
        os.closerange(0,resource.RLIMIT_NOFILE)
        os.open("/dev/null", os.O_RDONLY)
        os.open(out_file, os.O_RDWR|os.O_CREAT|os.O_APPEND )
        os.dup2(1,2)
    except:
        pass



if __name__ == "__main__":
    init_daemon(".","out.txt")
    print "{} {} {}".format( os.getcwd(), os.getpid(), os.getpgrp() )
    import time
    while True:
        time.sleep(30)
