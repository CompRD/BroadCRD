#!/usr/bin/env python
#
# qgrind - grind through a queue of jobs
#
# neilw@broadinstitute.org
#

import os
import sys
import argparse
import glob
from QGrindDefs import *
from QGrindDaemon import *
from daemonize import init_daemon

# this should move into the worklist class next week
def set_next_file(next_file, value):
        print('setting next job file {} to {}'.format(next_file, value))
        try:
            with FileLock(next_file,60,0.1,False) as lock:
                with open( next_file, 'w') as fd:
                    fd.write('{}\n'.format(value))
        except IOError as e:
            print('I/O error writing file {}:\n{}'.format(next_file,e))
            raise
        except FileLockException:
            print('failed to gain a lock on the file; exiting...')
            raise

def enable_disable(enable, disable, hosts):
    if not hosts: hosts = QGrindDefs._all_hosts
    for host in hosts:
        enable_file = QGrindDefs.enable_file(host)
        if enable:
            print 'enabling {}'.format(host),
            os.system("/wga/dev/local/bin/crds -s -R " + host + " qgrind  >& /dev/null")
            if not os.path.exists(enable_file):
                with open(enable_file, 'w'): pass
                print ''
            else:
                print '(already enabled)'
        elif disable:
            print'disabling {}'.format(host),
            os.system("/wga/dev/local/bin/crds -s -r " + host + "  >& /dev/null")
            if os.path.exists(enable_file):
                os.unlink(enable_file)
                print ''
            else:
                print '(already disabled)'

def stop(hosts):
    if not hosts: hosts = QGrindDefs._all_hosts
    for host in hosts:
        stop_file = QGrindDefs.stop_file(host)
	print 'stopping {}'.format(host),
	os.system("/wga/dev/local/bin/crds -s -r " + host + "  >& /dev/null")
	if not os.path.exists(stop_file):
		with open(stop_file, 'w'): pass
                print ''
	else:
                print '(already stopped/stopping)'

class QGStatus:
    def __init__(self,line):
        import re
        m = re.match(r'(\d+)-(\d+)-(\d+) (\d+):(\d+):(\d+).(\d+) (\S+) (\S+) *(.*)$', line)
        if not m: raise Exception("line doesn't match a status line: "+line)
        (m_day, m_month, m_year, m_hour, m_min, m_sec, m_micro, m_hostpid, m_state, m_info) = m.groups()
        self.time = datetime.datetime(int(m_day), int(m_month), int(m_year), int(m_hour), int(m_min), int(m_sec), int(m_micro))
        self.hostpid = m_hostpid
        self.state = m_state
        self.info = m_info

    def __repr__(self):
        return "{} {} {} {}".format( self.time, self.hostpid, self.state, self.info)

def sort_len_order(x,y):
    if len(x) < len(y): return -1
    elif len(x) > len(y): return +1
    elif x < y: return -1
    elif x > y: return +1
    else: return 0


def print_statuses(all=False):
    import os
    try:
        (rows, columns) = map(int, os.popen('stty size', 'r').read().split())
    except:
        columns = 80
    statuses = glob.glob(os.path.join(QGrindDefs._dir_base, "*.status"))
    statuses.sort(sort_len_order)
    saw_pending = False
    for status_file in statuses:
        with open(status_file,'r') as f:
            stat=QGStatus(f.readline().strip())
            host = "_".join( stat.hostpid.split('_')[:-1] )     # all but last underscore
            host_enable = os.path.join( QGrindDefs._dir_base, "{}.enable".format( host) )
            if stat.state == "DISABLED" or os.path.exists( host_enable ): pending = " "
            else:
                pending = "*"
                saw_pending = True
            line="{:12} ({:3.0f} min) {:1}{:8} {}".format(stat.hostpid, \
                    (datetime.datetime.now()-stat.time).total_seconds()/60, \
                    pending,stat.state, stat.info)
            if all:
                print line
            else:
                print line[:columns]

    if saw_pending: print "* indicates host disabled pending completion of current job"


def main( argv = [__name__]  ):
    default_log_base = "{}".format(QGrindDefs._hostname)

    parser=argparse.ArgumentParser(description='Overall queue list manager program and daemon',\
                                   epilog='Default host list for enable and disable: {}'.format(QGrindDefs._all_hosts) )
    parser.add_argument('--daemon','-d',help='cause the program to start a daemon on the current machine to process jobs',action='store_true')
    parser.add_argument('--next', '-n', help='(re)set the next job file to this value',type=int,default=-1)
    parser.add_argument('--enable', help='enable a list of hosts and exit', action='store_true')
    parser.add_argument('--disable', help='disable a list of hosts and exit', action='store_true')
    parser.add_argument('--stop', help='stop a list of hosts and exit', action='store_true')
    parser.add_argument('--dump', help='dump the configuration and exit', action='store_true')
    parser.add_argument('--status', help='output status of each daemon from .status files (full line)', action='store_true')
    parser.add_argument('-s', help='output status of each daemon from .status files', action='store_true')
    parser.add_argument('hostname', nargs='*', default=[])
    args=parser.parse_args(argv[1:])

    if args.hostname and not args.enable and not args.disable and not args.stop:
        print("you must specify --enable or --disable if you list hosts")
        return 1

    if args.enable and args.disable:
        print("you must specify EITHER --enable or --disable, not both")
        return 1

    if args.next >= 0:
        set_next_file(QGrindDefs._next_file, args.next)

    if args.enable or args.disable:
        return enable_disable(args.enable, args.disable, args.hostname)
    elif args.dump:
        QGrindDefs.dump()
    elif args.status or args.s:
        print_statuses(args.status)
    elif args.stop:
        return stop(args.hostname)
    elif args.daemon:
        print "about to start daemon -- process will run in the background"
        init_daemon(QGrindDefs._dir_base,default_log_base)
        # this is a kludge to update the pid after we fork -- awful
        QGrindDefs._status_file = QGrindDefs.status_file( QGrindDefs._hostname, os.getpid() )
        QGrindDaemon(QGrindDefs)
    else:
        print 'use -h for help'

if __name__ == "__main__":
    sys.exit(main(sys.argv))
