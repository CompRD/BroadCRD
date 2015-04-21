import socket
import datetime
import os

#class QGrindDefsError(Exception):
#    def __init__(self,msg):
#        print msg
#        Exception.__init__(self,msg)
#    def __new__(cls, *args, **kw):
#        raise QGrindDefsError("QGrindDefs is a static options class")

class QGrindDefs:

    @staticmethod
    def dump():
        for name in QGrindDefs.__dict__.keys():
           if name[0] == '_' and name[:2] != '__':
               print "\t{}={}".format(name,QGrindDefs.__dict__[name])

    @staticmethod
    def log(msg, level=0):
        if level <= QGrindDefs._debug_level:
            print "{} {}: {}".format(datetime.datetime.now(), QGrindDefs._hostname, msg)

    @staticmethod
    def status_file(hostname,pid):
        return '{}/{}_{}.status'.format(QGrindDefs._dir_base,hostname,pid)

    @staticmethod
    def enable_file(hostname):
        return '{}/{}.enable'.format(QGrindDefs._dir_base,hostname)

    @staticmethod
    def stop_file(hostname):
        return '{}/{}.stop'.format(QGrindDefs._dir_base,hostname)


    # files/dirs shared between instances
    _dir_base = '/wga/scr4/qgrind'

    _job_file = '{}/jobfile.txt'.format(_dir_base)
    _next_file = '{}/nextjob.txt'.format(_dir_base)
    _tmp_dir = '{}/tmp'.format(_dir_base)

    _admin_email = 'neilw@broadinstitute.org'
    _from_email = 'dexter@broadinstitute.org'

    # host-specific files (ignoring .broadinstitute.org part of hostname)
    _hostname = socket.gethostname().partition(".")[0]
    _pid = os.getpid()
    _status_file = '{}/{}_{}.status'.format(_dir_base,_hostname,_pid)
    _enable_file = '{}/{}.enable'.format(_dir_base,_hostname)
    _stop_file = '{}/{}.stop'.format(_dir_base,_hostname)

    # polling interval for a running job
    _job_poll_interval = 2
    # sleep time when daemon is dormant
    _sleep_time = 10
    # timeout length for a running job
    _job_timeout = 1800
    # debugging
    _debug_level = 2


    # all hosts -- not definitive; could run on other hosts
    _all_hosts = [ 'crd4', 'crd6', 'crd8', 'crd9', 'crd10', 'crd11', 'crd14', 'crd15' ]
