from FileLock import *
import time

class QGrindWorklistError(Exception):
    def __init__(self, msg):
        Exception.__init__(self,msg)

class QGrindWorklistMissingNext(Exception):
    def __init__(self,filename):
        self.filename = filename

class QGrindWorklist:
    def __init__(self, jobfile, nextfile, debuglevel = 0 ):
        self.jobfile = jobfile
        self.mtime = 0.0
        self.nextfile = nextfile
        self.debuglevel = debuglevel

    # conditionally refresh jobs -- if the jobs file modification time is greater than
    # the last time we read it, we read it.  There's no protection here, so we
    # intend for this to really only happen when "get_job" is called and locking occurs.
    def cond_refresh_jobs(self):
        try:
            with open(self.jobfile,'r') as fd:
               mtime = os.fstat( fd.fileno() ).st_mtime 
               if mtime > self.mtime:
                   self.jobs = map( str.strip, fd.readlines() )
                   self.mtime = mtime
                   if self.debuglevel > 0: print "we re-read the jobs file!"
        except IOError:
            raise QGrindWorklistError('I/O failure reading worklist {}'.format(self.jobfile) )


    def dump_jobs(self):
        self.cond_refresh_jobs()          # no locking here... could be a race if you modify the file
        for job in self.jobs:
            print "\t",job

    # return (jobno, cmdline) with jobno = -1 if no work is available
    # otherwise it "claims" a job in that the _next_file gets locked and incremented
    def get_job(self):
        with FileLock(self.nextfile,60,1,True) as lock:
            try:
                with open( self.nextfile, 'r') as job:
                    nextjob = int( job.readline().strip() )
            except IOError:
                raise QGrindWorklistMissingNext(self.jobfile)

            self.cond_refresh_jobs()

            if nextjob < len(self.jobs):
                try:
                    with open( self.nextfile, 'w' ) as job:
                        job.write('{}'.format(nextjob+1))
                except IOError:
                    raise QGrindWorklistError('I/O failure writing {} to next job list {}'.format( nextjob+1, self.nextfile) )
                return (nextjob, self.jobs[nextjob] )

        return (-1,"")
