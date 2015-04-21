from QGrindWorklist import *
import os
import time
import datetime
import shlex, subprocess
import signal
import smtplib
from email.mime.text import MIMEText


class QGrindDaemon:
    def __init__(self, defs):
        if defs._debug_level > 0:
             print "daemon startup pid={}, cwd={}".format(os.getpid(), os.getcwd())
             print "configuration:"
             defs.dump()
        self.defs = defs
        self.log_status('START pid={}, user={}, dir_base={}'.format(os.getpid(), os.getenv('USER'), defs._dir_base))
        self.worklist = QGrindWorklist(defs._job_file, defs._next_file, defs._debug_level)
        if defs._debug_level > 2:
            print "joblist:"
            self.worklist.dump_jobs()
        if os.path.exists(self.defs._stop_file):
            os.unlink(self.defs._stop_file)


        self.looper()

    def __del__(self):
        print "daemon dying"


    # write a single line to a status file
    # anything previously there is wiped out -- this is supposed to be a one line file
    def log_status(self,msg):
        try:
            open(self.defs._status_file,'w').write('{} {}_{} {}\n'.format(\
                datetime.datetime.now(), self.defs._hostname, os.getpid(), msg))
        except Exception as e:
            self.defs.log("can't update status file {} exception:\n{}\n".format(self.defs._status_file, e))

    def email_admin( self, text, subj = "qgrind alert" ):
        msg = MIMEText(text)
        from_addr = self.defs._from_email
        to_addr = self.defs._admin_email
        msg['Subject'] = subj
        msg['From'] = from_addr
        msg['To'] = to_addr

        s = smtplib.SMTP('localhost')
        s.sendmail(from_addr, [to_addr], msg.as_string())
        s.quit()

    def looper(self):
        while 1:
            if os.path.exists( self.defs._enable_file ):
                self.defs.log("scanning",2)
                try:
                    (jobno, jobcmd ) = self.worklist.get_job()
                    if jobno >= 0:
                        self.defs.log("running job {}: {}".format(jobno,jobcmd),2)
                        self.run_job( jobno, jobcmd )
                        continue;
                    else:
                        self.defs.log("no work available", 2)
                except QGrindWorklistMissingNext as e:
                    self.defs.log("worklist next job file is missing: {}".format(e.filename))
                except Exception as e:
                    status="FAIL exception: {}".format(e)
                    self.defs.log(status)
                    self.log_status(status)
                    return

                self.log_status('SLEEP')

            elif os.path.exists( self.defs._stop_file ):
                    status="STOPPED"
                    self.defs.log(status)
                    self.log_status(status)
                    os.rename(self.defs._status_file, self.defs._status_file + '.stopped')
                    return
            else:
                self.defs.log("this host is not enabled", 2)
                self.log_status('DISABLED')

            self.defs.log("sleeping for {} seconds".format(self.defs._sleep_time),2)
            time.sleep(self.defs._sleep_time)

    def ret_words(self, ret):
        if ret == None:
            return "bug -- no status!"
        elif ret == 0:
            return "normal exit"
        elif ret > 0:
            return "error exit code {}".format(ret)
        elif ret == -999:
            return "the process was killed for exceeding the time limit of {} seconds".format(self.defs._job_timeout)
        elif ret < 0:
            return "the process was killed by signal {}".format(-ret)

    def run_job(self, jobno, jobcmd):
        self.defs.log("job {} command {}".format(jobno, jobcmd));

        # make a temp directory for the job to run in
        job_dir=os.path.join(self.defs._tmp_dir, "job.{}".format(jobno))

        ## move a pre-existing directory out of the way
        if os.path.exists(job_dir):
            bak=0
            dst_dir = job_dir
            while os.path.exists(dst_dir):
                dst_dir=os.path.join(self.defs._tmp_dir, "job.{}.old{}".format(jobno,bak))
                bak=bak+1
            os.rename(job_dir,dst_dir)

        ## redundant check in case the moving failed; make new directory
        if not os.path.exists(job_dir):
            os.mkdir(job_dir)

        self.log_status("START job {} command {} directory {}".format(jobno, jobcmd, job_dir))

        job_dir_stderr=os.path.join(job_dir,'stderr')
        job_dir_stdout=os.path.join(job_dir,'stdout')
        job_dir_status=os.path.join(job_dir, 'status')
        job_dir_info=os.path.join(job_dir, 'info')

        # run the job
        start_time = time.time()
        argv=shlex.split(jobcmd)
        with open(job_dir_stdout,'w') as stdout:
            with open(job_dir_stderr,'w') as stderr:
                with open(job_dir_status,'w') as status_file:
                        try:
                            with open(job_dir_info,'w') as info_file:
                                info_file.write('host: {}\n'.format(self.defs._hostname))
                                info_file.write('command-line: {}\n'.format(argv))
                                info_file.write('pid: {}\n'.format(os.getpid()))
                            self.defs.log('about to start process {}'.format(argv), 2)
                            proc = subprocess.Popen(argv, stdout=stdout, stderr=stderr, \
                                                  cwd=job_dir, preexec_fn=os.setsid)
                            self.log_status('RUN job {} command {} directory {} pid {}'.format(jobno, jobcmd, job_dir, proc.pid))
                            ret = None
                            while ret == None:
                                self.defs.log('about to poll process', 2)
                                ret = proc.poll()
                                self.defs.log('polled process and found ret={}'.format(ret), 2)
                                if ret == None and time.time() - start_time > self.defs._job_timeout:
                                    self.defs.log('process overdue... sending SIGTERM', 2)
                                    os.killpg( proc.pid, signal.SIGTERM )
                                    time.sleep(15)
                                    self.defs.log('process overdue... sending SIGKILL', 2)
                                    os.killpg( proc.pid, signal.SIGKILL )
                                    ret = -999
                                elif ret == None:
                                    self.defs.log('process not finished...sleeping for {} seconds'.format(self.defs._job_poll_interval), 2)
                                    time.sleep(self.defs._job_poll_interval)
                            status_file.write("{}\n".format(ret))
                            if ret >= 0:
                                status_file.write("process exited\n")
                            elif ret == -999:
                                status_file.write("we killed the process for exceeding the time limit of {} seconds\n".format(self.defs._job_timeout))
                            else:
                                status_file.write("the job was killed by signal {}\n".format(-ret))

                            self.log_status("COMPLETE job {} status {}".format(jobno, self.ret_words(ret)))
                        except Exception as e:
                            status_file.write("-888\nAn exception was thrown:\n{}\n".format(e))
                            raise       # re-throw the exception so we can exit
