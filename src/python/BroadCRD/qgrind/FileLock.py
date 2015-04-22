#
# modified from:
# http://www.evanfosmark.com/2009/01/cross-platform-file-locking-support-in-python/
#
import os
import time
import errno
import socket
 
class FileLockException(Exception):
    pass
 
class FileLock(object):
    """ A file locking mechanism that has context-manager support so 
        you can use it in a with statement. This should be relatively cross
        compatible as it doesn't rely on msvcrt or fcntl for the locking.
    """
 
    def __init__(self, file_name, timeout=60, delay=.05, debug = False):
        """ Prepare the file locker. Specify the file to lock and optionally
            the maximum timeout and the delay between each attempt to lock.
        """
        self.is_locked = False
        if os.path.dirname(file_name):
            self.lockfile = "%s.lock" % file_name
        else:
            self.lockfile = os.path.join(os.getcwd(), "%s.lock" % file_name)
        self.file_name = file_name
        self.timeout = timeout
        self.delay = delay
        self.hostname = socket.gethostname()
        self.debug = debug
 
 
    def acquire(self):
        """ Acquire the lock, if possible. If the lock is in use, it check again
            every `wait` seconds. It does this until it either gets the lock or
            exceeds `timeout` number of seconds, in which case it throws 
            an exception.
        """
        start_time = time.time()
        while True:
            try:
                self.fd = os.open(self.lockfile, os.O_CREAT|os.O_EXCL|os.O_RDWR)
                os.write(self.fd, "{}\n".format(self.hostname))
                os.fsync(self.fd)
                if self.debug: print "{} acquired lock {}".format( self.hostname, self.lockfile )
                break;
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise 
                if (time.time() - start_time) >= self.timeout:
                    raise FileLockException("Timeout occured.")
                if self.debug: print "{} failed to acquire lock {}; sleeping".format( self.hostname, self.lockfile )
                time.sleep(self.delay)
        self.is_locked = True
 
 
    def release(self):
        """ Get rid of the lock by deleting the lockfile. 
            When working in a `with` statement, this gets automatically 
            called at the end.
        """
        if self.is_locked:
            os.close(self.fd)
            os.unlink(self.lockfile)
            self.is_locked = False
            if self.debug: print "{} released lock {}".format( self.hostname, self.lockfile )
 
 
    def __enter__(self):
        """ Activated when used in the with statement. 
            Should automatically acquire a lock to be used in the with block.
        """
        if not self.is_locked:
            self.acquire()
        return self
 
 
    def __exit__(self, type, value, traceback):
        """ Activated at the end of the with statement.
            It automatically releases the lock if it isn't locked.
        """
        if self.is_locked:
            self.release()
 
 
    def __del__(self):
        """ Make sure that the FileLock instance doesn't leave a lockfile
            lying around.
        """
        self.release()
