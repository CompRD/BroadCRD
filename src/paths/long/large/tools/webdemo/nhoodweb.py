#!/usr/bin/env python

import BaseHTTPServer as http
from SocketServer import ThreadingMixIn
import time
import re
import urlparse
import threading
import os
import argparse
import ssl
import sys
from BroadCRD.util.NAML import *
from BroadCRD.util.daemonize import *

class MultiThreadedHTTPServer(ThreadingMixIn, http.HTTPServer):
    pass

class SerialNo:
    def __init__(self, serial=0):
        self.serial=serial
        self.lock=threading.Lock()

    def click(self):
        self.lock.acquire(True)
        out=self.serial
        self.serial+=1
        self.lock.release()
        return out


class Server(MultiThreadedHTTPServer):
    def __init__(self, server_dir, handler, config, \
            noticefd = sys.stdout, logfd = sys.stderr ):
        print "init called for Server() with process id {}".format(os.getpid())
        # get config options with defaults
        host=config.get('host', '')
        port=config.get('port', 8000)
        ssl_cert=config.get('ssl_cert', None)
        self.proxy_path=config.get('proxy_path','/')
        # get required config options
        self.input_form=config.get('input_form')
        self.error_form=config.get('error_form')
        self.assets_dir=config.get('assets_dir')
        self.include_dir=config.get('include_dir')
        self.default_image=config.get('default_image', 'assets/broad-logo.jpg')
        self.error_image=config.get('error_image', self.default_image)
        self.javascript=config.get('javascript')
        self.noticefd = noticefd
        self.logfd = logfd
        self.server_dir=server_dir
        self.nthreads = self.count_thread_dirs()

        print "running with {} threads".format(self.nthreads)

        MultiThreadedHTTPServer.__init__(self, (host,port), handler )
        if ssl_cert:
            self.socket=ssl.wrap_socket( self.socket, certfile=ssl_cert, \
                    server_side=True )
        self.template=CachedTemplate(self.input_form, self.include_dir)
        self.error_template=CachedTemplate(self.error_form, self.include_dir)
        self.serialno=SerialNo()
        try:
            self.serve_forever()
        except KeyboardInterrupt:       # change to something more sensible
            pass
        self.server_close()

    def count_thread_dirs(self, retry=True):
        while True:
            x=0
            while True:
                path=os.path.join(self.server_dir, str(x+1))
                if not os.path.isdir( path ): break
                x=x+1

            if x != 0: return x
            elif not retry:
                raise Exception(
                        "no per-thread directories found in {}".format(self.server_dir) )

            # if no directories were found and we're going to retry, 
            # then we sleep first
            print "no per-thread directories yet, sleeping..."
            time.sleep(60)


# someday this should cached and check modification times
# AND do the substitution in a more efficient manner
class CachedTemplate:
    def __init__(self, filename, include_dir):
        self.filename=filename
        self.include_dir=include_dir
        self.reload()

    def reload(self):
        print "loading template from " + self.filename
        self.mtime=os.path.getmtime(self.filename)
        self.templ=open(self.filename,'r').read()

    def maybe_reload(self):
        mtime=os.path.getmtime(self.filename)
        if mtime != self.mtime: self.reload()

    def process_include(self, match):
        filename=match.groups()[0]
        filename=os.path.join(self.include_dir, filename)
        return open(filename,'r').read()

    def trans(self, **kwargs):
        self.maybe_reload()
        output=self.templ
        process = lambda match: self.process_include(match)
        output=re.sub('@include:(\S+)', process, output)
        for key,val in kwargs.items():
            output=re.sub( '@{}'.format(key), val, output )
        return output

class Handler(http.BaseHTTPRequestHandler):
    def __init__( self, request, client_addr, server ):
        http.BaseHTTPRequestHandler.__init__(self, request, client_addr, server)

    def log_write(self, msg): self.server.logfd.write(msg)
    def notice_write(self, msg): self.server.notifyfd.write(msg)

    def log_message(self, fmt, *args):
        fwd=self.headers.getheader('X-Forwarded-For')
        if fwd:
            client=fwd
            proxy='P'
        else:
            client=self.address_string()
            proxy='-'

        start="{} {} - [{}] ".format( client, proxy, self.log_date_time_string() )
        self.log_write(start)
        self.log_write( fmt % args )
        self.log_write("\n")

    def send_html(self, html, code = 200, msg = "OK" ):
        self.send_response(code, msg)
        self.send_header("Content-type", "text/html")
        self.end_headers()
        self.wfile.write(html)
        self.wfile.flush()

    def send_error(self, code, msg):
        self.send_response(code, msg)
        self.send_header("Content-type", "text/plain")
        self.end_headers()
        self.wfile.write("a server failure has occurred: \n" + msg)
        self.wfile.flush()

    def guess_mime(self, filename):
        table={".png":"image/png", ".svg":"image/svg+xml", \
                ".jpg":"image/jpeg", ".js":"text/javascript",
                ".ico":"image/x-icon" }
        (root,ext)=os.path.splitext(filename)
        if ext in table:
            return table[ext]
        else:
            raise Exception("no known mime type for file "+filename)


    def send_raw(self, rawfn):
        mimetype=self.guess_mime(rawfn)
        print "send_raw for rawfile {} mime-type {}".format(rawfn, mimetype)
        try:
            with open(os.path.join(self.server.server_dir,rawfn),'r') as rawfile:
                self.send_response(200, "OK")
                self.send_header("Content-type", mimetype )
                self.end_headers()
                self.wfile.write( rawfile.read() )
            self.wfile.flush()
        except:
            self.send_error(501,'failed to read png file')
            raise

    def turn_around(self, query):
        serial=self.server.serialno.click()
        handler=(serial % self.server.nthreads) + 1
        abbv="{}/request-{}".format(handler,serial)
        head=os.path.join(self.server.server_dir, abbv)
        with open(head+".req",'w') as fd:
            fd.write(query['n'][0])
        pngfile=head+".svg"
        txtfile=head+".txt"
        tries=0
        while tries < 20 and not os.path.isfile(txtfile):
            time.sleep(1)
            tries += 1
        if not os.path.isfile(txtfile): txtfile=None
        if not os.path.isfile(pngfile): pngfile=None
        print "turn around returning... {}".format((pngfile,txtfile))
        return (pngfile,txtfile)

    def response_503fail(self):
        output=self.server.error_template.trans(
                textoutput='',
                imagefile=self.server.error_image,
                proxypath=self.server.proxy_path )
        self.send_html(output, 503, "Service Unavailable")

    def check_server(self):
        pidfile=os.path.join(self.server.server_dir,'pid')
        try:
            pid=open(pidfile,'r').read().rstrip()
            cmdline=open(os.path.join('/proc',pid,'cmdline'),'r').read().split('\x00')
            if os.path.basename(cmdline[0]) != 'NhoodInfo':
                raise Exception('process id {} had cmdline {}'.format(
                    pid, " ".join(cmdline) ) )
        except Exception as e:
            print "SERVERCHECK: {}".format(e.strerror)
            return False

        return True


    def do_GET(self):
        reqc=self.path[1:].split('/')
        request=reqc[0]
        if len(reqc) > 1: subrequest=reqc[1]
        else: subrequest = None
        print "request is: {}".format(request)
        if request == '':
            print "we received a GET for root"
            try:
                # horribly inefficent... all of this.. but it's a start
                output=self.server.template.trans( textoutput='', \
                        value='S=11:90.26M D=1', \
                        imagefile=self.server.default_image,
                        proxypath=self.server.proxy_path,
                        javascript=self.server.javascript)
                self.send_html(output)
                print "we sent the output!"
            except:
                self.send_error(501, "missing template file")
                raise
        elif request[:4]=='req?':
            if not self.check_server():
                self.response_503fail()
                return
            query=urlparse.parse_qs( request[4:] )
            print "query={}".format(query)
            answer=self.turn_around(query)
            if not answer[1]:
                output=self.server.template.trans(textoutput=r"""
                The request to the server timed out.  This is rather
                unusual and has been logged.  Please try again in a
                short while.""", value='S=11:90.26M D=1', \
                imagefile=self.server.default_image,
                proxypath=self.server.proxy_path)
                self.send_html(output)
                print "we told someone to go away"
            else:
                try:
                    text=open(answer[1],'r').read()
                    if answer[0]:
                        pathc=answer[0].split(os.path.sep)
                        imagefile=os.path.join( pathc[-2], pathc[-1] )
                    else: imagefile=self.server.default_image
                    text="<br/>\n".join( text.split("\n") )
                    output=self.server.template.trans( textoutput=text, \
                            value=query['n'][0], \
                            imagefile=imagefile,
                            proxypath=self.server.proxy_path,
                            javascript=self.server.javascript)
                    self.send_html(output)
                except:
                    self.send_error(501, "failed to read file")
                    raise
        elif re.match("\d+", request) and subrequest and \
            re.match("request-(\d+)\.(...)", subrequest ):
            self.send_raw(os.path.join(request,subrequest))
        elif request == 'assets':
            print self.server.assets_dir
            self.send_raw( os.path.join( self.server.assets_dir,  \
                    self.path[len('/assets/'):] ) )
        elif request == 'failwhale':
            # test the 503 response
            self.response_503failure()
        else:
            print "request received: {}".format(request)
            self.send_response(404,"Unknown request")




if __name__=='__main__':
    parser=argparse.ArgumentParser(description='web server bridge for NhoodInfo server mode')
    parser.add_argument('--dry-run', '-n', help='dry-run; don\'t run server', \
            action='store_true' )
    parser.add_argument('--no-daemon', '-D', \
            help='run in the foreground with daemonizing; for testing', \
            action='store_true')
    parser.add_argument('--server_dir', help='NhoodInfo SERVER_DIR', required=True )
    parser.add_argument('--debug_file' )
    parser.add_argument('--log_file' )
    parser.add_argument('--pid_file' )
    parser.add_argument('configfile', help=".naml format configuration file")
    args= parser.parse_args()

    config=NAML(args.configfile)
    # get absolute base directory for the config file; may be used to
    # normalize any other relative directories because we'll daemonize
    # and current working directory will change
    basedir=os.path.dirname(os.path.abspath(args.configfile))
    if args.debug_file: args.debug_file=os.path.abspath(args.debug_file)
    if args.log_file: args.log_file=os.path.abspath(args.log_file)
    if args.pid_file: args.pid_file=os.path.abspath(args.pid_file)

    # normalize all relative paths to be absolute relative the parent
    # directory containing the config file.  Also normalize any other
    # input argument types.
    for f in [ 'ssl_cert', 'assets_dir', 'include_dir', \
            'input_form', 'error_form' ]:
        if not os.path.isabs( config[f] ):
            config[f] = os.path.join( basedir, config[f] )
    config['port']=int(config['port'])

    if not os.path.isfile( config['input_form'] ):
        raise Exception("\nyou specified a non-existent file: {}".format(\
                        config['input_form'] ))

    if not os.path.isdir( args.server_dir ):
        raise Exception("\nyou specified a non-existent directory: {}".format(\
                args.server_dir ))

    # files are a little backwards for historical reasons
    # the HTTP server logs connections to stderr
    # we therefore use stdout for debugging
    if not args.no_daemon:
        init_daemon( args.debug_file, args.log_file, args.pid_file )

    print "arguments as read from configuration file with pathnames normalized:"
    for k,v in config.iteritems(): print "{}: {}".format(k,v)

    if not os.path.isdir( config['assets_dir'] ):
        raise Exception("\nyou specified a non-existent directory: {}".format(\
                config['assets_dir'] ))
    if not os.path.isfile( config['ssl_cert'] ):
        raise Exception("\nyou specified a non-existent ssl_cert: {}".format(\
                config['ssl_cert'] ))
    if not args.dry_run:
        print "starting server..."
        Server(args.server_dir, Handler, config )
