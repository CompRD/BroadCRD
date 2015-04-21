
# starting and stopping the web service

Starting the web service:

1. log in as dexter
2. execute:

        ./start_web_service.sh

Stopping the web service:

1. log in as dexter
2. execute:

        kill -9 `cat nhoodinfo.pid` `cat nhoodweb.pid`

# output files

`nhoodweb.<process id>.access` - access logs; essentially stderr from
  the nhoodweb.py process (for historical reasons having to do with the
  Python web service).

`nhoodweb.<process id>.stdout` - stdout, but look here for tracebacks
  from the web service itself.  Again, you would think that this would
  be stderr, but no.

`nhoodinfo.<process id>-<restart>.log` - this is the log file for the
  NhoodInfo server.

Note that <process id> is the process id of the parent script which no
longer exists.  For the process ids of the daemons, check the `*.pid`
files.

<restart> indicates which instance of NhoodInfo this is the log for.
The shell script will loop and restart NhoodInfo if it crashes, causing
-2, -3, ... instances.

# configuration file

`config.naml`   - sets various parameters for a production setting,
including the outward facing web address (listed as `proxy_path`).  Also
indicates things such as port number, ssl certificate file, and template
file.  NOTE: the python script uses the location of this configuration
file as the assumed location of any relative pathnames specified for the
SSL cert, template file, etc.  So when `input_form` is specified as just
"filename," it's assumed to be in the same directory as this
configuration file.  The name of this configuration file is specified on
the command line to nhoodweb.py.

`testing.naml`  - same as `config.naml,` but altered for local testing
on a machine (rather than outward facing).  The main issue here is to
change the port not to conflict with the production port, if you're on
the same machine and to change the proxypath not to point to the
external web server.


# template file

The template file is essentially HTML (with embedded javascript, as you
wish) that will be relayed to the end user as the root document.  There
are a few @tokens that will be replaced (e.g. `@proxypath`) by some
other text.  They are as follows:

@proxypath - the value of the proxy_path (sorry for the underscore)
configuration parameter.

@imagefile - the graph image file presented as output (.png)

@textoutput - any text output from the NhoodInfo process (.txt)

@value - the value chosen by the user on the previous round (or the
default value)


