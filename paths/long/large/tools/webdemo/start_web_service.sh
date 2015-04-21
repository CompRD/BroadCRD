#!/bin/bash

if (( $# < 2 ))
then
    echo 1>&2 "usage: $0 webdemo_code_dir assembly_fin_dir"
    exit 1
fi

#toolsdir=/wga/dev/neilw/BroadCRD/paths/long/large/tools/webdemo
toolsdir=$1
if [[ ! -e $toolsdir/nhoodweb.py ]]
then
    echo 1>&2 "there is no nhoodweb.py in $toolsdir"
    exit 2
fi

if [[ -e ./nhoodweb.py ]]
then
    echo 1>&2 "do not run from the source code directory"
    echo 1>&2 "maybe run from /wga/scr4/webdemo/deployment"
    exit 3
fi

#assemblydir=/wga/scr4/jaffe/GapToy/49962/a.fin2
assemblydir=$2
if [[ ! -e $assemblydir/a.lines.stats ]]
then
    echo 1>&2 "$assemblydir does not appear to be an a.fin directory"
    exit 1
fi

eval `/broad/software/dotkit/init -b`
use .gcc-4.7.2
use Graphviz

export PATH=/wga/dev/local/bin:$PATH

bindir=/wga/scr4/webdemo/BroadCRD/bin
serverdir=/tmp/nhood.dexter
mailto=neilw@broadinstitute.org


if [[ "$PYTHONPATH" == "" ]]
then
    export PYTHONPATH=/wga/dev/neilw/BroadCRD/python
else
    export PYTHONPATH=/wga/dev/neilw/BroadCRD/python:$PYTHONPATH
fi

[[ -d $serverdir ]] || mkdir $serverdir

rm -f stdout access
ln -s nhoodweb.$$.access access
ln -s nhoodweb.$$.stdout stdout
$toolsdir/nhoodweb.py --debug_file nhoodweb.$$.stdout --log_file nhoodweb.$$.access \
    --pid_file nhoodweb.pid --server_dir $serverdir config.naml

basedir=`pwd`
(
cd /
trap "" HUP TERM INT
let runs=0
let maxruns=5
let tsleep=30
while (( $runs < $maxruns ))
do
    let runs=runs+1
    log=$basedir/nhoodinfo.$$-$runs.log
    $bindir/NhoodInfo DIR=$assemblydir SERVER_DIR=$serverdir 1>$log 2>&1
    mail -s "NhoodInfo failed" $mailto <<EOF
NhoodInfo has failed and is being restarted.
The log file is in $log
We restart $maxruns times
This is restart #$runs
EOF
    sleep $tsleep
done &
) 0<&- 1>/dev/null 2>&1 &
# This (below) is a KLUDGE -- it's hard to get the PID associated with a
# subshell.  It *seems* to always be $! plus one, but that may be a race
# condition.  I'm not sure that it's guaranteed to be.
let pid=$!+1
echo $pid > nhoodinfo.pid
