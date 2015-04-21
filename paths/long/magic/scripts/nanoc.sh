#!/bin/bash
#eval `/broad/software/dotkit/init -b`
#use .hdf5-1.8.9
if [[ "$PYTHONPATH"=="" ]]
then
    export PYTHONPATH=/wga/dev/local/lib/python2.7/site-packages:/wga/dev/neilw/BroadCRD/python
else
    export PYTHONPATH=/wga/dev/local/lib/python2.7/site-packages:/wga/dev/neilw/BroadCRD/python:$PYTHONPATH
fi

export PATH=/wga/dev/local/bin:$PATH
exec `dirname $0`/nanoc.py $*
