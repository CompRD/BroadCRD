#!/bin/bash
# nanoo
#
# A noo version of nanoc.  Really a complementary program for bulk converting a lot of reads, hopefully a
# bit faster.
#
# neilw - Nov 14, 2014

if [[ "$PYTHONPATH"=="" ]]
then
    export PYTHONPATH=/wga/dev/local/lib/python2.7/site-packages:/wga/dev/neilw/BroadCRD/python
else
    export PYTHONPATH=/wga/dev/local/lib/python2.7/site-packages:/wga/dev/neilw/BroadCRD/python:$PYTHONPATH
fi

export PATH=/wga/dev/local/bin:$PATH

exec `dirname $0`/nanoo.py $*

