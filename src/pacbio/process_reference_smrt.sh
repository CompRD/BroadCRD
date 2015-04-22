#!/usr/bin/env bash

#import a reference to a directory for smrt usage

# the following message can be ignored according to PacBio documentation
#SLF4J: Failed to load class "org.slf4j.impl.StaticLoggerBinder".
#SLF4J: Defaulting to no-operation (NOP) logger implementation
#SLF4J: See http://www.slf4j.org/codes.html#StaticLoggerBinder for further details.

set -e

if [ $# -ne 3 ] ; then
    echo "$0 smrt_ref_dir name fasta" >& 2
    exit 1
fi

tgt=$1
name=$2
ref=$3

setting=/broad/software/free/Linux/redhat_5_x86_64/pkgs/smrtanalysis_2.1.1-centos/install/smrtanalysis-2.1.1.128549/etc/setup.sh

if [ ! -d $tgt ] ; then
    echo "$tgt does not exist" >& 2
    exit 1
fi

if [ ! -f ${ref} ] ; then
    echo "${ref} does not exist" >& 2
    exit 1
fi

if [ ! -f ${setting} ] ; then
    echo "${setting} does not exist" >& 2
    exit 1
fi

source $setting
referenceUploader -c -p ${tgt} -n ${name} -f ${ref}
