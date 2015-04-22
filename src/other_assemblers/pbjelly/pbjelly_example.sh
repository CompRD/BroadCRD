#! /bin/bash -x


if [ $# -ne 1 ]; then
  echo "Function: run PBJelly installed in a specified directory." 1>&2
  echo "Usage: `basename $0` pbj_dir " 1>&2
  exit 1
fi


#starting_dir=${PWD}
#base_dir=${starting_dir}/sga_standalone

source /broad/software/scripts/useuse
use .smrtanalysis-1.3.1-centos

pbj_dir=$1

export JELLYPATH=${pbj_dir}
export SEYMOUR_HOME=/broad/software/free/Linux/redhat_5_x86_64/pkgs/smrtanalysis_1.3.1-centos

source ${SEYMOUR_HOME}/etc/setup.sh

#source ${pbj_dir}/exportPaths.sh

export PATH=${JELLYPATH}:$PATH

export test_dir=`pwd`/bayo.test

rm -f ${test_dir}

cp -r ${pbj_dir}/lambdaExample ${test_dir}

cd ${test_dir}


sed -e "s|/__PATH__/_TO_/PBJelly/lambdaExample/|${PWD}/|" Protocol.xml > tmp
mv tmp Protocol.xml

Jelly.py setup Protocol.xml
Jelly.py mapping Protocol.xml -x " -nproc 48 "
Jelly.py support Protocol.xml
Jelly.py extraction Protocol.xml
Jelly.py assembly Protocol.xml -x " --nproc=48 "
Jelly.py output Protocol.xml


