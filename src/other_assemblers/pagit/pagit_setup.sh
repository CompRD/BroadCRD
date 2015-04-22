#! /bin/bash -x


if [ $# -ne 1 ]; then
  echo "Function: to download, compile, and install Cortex in a target directory." 1>&2
  echo "Usage: `basename $0` target_directory" 1>&2
  exit 1
fi

#starting_dir=${PWD}
#base_dir=${starting_dir}/sga_standalone

base_dir=$1
if [ -e $base_dir ]; then
  echo "$base_dir exists" 1>&2
  exit 1
fi

mkdir $base_dir
if [ $? -ne 0 ]; then
  echo "Error making directory $base_dir" 1>&2
fi




cd ${base_dir}


pagit_package=PAGIT.V1.64bit.tgz

wget ftp://ftp.sanger.ac.uk/pub4/resources/software/pagit/${pagit_package}

tar xzf ${pagit_package}

bash installme.sh

source ${base_dir}/PAGIT/sourceme.pagit

cd ${base_dir}/PAGIT/exampleTestset

bash ./dotestrun.sh

cd ${base_dir}




