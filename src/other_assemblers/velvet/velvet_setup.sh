#! /bin/bash -x


if [ $# -ne 1 ]; then
  echo "Function: to download, compile, and install velvet in a target directory." 1>&2
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


#source /broad/software/scripts/useuse
#use .git-1.7.10.1
#use BamTools
#use GCC-4.4

#use


cd ${base_dir}

vo_base=VelvetOptimiser-2.2.5
vo_package=${vo_base}.tar.gz
wget http://www.vicbioinformatics.com/${vo_package}

tar xzf ${vo_package}


cd ${base_dir}

velvet_base=velvet_1.2.10
velvet_package=${velvet_base}.tgz

wget http://www.ebi.ac.uk/~zerbino/velvet/${velvet_package}

tar zxf ${velvet_package}

cd ${velvet_base}

make -j32 'LONGSEQUENCES=1' 'BUNDLEDZLIB=1' 'OPENMP=1' 'BIGASSEMBLY=1' 'MAXKMERLENGTH=161' 'CATEGORIES=4'

cd tests
pwd 
./run-tests.sh



