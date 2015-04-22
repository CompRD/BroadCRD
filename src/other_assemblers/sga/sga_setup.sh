#! /bin/bash -x


if [ $# -ne 1 ]; then
  echo "Function: to download, compile, and install SGA, Google sparsehash, Hoard, Abyss with Boost in a target directory." 1>&2
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


source /broad/software/scripts/useuse
use .git-1.7.10.1
use BamTools
use GCC-4.3

sga_dir=${base_dir}/sga
stuff_dir=${base_dir}/stuff
mkdir ${stuff_dir}
test_dir=${base_dir}/test
mkdir ${test_dir}

sparsehash_head="sparsehash-2.0.2"
sparsehash_file=${sparsehash_head}.tar.gz


#google hash
cd ${base_dir}
wget http://sparsehash.googlecode.com/files/${sparsehash_file}
tar zxf ${sparsehash_file}
rm ${sparsehash_file}
cd ${sparsehash_head}
./configure --prefix=${stuff_dir}
make -j16
make install

#Hoard
cd ${base_dir}
env GIT_SSL_NO_VERIFY=true git clone --recursive http://github.com/emeryberger/Hoard.git
cd Hoard/src
make linux-gcc-x86-64


#abyss
cd ${base_dir}
wget http://www.bcgsc.ca/platform/bioinfo/software/abyss/releases/1.3.5/abyss-1.3.5.tar.gz
tar xzf abyss-1.3.5.tar.gz
rm abyss-1.3.5.tar.gz
cd abyss-1.3.5
wget http://downloads.sourceforge.net/project/boost/boost/1.50.0/boost_1_50_0.tar.bz2
tar jxf boost_1_50_0.tar.bz2
ln -s boost_1_50_0/boost boost
./configure --prefix=${stuff_dir} CPPFLAGS=-I${stuff_dir}/include
make
make install

#sga
cd ${base_dir}
git clone git://github.com/jts/sga.git
cd ${sga_dir}/src
./autogen.sh
./configure --prefix=${stuff_dir} --with-sparsehash=${stuff_dir} --with-hoard=${base_dir}/Hoard/src
make -j16
make install


cd ${test_dir}
wget https://raw.github.com/jts/sga/master/src/examples/sga-ecoli-miseq.sh
wget https://raw.github.com/jts/sga/master/src/examples/sga-human-NA12878.txt
wget https://raw.github.com/jts/sga/master/src/examples/sga-celegans.sh

echo
echo "binary directory:" ${stuff_dir}/bin
echo
echo "official examples directory:" ${test_dir}
echo




