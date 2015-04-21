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


source /broad/software/scripts/useuse
use .git-1.7.10.1
use BamTools
use GCC-4.4

use

stuff_dir=${base_dir}/stuff
mkdir ${stuff_dir}
#test_dir=${base_dir}/test
#mkdir ${test_dir}



cd ${base_dir}
zlib_head=zlib-1.2.8
zlib_package=${zlib_head}.tar.gz
wget http://zlib.net/${zlib_package}
tar xzf ${zlib_package}
cd ${zlib_head}
./configure --prefix=${stuff_dir}
make
make install

export C_INCLUDE_PATH=${stuff_dir}/include:${C_INCLUDE_PATH}
export CPLUS_INCLUDE_PATH=${stuff_dir}/include:${CPLUS_INCLUDE_PATH}
export LD_LIBRARY_PATH=${stuff_dir}/lib:${LD_LIBRARY_PATH}

echo ${C_INCLUDE_PATH}
echo ${CPATH_INCLUDE_PATH}
echo ${LD_LIBRARY_PATH}

cd ${base_dir}
cortex_var_head=CORTEX_release_v1.0.5.15
cortex_var_head=CORTEX_release_v1.0.5.19
cv_package=${cortex_var_head}.tgz

wget "http://downloads.sourceforge.net/project/cortexassembler/cortex_var/latest/CORTEX_release_v1.0.5.19.tgz?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fcortexassembler%2Ffiles%2F&ts=1375119994&use_mirror=superb-dca2" -O ${cv_package}

#wget "http://downloads.sourceforge.net/project/cortexassembler/cortex_var/latest/CORTEX_release_v1.0.5.15.tgz?r=http%3A%2F%2Fcortexassembler.sourceforge.net%2Findex_cortex_var.html&ts=1367873131&use_mirror=iweb" -O ${cv_package}
#wget "http://downloads.sourceforge.net/project/cortexassembler/cortex_var/previous_releases/CORTEX_release_v1.0.5.15.tgz?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fcortexassembler%2Ffiles%2Fcortex_var%2Fprevious_releases%2F&ts=1375123470&use_mirror=iweb" -O ${cv_package}


tar xzf ${cv_package}

cd ${cortex_var_head}

bash install.sh

make cortex_var
make MAXK=63 cortex_var

cd  ${base_dir}

wget http://www.well.ox.ac.uk/bioinformatics/Software/Stampy-latest.tgz
tar zxf Stampy-latest.tgz


wget "http://downloads.sourceforge.net/project/vcftools/vcftools_0.1.9.tar.gz?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fvcftools%2Ffiles%2F&ts=1367951825&use_mirror=iweb" -O vcftools_0.1.9.tar.gz
tar xzf vcftools_0.1.9.tar.gz

wget "http://downloads.sourceforge.net/project/vcftools/vcftools_0.1.8a.tar.gz?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fvcftools%2Ffiles%2F%3Fsource%3Dnavbar&ts=1368026057&use_mirror=hivelocity" -O vcftools_0.1.8a.tar.gz
tar xzf vcftools_0.1.8a.tar.gz




#sparsehash_head="sparsehash-2.0.2"
#sparsehash_file=${sparsehash_head}.tar.gz


#google hash
#cd ${base_dir}
#wget http://sparsehash.googlecode.com/files/${sparsehash_file}
#tar zxf ${sparsehash_file}
#rm ${sparsehash_file}
#cd ${sparsehash_head}
#./configure --prefix=${stuff_dir}
#make -j16
#make install

#Hoard
#cd ${base_dir}
#env GIT_SSL_NO_VERIFY=true git clone --recursive http://github.com/emeryberger/Hoard.git
#cd Hoard/src
#make linux-gcc-x86-64


#abyss
#cd ${base_dir}
#wget http://www.bcgsc.ca/platform/bioinfo/software/abyss/releases/1.3.5/abyss-1.3.5.tar.gz
#tar xzf abyss-1.3.5.tar.gz
#rm abyss-1.3.5.tar.gz
#cd abyss-1.3.5
#wget http://downloads.sourceforge.net/project/boost/boost/1.50.0/boost_1_50_0.tar.bz2
#tar jxf boost_1_50_0.tar.bz2
#ln -s boost_1_50_0/boost boost
#./configure --prefix=${stuff_dir} CPPFLAGS=-I${stuff_dir}/include
#make
#make install

#sga
#cd ${base_dir}
#git clone git://github.com/jts/sga.git
#cd ${sga_dir}/src
#./autogen.sh
#./configure --prefix=${stuff_dir} --with-sparsehash=${stuff_dir} --with-hoard=${base_dir}/Hoard/src
#make -j16
#make install


#cd ${test_dir}
#wget https://raw.github.com/jts/sga/master/src/examples/sga-ecoli-miseq.sh
#wget https://raw.github.com/jts/sga/master/src/examples/sga-human-NA12878.txt
#wget https://raw.github.com/jts/sga/master/src/examples/sga-celegans.sh

#echo
#echo "binary directory:" ${stuff_dir}/bin
#echo
#echo "official examples directory:" ${test_dir}
#echo




