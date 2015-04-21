#! /bin/bash -x


if [ $# -ne 1 ]; then
  echo "Function: to download, compile, and install PBJelly in a target directory." 1>&2
  echo "Usage: `basename $0` target_directory" 1>&2
  exit 1
fi


#starting_dir=${PWD}
#base_dir=${starting_dir}/sga_standalone

base_dir=$1
#if [ -e $base_dir ]; then
#  echo "$base_dir exists" 1>&2
#  exit 1
#fi

#mkdir $base_dir
#if [ $? -ne 0 ]; then
#  echo "Error making directory $base_dir" 1>&2
#fi


#source /broad/software/scripts/useuse
#use .git-1.7.10.1
#use BamTools
#use GCC-4.4

#use

stuff_dir=${base_dir}/stuff
mkdir ${stuff_dir}


cd ${base_dir}

smrt_base=smrtanalysis-1.3.1
smrt_package=${smrt_base}-centos.tgz

wget http://files.pacb.com/software/smrtanalysis/1.3.1/${smrt_package}
tar -C ${base_dir} -xvvzf ${smrt_package}


smrt_dir=${base_dir}/smrtanalysis
ln -s ${smrt_base} ${smrt_dir}
chown -R `id -u -n`:`id -g -n` ${smrt_base}

sed -e "s|/opt/smrtanalysis|${smrt_dir}|"  ${smrt_dir}/etc/setup.sh > bayo.tmp
mv bayo.tmp ${smrt_dir}/etc/setup.sh

export SEYMOUR_HOME=${smrt_dir}

#${smrt_dir}/etc/scripts/postinstall/configure_smrtanalysis.sh

cd ${base_dir}

pbj_base=PBJelly_12.9.14
pbj_package=${pbj_base}.tgz
wget "http://downloads.sourceforge.net/project/pb-jelly/PBJelly_12.9.14.tgz?r=&ts=1376667478&use_mirror=softlayer-dal" -O ${pbj_package}

tar zxf ${pbj_package}



