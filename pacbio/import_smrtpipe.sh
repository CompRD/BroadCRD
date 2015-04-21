#!/usr/bin/env bash

#clone the setting of pacbio job in $1 to the directory $2, invokable by job.sh
#some, but not all, algorithmic settings, ie quiver coverage can be edited in settings.xml

set -e

jid=$1
tgt=$2
pacbio_root=/seq/pacbio_results/userdata/jobs/

if [ $# -ne 2 ] ; then
    echo "$0 jid tgt_dir" >& 2
fi

if [ -e $tgt ] ; then
    echo "$tgt exists" >& 2
fi


first3=${jid:0:3}

pacbio_loc=${pacbio_root}/${first3}/${jid}

echo $pacbio_loc

if [ ! -d $pacbio_loc ] ; then
    echo "$pacbio_loc does not exist" >& 2
fi

src_files=("job.sh" "input.xml" "settings.xml")

for file in "${src_files[@]}" ; do
    full_path=${pacbio_loc}/$file
    if [ ! -f $full_path ] ; then
        echo "$full_path does not exist" >& 2
    fi
    unset full_path
done



mkdir $tgt
tgt=`readlink -e $tgt`

for file in "${src_files[@]}" ; do
    cp ${pacbio_loc}/$file $tgt
    chmod 666 ${tgt}/${file}
done

tmp_dir=${tgt}/tmp_crd
mkdir ${tmp_dir}


setup_script=`cat ${tgt}/job.sh | awk '{print $2}'`
echo ". $setup_script && export MPLCONFIGDIR=${tgt}/.m && ([[ ! -e \$MPLCONFIGDIR ]] && mkdir \$MPLCONFIGDIR ) || true && smrtpipe.py -D TMP=${tmp_dir} -D NPROC=48 --params=${tgt}/settings.xml --output=${tgt} xml:${tgt}/input.xml" > ${tgt}/job.sh
chmod 700 ${tgt}/job.sh


sed -e s/${jid}/${jid}.mod/g < ${tgt}/input.xml > ${tgt}/tmp
mv ${tgt}/tmp ${tgt}/input.xml
chmod 666 ${tgt}/input.xml
