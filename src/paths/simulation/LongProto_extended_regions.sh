#/bin/bash -f
description()
{
  echo 1>& 2
  echo 'Setup extend-fosmid test. Creates subfolder ${region_list}.extended_fosmid,' 1>& 2
  echo ' goes through each entry of ${region_list} and run LongProto,' 1>& 2
  echo ' with the option of extending the regions by a number of basis' 1>& 2
  echo 1>& 2
  echo "Usage: `basename $0` region_list [length_of_extension_to_fosmid]" 1>& 2
  echo 1>& 2
}



if [ $# -lt 1 ]; then
  description
  exit 1
fi

# optional arguments to modify the range
if [ $# -gt 1 ]; then
  extension=$2
  extension_suffix="_$2"
fi

#region file, 1st column is fosmid id, 2nd column is C:S-E, C is chromosome name, S is start, E is end
region_file=$1

# SAMPLE argument in long-proto
out_header=output

# SAMPLE argument in long-proto
sample_spec=human

# make a sub-directory in $PWD
target_root=${PWD}/`basename ${region_file}`.extended_fosmid

#make sure we don't overwrite files
if [ ! -e ${target_root} ]; then
  mkdir ${target_root}
elif [ ! -d ${target_root} ]; then
  echo "ERROR: ${target_root} exists and is not a directory" 1>& 2
  exit 1
fi

#make a copy
cp ${region_file} ${target_root}
region_file=${target_root}/`basename ${region_file}`

#takes log
log_file=${target_root}/`basename ${region_file}`.run_log


#loop over each line of region file
nLines=`wc -l < ${region_file}`
ll=0

while [ $ll -lt $nLines ]; do
  let ++ll
  line=`head -n $ll $region_file | tail -1`

  fosmid_id=`echo $line | awk '{print $1}'`
  range_spec=`echo $line | awk '{print $2}'`


  # extend range if specified
  unset reference_genome
  if [[ $extension && ${extension-x} ]]; then
    range_chromo=`echo $range_spec | sed -e 's/:.*//'`
    range_start=`echo $range_spec | sed -e 's/.*://' | sed -e 's/-.*//'`
    range_end=`echo $range_spec | sed -e 's/.*-//'`
    range_spec=$range_chromo:`echo $range_start'-'$extension|bc`-`echo $range_end'+'$extension|bc`

    reference_genome=${target_root}/${fosmid_id}/${out_header}.hbv
  fi

  child_directory=${target_root}/${fosmid_id}${extension_suffix}
  child_log=${target_root}/${fosmid_id}${extension_suffix}.log

  if [ -e ${child_directory} ]; then
    echo "skipping region $fomid_id: ${child_directory} record exists" >> ${log_file}
    continue
  fi

  mkdir $child_directory
  cd $child_directory


  if [[ $reference_genome && ${reference_genome-x} ]]; then
    cmd="LongProto SAMPLE=${sample_spec} X=${range_spec} OUT_HEAD=${child_directory}/${out_header} TMP=${child_directory}/tmp REF_TRACE=True IN_GENOME=${reference_genome}"
  else
    cmd="LongProto SAMPLE=${sample_spec} X=${range_spec} OUT_HEAD=${child_directory}/${out_header} TMP=${child_directory}/tmp REF_TRACE=True"
  fi

  echo "running: $cmd" >> ${log_file}
  $cmd >> ${child_log}

  echo >> ${log_file}
  echo "----------------------------------------------------------------------------------------------" >> ${log_file}
  echo  >> ${log_file}

  cd ${target_root}
done

echo "directory A_B stores LongProto run for fosmid A, extended by B basis on both ends" > ${target_root}/README
echo "if B=0, direcotry A is used ("_B" ommited)" > ${target_root}/README

