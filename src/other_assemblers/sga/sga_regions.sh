#/bin/bash -f

file_name=$1
sga_script=~/wga/src/BroadCRD/other_assemblers/sga/fromLongProto.sh
sga_summary=/home/unix/blau/wga/src/BroadCRD/bin/SummarizeSGA

org_dir=${PWD}


cat ${file_name} | while read -r LINE;do

  index=`echo $LINE | awk '{print $1}'`
  chrom=`echo $LINE | awk '{print $2}'`
  coord=`echo $LINE | awk '{print $3}'`

  chrom_shifted=`echo $chrom-1| bc -l` 

  range="${chrom_shifted}:${coord}"

  echo $index $chrom $coord $chrom_shifted ${range}

  work_dir=${org_dir}/${index}/sga

  lp_head=${org_dir}/${index}/LP/tmp.xxx/frag_reads_orig

  rm -rf ${work_dir}
  mkdir -p ${work_dir}
  cd ${work_dir}

  ${sga_script} ${lp_head} | tee run.log

  cp ${org_dir}/${index}/LP/ref.fastb genome.fastb

#  ${sga_summary} IN_DIR=. | tee sum.log


  cd ${org_dir}
done
