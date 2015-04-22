#/bin/bash -f

file_name=$1
cortex_script=~/wga/src/BroadCRD/other_assemblers/cortex/fromLongProto.sh
cortex_loc=/wga/scr4/assemblers/cortex/CORTEX_release_v1.0.5.15
#sga_summary=/home/unix/blau/wga/src/BroadCRD/bin/SummarizeSGA

org_dir=${PWD}


cat ${file_name} | while read -r LINE;do

  index=`echo $LINE | awk '{print $1}'`
  chrom=`echo $LINE | awk '{print $2}'`
  coord=`echo $LINE | awk '{print $3}'`

  chrom_shifted=`echo $chrom-1| bc -l` 

  range="${chrom_shifted}:${coord}"

  echo $index $chrom $coord $chrom_shifted ${range}

  length=`echo \($coord\)*-1 | bc -l`
  echo $length

  work_dir=${org_dir}/${index}/cortex

  lp_head=${org_dir}/${index}/LP/tmp.xxx/frag_reads_orig

  rm -rf ${work_dir}
  mkdir -p ${work_dir}
  cd ${work_dir}

  ${cortex_script} ${cortex_loc} ${lp_head} ${org_dir}/${index}/LP/ref.fastb ${length} | tee run.log

#  cp ${org_dir}/${index}/LP/ref.fastb genome.fastb

#  ${sga_summary} IN_DIR=. | tee sum.log


  cd ${org_dir}
done
