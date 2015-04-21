#/bin/bash -f

file_name=$1
LP_EXE=/wga/scr4/blau/src/BroadCRD/bin/LongProto

org_dir=${PWD}


cat ${file_name} | while read -r LINE;do
  index=`echo $LINE | awk '{print $1}'`
  chrom=`echo $LINE | awk '{print $2}'`
  coord=`echo $LINE | awk '{print $3}'`

  chrom_shifted=`echo $chrom-1| bc -l` 

  range="${chrom_shifted}:${coord}"

  echo $index $chrom $coord $chrom_shifted ${range}


  work_dir=${org_dir}/${index}/LP
  rm -rf ${work_dir}
  mkdir -p ${org_dir}/${index}/LP
  cd ${org_dir}/${index}/LP

  #${LP_EXE} SAMPLE=human READS=#picard X=${range} TMP=${PWD}/tmp.xxx OUT_INT_HEAD=${PWD}/assembly OUT_GENOME=${PWD}/ref.fastb  LOGGING="{REFTRACE_VARIANTS=True}"  LOGGING=PRINT_TIME_USED=True | tee run.log

  ${LP_EXE} SAMPLE=human READS=#picard X=${range} TMP=${PWD}/tmp.xxx OUT_INT_HEAD=${PWD}/assembly OUT_GENOME=${PWD}/ref.fastb  LOGGING="{REFTRACE_VARIANTS=True,PRINT_TIME_USED=True}"  | tee run.log

  cd ${org_dir}
done
