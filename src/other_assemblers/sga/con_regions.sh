#/bin/bash -f

file_name=$1
sga_script=~/wga/src/BroadCRD/other_assemblers/sga/fromLongProto.sh
sga_summary=/home/unix/blau/wga/src/BroadCRD/bin/SummarizeSGA
sga_postprocess=/home/unix/blau/wga/src/BroadCRD/other_assemblers/sga/post_process.sh

org_dir=${PWD}

rm -f sequence_graph*


cat ${file_name} | while read -r LINE;do
  index=`echo $LINE | awk '{print $1}'`
  chrom=`echo $LINE | awk '{print $2}'`
  coord=`echo $LINE | awk '{print $3}'`

  chrom_shifted=`echo $chrom-1| bc -l` 

  range="${chrom_shifted}:${coord}"

  echo $index $chrom $coord $chrom_shifted ${range}

#  index=`echo $LINE | awk '{print $1}'`
#  range=`echo $LINE | awk '{print $2}'`

#  echo ">> " $index $range

  work_dir=${org_dir}/${index}/sga

  lp_head=${org_dir}/${index}/LP/tmp.xxx/frag_reads_orig

#  rm -rf ${work_dir}
#  mkdir -p ${work_dir}
  cd ${work_dir}

#  ${sga_script} ${lp_head} | tee run.log

  cp ${org_dir}/${index}/LP/ref.fastb genome.fastb

  rm -f sequence_graph.dot sequence_graph.fasta

  ${sga_summary} IN_DIR=. VARIANTS=False | tee sum.log
  cp sequence_graph.dot ${org_dir}/sequence_graph_${index}.dot
  cp sequence_graph.fasta ${org_dir}/sequence_graph_${index}.fasta

#  ${sga_postprocess}

  ${sga_summary} IN_DIR=. VARIANTS=True | tee sum.log

  ${sga_summary} IN_DIR=. | tee sum.log
  cp sequence_graph.dot ${org_dir}/sequence_graph_${index}v.dot
  cp sequence_graph.fasta ${org_dir}/sequence_graph_${index}v.fasta


  cd ${org_dir}

  neato -Tpng -osequence_graph_${index}v.png sequence_graph_${index}v.dot 
  neato -Tpng -osequence_graph_${index}.png sequence_graph_${index}.dot 

  chmod -x *dot
  chmod -x *png
  chmod -x *fasta
done
