#!/bin/bash 

if [ $# -ne 5 ]; then
  echo "Usage: `basename $0` velvet_dir long_proto_run" 1>&2
  exit 1
fi

velvet_dir=$1
fosmid_id=$2

K=$3
INS=$4
INS_SD=$5

CUTOFF_COV=250

base_dir=`pwd`

#long_proto_run=lp

#if [[ "${fosmid_id}" -ge 56 ]] ; then
#    SAMPLE=hpool3
#else
#    SAMPLE=hpool2
#fi

#mkdir ${long_proto_run}

#cd ${long_proto_run}
#export MALLOC_PER_THREAD=1
#LongProto READS=#picard SAMPLE=${SAMPLE} TMP=tmp X=${fosmid_id} OUT_INT_HEAD=oih OUT_GENOME=ref.fastb EXIT=CORRECT
#cd ${base_dir}

#ref_length=`FastbStats FASTB=lp/ref.fastb | grep 'total length' | awk '{print $3}'`
#n_reads=`FastbStats FASTB=lp/tmp/frag_reads_orig.fastb | grep 'total objects count' | awk '{print $4}'`
#echo $ref_length $n_reads

#read_length=250


#EXP_COV=`echo " $n_reads * ( $read_length - $K ) * ( $INS -$K) / 2 / ($read_length-$K)  / ${ref_length} " | bc -l`

EXP_COV=`cat output/stats.txt | awk '(NR>1){if($6<250){nbases+=$2;cum+=$6*$2}}END{if(nbases>0){print cum/nbases}else{print 150}}' `
COV_CUT=`echo "$EXP_COV * 0.5" | bc -l`

#echo "cov,cov_cutoff: $EXP_COV,$COV_CUT"


#Fastb2Fasta IN=${long_proto_run}/ref.fastb OUT=ref.fasta
#MakeLookupTable SOURCE=ref.fasta OUT_HEAD=ref LO=True
#FastbQualbToFastq HEAD_IN=${long_proto_run}/tmp/frag_reads_orig HEAD_OUT=f PAIRED=True PHRED_OFFSET=33
#python ~/src/python/stagger_pair_end_fastq.py f.A.fastq f.B.fastq > input.fastq

#Fastb2Fasta IN=${long_proto_run}/tmp/frag_reads_orig.fastb OUT=reads.fasta


source /broad/software/scripts/useuse
use -q Perl-5.8.8
export PERL5LIB=/broad/tools/Linux/x86_64/lib/perl5/site_perl/5.8.8:${PERL5LIB}

export PATH=${velvet_dir}/VelvetOptimiser-2.2.5:${velvet_dir}/velvet_1.2.10:${PATH}



old_dir=output.`date +%D-%T | tr '/:' '-'`

mv output ${old_dir}



echo ----------------------------------------------
echo starting velvet
echo ----------------------------------------------

velveth ${PWD}/output ${K} -fastq -shortPaired ${PWD}/input.fastq

velvetg ${PWD}/output -ins_length ${INS} -ins_length_sd ${INS_SD}  -exp_cov ${EXP_COV} -cov_cutoff ${COV_CUT}


BreakAtNs IN=output/contigs.fa OUT=output/contigs.noN.fa MIN_TO_BREAK=1

QueryLookupTable K=12 MM=12 MC=0.15 SEQS=output/contigs.noN.fa L=ref.lookup PARSEABLE=False SMITH_WAT=True VISUAL=False QUERY_NAMING=from_record

Fasta2Fastb IN=output/contigs.noN.fa
MakeLookupTable SOURCE=output/contigs.noN.fa OUT_HEAD=output/contigs.noN LO=True

FastbStats FASTB=output/contigs.noN.fastb
echo "cov,cov_cutoff: $EXP_COV,$COV_CUT"

