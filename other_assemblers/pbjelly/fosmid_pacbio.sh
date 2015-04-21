#! /bin/bash 


source /broad/software/scripts/useuse
#use .git-1.7.10.1
use Samtools

pb_loc=/home/unix/blau/wga/pbjelly/software/PBJelly_12.9.14/

pb_base=/wga/scr4/picard/pacbio


pb_runs=( '019892' '019891' '019899' '019900' '019914' '019915' '019916' )

for run in "${pb_runs[@]}"
do
    samtools view -h  ${pb_base}/${run}/comprd-data/filtered_subreads.bam | SAM2CRDDump OUT_HEAD=${run} \
        SEP=-15 DEV=12 NOMINAL_READ_LEN=251 USE_OQ=True \
        NH=True LOG_TO_CERR=False WRITE_ALIGNS=False WRITE_NAMES=False \
        LIBINFO=/wga/scr4/dexter/libinfo/dexter_libs
    FastbQualbToFastq HEAD_IN=${run} HEAD_OUT=${run} PAIRED=False PHRED_OFFSET=33
    ${pb_loc}/fastqSplit.py ${run}.fastq -o ${run}
done



