#!/bin/bash
# vim: tw=1000

eval `/broad/software/dotkit/init -b`
use .bwa-0.6.2

PICARD=/wga/scr4/picard/A2925/C1-508_2012-11-12_2012-11-14/1

for lib in 127359 127365;
do
    SOL=Solexa-${lib}
    BAM=$PICARD/$SOL/A2925.1.aligned.duplicates_marked.bam
    GENOME=genome.${SOL}-fosmids.fasta

    FloodNFasta OUTPUT_FILE=$GENOME REGION_FILE=${SOL}.regions FASTA_FILE=/wga/scr4/bigrefs/human19/genome.fasta PAD_LENGTH=4000 || { echo 1>&2 "error FloodNFasta"; exit 1; }

    bwa index $GENOME

    picard SamToFastq INPUT=$BAM FASTQ=${SOL}.1.fastq.gz SECOND_END_FASTQ=${SOL}.2.fastq.gz

    bwa aln $GENOME -q 5 -l 32 -k 2 -t 48 -o 1 -f ${SOL}.1.sai ${SOL}.1.fastq.gz
    bwa aln $GENOME -q 5 -l 32 -k 2 -t 48 -o 1 -f ${SOL}.2.sai ${SOL}.2.fastq.gz
    bwa sampe -P -a 900 -f ${SOL}_aligned.sam $GENOME ${SOL}.1.sai ${SOL}.2.sai ${SOL}.1.fastq.gz ${SOL}.2.fastq.gz

    samtools view -b -S ${SOL}_aligned.sam -o ${SOL}_aligned.bam 
    samtools sort  -m 2000000000 ${SOL}_aligned.bam ${SOL}_aligned.sorted
    samtools index ${SOL}_aligned.sorted.bam ${SOL}_aligned.sorted.bai
done

