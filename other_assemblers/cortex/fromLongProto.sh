
if [ $# -ne 4 ]; then
  echo "Usage: `basename $0` cortex_position long_proto_read_file_head long_proto_reference length" 1>&2
  exit 1
fi

cortex_loc=$1
long_proto_read_head=$2
long_proto_reference=$3
length=$4

echo running with $length

cortex_bin_31=${cortex_loc}/bin/cortex_var_31_c1
cortex_bin_63=${cortex_loc}/bin/cortex_var_63_c1
cortex_bin_95=${cortex_loc}/bin/cortex_var_95_c1
stampy_bin=/wga/scr4/assemblers/cortex/stampy-1.0.21/stampy.py
vcf_loc=/wga/scr4/assemblers/cortex/vcftools_0.1.8a


#environment variables
if [ ! -d $cortex_loc ]; then
  echo "$cortex_loc does not exist as a directory" 1>&2
  exit 1
fi
d1=${cortex_loc}/scripts/analyse_variants/bioinf-perl/lib
d2=${cortex_loc}/scripts/calling
d3=${cortex_loc}/scripts/analyse_variants/needleman_wunsch
if [ ! -d $d1 ]; then
  echo "${d1} does not exist as a directory" 1>&2; exit 1
fi

if [ ! -d $d2 ]; then
  echo "${d2} does not exist as a directory" 1>&2; exit 1
fi

if [ ! -d $d3 ]; then
  echo "${d3} does not exist as a directory" 1>&2; exit 1
fi
export PERL5LIB=${d1}:${d2}
export PATH=${d3}:${PATH}


#import reads from long proto
read_head=${PWD}/`basename ${long_proto_read_head}`

FastbQualbToFastq HEAD_IN=${long_proto_read_head} HEAD_OUT=${read_head} PAIRED=True PHRED_OFFSET=33

echo ${read_head}.A.fastq > pe_filelist1
echo ${read_head}.B.fastq > pe_filelist2
printf 'fromLP\t.\tpe_filelist1\tpe_filelist2\n' > INDEX

#import reference  from long proto
Fastb2Fasta IN=${long_proto_reference} OUT=ref.fasta
echo ref.fasta > list_ref_fasta
${stampy_bin} -G 0 ref.fasta
${stampy_bin} -g 0 -H 0 

#gzip ${read_head}.A.fastq
#gzip ${read_head}.B.fastq


#make reference, page 19-20 in manual
echo ref.fasta > ref_filelist
${cortex_bin_31} --kmer_size 31 \
                 --mem_height 17 \
                 --mem_width 100 \
                 --se_list ref_filelist \
                 --max_read_len 10000 \
                 --dump_binary ref.k31.ctx \
                 --sample_id REF

${cortex_bin_63} --kmer_size 61 \
                 --mem_height 17 \
                 --mem_width 100 \
                 --se_list ref_filelist \
                 --max_read_len 10000 \
                 --dump_binary ref.k61.ctx \
                 --sample_id REF

${cortex_bin_95} --kmer_size 91 \
                 --mem_height 17 \
                 --mem_width 100 \
                 --se_list ref_filelist \
                 --max_read_len 10000 \
                 --dump_binary ref.k91.ctx \
                 --sample_id REF




perl ${cortex_loc}/scripts/calling/run_calls.pl --first_kmer 31    \
                  --last_kmer 91    \
                  --kmer_step 30    \
                  --fastaq_index INDEX --auto_cleaning yes    \
                  --bc yes --pd no    \
                  --outdir ${PWD}/out    \
                  --outvcf vcf_out    \
                  --ploidy 2    \
                  --stampy_hash 0    \
                  --stampy_bin ${stampy_bin}    \
                  --list_ref_fasta list_ref_fasta    \
                  --refbindir ${PWD}    \
                  --genome_size ${length}    \
                  --qthresh 10    \
                  --mem_height 17  --mem_width 100    \
                  --vcftools_dir ${vcf_loc}/ \
                  --do_union yes    \
                  --ref CoordinatesAndInCalling    \
                  --workflow independent    \
                  --logfile logfile log.txt    \
                  --apply_pop_classifier




date

