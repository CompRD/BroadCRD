#! /bin/bash -x


if [ $# -ne 1 ]; then
  echo "Function: run velvet installed in a specified directory." 1>&2
  echo "Usage: `basename $0` velvet_dir" 1>&2
  exit 1
fi


#starting_dir=${PWD}
#base_dir=${starting_dir}/sga_standalone

velvet_dir=$1

source /broad/software/scripts/useuse
use -q Perl-5.8.8
export PERL5LIB=/broad/tools/Linux/x86_64/lib/perl5/site_perl/5.8.8:${PERL5LIB}

export PATH=${velvet_dir}/VelvetOptimiser-2.2.5:${velvet_dir}/velvet_1.2.10:${PATH}



vdata=${velvet_dir}/velvet_1.2.10/data

#VelvetOptimiser.pl -v -s 25 -e 95 -x 10 -t 12  -f "-fastq -shortPaired ${PWD}/input.fastq"

velveth ${PWD}/output 61 -fastq -shortPaired ${PWD}/input.fastq

velvetg ${PWD}/output -ins_length 300 -ins_length_sd 50 -exp_cov auto -cov_cutoff auto
