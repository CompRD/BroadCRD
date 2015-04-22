#! /bin/bash -x


if [ $# -ne 2 ]; then
  echo "Function: run PBJelly installed in a specified directory." 1>&2
  echo "Usage: `basename $0` velvet_dir bam_file" 1>&2
  exit 1
fi


#starting_dir=${PWD}
#base_dir=${starting_dir}/sga_standalone

velvet_dir=$1
bam_file=$2

source /broad/software/scripts/useuse
use -q Perl-5.8.8
export PERL5LIB=/broad/tools/Linux/x86_64/lib/perl5/site_perl/5.8.8:${PERL5LIB}

export PATH=${velvet_dir}/VelvetOptimiser-2.2.5:${velvet_dir}/velvet_1.2.10:${PATH}



vdata=${velvet_dir}/velvet_1.2.10/data

#export OMP_NUM_THREADS=4

#VelvetOptimiser.pl -c -s 25 -x 10 -e 95 -t 1 -g 2000000 -f "-sam -shortPaired ${sam_file}"

#VelvetOptimiser.pl -s 27 -e 61 -f "-short -sam ${vdata}/test_reads.sam"
#VelvetOptimiser.pl -c -s 25 -x 10 -e 95 -f "-short -sam ${vdata}/test_reads.sam"

export OMP_NUM_THREADS=2

VelvetOptimiser.pl -s 25 -e 95 -x 10 -t 16  -f "-bam -shortPaired ${bam_file}"
