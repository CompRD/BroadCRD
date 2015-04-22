#!/usr/bin/perl -w
###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2012) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# ErrorCorrectReads.pl
#
# Script to run the ALLPATHS-LG error correction code, and optionally also the
# fragment filling code if given overlapping paired reads.
#
# This script requires the full installation of ALLPATHS-LG. You should first
# test your installation using the example data before attempting to use this
# script. It will run the following key modules from the ALLPATHS-LG pipeline:
# RemoveDodgyReads, FindErrors, CleanCorrectedReads, FillFragments (optional).
#
# Note that the error correction code in ALLPATHS-LG was designed specifically
# for the needs of assembly and may not be appropriate for your application.
# It operates best given a coverage level of between 20x and 200x. You may
# encounter problems if used outside this range.
#
# The code expects Illumina unpaired or paired end (fragment) reads of size
# 50 bases to 300 bases. It is generally not suitable for error correcting mate
# pair libraries unless they have been appropriately trimmed. It you attempt to
# use it with other types of data you may experience problems.
#
# Please supply your reads in FASTQ format. Note that ALLPATHS-LG discards the
# read names, but does generate a table that allows you to map the error
# corrected read IDs back to the original ones. This is important, as during
# error correction some reads may be discarded as uncorrectable and so there
# is not a one to one correspondence between input and output reads.
#
# If your reads are unpaired, or you don't care about the pairing information
# then use the UNPAIRED_READS_IN argument to specify your input FASTQ file.
# You should also use this if you want to use the ID table to map the original
# read names to the error corrected reads, even for paired reads.
#
# If you have paired reads and would like to retain the pairing information
# after error correction then use the combination of PAIRED_READS_A_IN and
# PAIRED_READS_B_IN to specify the FASTQ files containing the 1st and 2nd reads
# of each pair respectively. You can combine this with UNPAIRED_READS_IN if
# you would also like to include some unpaired reads.
#
# The script will produce a FASTQ file containing all the error corrected reads:
#   READS_OUT.fastq
# Plus the associated ID mapping table
#   READS_OUT.fastq.ids
#
# If you have paired reads then the script will also generate files containing
# the error corrected 1st and 2nd reads of each pair:
#   READS_OUT.paired.A.fastq
#   READS_OUT.paired.B.fastq
# Plus unpaired reads, which includes those reads that were originally paired
# but whose partners have been discared during error correction:
#   READS_OUT.unpaired.fastq
#
# Quality scores are passed through unchanged, except for those bases which were
# modified - they are set to quality 0.
#
# The error correction process also produces kmer specta. You can find them in
# the directory:
#   READS_OUT.fastq.kspec
#
# Fragment Filling
#
# If you have paired reads from a library that overlaps, or nearly overlaps,
# then you can also run the extra fragment filling stage. This tries to merge
# the overlapping pairs to produce 'filled fragments' or 'super reads'. This
# only works if your fragment is short enough that the reads overlap.
# For example, if your reads are 100 bases long and from a fragment that is 200
# bases long, then you will be able to merge them together. However, 100 base
# reads from a more standard Illumina library with a fragment size of ~400 bases
# will not work and an error will occur. Please see the ALLPATHS-LG manual and
# blog for more information.
#
# To use your paired end reads to generate filled fragments use the option:
#   FILL_FRAGMENTS=True
#
# You should also specify the geometry of your library with the options:
#   PAIRED_SEP
#   PAIRED_STDEV
#
# Where:
#                    <---- PAIRED SEP ---->
#   -----read A----->                      <-----read B-----
#
# Note that for overlapping reads, PAIRED_SEP will be negative.
#
# The filled fragments produced can be found in the file:
#   READS_OUT.filled.fastq
#
# Filled fragments do not have real quality scores. In the FASTQ file they are set to Q50.
#
# Additional arguments may be passed to RunAllPathsLG. Simply add any valid RunAllPathsLG
# argument to your command line when calling ErrorCorrectReads.pl.

use strict;
use FindBin;
use File::Basename;
use File::Path;

# ---- Local libraries, from elsewhere in the BroadCRD repository.
use lib "$FindBin::Bin/../";
use lib "$FindBin::Bin/";

use PerlRunTime;    # run_or_die
use PerlUtils;      # ISO_date, abort
use ArachneArgs;    # Command-line ARG=VALUE parser


# ---- CONTROL BEGINS HERE
#      Parse command-line options of the form KEY=value.
#      This function comes from the ArachneArgs.pm module.

my %args = getCommandArguments
  (UNPAIRED_READS_IN     => { value => "",
			      help  => "Fastq file containing unpaired reads to error correct." },
   PAIRED_READS_A_IN     => { value => "",
			      help  => "Fastq file containing the 1st (-->) paired reads to error correct." },
   PAIRED_READS_B_IN     => { value => "",
			      help  => "Fastq file containing the 2nd (<--) paired reads to error correct." },
   READS_OUT             => { value => undef,
			      help  => "Error corrected output fastq file." },
   PHRED_ENCODING        => { value => undef,
			      help  => "Fastq PHRED quality scores encoding (offset). Valid values are 33 or 64." },
   THREADS               => { value => "max",
			      help  => "Number of threads to use in error correction." },
   MAX_MEMORY_GB         => { value => "0",
			      help  => "Restrict maximum memory usage - default is to use it all." },
   PLOIDY                => { value => "2",
			      help  => "Organism's ploidy. Used for evaluation purposes only. Valid values are 1 or 2." },
   FORCE_PHRED           => { value => "False",
			      help  => "If True, ignore PHRED encoding warnings." },
   REMOVE_DODGY_READS    => { value => "True",
			      help  => "Run RemoveDodgyReads as part of the ALLPATHS-LG error correction pipeline." },
   FILL_FRAGMENTS        => { value => "False",
			      help  => "Perform fragment filling step on error corrected paired reads."},
   PAIRED_SEP            => { value => "-20",
			      help  => "Paired read library mean separation. ( --A-->  sep  <--B-- )" },
   PAIRED_STDEV          => { value => "30",
			      help  => "Paired read library standard deviation." },
   KEEP_KMER_SPECTRA     => { value => "False",
			      help  => "Save the kmer spectra generated during error correction" },
   KEEP_INTERMEDIATES    => { value => "False",
			      help  => "Keep intermediate files generated by ALLPATHS-LG pipeline." });


# ---- Validate arguments

my $phred_encoding = $args{PHRED_ENCODING};
if ($phred_encoding != 33 && $phred_encoding != 64) {
  die "Invalid value for PHRED_ENCODING. Must be either 33 or 64\n";
}

my $ploidy = $args{PLOIDY};
if ($ploidy != 1 && $ploidy != 2) {
  die "Invalid value for PLOIDY. Must be either 1 or 2\n";
}

if (!$args{UNPAIRED_READS_IN} && !($args{PAIRED_READS_A_IN} && $args{PAIRED_READS_B_IN})) {
  die "You must provide unpaired reads and/or paired reads (A & B).\n";
}

if (($args{PAIRED_READS_A_IN} xor $args{PAIRED_READS_B_IN})) {
  die "You must provide both A and B paired reads.\n";
}

if (($args{FILL_FRAGMENTS} && !($args{PAIRED_READS_A_IN} && $args{PAIRED_READS_B_IN}))) {
  die "You must provide paired reads in order to generate filled fragments.\n";
}

# ---- Parse arguments

my $phred_64 = ($args{PHRED_ENCODING} == 64);
my $force_phred = $args{FORCE_PHRED};
my $lg_target = ($args{FILL_FRAGMENTS} ? "filled_fragments" : "error_correction");

# ---- Determine filenames and directories to use

my($fastq_head_out, $output_dir, $suffix) = fileparse($args{READS_OUT}, ".fastq");
my $unpaired_fastq_in = $args{UNPAIRED_READS_IN};
my $paired_a_fastq_in = $args{PAIRED_READS_A_IN};
my $paired_b_fastq_in = $args{PAIRED_READS_B_IN};

my $fastq_out = "$output_dir$fastq_head_out.fastq";
my $lg_ref_name = "$fastq_head_out.allpaths-lg";
my $lg_ref_dir = "$output_dir$lg_ref_name";
my $lg_data_dir = "$lg_ref_dir/data";
my $lg_run_dir = "$lg_data_dir/run";

# ---- Make sure input FASTQ files exist

if ($args{UNPAIRED_READS_IN} && -e $args{UNPAIRED_READS_IN} == 0) {
  abort("Fastq file '$args{UNPAIRED_READS_IN}' does not exist.");
}

if ($args{PAIRED_READS_A_IN} && -e $args{PAIRED_READS_A_IN} == 0) {
  abort("Fastq file '$args{PAIRED_READS_A_IN}' does not exist.");
}

if ($args{PAIRED_READS_B_IN} && -e $args{PAIRED_READS_B_IN} == 0) {
  abort("Fastq file '$args{PAIRED_READS_B_IN}' does not exist.");
}

# ---- Create allpaths-lg pipeline directory

if (-e $lg_ref_dir) {
  die "Found existing allpaths-lg pipeline directory: " .
  $lg_ref_dir . "\nPlease erase to continue.\n";
}

mkdir $lg_ref_dir;
mkdir $lg_data_dir;

my $cmd;

# ---- Setup Ploidy

run_or_die("echo $args{PLOIDY} > $lg_data_dir/ploidy");

# ---- Create dummy jump files to fool ALLPATHS-LG

run_or_die("touch $lg_data_dir/jump_reads_orig.fastb");
run_or_die("touch $lg_data_dir/jump_reads_orig.qualb");
run_or_die("touch $lg_data_dir/jump_reads_orig.pairs");


# ---- Run FastqToFastbQualb for unpaired reads

if ($args{UNPAIRED_READS_IN}) {
  fastq_to_fastb($unpaired_fastq_in, "$lg_data_dir/unpaired", $phred_64, $force_phred );

# ---- Run PairsFake

  print "\nRunning PairsFake.\n";
  $cmd = ("$FindBin::Bin/PairsFake" .
	" HEAD=$lg_data_dir/unpaired" );
  run_or_die($cmd);

}

# ---- Run FastqToFastbQualb for paired reads

if ($args{PAIRED_READS_A_IN} && $args{PAIRED_READS_B_IN}) {
  fastq_to_fastb($paired_a_fastq_in, "$lg_data_dir/paired_a", $phred_64, $force_phred );
  fastq_to_fastb($paired_b_fastq_in, "$lg_data_dir/paired_b", $phred_64, $force_phred );
}

# ---- Merge reads

if ($args{PAIRED_READS_A_IN} && $args{PAIRED_READS_B_IN}) {

  print "\nRunning MergePairedFastbs.\n";
  $cmd = ("$FindBin::Bin/MergePairedFastbs" .
	  " HEAD1_IN=$lg_data_dir/paired_a" .
	  " HEAD2_IN=$lg_data_dir/paired_b" .
	  " HEAD_OUT=$lg_data_dir/paired" .
	  " LIB_NAME=paired" .
	  " LIB_SEP=$args{PAIRED_SEP}" .
	  " LIB_STDEV=$args{PAIRED_STDEV}" .
	  " WRITE_PAIRS=True");
  run_or_die($cmd);

  if ($args{UNPAIRED_READS_IN}) {
    
    print "\nRunning MergeReadSets.\n";
    $cmd = ("$FindBin::Bin/MergeReadSets" .
	    " DIR=$lg_data_dir" .
	    " READS_IN=\"{paired,unpaired}\"" .
	    " READS_OUT=frag_reads_orig");
    run_or_die($cmd);

  } else {
    symlink("paired.fastb", "$lg_data_dir/frag_reads_orig.fastb");
    symlink("paired.qualb", "$lg_data_dir/frag_reads_orig.qualb");
    symlink("paired.pairs", "$lg_data_dir/frag_reads_orig.pairs");
  }

} else { # Unpaired only
  symlink("unpaired.fastb", "$lg_data_dir/frag_reads_orig.fastb");
  symlink("unpaired.qualb", "$lg_data_dir/frag_reads_orig.qualb");
  symlink("unpaired.pairs", "$lg_data_dir/frag_reads_orig.pairs");
}


# ---- Run ALLPATHS-LG
print "\nRunning ALLPATHS-LG.\n";
$cmd = ("$FindBin::Bin/RunAllPathsLG" .
	" PRE=$output_dir" .
	" REFERENCE_NAME=$lg_ref_name" .
	" DATA_SUBDIR=data" .
	" RUN=run" .
	" THREADS=$args{THREADS}" .
	" MAX_MEMORY_GB=$args{MAX_MEMORY_GB}" .
	" VALIDATE_INPUTS=False" .
	" TARGETS=$lg_target" .
	" REMOVE_DODGY_READS_FRAG=" . ($args{REMOVE_DODGY_READS} != 0 ? "True" : "False") .
	" $args{'_extra'}");
run_or_die($cmd);

# ---- Convert results to fastq

if ($args{PAIRED_READS_A_IN} && $args{PAIRED_READS_B_IN}) {

# ---- Run SplitReadsByLibrary

  print "\nRunning SplitReadsByLibrary.\n";
  $cmd = ("$FindBin::Bin/SplitReadsByLibrary" .
	  " READS_IN=$lg_run_dir/frag_reads_corr" .
	  " SPLIT_PAIRS=True" .
	  " QUALS=True" );
  run_or_die($cmd);

  fastb_to_fastq("$lg_run_dir/frag_reads_corr.paired.A",
		 "$output_dir$fastq_head_out.paired.A", $phred_encoding);
  fastb_to_fastq("$lg_run_dir/frag_reads_corr.paired.B",
		 "$output_dir$fastq_head_out.paired.B", $phred_encoding);
  if (-e "$lg_run_dir/frag_reads_corr.unpaired.fastb") {
    fastb_to_fastq("$lg_run_dir/frag_reads_corr.unpaired",
		   "$output_dir$fastq_head_out.unpaired", $phred_encoding);
  }

  if ($args{FILL_FRAGMENTS}) {
    fastb_to_fastq("$lg_run_dir/filled_reads",
		   "$output_dir$fastq_head_out.filled", $phred_encoding);
  }

}

fastb_to_fastq("$lg_run_dir/frag_reads_corr",
	       "$output_dir$fastq_head_out", $phred_encoding);


# ---- Run ReadTrack

print "\nRunning ReadTrack. (silent)\n";
$cmd = ("$FindBin::Bin/ReadTrack" .
	" NH=True" .
	" READS=$lg_run_dir/frag_reads_corr" .
	" TRACK_BACK=True" .
	" > $lg_run_dir/frag_reads_corr.tracked");
run_or_die($cmd);

# ---- Extract read ID tracking information

print "\nParsing read tracking information. (silent)\n";
parse_readtrack("$lg_run_dir/frag_reads_corr.tracked", "$fastq_out.ids");

# ---- Keep kmer spectra

if ($args{KEEP_KMER_SPECTRA} != 0) {
  mkdir "$fastq_out.kspec";
  run_or_die("cp $lg_run_dir/*.kspec $fastq_out.kspec");
}

# ---- Clean up

if ($args{KEEP_INTERMEDIATES} == 0) {
  rmtree($lg_ref_dir);
}

print "\nDone.\n";




sub parse_readtrack {

  my ($file_in, $file_out) = @_;

  my $input;
  open ($input, "<", $file_in) or die "Unable to open: $file_in\n";

  my $output;
  open ($output, ">", $file_out) or die "Unable to write to: $file_out\n";

  print $output "#New_ID,  Original_ID\n";

  my $line;
  while ($line = <$input>) {
    chomp $line;
    if ($line =~ /^frag_reads_corr:(\d+)\s+.*:(\d+)\s*$/ ) {
      print $output "$1,  $2\n";
    }
  }

  close $output;
  close $input;
  
}


sub fastq_to_fastb {

  my ($file_in, $head_out, $phred_64, $force_phred) = @_;

  print "\nRunning FastqToFastbQualb.\n";
  $cmd = ("$FindBin::Bin/FastqToFastbQualb" .
	  " FASTQ=$file_in" .
	  " OUT_HEAD=$head_out" .
	  " PHRED_64=" . bool_str($phred_64) .
	  " FORCE_PHRED=" . bool_str($force_phred) );
  run_or_die($cmd);
  
}


sub fastb_to_fastq {

  my ($head_in, $head_out, $phred_encoding) = @_;

  my $quals = (-e "$head_in.qualb");

  print "\nRunning FastbQualbToFastq.\n";
  $cmd = ("$FindBin::Bin/FastbQualbToFastq" .
	  " HEAD_IN=$head_in" .
	  " HEAD_OUT=$head_out" .
	  " NAMING_PREFIX=read" .
	  " PAIRED=false" .
	  " PHRED_OFFSET=$phred_encoding" .
	  " SIMULATE_QUALS=" . bool_str(!$quals) );

  run_or_die($cmd);
}

