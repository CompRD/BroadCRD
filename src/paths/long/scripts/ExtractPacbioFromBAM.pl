#!/usr/bin/perl -w
###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2012) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################

# ExtractPacbioFromBAM.pl
#
# Script to extract reads for a region from one or more pacbio runs.  
#

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

use File::Temp qw/ tempfile /;          # generate tempfile names

# ---- Executables
my %execs = ( SAMTOOLS => "/broad/software/free/Linux/redhat_5_x86_64/pkgs/samtools/samtools_0.1.16/bin/samtools",
              MERGESAMFILES => "java -jar /seq/software/picard/current/bin/MergeSamFiles.jar TMP_DIR=. USE_THREADING=true",
              SAMTOFASTQ => "java -jar /seq/software/picard/current/bin/SamToFastq.jar TMP_DIR=.",
              FASTQ2FASTA => "Fastq2Fasta",
              FASTQTOFASTBQUALB => "FastqToFastbQualb"
            );

# ---- CONTROL BEGINS HERE
#      Parse command-line options of the form KEY=value.
#      This function comes from the ArachneArgs.pm module.

my %args = getCommandArguments(
    PACBIO_CACHE         => { value => "/wga/scr4/picard/pacbio",
                                help => "Picard cache directory with 1) PacBio runs and " . 
                                "2) comprd-data directory in each run with aligned reads." },
    RUNS                 => { value => undef,
                                help => "PacBio run numbers corresponding to directory names under PACBIO_CACHE." },
    HEAD_OUT             => { value => "pacbio_reads",
                                help => "HEAD for output FASTA filename."},
    RANGE                => { value => "",
                                help => "Range to extract for aligned reads.  Should be contig:0-10000."}
  );


# ---- Validate arguments

my $pacbio_cache = $args{PACBIO_CACHE};
-d $pacbio_cache || die "PACBIO_CACHE directory $pacbio_cache does not appear to exist\n";

my $head_out = $args{HEAD_OUT};

# ---- Validate and parse out RANGE variable, if specified
my $range = $args{RANGE};
$range eq "" || $range =~ /([^:]+):(\d+)-(\d+)/ || die "RANGE value seems to be invalid: \"$range\"\n";

my %range_exp;
if ( $range ) {
    %range_exp = ( CHR => $1, START => $2, END => $3 );

    print "range parsed as chr=$range_exp{CHR}, start=$range_exp{START}, end=$range_exp{END}\n";
}


# ---- Parse out comma-separated runs and make six-digit numbers
my @runs = map { sprintf("%06d", $_); } split( /,/, $args{RUNS} );


# ---- Check for input files
foreach my $run (@runs) {
    my $bamfile="$pacbio_cache/$run/comprd-data/filtered_subreads.bam";
    -e $bamfile || die "can't find BAM file $bamfile\n";
}


# ---- Delete temp files (names in @tempfiles) upon normal or abnormal exit
my @tempfiles;                           
sub cleanup {                           
    map( unlink, @tempfiles );
}
$SIG{__DIE__}=\&cleanup;                 # delete the temp files on abnormal exit via die()
END { map( unlink, @tempfiles) };        # delete temp files on normal exit


# ---- Extract Region of Interest from each BAM using "samtools view"
my @bamfiles;
foreach my $run ( @runs ) {
    my $inputfile="$pacbio_cache/$run/comprd-data/filtered_subreads.bam";
    (undef, my $tmpfile) = tempfile("tempfileXXXXXXX", OPEN => 1, DIR => ".", SUFFIX => ".bam"  );

    my $chr = $range_exp{CHR};
    my $start1 = $range_exp{START}+1;
    my $end1 = $range_exp{END}+1;
    my $cmd = "$execs{SAMTOOLS} view -b -o $tmpfile $inputfile $chr:$start1-$end1";
    print "CMD: $cmd\n";
    push(@tempfiles, $tmpfile);
    push(@bamfiles, $tmpfile);
    run_or_die($cmd);
}


# ---- Merge temp files
(undef, my $mergefile) = tempfile("tempfileXXXXXXX", OPEN => 1, DIR => ".", SUFFIX => ".bam"  );
my $inputs = "I=" . join( " I=", @bamfiles);
push(@tempfiles, $mergefile);
{ my $cmd = "$execs{MERGESAMFILES} $inputs O=$mergefile"; print "CMD: $cmd\n"; run_or_die($cmd); }

# ---- Convert to Fastq
(undef, my $fastqfile) = tempfile("tempfileXXXXXXX", OPEN => 1, DIR => ".", SUFFIX => ".fastq"  );
push(@tempfiles, $fastqfile);
{ my $cmd = "$execs{SAMTOFASTQ} I=$mergefile F=$fastqfile"; print "CMD: $cmd\n"; run_or_die($cmd); }


# ---- Convert to Fastb
{ my $cmd = "$execs{FASTQTOFASTBQUALB} FASTQ=$fastqfile OUT_HEAD=$head_out"; print "CMD: $cmd\n"; run_or_die($cmd); }


print "\nDone.\n";
