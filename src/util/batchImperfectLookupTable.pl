#!/usr/bin/perl -w

use POSIX qw(floor);

###########################################################################
#                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   #
#       This software and its documentation are copyright (2007) by the   #
#   Broad Institute/Massachusetts Institute of Technology.  All rights    #
#   are reserved.  This software is supplied without any warranty or      #
#   guaranteed support whatsoever. Neither the Broad Institute nor MIT    #
#   can be responsible for its use, misuse, or functionality.             #
###########################################################################


#Create solexa files from a Bustard directory.

# @cmdline = @ARGV;
  
use Getopt::Long;
use IO::Handle;
use File::Basename;
use Cwd 'abs_path';
  


my $START=0;
my $END=-1;
my $NUMPROCS=1;
my $BATCHQUEUE="broad";
my $LOOKUP_TABLE="";
my $IN_SUFFIX=".fastb";
my $OUT_PREFIX="";
my $OUT_SUFFIX=".ilt.qltout";
my $CHUNK_SIZE="";
my $K="";
my $ERR_DIFF="";
my $MAX_ERRS="";
my $WRITE_MINALIGNS="";
my $PRINT="";
my $GENOME="";
my $FWRC="";
my $QUALS="";
my $MYSEQS="";
my $KEEPLOGS="";
my $TIMER="";
my $STOP_CAT=0;
my $MODE="";
my $DRY_RUN = 0;
my $MAX_FREQ = "";

# STOP_CAT stops the concatentation.  For large jobs (> few hundred), the 
# date bsub will fail because the command line argument will be too long.
# So, it would have to be done manually.

# Add , and other ShortQueryLookup options as needed.
GetOptions ('NUMPROCS=i' => \$NUMPROCS, 
	    'L:s' => \$LOOKUP_TABLE,
            'BATCHQUEUE:s' => \$BATCHQUEUE,
            'CHUNK:s' => \$CHUNK_SIZE, 
	    'K:s' => \$K,
            'START:i' => \$START, 
	    'END:i' => \$END,
            'ERR_DIFF:s' => \$ERR_DIFF,
            'FWRC:s' => \$FWRC,
            'WRITE_MINALIGNS:s' => \$WRITE_MINALIGNS,
            'MAX_ERRS:s' => \$MAX_ERRS,
            'PRINT:s' => \$PRINT,
            'GENOME:s' => \$GENOME,
            'QUALS:s' => \$QUALS,
            'SEQS:s' => \$MYSEQS,
            'TIMER:s' => \$TIMER,
            'STOP_CAT' => \$STOP_CAT,
            'KEEPLOGS:s' => \$KEEPLOGS,
	    'MODE:s' => \$MODE,
            'O:s' => \$OUT_PREFIX, 
	    'OUT_SUFFIX:s' => \$OUT_SUFFIX,
	    'MAX_FREQ:i' => \$MAX_FREQ,
	    'DRY_RUN' => \$DRY_RUN) ||
die "Could not process options correctly";

# STOP_CAT exits the program after submitting the jobs.  The point is that
# for NUMPROCS that are very large, the synchronization barrier "date" batch
# job will fail since there are too many stop conditions (too long of a
# command line).  Thus, we stop and do it manually later.

sub usage {
  print "\n Usage: batchImperfectLookupTable.pl HEAD --NUMPROCS=n --L=lookup_table --BATCHQUEUE=queuename --CHUNK=chunksize --START=first --END=last --MODE=mode --MAX_FREQ=maxfreq --O=out_prefix --OUT_SUFFIX=outsuffixinclleadingdot --K=number --ERR_DIFF=number --FWRC=number --WRITE_MINALIGNS=bool --MAX_ERRS=number --PRINT=bool --GENOME=string --QUALS=string --KEEPLOGS=bool --STOP_CAT\n\n";

  print "Defaults explicitly set in this script:\n";
  print "  L (lookup table) = HEAD.reference.lookup\n";
  print "  OUT_PREFIX       = HEAD\n";
  print "For other parameters, native defaults of ImperfectLookupTable will be used\n\n";

  print " example: batchImperfectLookupTable.pl /wga/scr2/ohsumit/fc/070727_7582_batch/7582.5 --NUMPROCS=20 --K=12\n";
  print "   will do what ImperfectLookupTable does except it will split it into 20 jobs on the gpfs queue.\n";

}

if (@ARGV < 1) {
  usage();
  exit;
}

my $HEAD = $ARGV[0];

my $metricfile = $HEAD . ".metrics";
my $numReads;
my $seqfile;
if ($MYSEQS ne "") {
  $seqfile = $MYSEQS;
} else {
  $seqfile = $HEAD . $IN_SUFFIX;
}

if ($START != 0) {
   print "Option not implemented yet.  Only start = 0.\n";
   exit(1); 
}


if ($END == -1) {
   my $feudallines = `FeudalSize FILE=$seqfile NH=True`;
   @myarray = ($feudallines =~ m/(\d+)/g);
   $numReads = $myarray[0];
} else {
   $numReads = $END;
}



my $readsize = floor($numReads/$NUMPROCS);
my $finalindex = $NUMPROCS*$readsize;
my $readrem = $numReads % $NUMPROCS;


$|=1; #Flush all output immediately.


my $i; 
my $k = 0;
my $command;

if ($LOOKUP_TABLE eq "") {
  $LOOKUP_TABLE = $HEAD . ".reference.lookup";
}

if ($OUT_PREFIX eq "") {
  die ("At least one of OUT_PREFIX and HEAD arguments must be specified" ) if ( $HEAD eq "" );
  $OUT_PREFIX = $HEAD;
}

my $midcomm = "";
if ($K ne "") {
  $midcomm = $midcomm . " K=$K";
}
if ($CHUNK_SIZE ne "") {
  $midcomm = $midcomm . " CHUNK=$CHUNK_SIZE";
}
if ($FWRC ne "") {
  $midcomm = $midcomm . " FWRC=$FWRC";
}
if ($WRITE_MINALIGNS ne "") {
  $midcomm = $midcomm . " WRITE_MINALIGNS=$WRITE_MINALIGNS";
}
if ($MAX_ERRS ne "") {
  $midcomm = $midcomm . " MAX_ERRS=$MAX_ERRS";
}
if ($PRINT ne "") {
  $midcomm = $midcomm . " PRINT=$PRINT";
}
if ($GENOME ne "") {
  $midcomm = $midcomm . " GENOME=$GENOME";
}
if ($TIMER ne "") {
  $midcomm = $midcomm . " TIMER=$TIMER";
}
if ($MAX_FREQ ne "") {
  $midcomm = $midcomm . " MAX_FREQ=$MAX_FREQ";
}
if ($QUALS ne "") {
  $midcomm = $midcomm . " QUALS=$QUALS";
}
if ($MODE ne "") {
  $midcomm = $midcomm . " MODE=$MODE";
}

my $commout;
my @jobids;



for ($i = 0; $i < $finalindex; $i = $i + $readsize, $k++) {
  my $j = $i + $readsize;
  $command = "bsub -R \"rusage[mem=4096]\" -q $BATCHQUEUE -o $OUT_PREFIX.tmp$k.log ImperfectLookupTable SEQS=$seqfile L=$LOOKUP_TABLE O=$OUT_PREFIX.$k OUT_SUFFIX=$OUT_SUFFIX";
  $command = $command . $midcomm . " NH=False START=$i END=$j";
  if ( $DRY_RUN ) {
    print($command . "\n");
  } else {
    $commout = `$command`;
    @temparray = ($commout =~ m/(\d+)/g);
    push(@jobids, $temparray[0]);
  }
}

if ($readrem != 0) {
  $command = "bsub -R \"rusage[mem=4096]\" -q $BATCHQUEUE -o $OUT_PREFIX.tmp$k.log ImperfectLookupTable SEQS=$seqfile L=$LOOKUP_TABLE O=$OUT_PREFIX.$k OUT_SUFFIX=$OUT_SUFFIX";
  $command = $command . $midcomm . " NH=False START=$i END=$END";
  if ( $DRY_RUN ) {
    print($command."\n");
  } else {
    $commout = `$command`;
    @temparray = ($commout =~ m/(\d+)/g);
    push(@jobids, $temparray[0]);
  }
  $k++;
}

if ( $DRY_RUN ) {
  print "\nThis was a dry run, no jobs were submitted to the queue. Good bye.\n";
  exit(1);
}

my $bcommid = "";
my $firsttime = 1;
foreach $jobid (@jobids) {
  if ($firsttime == 1) {
    $firsttime = 0;
  } else {
    $bcommid = $bcommid . " && ";
  }
  # We do ended below instead of done since if there is a problem with the
  # job (either the command or the queue) this script will hang indefinitely
  # on the next bsub command.
  $bcommid = $bcommid . " ended($jobid) ";
}


$command = "bsub -q $BATCHQUEUE -o $OUT_PREFIX.tmp$k.log -K -w \"$bcommid\" date";
#print($command . "\n");
system($command);


if ($STOP_CAT != 0) {
  exit(0);
}


my $comm1 = "cat ";
my $comm2 = "cat ";
for ($i = 0; $i < $k; $i++) {
  $comm1 = $comm1 . " $OUT_PREFIX.$i.minAlignErrors.txt";
  $comm2 = $comm2 . " $OUT_PREFIX.$i$OUT_SUFFIX";
}
$comm1 = $comm1 . " > $OUT_PREFIX.minAlignErrors.txt";
$comm2 = $comm2 . " > $OUT_PREFIX$OUT_SUFFIX";

if ( $WRITE_MINALIGNS eq "True" ) {
  system($comm1);
}
system($comm2);

if ($KEEPLOGS ne "True") {
  $command = "/bin/rm $OUT_PREFIX.tmp*.log $OUT_PREFIX.*.minAlignErrors.txt $OUT_PREFIX.*$OUT_SUFFIX";
  system($command);
}

$command = `date`;
$command =~ s/(^\s+)|(\s+$)//g;

# Need this for David's scripts
print $command . ": done.";

