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
my $OUT_SUFFIX=".sql.qltout";
my $CHUNK_SIZE="";
my $ERR_DIFF="";
my $MAX_ERRS="";
my $MAX_ERR_PERCENT="";
my $MAX_FREQ="";
my $UNIQUE="";
my $UNIQUE_ONLY="";
my $MAX_INDEL_LEN="";
my $MYSEQS="";
my $KEEPLOGS="";


# Add READABLE, UNIQUE, and other ShortQueryLookup options as needed.
GetOptions ('NUMPROCS=i' => \$NUMPROCS, 'L:s' => \$LOOKUP_TABLE,
            'BATCHQUEUE:s' => \$BATCHQUEUE,
            'CHUNK:s' => \$CHUNK_SIZE,
            'START:i' => \$START, 'END:i' => \$END,
            'ERR_DIFF:s' => \$ERR_DIFF,
            'MAX_ERRS:s' => \$MAX_ERRS,
            'MAX_ERR_PERCENT:s' => \$MAX_ERR_PERCENT,
            'MAX_FREQ:s' => \$MAX_FREQ,
            'UNIQUE_ONLY:s' => \$UNIQUE_ONLY,
            'UNIQUE:s' => \$UNIQUE,
            'MAX_INDEL_LEN:s' => \$MAX_INDEL_LEN,
            'SEQS:s' => \$MYSEQS,
            'KEEPLOGS:s' => \$KEEPLOGS,
            'O:s' => \$OUT_PREFIX, 'OUT_SUFFIX:s' => \$OUT_SUFFIX) ||
die "Could not process options correctly";

sub usage {
  print "\n Usage: batchShortQueryLookup.pl HEADDIR --NUMPROCS=n --L=lookup_table --BATCHQUEUE=queuename --CHUNK=chunksize --START=first --END=last --O=out_prefix --OUT_SUFFIX=outsuffixinclleadingdot --MAX_ERRS=num --MAX_ERR_PERCENT=num --MAX_FREQ=num --UNIQUE_ONLY=bool --MAX_INDEL_LEN=num --SEQS=fullfilename --KEEPLOGS=bool\n\n";

  print " example: batchShortQueryLookup.pl /wga/scr2/ohsumit/fc/070727_7582_batch/7582.5 --NUMPROCS=20 \n";
  print "   will do what ShortQueryLookUp does except it will split it into 20 jobs on the gpfs queue.\n";

}

if (@ARGV < 1) {
  usage();
  exit;
}

my $HEAD = $ARGV[0];

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
   $OUT_PREFIX = $HEAD;
}

my $midcomm = "";
if ($ERR_DIFF ne "") {
  $midcomm = $midcomm . " ERR_DIFF=$ERR_DIFF";
}
if ($CHUNK_SIZE ne "") {
  $midcomm = $midcomm . " CHUNK=$CHUNK_SIZE";
}
if ($MAX_ERRS ne "") {
  $midcomm = $midcomm . " MAX_ERRS=$MAX_ERRS";
}
if ($MAX_ERR_PERCENT ne "") {
  $midcomm = $midcomm . " MAX_ERR_PERCENT=$MAX_ERR_PERCENT";
}
if ($MAX_FREQ ne "") {
  $midcomm = $midcomm . " MAX_FREQ=$MAX_FREQ";
}
if ($UNIQUE ne "") {
  $midcomm = $midcomm . " UNIQUE=$UNIQUE";
}
if ($UNIQUE_ONLY ne "") {
  $midcomm = $midcomm . " UNIQUE_ONLY=$UNIQUE_ONLY";
}
if ($MAX_INDEL_LEN ne "") {
  $midcomm = $midcomm . " MAX_INDEL_LEN=$MAX_INDEL_LEN";
}

my $commout;
my @jobids;

for ($i = 0; $i < $finalindex; $i = $i + $readsize, $k++) {
  my $j = $i + $readsize;
  $command = "bsub -R\"rusage[mem=4096]\" -q $BATCHQUEUE -o $OUT_PREFIX.tmp$k.log ShortQueryLookup SEQS=$seqfile LOOKUP_TABLE=$LOOKUP_TABLE OUT_PREFIX=$OUT_PREFIX.$k OUT_SUFFIX=$OUT_SUFFIX";
  $command = $command . $midcomm . " NH=False START=$i END=$j";
   # print($command . "\n");
   $commout = `$command`;
   @temparray = ($commout =~ m/(\d+)/g);
   push(@jobids, $temparray[0]);
}
if ($readrem != 0) {
  $command = "bsub -R\"rusage[mem=4096]\" -q $BATCHQUEUE -o $OUT_PREFIX.tmp$k.log ShortQueryLookup SEQS=$seqfile LOOKUP_TABLE=$LOOKUP_TABLE OUT_PREFIX=$OUT_PREFIX.$k OUT_SUFFIX=$OUT_SUFFIX";
  $command = $command . $midcomm . " NH=False START=$i END=$END";
  # print($command . "\n");
  $commout = `$command`;
  @temparray = ($commout =~ m/(\d+)/g);
  push(@jobids, $temparray[0]);
  $k++;
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
# print($command . "\n");
system($command);

# my $comm1 = "cat ";
my $comm2 = "cat ";
for ($i = 0; $i < $k; $i++) {
# Will suggested this is probably no longer needed
#  $comm1 = $comm1 . " 7582.5.$i.minAlignErrors.txt";
  $comm2 = $comm2 . " $OUT_PREFIX.$i$OUT_SUFFIX";
}
# $comm1 = $comm1 . " > 7582.5.minAlignErrors.txt.new";
$comm2 = $comm2 . " > $OUT_PREFIX$OUT_SUFFIX";
# system($comm1);
system($comm2);

if ($KEEPLOGS ne "True") {
  $command = "/bin/rm $OUT_PREFIX.tmp*.log $OUT_PREFIX.*.minAlignErrors.txt $OUT_PREFIX.*$OUT_SUFFIX";
  system($command);
}


