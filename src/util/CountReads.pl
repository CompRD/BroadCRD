#!/util/bin/perl -w
use strict;



##
#
# CountReads.pl
#
#
# It counts the number of reads for a given project in a directory.
# The idea is that usually reads are stored inside plate-named
# directories.
#
#
# How to use (examples):
#  CountReads.pl G592 /seq/comparative11/G592/chromat_dir/hrs
#  CountReads.pl G592 /seq/comparative11/G592/seq_dir/hrs
#  CountReads.pl G592 /seq/comparative11/G592/seq_dir/hrs
#  etc.
#
#



# Parse arguments.
my $project = $ARGV[0];
my $data_dir = $ARGV[1];

my $plates_count = 0;
my $reads_count = 0;

my @all_plates;
my @all_reads;
my $plate;
my $read;

# Project directory.
opendir( DIR, $data_dir ) || die "can't opendir $data_dir: $!";
@all_plates = readdir( DIR );

# Each entry is a plate directory.
foreach $plate ( @all_plates ) {
  if ( $plate =~ /$project/ ) {
    $plates_count += 1;

    my $plate_dir = "$data_dir" . "/$plate";
    opendir( PLATE_DIR, $plate_dir ) || die "can't opendir $plate_dir: $!";
    @all_reads = readdir( PLATE_DIR );
    foreach $read ( @all_reads ) {
      if ( $read =~ /$project/ ) {
	$reads_count += 1;
      }
    }
    
    closedir PLATE_DIR;
  }
}

closedir DIR;  

# Print out.
print "There are $plates_count plates, and $reads_count reads.\n";
