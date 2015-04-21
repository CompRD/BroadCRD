#!/usr/bin/perl -w
use strict;
use lib '/wga/scr2/jburton/Arachne';
use PerlRunTime;
use LSF::LSFInterface; # allows parallelization

##
# RunDetectUnassembleable
#
# A wrapper script for DetectUnassembleable.  It runs DetectUnassembleable
# on a swath of ALLPATHS data directories and gathers data on the percent
# unassembleable; then it produces charts in the file charts.txt.
#
# There are three steps:
#
# Pre-processing: Create the unipaths for DetectUnassembleable to use.
#     Depending on which files already exist in DATA, thia might
#     entail running CommonPather, MakeRcDb, and/or Unipathser.
#
# Run: Actually call DetectUnassembleable.  Loop over all species as
#     well as all values of K, ANCHOR, and BRIDGE; the values to loop
#     over can be edited below.
#
# Chart creation: Write to the file charts.txt, in $PRE/output.
#
# Josh Burton
# August 2008
#


##################################
# PARAMETRIC RANGES (EDIT THESE) #
##################################

my @Ks = ( 20, 40, 48, 64, 96, 100, 200 );
my @BRIDGEs = ( 6000, 40000 );
my @ANCHORs = ( 250, 500, 2000 );
my @species = ( 'E.coli',
		'C.jejuni',
		'R.sphaeroides',
		'C.albicans',
		'N.crassa',
		'G.aculeatus',
		'G.aculeatus2',
		'Human',
		'Chimp',
		'Mouse',
		'Dog',
		'Horse',
		'Monodelphis' );
my %big_species = ( 'E.coli' => 0,
		    'C.jejuni' => 0,
		    'R.sphaeroides' => 0,
		    'C.albicans' => 0,
		    'N.crassa' => 0,
		    'G.aculeatus' => 1,
		    'G.aculeatus2' => 1,
		    'Human' => 1,
		    'Chimp' => 1,
		    'Mouse' => 1,
		    'Dog' => 1,
		    'Horse' => 1,
		    'Monodelphis' => 1,
		    'Cavia' => 1 );

# Run directory
my $PRE = "/wga/scr2/jburton/DetectUnassembleable";
my $output_dir = "$PRE/output";
run_or_die ("mkdir -p $output_dir");

# Number of CPUs to use
my $n_cpus = 4;

# Shell command aliases
my $rm = "rm -rf --";
my $cp = "cp -pr --";

# Options
my ( $preprocess, $DU, $make_charts ) = ( 0, 0, 1 );


######################
# SCRIPT BEGINS HERE #
######################

# Set up LSF interface
my $lsf_dir = "$PRE/LSF";
system("$rm $lsf_dir/*");
my $lsf = new LSFInterface("queue" => "localhost", "outputDir" => "$lsf_dir");
my @IJids = ( );


##################
# PRE-PROCESSING #
##################

if ( $preprocess ) {
  # Set up input data for DetectUnassembleable (by running CommonPather, etc.)
  foreach my $species (@species) {
    foreach my $K (@Ks) {

      my $big = $big_species{$species};

      # Input directories
      my $DATA = $species;
      my $PDR = "PRE=$PRE DATA=$DATA RUN=";
      my $genome_size = `cat $PRE/$DATA/genome.size`;
      chomp $genome_size;

      # Pre-processing commands
      my $S_command = "cd $species; ./source; cd ..";
      my $CP_command = "CommonPather K=$K DIR_IN=$PRE/$DATA READS_IN=genome GENOME_SIZE=$genome_size";
      my $M_command = "MakeRcDb K=$K $PDR READS=genome";
      my $U_command = "Unipather K=$K $PDR READS=genome SIM=True SHOW_PLACEMENTS=True";

      # Account for bigness
      my $suffix = "";
      if ($big_species{$species}) {
	$M_command .= " FORCE_BIG=True";
	$U_command =~ s/Unipather/UnipatherBig/;
	$suffix = "_big";
      }

      # Determine which pre-processing commands are necessary
      my $command;
      if (-e "$PRE/$DATA/genome.unipaths.k$K" &&
	  -e "$PRE/$DATA/genome.unipaths.k$K.locs") {
	# i.e., no command
	$command = 'echo';
      } elsif (-e "$PRE/$DATA/genome.paths_rc.k$K" &&
	       -e "$PRE/$DATA/genome.pathsdb$suffix.k$K") {
	# run Unipather
	$command = $U_command;
      } elsif (-e "$PRE/$DATA/genome.paths.k$K") {
	# run MakeRcDb & Unipather
	$command = "$M_command; $U_command";
      } elsif (-e "$PRE/$DATA/genome.fastb") {
	# run CommonPather, MakeRcDb, Unipather
	$command = "$CP_command; $M_command; $U_command";
      } else {
	# setup source genome, then run CommonPather, MakeRcDb, Unipather
	$command = "$S_command; $CP_command; $M_command; $U_command";
      }

      # Run pre-processing commands
      print "Running: $command\n\n" unless $command eq 'echo';
      my $IJid = $lsf->submit("$command");
      push @IJids, $IJid;

      # Before submitting another job, wait until the number of active
      # jobs is less than the number of CPUs
      while ( 1 ) {
	my $n_active_jobs = 0;
	map {$n_active_jobs++ if ($lsf->getIsRunning($_)) } @IJids;
	last if ($n_active_jobs < $n_cpus);
	sleep 30;
      }
    }
  }

  # Wait for all pre-processing commands to finish
  $lsf->wait( );
}

@IJids = ( );


#######
# RUN #
#######

if ( $DU ) {

  # Run DetectUnassembleable for all cases, and gather data
  foreach my $species (@species) {
    foreach my $K (@Ks) {
      foreach my $BRIDGE (@BRIDGEs) {
	foreach my $ANCHOR (@ANCHORs) {

	  # Run DetectUnassembleable
	  my $output_file = "$output_dir/$species.$K.$BRIDGE.$ANCHOR.txt";
	  my $DU_command = "DetectUnassembleable K=$K PRE=$PRE DATA=$species ANCHOR=$ANCHOR BRIDGE=$BRIDGE OUTFILE=$output_file";
	  print "Running: $DU_command\n";

	  my $IJid = $lsf->submit("$DU_command");
	  push @IJids, $IJid;

	  # Before submitting another job, wait until the number of
	  # active jobs is less than the number of CPUs
	  while ( 1 ) {
	    my $n_active_jobs = 0;
	    map {$n_active_jobs++ if ($lsf->getIsRunning($_)) } @IJids;
	    last if ($n_active_jobs < $n_cpus);
	    sleep 5;
	  }
	}
      }
    }
  }

  # Wait for all DetectUnassembleable commands to finish
  $lsf->wait();
}


##################
# CHART CREATION #
##################

if ( $make_charts ) {

  # Grep the DetectUnassembleable output files and fill %data
  my %data;
  print "Gathering data from DetectUnassembleable output...\n";

  foreach my $species (@species) {
    foreach my $K (@Ks) {
      foreach my $BRIDGE (@BRIDGEs) {
	foreach my $ANCHOR (@ANCHORs) {

	  # Grep the output file and get the % unassembleable
	  my $output_file = "$output_dir/$species.$K.$BRIDGE.$ANCHOR.txt";
	  my @lines = `cat $output_file`;
	  my $pct = '-';
	  foreach (@lines) {
	    if (/\% bad vs total genome size\: (.*)\n/) {
	      $pct = $1;
	    }
	  }

	  $data{$species}{$K}{$BRIDGE}{$ANCHOR} = $pct;
	}
      }
    }
  }

  # Print data in chart file
  my $chart_file = "$output_dir/charts.txt";
  print "Writing to chart file $chart_file\n\n";
  open CHART, '>', $chart_file
    or die "Can't open $chart_file: $!";
  print CHART "PERCENT UNASSEMBLEABLE\n\n\n\n";

  foreach my $BRIDGE (@BRIDGEs) {
    foreach my $ANCHOR (@ANCHORs) {

      print CHART "BRIDGE: $BRIDGE\nANCHOR: $ANCHOR\n";

      # Header row
      print CHART "\tK=";
      map {print CHART "\t$_"} (@Ks);
      print CHART "\n";

      # Chart rows
      foreach my $species (@species) {
	printf CHART "%-16s", $species;
	foreach my $K (@Ks) {
	  print CHART $data{$species}{$K}{$BRIDGE}{$ANCHOR}, "\t";
	}
	print CHART "\n";
      }
      print CHART "\n\n\n";
    }
  }

  close CHART;
}
