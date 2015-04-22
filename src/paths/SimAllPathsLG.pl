#!/usr/bin/perl -w
use strict;
use lib "/wga/scr2/jburton/Arachne"; # hard-wired connection for local libraries
use PerlRunTime; # run_or_die

# SimAllPathsLG.pl
#
# SimAllPathsLG: A tool to assemble a complete genome using simulated data.
# Mostly this consists of a wrapper to RunAllPathsLG.
#
# For syntax info, run SimAllPathsLG.pl with no arguments.
#
# Josh Burton
# Work begun 2009-07-28



# Hash of scenarios.  Each scenario name is keyed to a DATA directory and a SNP
# rate (0 if haploid.)
# DATA must be an existing subdirectory of ALLPATHS_BASE=/wga/scr4/ALLPATHS and
# must contain files genome.fastb and genome.size.
# Feel free to add to this list as necessary!
my %scenarios = (
    'E.coli' => [ "E.coli", 0 ],
    'Neurospora' => [ "N.crassa", 0 ],
    'Stickleback' => [ "G.aculeatus.bearpaw", 0.0025 ],
    'Bushbaby' => [ "O.garnettii", 0.001 ],
);


# Forward definitions.
sub print_help();

# PROCESSING BEGINS HERE

# If not enough arguments are given, explain syntax and quit.
if (scalar @ARGV < 2) {
    print_help();
    exit 0;
}


# Parse command-line arguments and load scenario.
my ($scenario, $label) = @ARGV;
die "ERROR: Scenario '$scenario' unknown.  If you want to use this scenario,\nyou can add it to the hash '%scenarios' in $0."
    unless exists $scenarios{$scenario};
my $data_rel = $scenarios{$scenario}[0]; # data dir, relative to ALLPATHS_BASE
my $SNP_rate = $scenarios{$scenario}[1];
my $ALLPATHS_BASE = "$ENV{'ARACHNE_PRE'}/projects/ALLPATHS";
my $data_dir = "$ALLPATHS_BASE/$data_rel";

# Verify that this scenario makes sense.
die "ERROR: invalid scenario!  Data directory '$data_dir' must already exist\n"
    unless (-d $data_dir);
die "ERROR: invalid scenario!  Data directory '$data_dir' must contain a file genome.size\n"
    unless (-f "$data_dir/genome.size");
die "ERROR: invalid scenario!  Data directory '$data_dir' must contain a file genome.fastb\n"
    unless (-f "$data_dir/genome.fastb");
die "ERROR: invalid scenario!  Can't parse SNP rate '$SNP_rate' as a number\n"
    if ($SNP_rate =~ /[^0-9\.]/);
die "ERROR: invalid scenario!  SNP rate must be less than 1\n"
    if ($SNP_rate >= 1);

# Find the 'label' subdirectory.
my $data_subdir = "$data_dir/$label";
die "ERROR: Data subdirectory '$data_dir/$label' must NOT already exist\n"
    if (-e "$data_subdir");


# Deduce the ploidy from the SNP rate.
my $ploidy = $SNP_rate == 0 ? 1 : 2;


# Announce the beginning of this run.
print "$0 @ARGV\n\n";
print localtime() . ": Begun SimAllPathsLG!\n\n";
print "Using scenario '$scenario':\n";
print "DATA = $data_dir\n";
print "Ploidy = $ploidy\n";
if ($ploidy == 2) {
    my $inverse_SNP_rate = 1 / $SNP_rate;
    print "SNP rate = $SNP_rate (1 in $inverse_SNP_rate)\n";
}
print "\n";
print "Input is genome.fastb file in $data_dir.\n";
print "Output is going to new directory $data_subdir.\n";


# Create the data subdirectory and put a simple README into it.
run_or_die( "mkdir -p $data_subdir" );
my $readme = "This directory was created by SimAllPathsLG.pl at " . localtime() .
    ".\nThe genome.fastb file is from $data_subdir.";
run_or_die( "echo '$readme' > $data_subdir/README" );

# Copy genome.fastb into the data directory.
run_or_die( "cp -rf $data_dir/genome.fastb $data_subdir" );

# Get genome size.  This may not match up with the size of genome.fastb.
my $genome_size = `cat $data_dir/genome.size`;
chomp $genome_size;

# Create auxiliary data structures for the reference genome.
run_or_die( "echo '$ploidy' > $data_subdir/ploidy" );
run_or_die( "MakeLookupTable SOURCE=$data_subdir/genome.fastb OUT_HEAD=$data_subdir/genome LOOKUP_ONLY=True" );
run_or_die( "GenomeUnipathStats K=96 GENOME=$data_subdir/genome DOT=$data_subdir/genome.dot" );

# Add SNPs to the reference genome if necessary.
if ($ploidy > 1) {
    run_or_die( "MutateReference FASTB_IN=$data_dir/genome.fastb FASTB_OUT=$data_dir/genome.diploid.fastb SUB_RATE=$SNP_rate" );
}


# Get date in format YYYY-MM-DD.  This is used as the name of the RUN directory.
my $today = `date +"%Y-%m-%d"`;
chomp $today;


# Set up fragment/jumping libraries and error generator for simulation.
# Libraries (including coverages, etc.) are hard-wired.
my $jlib = "dev=10%,jmean=205,jdev=0%,fc_freq=0,ic_freq=0.02,rc_avoid=True";
my $frag_lib  = "n=100,N=180,dev=10%,C=30";
my $jump_lib  = "n=100,N=4000,C=10,$jlib:n=100,N=40000,C=0.2,$jlib";
my $error_gen = "error_templates/template_30BCP.8_100";


# Call RunAllPathsLG.  This is the main task, and it will take a while!
my $RAPLG = "RunAllPathsLG REFERENCE_NAME=$data_rel DATA_SUBDIR=$label RUN=$today K=80 FRAG_K=28 OVERWRITE=True MAXPAR=3 THREADS=crd-32 EVALUATION=FULL FRAG_LIBS=$frag_lib JUMP_LIBS=$jump_lib ERROR_GENERATOR_NAME=$error_gen TARGETS=full_eval GENOME_SIZE=$genome_size DO_UNIPATH_EVAL=False";

run_or_die( "$RAPLG" );


# Done!
print "\n", localtime() . ": SimAllPathsLG is done!\n\n";













# Subroutine: print usage syntax and options.
sub print_help() {
    my @scenarios = sort keys %scenarios;
    print <<"SYNTAX_END";
    
SimAllPathsLG: A tool to assemble a complete genome using simulated data.

[01mSYNTAX[0m
\tSimAllPathsLG.pl <scenario> <label>

[01mOPTIONS[0m

'scenario' must be one of:
\t[01m@scenarios[0m
\tIf you want to run a scenario that is not on this list, you will need to
\tadd it to the hash '%scenarios' in $0.
'label' can be any name you want.  The directory 'label' will be created
\tin the data directory of this scenario (it must not already exist) and
\tall the output of SimAllPathsLG.pl will go there.

SYNTAX_END
}
