#!/usr/bin/perl -w
use strict;

# ExtractFromBAM.pl
#
# A script to extract reads from the BAM files in the Picard pipeline and place
# them in new {fastb,qualb,pairs} files.  Uses samtools view, SAM2CRDDump, 
# FastbStats, and MergeReadSets.
#
# Arguments
#
# LANES: The set of flowcell/lane combinations to extract from.  Nesting is
#        supported, e.g., "{42EWK.{2,3},42EYC.{6,7,8}}" expands to the set
#        "42EWK.2, 42EWK.3, 42EYC.6, 42EYC.7, 42EYC.8".
# VERSION: If "oldest", use the oldest Picard directory; otherwise the newest.
# RANGE: If specified, only extract reads that match the reference in this
#        range - e.g., "scaffold_9.1-166029:0-150000".
# SEP, DEV: Passed to SAM2CRDDump.  At present, all read pairs get this sep/dev.
# TARGET_HEAD: The output files will be at <TARGET_HEAD>.{fastb,qualb,pairs}.
# DRY_RUN: If you this to 1 (true), ExtractFromBAM will run through its commands
# but not execute any of them.
# PICARD_HEAD: should be /seq or /wga/scr3.  The default is /seq.
# TMP: Temporary directory.
#
# RUN EXAMPLE:
# ExtractFromBAM.pl LANES="{42EWK.{2,3},42EYC.{6,7,8}}" RANGE=scaffold_9.1-166029:0-150000 TARGET_HEAD=ECOP15
#
# Josh Burton
# March 2010

use FindBin;

# Local modules
use lib "$FindBin::Bin/../";
use lib "$FindBin::Bin/";

use ArachneArgs; # allows parsing of command-line arguments
use PerlRunTime; # allows function 'run_or_die'

# Subroutine declarations.
sub get_newest_dir(@);
sub get_oldest_dir(@);

# This convoluted syntax is how we document arguments.
my %lanes_arg = ('value' => undef,
		 'help' => 'e.g., "{42EWK.{2,3},42EYC.{6,7,8}}"');
my %version_arg = ('value' => 'newest',
                   'help' => 'oldest or newest');
my %range_arg = ('value' => '',
		 'help' => 'If specified, only take reads in this range. Format: <scaffold name>:<begin>-<end>');
my %target_arg = ('value' => 'reads',
		  'help' => 'The output files will be at <TARGET_HEAD>.{fastb,qualb,pairs}');
my %readlen_arg = ('value' => 0,
		  'help' => 'use read length to obtain paired reads separation instead of insert sizes');
my %picard_head_arg = ('value' => '/seq',
		  'help' => 'set to either /seq or /wga/scr3');
my %tmp_arg = ('value' => undef,
		  'help' => 'set to location of directory for temporary files');
my %require_libinfo_arg = ('value' => 0,
		  'help' => 'Library names must be present in LIBINFO file - unknown librarys not allowed');
my %unmapped_arg = ('value' => 0,
		  'help' => 'Use unmapped bam file, not the mapped one');
my %libinfo_arg = ('value' => '',
                   'help' => 'Library information file, used by SAM2CRDDump to detemine library stats');
my %write_aligns_arg = ('value' => 1,
		    'help'  => 'Write alignments as .qltout and .mapq files ');
my %write_names_arg = ('value' => 1,
		    'help'  => 'Write original reads names as .names file ');
my %write_pairs_arg = ('value' => 1,
		    'help'  => 'Write pairing information as .pairs file ');
my %mapped_pairs_only_arg = ('value' => 0,
		    'help'  => 'Use only pairs with both read mapped');


# Get command-line arguments.
my %args = getCommandArguments('LANES' => \%lanes_arg,
                               'VERSION' => \%version_arg,
			       'RANGE' => \%range_arg,
			       'TARGET_HEAD' => \%target_arg,
                               'NOMINAL_READ_LEN' => \%readlen_arg,
                               'PICARD_HEAD' => \%picard_head_arg,
			       'LIBINFO' => \%libinfo_arg,
			       'REQUIRE_LIBINFO' => \%require_libinfo_arg,
			       'UNMAPPED' => \%unmapped_arg,
			       'WRITE_ALIGNS' => \%write_aligns_arg,
			       'WRITE_NAMES' => \%write_names_arg,
			       'WRITE_PAIRS' => \%write_pairs_arg,
                               'MAPPED_PAIRS_ONLY' => \%mapped_pairs_only_arg,
                               'TMP' => \%tmp_arg,
			       'SEP' => -15,
			       'DEV' => 12,
			       'DRY_RUN' => 0 );


my $require_libinfo_str = $args{REQUIRE_LIBINFO} ? "True" : "False";
my $write_aligns_str = $args{WRITE_ALIGNS} ? "True" : "False";
my $write_names_str = $args{WRITE_NAMES} ? "True" : "False";
my $write_pairs_str = $args{WRITE_PAIRS} ? "True" : "False";
my $mapped_pairs_only_str = $args{MAPPED_PAIRS_ONLY} ? "True" : "False";

my $bam_name = $args{UNMAPPED} ? "unmapped.bam" : "aligned.duplicates_marked.bam";

# Parse the LANES argument into a set of flowcell-lane combinations.
my (@flowcells, @lanes, @libraries);

# LANES may consist of several arguments, split up by the command-line parser,
# so we must parse them separately.
foreach my $lane ( @{$args{'LANES'}} ) {
    die "ERROR: Nonsensical flowcell/lane name \"$lane\" in LANES braced set $_\n"
        unless ( $lane =~ /^[\'\"]?(\w{5,9})\.(\d)(?:\.([\.\w\d-]+))*[\'\"]?$/ );
    push @flowcells, $1;
    push @lanes, $2;
    if ($3) {
      push @libraries, $3;
    } else {
      push @libraries, "*";
    }
}

my $n_lanes = scalar (@flowcells);
die "ERROR: No lanes specified\n" if ( $n_lanes == 0 );


# Print information about the run we are about to do.
print "THIS IS A DRY RUN\n" if ( $args{'DRY_RUN'} );
print "USING $n_lanes LANES:\n\n";
map { print "  $flowcells[$_].$lanes[$_].$libraries[$_]\n" } (0..$n_lanes-1);
print "\n";
if ($args{'RANGE'} eq '') {
    print "Extracting every read in these lanes.  WARNING: this will take a while.\n";
}
else {
    print "Extracting only reads in the range $args{'RANGE'}.\n";
}
print "\n";

my $samtools = 'samtools';

# Prepare temp directory.
my $temp_dir = $args{'TMP'};
mkdir( $temp_dir );


my $n_reads = 0;

# Loop over all lanes.
foreach ( 0..$n_lanes-1 ) 
{
    my ($fc, $lane, $library) = ($flowcells[$_], $lanes[$_], $libraries[$_]);
    
    # Find the Picard directory corresponding to this lane.
    my @Picard_dirs = glob "$args{'PICARD_HEAD'}/picard/${fc}/*/$lane/$library";
    die "ERROR: Can't find dir: $args{'PICARD_HEAD'}/picard/${fc}/*/$lane/$library\n" .
        "       for library $fc.$lane.$library\n"
        unless ( @Picard_dirs );
    my $Picard_dir;
    if ( $args{'VERSION'} eq 'newest' ) 
    {    $Picard_dir = get_newest_dir( @Picard_dirs );    }
    else
    {    $Picard_dir = get_oldest_dir( @Picard_dirs );    }
    
    # Determine the library name, if we don't already know it
    die "Unable to determine library name from $Picard_dir\n"
        unless ($Picard_dir =~ /^.*\/(.+)$/ );
    $library = $1;
    $libraries[$_] = $library;
    my $fullname = "$fc.$lane.$library";
    
    print localtime() . ": Retrieving lane $fullname from $Picard_dir\n";
    
    # Find the BAM file name, ignoring any other BAM files present. See arg UNMAPPED.
    my $BAM = "$Picard_dir/${fc}.$lane.$bam_name";
    die "ERROR: Can't find BAM file at $BAM\n" unless ( -e $BAM );
    print "Using BAM: $bam_name\n";
    
    # Use samtools view to grab the wanted range of data out of the BAM file.
    # Then, use SAM2CRDDump to parse these reads into CRD format (fastb, etc.)
    my $out_head = "$temp_dir/$$.$fullname";
    my $samtools_cmd = "$samtools view -h $BAM $args{'RANGE'}";
    my $SAM2CRDDump_cmd = "SAM2CRDDump OUT_HEAD=$out_head SEP=$args{'SEP'} DEV=$args{'DEV'} NOMINAL_READ_LEN=$args{'NOMINAL_READ_LEN'} USE_OQ=True NH=True LOG_TO_CERR=False LIBINFO=$args{'LIBINFO'} REQUIRE_LIBINFO=$require_libinfo_str WRITE_ALIGNS=$write_aligns_str WRITE_NAMES=$write_names_str WRITE_PAIRS=$write_pairs_str MAPPED_PAIRS_ONLY=$mapped_pairs_only_str";
    if ( $args{'DRY_RUN'} ) {
	print "DRY RUN: $samtools_cmd | $SAM2CRDDump_cmd\n";
	next;
    }
    run_or_die( "$samtools_cmd | $SAM2CRDDump_cmd" );
    
    
    # Report on the number of reads we've extracted.
    my $FastbStats = run_or_die( "FastbStats FASTB=$out_head.fastb QUICK=True" );
    if ($FastbStats =~ /total objects count\: (\d+)/) {
        $n_reads += $1;
        print "\tReads extracted: $1\n";
    }
    else {
        print "$FastbStats\n";
        die "ERROR: No reads found in lane $fullname (did you specify a valid range?)\n";
    }
}


# Prepare an argument to MergeReadSets.
my @all_lanes;
map { push @all_lanes, "$temp_dir/$$.$flowcells[$_].$lanes[$_].$libraries[$_]" } (0..$n_lanes-1);
my $all_lanes = join( ',', @all_lanes );


# If there is only one lane, we simply copy its files to the appropriate place.
if ( $n_lanes == 1 ) {
    
    # Set up the copy commands.
    my @cp_files = ("fastb", "qualb" );
    if ($args{'WRITE_PAIRS'}) { push @cp_files, "pairs"; }
    if ($args{'WRITE_NAMES'}) { push @cp_files, "names"; }
    if ($args{'WRITE_ALIGNS'}) { push @cp_files, "qltout", "mapq"; }
    my @cp_cmds;
    foreach ( @cp_files ) {
	push @cp_cmds, "cp $temp_dir/$$.$flowcells[0].$lanes[0].$libraries[0].$_ $args{'TARGET_HEAD'}.$_";
    }
    foreach ( @cp_cmds ) {
	if ( $args{'DRY_RUN'} ) { print "DRY RUN: $_\n" }
	else { run_or_die( "$_" ) }
    }
}

# If there are multiple lanes, we must call MergeReadSets to merge them.
else {
    
    # Call MergeReadSets to merge all of the datasets we have extracted.
    print localtime() . ": Using MergeReadSets to merge the extracted data files\n";
    my $MergeReadSets_cmd = "MergeReadSets NO_HEADER=True READS_IN=\"{$all_lanes}\" READS_OUT=$args{'TARGET_HEAD'} TRACK_READS=False";
    if ( $args{'DRY_RUN'} ) {
	print "DRY RUN: $MergeReadSets_cmd\n";
    }
    else {
	run_or_die( "$MergeReadSets_cmd" );
    }
}


# Clean up the temp files.
run_or_die( "rm -f $temp_dir/$$.*" ) unless $args{'DRY_RUN'};
map { run_or_die( "rm -f $_.*" ) } @all_lanes unless $args{'DRY_RUN'};



# Done!
print "\n";
if ( $args{'DRY_RUN'} ) {
    print localtime() . ": DRY RUN COMPLETE\n\n";
}
else {
    print localtime() . ": Done!\n";
    print "ExtractFromBAM has extracted $n_reads reads out of $n_lanes lanes into files at $args{'TARGET_HEAD'}.\{fastb,qualb,pairs,qltout,mpaq\}\n\n";
}



# SUBROUTINES

# get_newest_dir: From the input directory names, choose the one with the newest
# time stamp.
sub get_newest_dir(@)
{
    my ($min_M, $dir) = (1e8, '');
    # -M <pathname> = number of days since modification time of <pathname>
    foreach (@_ ) {
	if ( $min_M > -M $_ ) {
	    $min_M = -M $_;
	    $dir = $_;
	}
    }
    return $dir;
}

# get_oldest_dir: From the input directory names, choose the one with the oldest
# time stamp.
sub get_oldest_dir(@)
{
    my ($max_M, $dir) = (0, '');
    # -M <pathname> = number of days since modification time of <pathname>
    foreach (@_ ) {
	if ( $max_M < -M $_ ) {
	    $max_M = -M $_;
	    $dir = $_;
	}
    }
    return $dir;
}
