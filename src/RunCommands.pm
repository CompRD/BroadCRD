# RunCommands.pm
#
# A Perl module to allow for large-scale scripting of Arachne modules,
# through the parsing of START and STOP commands.
#
# prepare_commands: Prepares a set of commands (e.g., Arachne modules)
#     for execution, parsing the arguments START and STOP.  Returns the same
#     list, but shortened in accordance with START/STOP.
# explicitize: Makes each command in the list "explicit" by appending it
#     with a number indicating its appearance in the list (i.e., the first
#     instance of ModuleName is renamed ModuleName.1, the second instance is
#     ModuleName.2, etc.)  These names are the ones used by START and STOP.
#     If $longform = 1, the commands keep their arguments; otherwise
#     (by default) they are pruned down to the explicitized command itself.
#
# The effect of START and STOP is to denote a section of the list to be run.
# START and/or STOP (you can give either, or both) are command names.  The
# names are explicitized, so you can specify exactly which command you want to
# start with with (i.e., "START=ModuleName.1").  run_commands will run
# everything from START to STOP, inclusive.
#
# Syntax:
#
# @cmds = prepare_commands(\@cmds, \%args);
#
# @cmds = explicitize(\@cmds, \%args, $longform);
# @cmds = explicitize(\@cmds, \%args); # $longform defaults to 0
# @cmds = explicitize(\@cmds, (), $longform); # to explicitize the whole list
#
# Be sure to pass @cmds and %args by reference!  If you write:
# prepare_commands(@cmds);
# you will get an error: "Can't use string as an ARRAY ref..."
#
#
# Josh Burton
# May 2008


package RunCommands;
use strict;
use Cwd; # for 'getcwd'

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);
require Exporter;
@ISA = qw(Exporter AutoLoader);
@EXPORT = qw(prepare_commands explicitize relative_link);




# See top of module for documentation
sub prepare_commands {
    my ($commands, $args) = @_;
    
    # Get the range from START to STOP
    my ($start_id, $stop_id) = find_start_stop($commands, $args);
    return @$commands[$start_id..$stop_id];
}



# See top of module for documentation
sub explicitize {
    my ($commands, $args, $longform) = @_; # $longform defaults to 0
    
    # Get the range from START to STOP
    my ($start_id, $stop_id) = find_start_stop($commands, $args);
    
    # Find the explicitized form of the commands in this range
    my @explicit_cmds = explicitize_array(@$commands);
    
    # If not using longform, prune the commands to remove their arguments
    if (!$longform) {
	map {s/\s.+$//} @explicit_cmds;
    }
    
    return @explicit_cmds[$start_id..$stop_id];
}



# Determine the range of commands in @commands to run, according to
# the given values (if any) of START and STOP in %args
sub find_start_stop {
    my ($commands, $args) = @_;
    my @commands = $commands ? @$commands : ();
    my %args     = $args     ? %$args     : ();
    
    die "find_start_stop: Received an empty \@commands list.\n"
	if (scalar @commands == 0);
    
    # Explicitize these commands' names
    my $n_cmds = scalar @commands;
    my @explicit_cmds = explicitize_array(@commands);
    
    # Get START, STOP
    my ($start, $stop) = ($args{'START'}, $args{'STOP'});
    my ($start_id, $stop_id) = (0, $n_cmds-1);
    
    if ($start) {
	# Find the command whose name starts with START
	while ($explicit_cmds[$start_id] !~ /^$start/) {
	    $start_id++;
	    die "find_start_stop: Can't find START command ($start).\n"
		if ($start_id == $n_cmds);
	}
    }
    if ($stop) {
	# Find the command whose name starts with STOP
	while ($explicit_cmds[$stop_id] !~ /^$stop/) {
	    $stop_id--;
	    die "find_start_stop: Can't find STOP command ($stop).\n"
		if ($stop_id < 0);
	}
    }
    
    die "find_start_stop: Your STOP command ($stop) appears earlier in the command sequence than your START command ($start).\n"
	if ($stop_id < $start_id);
    
    # Return the range from START to STOP
    # (or the whole range, if no START/STOP were given)
    return ($start_id, $stop_id);
}



# Useful with the command-line arguments NOGO, START, STOP
sub explicitize_array {
    my @commands = @_;
    my %ids = ();
    
    foreach my $command (@commands) {
	# Find the first argument (space-delineated)
	$command =~ /^(\S+)\s/;
	my $first_arg = $1;
	
	# Change the value of $command (thereby changing @commands)
	$ids{$first_arg}++;
	$command =~ s/^$first_arg/$first_arg.$ids{$first_arg}/;
    }
    
    return @commands;
}





# TODO: These commands are ultimately going to be moved out of RunCommands and
# into a separate Perl module (UsefulShellLikeFunctions.pm? we'll see.)

# Input: $source, $target <, $root>
# Make a symbolic link from $source to $target, defined relatively
# (i.e. using ../../ instead of hardwiring the entire path)
# By default, uses the smallest number of `../'s, but if $root is given,
# uses that as the root
sub relative_link {
    my ($source, $target, $root) = @_;
    
    if ($source eq $target) {
	warn "relative_link: source and target files are the same (`$source')\n";
	return;
    }
    
    $root = '' unless $root;
    # First: Make the directories absolute, if they aren't already
    foreach ($source, $target, $root) {
	$_ = getcwd . "/$_" if /^[^\/]/;
    }
    
    # If $root was given: Use it
    if ($root) {
	die "relative_link: force_root $root must be a root directory for source and target\n"
	    unless ($source =~ /^$root/ && $target =~ /^$root/);
    }
    # Otherwise: Find the shared root of the source and target file
    else {
	my @s_dirs = split /\//, $source;
	my @t_dirs = split /\//, $target;
	my $n = 0;
	$n++ while ($s_dirs[$n] eq $t_dirs[$n]);
	$root = join('/', @s_dirs[0..$n-1]);
    }
    
    # Determine (by counting the number of backslashes)
    # the number of '../'s to prepend onto $source
    $root .= '/' unless ($root =~ /\/$/);
    my $n_back_dirs = ($target =~ tr{/}{/}) - ($root =~ tr{/}{/});
    
    # Re-route $source to be relative to $target
    $source =~ s/^$root//;
    $source = '../'x$n_back_dirs . $source;
    
    #print "root: $root\n";
    #print "$source -> $target\n";
    # TODO: if the link already exists and points correctly, don't overwrite it
    system "ln -sf -- $source $target";
}


# Prints a message to STDOUT with a timestamp
sub time_message ($ ) {
    my $time = localtime();
    print localtime() . ": $_[0]\n";
}

sub time_warning ($ ) {
    my $time = localtime();
    print STDERR "$time: $_[0]\n";
}

sub time_error ($ ) {
    my $time = localtime();
    print STDERR "$time: $_[0]\n";
    exit;
}






1;
