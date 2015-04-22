#!/usr/bin/perl -w
##########################################################################
#                  SOFTWARE COPYRIGHT NOTICE AGREEMENT                   #
#      This software and its documentation are copyright (2013) by the   #
#  Broad Institute.  All rights are reserved.  This software is supplied #
#  without any warranty or guaranteed support whatsoever. The Broad      #
#  Institute is not responsible for its use, misuse, or functionality.   #
##########################################################################
# ---------------------------------------------------------------------- #
#                                                                        #
#  Utility to merge a log file with memory usage statistics obtained     #
#  with MemMonitor                                                       #
#                                                                        #
# ---------------------------------------------------------------------- #

use strict;
use FindBin;

# Local libraries, from elsewhere in the BroadCRD repository.
use lib "$FindBin::Bin/../";
use lib "$FindBin::Bin/";

use PerlRunTime;    # run_or_die
use PerlUtils;      # ISO_date, abort
use ArachneArgs;    # Command-line ARG=VALUE parser
use Time::Local;


my %args = getCommandArguments
    (MM         => { value => undef,
                     help  => "MemMonitor file to merge."},
     LOG        => { value => undef,
                     help  => "Log file to merge."},
     LOG_ALL    => { value => "False",
		     help  => "If 'True' include log lines without timestamps." },
     MM_ALL     => { value => "False",
		     help  => "If 'True' then show all MM lines, not just the peak." },
     TABLE      => { value => "False",
		     help  => "If 'True' then generate in a table instead." },
     TIME_SEC   => { value => "False",
		     help  => "If 'True' then give the timestamp in elapsed seconds" },
     DUR_SEC    => { value => "True",
		     help  => "If 'False' then give the duration in hours, mins, secs" },
 );

# Parse Args

my $fn_mm_arg = $args{MM};

my $fn_log_arg = $args{LOG};

my $log_all_arg = $args{LOG_ALL};

my $mm_all_arg = $args{MM_ALL};

my $table_arg = $args{TABLE};

my $time_sec_arg = $args{TIME_SEC};

my $dur_sec_arg = $args{DUR_SEC};

if ($table_arg) {
  $mm_all_arg = 0;
  $log_all_arg = 0;
}

# Merge MM and LOG

mm_merge_logs ($fn_mm_arg, $fn_log_arg, $log_all_arg, $mm_all_arg, $table_arg, $time_sec_arg, $dur_sec_arg);

# Finished!


sub print_mem
{
  my ($time, $message, $memory) = @_;

  my $memory_gb = sprintf('%1.1f', $memory/(1024**2));
  print gmtime($time). ": " . $message . ": " . $memory_gb . " GB\n";
}

sub get_hms
{
  my ($duration) = @_;
  my $duration_h = int(($duration)/(60**2));
  my $duration_m = int(($duration - $duration_h*(60**2))/60);
  my $duration_s = int($duration - $duration_h*(60**2) - $duration_m*60);
  return sprintf('%3dh%02dm%02ds', $duration_h, $duration_m, $duration_s);
}

sub mm_merge_logs
{
  my ($fn_mm, $fn_log, $log_all, $mm_all, $table, $time_sec, $dur_sec) = @_;

  # Month name to month of year
  my %mon2num = qw( jan 0  feb 1  mar 2  apr 3  may 4  jun 5
		    jul 6  aug 7  sep 8  oct 9 nov 10 dec 11 );

  # local memory usage peak and the corresponding time
  my $peak_mem_time = 0;
  my $peak_mem = 0;

  # global memory usage peak and the corresponding time
  my $global_peak_mem_time = 0;
  my $global_peak_mem = 0;

  # global start and end time
  my $start_time = 0;
  my $end_time = 0;

  # log entry times
  my $log_time = 0;
  my $previous_log_time = 0;

  # log entry messages
  my $log_message = "";
  my $previous_log_message = "";

  # last MM log entry examined
  my $previous_time = 0;
  my $previous_mem = 0;

  # create table header
  if ($table) {
    my $time_str = ($time_sec ?  sprintf("%7s", "Time") : sprintf("%24s", "Time"));
    my $dur_str = ($dur_sec ? sprintf("%6s", "Dur.") : sprintf("%10s", "Duration"));
    printf("%s,  %s,  %5s,  %s\n", $time_str, $dur_str, "Mem.", "Log");
  }

  open(FILE_MM, "<$fn_mm") or die "\n**** can't open '$fn_mm' for reading.\n\n";
  open(FILE_LOG, "<$fn_log") or die "\n**** can't open '$fn_log' for reading.\n\n";

  # Start merging the logs together, starting with the LOG file

  while (my $log_line = <FILE_LOG>) {

    # Weekday Month Day HH:MM:SS Year : message
    if ($log_line =~ /^(\w{3})\s(\w{3})\s(\d+)\s(\d+):(\d+):(\d+)\s(\d{4})\s*:*\s*(.*)/) {

      # Save previous log entry details
      $previous_log_time = $log_time;
      $previous_log_message = $log_message;

      # parse log line
      my $day = $1; my $month = $mon2num{ lc $2 }; my $mday = $3;
      my $hour = $4; my $minute = $5; my $second = $6; my $year = $7;
      $log_message = $8;

      # compute unix time for LOG file entry
      $log_time = timelocal($second, $minute, $hour, $mday, $month, $year);

      # Look forward in the MM file only if need to
      if ($log_time > $previous_time) {

	# initialize local peak memory stats using cached values from previous search
	$peak_mem_time = $previous_time;
	$peak_mem = $previous_mem;

	# print cached memory usage if requested (MM_ALL=True) - special case
	print_mem($previous_time, "Memory Usage", $previous_mem) if ($mm_all && $previous_time != 0);

	# Examine MM log until we pass the current LOG file entry time
	while (my $mm_line = <FILE_MM>) {

	  # (MM,<PID>) time(sec) etime detime utime dutime vmsize(KB) vmrss vmdata vmstk vmexe vmlib
	  if ($mm_line =~ /^.MM\S+\s+(\S+)\s*\S+\s*(\S+)\s*\S+\s*(\S+)\s*\S+\s*(\S+)\s*\S+\s*\S+\s*\S+/) {
	    my $mm_time = $1;
	    my $mm_mem = $4;

	    # Update global statistics
	    $start_time = $mm_time if ($start_time == 0);
	    $end_time = $mm_time;
	    $global_peak_mem_time = $mm_time if ($mm_mem > $global_peak_mem);
	    $global_peak_mem = $mm_mem if ($mm_mem > $global_peak_mem);

	    # Update local statistics, or else cache value to start next search
	    if ($mm_time <= $log_time) {
	      if ($mm_mem > $peak_mem) {
		# New local peak found, record it
		$peak_mem_time = $mm_time;
		$peak_mem = $mm_mem;
	      }
	      # print memory usage if requested (MM_ALL=True)
	      print_mem($mm_time, "Memory Usage", $mm_mem) if ($mm_all);

	    } else {
	      # Looked too far ahead, caching values for next search
	      $previous_time = $mm_time;
	      $previous_mem = $4;
	    }
	    # Bail if we have looked too far ahead
	    last if ( $log_time < $mm_time)
	  }
	}
	# Print peak memory usage
	print_mem($peak_mem_time, "Peak Memory", $peak_mem) unless ($table);
      }
      # Print table row
      if ($table) {
	$previous_log_time = $log_time if ($previous_log_time == 0);
	my $duration = ($log_time - $previous_log_time);
	my $time_str = ($time_sec ?  sprintf("%7d", $duration) : sprintf("%24s", scalar(gmtime($log_time))));
	my $dur_str = ($dur_sec ? sprintf("%6d", $duration) : sprintf("%10s", get_hms($duration)));
	printf("%s,  %s,  % 5.1f,  %s \n", $time_str, $dur_str, $peak_mem/(1024**2), $log_message);
      } else {
      # Print LOG message
	print gmtime($log_time). ": " . $log_message . "\n";
      }
    } elsif ($log_all) {
      # Print LOG messages with no timestamps if requested (LOG_ALL=True)
      print $log_line;
    }
  }
  close(FILE_LOG);
  close(FILE_MM);

  # Write the summary
  
  my $runtime_hr = sprintf('%1.1f', ($end_time - $start_time)/(60**2));
  print "Total Runtime: " . $runtime_hr . " hours\n";
  my $peak_mem_gb = sprintf('%1.1f', $global_peak_mem/(1024**2));
  print "Peak memory usage of " . $peak_mem_gb . " GB occured on " . gmtime($global_peak_mem_time) . "\n";
  return;
}
