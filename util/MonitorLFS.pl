#!/util/bin/perl
use strict;
use Getopt::Long;

#
# usage: perl MonitorLFS [-l log_file] [-p poll_rate] [-r restarts] [-h] query start stop incr results_dir
#
# Monitor a program running on LFS and log errors
#
$Getopt::Long::passthrough = 1; # passthrough ignored args to @ARGV
$Getopt::Long::ignorecase = 0; # case sensitive args

# process command line arguments
my $options = {};
GetOptions($options, "--l=s", "--p=i", "--r=i", "--h");

my $log_file = undef;
my $poll_rate_secs = 60;
my $restarts = 2;
my $query = undef;
my $start = undef;
my $stop = undef;
my $incr = undef;
my $results_dir = undef;
if ( $options->{l} ){
  $log_file = $options->{l};
}
if ( $options->{p} ){
  $poll_rate_secs = $options->{p};
} 
if ( $options->{r} ){
  $restarts = $options->{r};
}

if ( scalar(@ARGV) != 5 ){
  print "@ARGV\n";
  Feedback();
  Usage();
} 

$query = $ARGV[0];
$start = $ARGV[1];
$stop = $ARGV[2];
$incr = $ARGV[3];
$results_dir = $ARGV[4];
if ( !defined($log_file) ) {
  $log_file= $query . "." . $start . "." . $stop . "." . $incr . ".log";
}
#Feedback();

open(LOGFILE, "> $log_file")
  or die "Couldn't open $log_file for writing: $!\n";

# don't buffer log file output
my  $old_fh = select(LOGFILE);
$| = 1;
select($old_fh);


my @job_names;
my @hits_files;
for ( my $i = $start; $i < $stop; $i += $incr ){
  my $job_stop = $i + $incr;
  if ( $job_stop > $stop ) { $job_stop = $stop; }
  my $job_name = $query. "." . $i . "." . $job_stop;
  my $hits_file = $results_dir . "/" . $query . ".hits." . $i . "." . $job_stop;
  push(@job_names, $job_name);
  push(@hits_files, $hits_file);
}

my %restarts_per_job;
my $job_name;
foreach $job_name ( @job_names ){
  $restarts_per_job{$job_name} = 0;
}

my $unfinished_jobs = 1; # unfinished jobs?
while($unfinished_jobs)
{
    my @finished_pos = ();
    $unfinished_jobs = 0; # reset for each iteration

    # get info about LSF jobs 
    my $pos = 0;
    foreach $job_name ( @job_names ){
        my @my_job_info = `bjobs -J $job_name`;

    	# see if job is finished 
	my $job_unfinished = 1;
	# Job <$job_name> is not found
    	if ( scalar(@my_job_info) == 1 ){
	  if ( $my_job_info[0] =~ /is not found/ ) {
	     $job_unfinished = 0;
	     print LOGFILE "$job_name $pos is not found\n";
	     }
    	} else {
    	  # remove title line
    	  shift(@my_job_info);

	  # process
	     # if job finished, remove from list
	  $job_unfinished = Process($my_job_info[0], $job_name);
        }
	if ( $job_unfinished ){
	  $unfinished_jobs++;
	} else {
	  print LOGFILE "$job_name is finished\n";
	  push(@finished_pos, $pos);
	}
	$pos++;
    }
    # sort numerically in descending order
    if ( scalar(@finished_pos) ){
       @finished_pos = sort { $b <=> $a } @finished_pos;
       }

    # remove finished jobs in descending order
    my $fin_pos;
    foreach $fin_pos (@finished_pos){
       my @return = splice(@job_names, $fin_pos, 1); 
       }

    sleep($poll_rate_secs);
}

# check to see if all jobs are done
if ( scalar(@job_names) == 0 ){
    print LOGFILE "ALL JOBS ARE DONE\n";
} else {
    print LOGFILE "ERROR: NOT ALL JOBS ARE DONE!!\n";
}

# now check output files for DONE! near the end
my $hits_file;
foreach $hits_file (@hits_files)
  {
    my $job_done = 0;
    open(HITSFILE, $hits_file)
      or die "Couldn't open $hits_file for reading: $!\n";
    my @lines = reverse <HITSFILE>;
    # only check last 5 lines at end of file
    for ( my $num = 0; $num < 5; $num++ ){
      if ( $lines[$num] =~ /DONE!/ ){
	 $job_done = 1;
      }
    }
    close(HITSFILE);
    if ( !$job_done ){
      print LOGFILE "$hits_file is not DONE\n";
    } else {
      #print "$hits_file is DONE\n";
    }
  }

close(LOGFILE);

##
## ----------------------------------------------------------------------------
##

sub Usage
  {
  print "Usage: perl MonitorLFS [-h] [-l log_file] [-p poll_rate_secs] [-r restarts]\n";
  print "                       query start stop incr results_dir\n";
  print "\t-h = help\n";
  print "\t-l = log file for errors, default is job_names.log\n";
  print "\t-p = check bjobs every poll_rate seconds, default is 60\n";
  print "\t-r = number of times to restart a suspended job, default is 2\n";
  print "\tquery = genome to query\n";
  print "\tstart = starting postion on genome\n";
  print "\tstopping = stopping position on genome\n";
  print "\tincr = processing increment on genome\n";
  print "\tresults_dir = results directory for QueryLookupTable output\n";
  exit(0);
  }

sub Feedback
  {
  print "\tlog_file = $log_file poll_rate_secs = $poll_rate_secs restarts = $restarts query = $query start = $start ";
  print "stop = $stop incr = $incr results_dir = $results_dir\n";
  }

sub Process
  {     
  my ($my_job_info, $job_name ) = @_;
  chomp($my_job_info);

  my $unfinished = 1;

  if ( $my_job_info =~ /\s+DONE\s+/ ){ 
    $unfinished=0;
    print LOGFILE "DONE:$job_name\n";
  } elsif ( $my_job_info =~ /\s+PEND\s+/ ){
    my @pend_reason = `bjobs -p -J $job_name`;
    print LOGFILE "PEND:$job_name\n$my_job_info\n@pend_reason\n";  
  } elsif ( $my_job_info =~ /\s+PSUSP\s+/ ){ 
    my @pend_reason = `bjobs -p -J $job_name`;
    if ( $restarts_per_job{$job_name} < $restarts ){
       `bresume -J $job_name`;
       $restarts_per_job{$job_name}++;
    }
    print LOGFILE "PSUSP:$job_name\n$my_job_info\n@pend_reason\n";  
  } elsif ( $my_job_info =~ /\s+USUSP\s+/ ){ 
    my @susp_reason = `bjobs -s -J $job_name`;
    if ( $restarts_per_job{$job_name} < $restarts ){
       `bresume -J $job_name`;
       $restarts_per_job{$job_name}++;
    }
    print LOGFILE "USUSP:$job_name\n$my_job_info\n@susp_reason\n";  
  } elsif ( $my_job_info =~ /\s+SSUSP\s+/ ){
    my @susp_reason = `bjobs -s -J $job_name`;
    print LOGFILE "SSUSP:$job_name\n$my_job_info\n@susp_reason\n";  
  } elsif ( $my_job_info =~ /\s+EXIT\s+/ ){ 
    print LOGFILE  "EXIT:$job_name\n$my_job_info\n"; 
  } elsif ( $my_job_info =~ /\s+ZOMBI\s+/ ){ 
    print LOGFILE "ZOMBI:$job_name\n$my_job_info\n";
  } elsif ( $my_job_info =~ /\s+UNKWN\s+/ ){ 
    print LOGFILE "UNKWN:$job_name\n$my_job_info\n"; 
  } elsif ( $my_job_info =~ /\s+RUN\s+/ ){ 
  } elsif ( $my_job_info =~ /\s+WAIT\s+/ ){ 
    print LOGFILE "WAIT:$job_name\n$my_job_info\n"; 
  } else { 
    $unfinished = 0;
    print LOGFILE "FINISHED:$job_name\n";
  }

  return($unfinished);
  }
