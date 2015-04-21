#! /usr/bin/perl -w

# This script is for convenient displaying of information stored in "summary.log" files in "make_log"
# directories of the AllPathsLG pipeline run. One summary.log file specified as argument displays times 
# in seconds and hours for each module ran. Two summary.log files result in module by module comparison 
# being displayed.

($fileIn1,$fileIn2) = @ARGV;

# $state must be {Starting,Finished}

sub hours {
    my ($secs) = @_;
    my $neg = $secs < 0 ? 1 : 0;
    $secs = -$secs if ( $secs < 0 );
    my $h = sprintf( "%1d", $secs / 3600 );
    my $m = sprintf( "%02d", ($secs - 3600 * $h) / 60 );
    my $s = sprintf( "%02d", $secs - 3600 * $h - 60 * $m );
    my $hours = $neg ? "-$h:$m:$s" : "$h:$m:$s";
    return $hours;
} 

open(FH,$fileIn1) || die "could not open fileIn1=$fileIn1";
while(<FH>){
    ($date, $module, $state ) = ($_ =~ /\[\S*\]\s(.*)\s+:\s+(\S+)\s+(\S+)/ );
    next if ($module =~ /^ln-{0,1}/);
    $time = `date +%s -d "$date"`;
    $hd1{$module}{$state} = $time;
    $hnames{$module} =1;
}
close FH;

foreach $module (sort {$a cmp $b } keys %hd1 ){
    $diff = $hd1{$module}{"Finished"} - $hd1{$module}{"Starting"};
    $hr1{$module}=$diff;
}

if ( ! defined $fileIn2 ){
    foreach $module (sort {$a cmp $b } keys %hr1 ){
	$time = $hr1{$module};
	$hours = hours($time);
	print "$module\t$time\t$hours\n";
    }

}else{
    open(FH,$fileIn2) || die "could not open fileIn2=$fileIn2";
    while(<FH>){
	($date, $module, $state ) = ($_ =~ /\[\S*\]\s(.*)\s+:\s+(\S+)\s+(\S+)/ );
	next if ($module =~ /^ln-{0,1}/);
	$time = `date +%s -d "$date"`;
	$hd2{$module}{$state} = $time;
	$hnames{$module}=1;
    }
    close FH;

    foreach $module (sort {$a cmp $b } keys %hd2 ){
	$diff = $hd2{$module}{"Finished"} - $hd2{$module}{"Starting"};
	$hr2{$module}=$diff;
    }
    
    

    print sprintf("%30s","module_name")."\t".sprintf("%7s","time2_s")."\t".sprintf("%10s","time2_h")."\t".sprintf("%7s","time2_s")."\t".sprintf( "%10s","time2_h" )."\t".sprintf( "%7s","diff_s" )."\t".sprintf( "%10s","diff_h" )."\n";
    $tot1 = $tot2 = 0;
    foreach $module (keys %hnames){
	if ( exists $hr1{$module} && exists $hr2{$module} ){
	    $diff = $hr2{$module} - $hr1{$module};
	    $tot1 += $hr1{$module};
	    $tot2 += $hr2{$module};
	    print sprintf("%30s",$module)."\t".sprintf("%7d",$hr1{$module})."\t".sprintf("%10s",hours($hr1{$module}))."\t".sprintf("%7d",$hr2{$module})."\t".sprintf("%10s",hours($hr2{$module}))."\t".sprintf("%7d",$diff)."\t".sprintf("%10s",hours($diff))."\n"; 
	}
    }
    print "------------------------------\n";
    print sprintf("%30s","total")."\t".sprintf("%7d",$tot1)."\t".sprintf("%10s",hours($tot1))."\t".sprintf("%7d",$tot2)."\t".sprintf( "%10s",hours($tot2) )."\t".sprintf( "%7d",$tot2-$tot1 )."\t".sprintf( "%10s",hours($tot2-$tot1) )."\n";
#    
}
