#!/usr/bin/perl
#///////////////////////////////////////////////////////////////////////////////
#//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
#//       This software and its documentation are copyright (2009) by the     //
#//   Broad Institute.  All rights are reserved.  This software is supplied   //
#//   without any warranty or guaranteed support whatsoever. The Broad        //
#//   Institute is not responsible for its use, misuse, or functionality.     //
#///////////////////////////////////////////////////////////////////////////////
#
# Parses .../makeinfo/<filename>.mm.<modulename> for memmonitor data:
#
# (MM,<pid>) fields: time(sec), etime, detime, utime, dutime, vmsize(KB), vmrss, vmdata, vmstk, vmexe, vmlib
#
# outputs report to '<SUBDIR>/mm.report' 
#

use strict;

die "usage: $0 SUBDIR=<full SUBDIR dir>\n" unless (@ARGV);

my $subdir = ("@ARGV" =~ /SUBDIR=\s*(\S+)\s*$/)[0];
$subdir =~ s/\/$//;

die "**** invalid SUBDIR directory '$subdir'.\n" unless (-d $subdir);

my $out_fn = "$subdir/mm.report";

# open extended report file

open(FILE, ">$out_fn") or die "**** Can't open '$out_fn'.\n";
print FILE "# $0 @ARGV\n#\n";


# ---- find all the '*.mm.*' files
my @mm_fns = (find_lower_makeinfo_files($subdir, "mm"), 
              find_upper_makeinfo_files($subdir, "mm"));


# ---- parse all the '*.mm.*' files

my @mods = ();

foreach my $mm_fn (@mm_fns) 
{
    my ($mod_t, $mod_n) = $mm_fn =~ /([^\/]+)\.mm\.(\S+)$/;

    my $t1 = (stat($mm_fn))[9];   # the file's modification time

    my @line = ();
    my $iter = 0;
    while (scalar @line == 0 && $iter < 10) {
        @line = split ' ', `grep ' Summary: ' $mm_fn`;  # get memmonitor data
        $iter++;
        sleep(30) if (scalar @line == 0);  # wait for MemMonitor to finish 
    }

    if (scalar @line > 0) {
	push @mods, {time0   => $line[3] - $line[4],
		     time1   => $line[3],
		     etime   => $line[4], 
		     detime  => $line[5], 
		     utime   => $line[6], 
		     dutime  => $line[7], 
		     vmsize  => $line[8] / 1024, # convert from KB to MB
		     vmrss   => $line[9] / 1024, # convert from KB to MB
		     vmdata  => $line[10],
		     vmstk   => $line[11],
		     vmexe   => $line[12],
		     vmlib   => $line[13],
		     name    => $mod_n,
		     target  => $mod_t};
    }
}



my $max_vmrss = 0;
my $max_vmsize = 0;

if (1) {
    foreach my $mod (@mods) {
        $max_vmrss = max ($mod->{vmrss}, $max_vmrss);
        $max_vmsize = max ($mod->{vmsize}, $max_vmsize);
    }     
}
else {
    # ---- find worst-case max memory usage 
    
    # create an array with all the start(0) and end(1) times pegged with 
    # the effects on memory usage (+ on start, - on end)

    my @ts = ((map { t => $_->{time0}, dvmrss =>   $_->{vmrss}, dvmsize =>   $_->{vmsize} }, @mods),
              (map { t => $_->{time1}, dvmrss => - $_->{vmrss}, dvmsize => - $_->{vmsize} }, @mods));

    my $vmrss = 0;
    my $vmsize = 0;
    
    my @ts_sorted = sort { $a->{t} <=> $b->{t} } @ts;
    
    foreach my $t (@ts_sorted) 
    {
        $vmrss += $t->{dvmrss};
        $vmsize += $t->{dvmsize};
        
        $max_vmrss = max($max_vmrss, $vmrss);
        $max_vmsize = max($max_vmsize, $vmsize);
    }
}

# ---- prepare for output

my @mods_sorted = sort { $a->{time0} <=> $b->{time0} } @mods;

my $etime = 0;
my $utime = 0;

# ---- condensed report to standard out

my $sep_c =   "--------------------------------------------------------------------------------\n";

my $out_c = " e_time(sec)  u_time(sec)  vmrss(MB)  vmsize(MB)   module name\n";
$out_c .= $sep_c;

foreach (@mods_sorted) 
{ 
    $out_c .= sprintf("%12d %12d %10d %11d   %-29s\n", 
                      $_->{etime}, 
                      $_->{utime}, 
                      $_->{vmrss}, $_->{vmsize}, 
                      substr($_->{name},0,29),);
    $etime += $_->{etime};
    $utime += $_->{utime};
}
$out_c .= $sep_c;
$out_c .= sprintf("%12d %12d %10d %11d   Total/Peak %d modules\n", 
		  $etime, $utime, $max_vmrss, $max_vmsize, scalar @mods_sorted);

print $out_c;


# ---- extended report to file

my $sep_e = "# ------------------------------------------------------------------------------------------------------------------------------------------------\n";

my $out_e = "#               start  e_time(sec)  u_time(sec)   vmrss(MB)  vmsize(MB)    module name                    module target\n";
$out_e .= $sep_e;

foreach (@mods_sorted) 
{ 
    $out_e .= sprintf("  %19s %12d %12d %11d %11d    %-30s %-40s\n", 
                      t_fmt($_->{time0}), 
                      $_->{etime}, 
                      $_->{utime}, 
                      $_->{vmrss}, $_->{vmsize}, 
                      $_->{name}, $_->{target});
}
$out_e .= $sep_e;
$out_e .= sprintf("# Total/Peak:         %12d %12d %11d %11d    %d modules\n", 
		  $etime, $utime, $max_vmrss, $max_vmsize, scalar @mods_sorted);


print FILE $out_e;
close(FILE);






















sub find_lower_makeinfo_files 
{
    my ($curdir, $label) = @_;
    $curdir =~ s/\/?$//;

    my @fns = (-d "$curdir/makeinfo") ? `find $curdir/makeinfo -name '*.$label.*'` : ();
    #print "lower $curdir\n";

    chomp @fns;

    foreach my $dir (`ls $curdir`) 
    {
        chomp $dir;
        if (-d "$curdir/$dir" and 
            $dir ne "seed" and 
            $dir ne "makeinfo") {
            push @fns, find_lower_makeinfo_files("$curdir/$dir", $label)
        }
    }
    return @fns;
}



sub find_upper_makeinfo_files 
{
    my ($curdir, $label) = @_;

    return () if (!$curdir or -d "$curdir/make_log");

    $curdir =~ s/\/[^\/]+\/?$//;   # get rid of top level directory

    my @fns = (-d "$curdir/makeinfo") ? `find $curdir/makeinfo -name '*.$label.*'` : ();

    #print "higher $curdir\n";
    chomp @fns;
    
    return (@fns, find_upper_makeinfo_files($curdir, $label));
}




sub t_fmt
{
    my ($t) = @_;
    my @aux = localtime($t);
    $aux[5] += 1900; # dates start in 1900
    $aux[4]++; # month would otherwise be in [0..11]
    return sprintf "%4d.%02d.%02d-%02d:%02d:%02d", @aux[5,4,3,2,1,0];
}


sub max
{
    my ($a, $b) = @_;
    return ($a > $b) ? $a : $b;
}
