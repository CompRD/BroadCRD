#!/usr/bin/perl -w
##########################################################################
#                  SOFTWARE COPYRIGHT NOTICE AGREEMENT                   #
#      This software and its documentation are copyright (2012) by the   #
#  Broad Institute.  All rights are reserved.  This software is supplied #
#  without any warranty or guaranteed support whatsoever. The Broad      #
#  Institute is not responsible for its use, misuse, or functionality.   #
##########################################################################
# ---------------------------------------------------------------------- #
#                                                                        #
#  Utility to plot the memory and cpu usage statistics                   #
#  obtained with MemMonitor                                              #
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


my %args = getCommandArguments
    (MM         => { value => undef,
                     help  => "MemMonitor file to plot."},
     WINDOW     => { value => "0:100",
                     help  => "Relative window to display."},
     GIF        => { value => 0,
                     help  => "Whether to convert .eps to .gif (better for powerpoint)."} );

my $fn_mm = $args{MM};


my $stats = mm_file_stats($fn_mm);

# ---- adjust time limits 
{
    my ($p0, $p1) = $args{WINDOW} =~ /^(.+):(.+)$/;
    my $dt = $stats->{t1} - $stats->{t0};
    $stats->{t0} += 0.01 * $p0 * $dt;
    $stats->{t1} += (0.01 * $p1 - 1) * $dt;

    $dt = $stats->{t1} - $stats->{t0};
    $stats->{t0} -= 0.1 * $dt;
    $stats->{t1} += 0.1 * $dt;

    die "**** no data points to plot!\n"
        unless ($stats->{t1} > $stats->{t0});
}

# ---- get time zone information 

my $t_file = file_time($fn_mm); 
my $dtz = (split " ", `date +"\%z" -d"\@$t_file"`)[0] / 100 * 3600; 

my $t0 = $stats->{t0} + $dtz;
my $t1 = $stats->{t1} + $dtz;


# ---- figuring out point and line type. 
#      cumbersome so that it might be extended to multiple files

my @pt = (5, 7, 9, 11, 13); # point types
my $lts = line_types($stats->{cmd});

my @titles = ($stats->{cmd});

my @ltypes = ();
my @ptypes = ();
{
    my $pt = $pt[0];
    my $lt = $lts->{$stats->{cmd}};
    $ltypes[0] = "lt rgb '#$lt' lw 3*scale";
    $ptypes[0] = "pt $pt ps 1.8*scale";
}

# ----------------------------------------
#  create gnuplot script and execute it
# ----------------------------------------

my $gp_fn = "./$fn_mm.gp"; # the gnuplot script
my $ps_fn = "$fn_mm.ps"; # the postscript

print "     generating gnuplot script\n";
my @gp = <DATA>; # get the gnuplot template from the __DATA__ section at the end 


foreach my $line (@gp)
{
    $line =~ s/__TITLE__/$fn_mm/;
    $line =~ s/__MAXMEMY__/$stats->{memGB}/;
    $line =~ s/__NCPUS__/$stats->{cpu_max}/;
    $line =~ s/__PS_FN__/$ps_fn/;
    $line =~ s/__T_MIN__/$t0/;
    $line =~ s/__T_MAX__/$t1/;
    
    if ($line =~ /__CPU_PLOT_LINES__/)
    {
        $line = "";
        $line .= "'$fn_mm' u (timecolumn(2) + $dtz):(\$6/\$4) not w lp $ltypes[0] $ptypes[0], \\\n"; 
        $line .= "'$fn_mm' u (timecolumn(2) + $dtz):(\$5/\$3) not w l $ltypes[0]\n"; 
    }
    elsif ($line =~ /__MEM_PLOT_LINES__/)
    {
        $line = "";
        $line .= "'$fn_mm' u (timecolumn(2) + $dtz):(fmem(\$7)) not w l $ltypes[0], \\\n"; 
        $line .= "'$fn_mm' u (timecolumn(2) + $dtz):(fmem(\$8)) t '$stats->{cmd}' w lp $ltypes[0] $ptypes[0]\n"; 
    }
}
file_write_or_die($gp_fn, @gp);

print "     generating plot\n";
system "chmod --quiet +x $gp_fn ; $gp_fn";

print "     showing plot\n";
exec "gv ./$ps_fn";



sub cryptchr_to_uint6
{
    my $n = ord($_[0]);
    
    return ($n - ord("."))      if ($n == ord(".") or $n == ord("/"));
    return ($n - ord("0") + 2)  if ($n >= ord("0") and $n <= ord("9"));
    return ($n - ord("A") + 12) if ($n >= ord("A") and $n <= ord("Z"));
    return ($n - ord("a") + 38) if ($n >= ord("a") and $n <= ord("z"));
    die;
}

sub rgb_from_string 
{
    my ($str) = @_;
    my $salt = "11";
    my $sc = crypt($str, $salt);
    
    my $r = cryptchr_to_uint6(substr($sc, -1, 1));
    my $g = cryptchr_to_uint6(substr($sc, -2, 1));
    my $b = cryptchr_to_uint6(substr($sc, -3, 1));

    # $r, $g, $b are 64 bit numbers
    
    my @hex = (0..9, "A".."F");
    
    my $rs = $hex[$r >> 2] . $hex[($r % 4) << 2]; 
    my $gs = $hex[$g >> 2] . $hex[($g % 4) << 2]; 
    my $bs = $hex[$b >> 2] . $hex[($b % 4) << 2]; 
    
    return "$rs$gs$bs";
}



sub line_types 
{
    my %lt = ();
    foreach my $l (@_) 
    {
        if (!exists $lt{$l}) {
            $lt{$l} = rgb_from_string($l);
        }
    }
    return \%lt;
}



# -------------------------------
#  subroutines i/o
# -------------------------------

sub mm_file_stats
{
    my ($fn) = @_;

    my $cmd = 0;
    my $t0 = 0;
    my $t1 = 0;
    my $mem_kB_max = 0;
    my $cpu_max = 0;
        
    open(FILE, "<$fn") or die "\n**** can't open '$fn' for reading.\n\n";
    while (<FILE>) {
        # (MM,<PID>) time(sec) etime detime utime dutime vmsize(KB) vmrss vmdata vmstk vmexe vmlib
        
        if (/^.MM\S+\s+(\S+)\s*\S+\s*(\S+)\s*\S+\s*(\S+)\s*\S+\s*(\S+)\s*\S+\s*\S+\s*\S+/) {
            $t0 = $1 if ($t0 == 0);
            $t1 = $1;
            $mem_kB_max = $4 if ($4 > $mem_kB_max);
            my $cpu = $3/$2;
            $cpu_max = $cpu if ($cpu > $cpu_max);
        }
        elsif (/cmdline:\s.(\S+)\s/) {
            $cmd = $1;
        }
    }
    close(FILE);

    die "\n**** invalid MemMonitor file '$fn'.\n\n" unless $cmd;

    return { cmd => $cmd, 
             t0 => $t0, 
             t1 => $t1, 
             memGB => 1.0 * $mem_kB_max / 1024**2, 
             cpu_max => 1.0 * $cpu_max};
}


sub file_write_or_die
{
    my ($fn, @l) = @_;
    open(FILE, ">$fn") or die "\n**** can't open '$fn' for writing.\n\n";
    print FILE @l;
    close(FILE);
}


sub file_time
{
    my ($fn) = @_;
    return (-e $fn) ? (stat($fn))[9] : 0;
}


__DATA__
#!/usr/bin/env gnuplot

scale = 0.43

pt_pl  =  1     # +
pt_cr  =  2     # x
pt_st  =  3     # *
pt_eu4 =  4     # empty square
pt_fu4 =  5     # full  
pt_ec  =  6     # empty circle 
pt_fc  =  7     # full
pt_eu3 =  8     # empty up triangle
pt_fu3 =  9     # full
pt_ed3 = 10     # empty down triangle
pt_fd3 = 11     # full
pt_ed4 = 12     # empty diamond
pt_fd4 = 13     # full
pt_e5 =  14     # empty pentagon
pt_f5 =  15     # full



set encoding iso_8859_1

#set terminal postscript {landscape | portrait | eps | default}
#                        {enhanced | noenhanced}
#                        {color | monochrome} {solid | dashed}
#                        {<duplexing>}
#                        {"<fontname>"} {<fontsize>}
#
#set terminal postscript eps enhanced monochrome dashed 'Helvetica' 22*scale
#black = 1; grey = 0
set terminal postscript landscape noenhanced color solid 'Helvetica' 22 * scale
black = 7; grey = 9



set border 31 lt black lw 2*scale
#set zeroaxis lt 7 lw 1
set grid  lt grey lw scale

size_h    = 1.40 * scale 
margin_l  = 0.23 * scale
margin_ih = 0.10 * scale
margin_r  = 0.75 * scale

size_v    = 0.90 * scale
margin_t  = 0.20 * scale
margin_iv = 0.17 * scale
margin_b  = 0.10 * scale


nh = 1
nv = 2

set size margin_l + nh*size_h +(nh-1)*margin_ih + margin_r, \
         margin_b + nv*size_v +(nv-1)*margin_iv + margin_t 

set lmargin  0; set rmargin  0
set tmargin  0; set bmargin  0

#set nokey


set tics scale 0.0 * scale

set output '__PS_FN__'

set multiplot

  # -------------------
  # cpu usage pcpus
  # -------------------

  ih = 0; iv = 0;

  #ymin = 0.0625;  ymax = 64 
  ymin = 0; ymax = ceil(__NCPUS__)

  set timefmt "%s"
  set xdata time
  unset log x
  unset log y
  
  set xlabel 'time' font 'Helvetica'
  set format x ''
  set xtics
  
  set ylabel 'fcpu' font 'Helvetica'
  set format y '%g'
  set ytics 4
  #set ytics ("1/16" 1./16.,\
  #           "1/8"  1./8.,\
  #           "1/4"  1./4.,\
  #           "1/2"   1./2.,\
  #             "1"   1,\
  #             "2"   2,\
  #             "4"   4,\
  #             "8"   8,\
  #            "16"   16,\
  #            "32"   32,\
  #            "64"   64)
  set mytics 1
  set my2tics 1
  
  fpcpu(x) = x/100.0
  
  set grid xtics ytics lt grey
  #set key left center at graph 1.03,1.10 vertical Left reverse height -40
  set key off
  set size size_h, size_v
  set origin margin_l + ih*(size_h + margin_ih), margin_b+ iv*(size_v + margin_iv)
   
  # time limits need to be in quotes to be interpreted as 1970 unix time
  # otherwise they are interpreted as seconds from 2000 (internal gnuplot representation)
  plot ['__T_MIN__':'__T_MAX__'] [ymin:ymax] \
       __CPU_PLOT_LINES__



  # -------------------
  # memory usage (vmrss)
  # -------------------

  set title '__TITLE__'

  ih = 0; iv = 1;

  ymin = 0; ymax = 2**(int(log( __MAXMEMY__) / log(2)) + 1 )

  set timefmt "%s"
  set xdata time
  unset log x
  unset log y
  unset log y2

  set xlabel ''  font 'Helvetica'
  set format x "%Y\n%b %d\n%02H:%02M" 
  set xtics
  set x2label '' font 'Helvetica'
  set format x2 '' 
  set x2tics

  set ylabel 'rss - vsize' font 'Helvetica'
  set format y '%g GB'
  #set ytics ymin, __MEMYTICS__, ymax
  set ytics ymin, ymax/8.0, ymax
  #set mytics 1

    # set ytics ("16 MB"     1./64.,\
    #         "64 MB"     1./16.,\
    #         "256 MB"    1./4.,\
    #         "1 GB"      1,\
    #         "4 GB"      4,\
    #         "16 GB"    16,\
    #         "64 GB"    64,\
    #         "256 GB"  256,\
    #         "1 TB"   1024,\
    #         "4 TB"   4096)
     

  set y2label '' font 'Helvetica'
  set format y2 ''
  unset y2tics 
  set my2tics 1

  fmem(x) = x/1024**2

  set grid xtics ytics lt grey
  set key left center at graph 1.03,-0.10 vertical Left reverse height -40  
  set size size_h, size_v
  set origin margin_l + ih*(size_h + margin_ih), margin_b + iv*(size_v + margin_iv)

  # time limits need to be in quotes to be interpreted as 1970 unix time
  # otherwise they are interpreted as seconds from 2000 (internal gnuplot representation)
  plot ['__T_MIN__':'__T_MAX__'] [0:ymax] \
       __MEM_PLOT_LINES__




unset multiplot



#guide for line and point styles:
#  0  ..............  .                    broken line
#  1  --------------  +                    red
#  2  -- -- -- -- --  x                    green
#  3  -  -  -  -  -   *                    blue
#  4  ..............  empty square         magenta
#  5  __.__.__.__.__  full  square         cyan
#  6  _ . _ . _ . _   empty circle         yellow
#  7  - -  - -  - -   full  circle         black
#  8  - - -  - - -    empty up triangle    brown
#  9  - - - -  - - -  full  up triangle    grey
# 10 (1)              empty down triangle
# 11 (2)              full  down triangle
# 12 (3)              empty diamond
# 13 (4)              full  diamond
# 14 (5)              empty pentagon
# 15 (6)              full  pentagon
# 16-31               watches



# ------------
#    done
# ------------



