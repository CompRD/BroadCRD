#!/usr/bin/perl -w
###############################################################################
#                    SOFTWARE COPYRIGHT NOTICE AGREEMENT                      #
#        This software and its documentation are copyright (2011) by the      #
#    Broad Institute.  All rights are reserved.  This software is supplied    #
#    without any warranty or guaranteed support whatsoever. The Broad         #
#    Institute is not responsible for its use, misuse, or functionality.      #
###############################################################################
use strict;
use FindBin;
use POSIX ":sys_wait_h";     # for the WNOHANG signal to wait pid

# Local libraries, from elsewhere in the BroadCRD repository.
use lib "$FindBin::Bin/../";
use lib "$FindBin::Bin/";

use PerlRunTime;    # run_or_die
use PerlUtils;      # ISO_date, abort
use ArachneArgs;    # Command-line ARG=VALUE parser


my %args = getCommandArguments
    (K        => { value => undef,
                   help  => "Kmer size."},
     SPECTRA  => { value => undef,
                   help  => "List of the files to plot."},
     OUT_HEAD  => { value => "",
                   help  => "Prefix for output files."},
     MIN_FREQ => { value => 1,
                   help  => "Minimum kmer frequency."},
     MAX_FREQ => { value => "",
                   help  => "Maximum kmer frequency."},
     CUMULATIVE => { value => 1,
                   help  => "Include cumulative lines in plots."},
     GIF      => { value => 0,
                   help  => "Whether to convert .eps to .gif (better for powerpoint)."} );

my $K = $args{K};
my $outhead = "";
$outhead = $args{OUT_HEAD} . "." if $args{OUT_HEAD};

my @fns = array_from_ref_or_value($args{SPECTRA});

foreach my $fn (@fns)
{
    die "can't find '$fn'.\n" if (! -e $fn || -d $fn);
}

my ($x_max, $y_max, $x2_max, $y2_max) = limits(@fns);


my $gp_fn = "$K-mer_spectra.gp";

my @lts = (1,3,4,2,5,8,6,9);

# columns in kmer spectrum files
# 1: kmer_frequency
# 2: num_distinct_kmers
# 3: frac_distinct_kmers = 2/total_num_distinct_kmers
# 4: cummulative(2)
# 5: cummulative(3) = 4/total_num_distinct_kmers
# 6: num_kmers = 1*2
# 7: frac_kmers = 6/total_num_kmers
# 8: cummulative(6)
# 9: cummulative(7) = 8/total_num_kmers


my $pl1 = "";  # num_distinct_kmers
for (my $i = 0; $i != @fns; $i++) 
{
    $pl1 .= ", " if ($pl1 ne "");
    $pl1 .= "'$fns[$i]' u 1:2 not w l lt $lts[$i] lw 4*scale";
    if ($args{CUMULATIVE}) {
	$pl1 .=", ";
	$pl1 .= "'$fns[$i]' u 1:4 not w l lt $lts[$i] lw 2*scale";
    }
}

my $pl2 = ""; # frac_distinct_kmers
for (my $i = 0; $i != @fns; $i++) 
{
    $pl2 .= ", " if ($pl2 ne "");
    $pl2 .= "'$fns[$i]' u 1:3 not w l lt $lts[$i] lw 4*scale";
    if ($args{CUMULATIVE}) {
	$pl2 .=", ";
	$pl2 .= "'$fns[$i]' u 1:5 not w l lt $lts[$i] lw 2*scale";
	$pl2 .=", ";
	$pl2 .= "'$fns[$i]' u 1:(1-\$5) not w l lt $lts[$i] lw 2*scale";
    }
}


my $pl3 = "";  # num_kmers
for (my $i = 0; $i != @fns; $i++) 
{
    $pl3 .= ", " if ($pl3 ne "");
    $pl3 .= "'$fns[$i]' u 1:6 not w l lt $lts[$i] lw 4*scale";
    if ($args{CUMULATIVE}) {
	$pl3 .=", ";
	$pl3 .= "'$fns[$i]' u 1:8 not w l lt $lts[$i] lw 2*scale";
    }
}

my $pl4 = ""; # frac_kmers
for (my $i = 0; $i != @fns; $i++) 
{
    $pl4 .= ", " if ($pl4 ne "");
    $pl4 .= "'$fns[$i]' u 1:7 not w l lt $lts[$i] lw 4*scale";
    if ($args{CUMULATIVE}) {
	$pl4 .=", ";
	$pl4 .= "'$fns[$i]' u 1:9 not w l lt $lts[$i] lw 2*scale";
	$pl4 .=", ";
	$pl4 .= "'$fns[$i]' u 1:(1-\$9) not w l lt $lts[$i] lw 2*scale";
    }
}




my @gp = <DATA>; # get the gnuplot template from the __DATA__ section at the end 


foreach my $line (@gp)
{
    $line =~ s/__HEAD__/$outhead/g;
    $line =~ s/__K__/$K/g;
    $line =~ s/__MIN_FREQ__/$args{MIN_FREQ}/g;
    $line =~ s/__MAX_FREQ__/$args{MAX_FREQ}/g;
    $line =~ s/__X_MAX__/$x_max/g;
    $line =~ s/__Y_MAX__/$y_max/g;
    $line =~ s/__X2_MAX__/$x2_max/g;
    $line =~ s/__Y2_MAX__/$y2_max/g;
    $line =~ s/__PLOT_LINES_1__/$pl1/g;
    $line =~ s/__PLOT_LINES_2__/$pl2/g;
    $line =~ s/__PLOT_LINES_3__/$pl3/g;
    $line =~ s/__PLOT_LINES_4__/$pl4/g;
}

file_write_or_die($gp_fn, @gp);

print "     generating plots.\n";

#system "cd $dn ; chmod +x $gp_fn ; $gp_fn";
system "chmod --quiet +x $gp_fn ; ./$gp_fn";

foreach my $eps_fn (`ls *$K*.eps`)
{
    chomp $eps_fn;
    print "     generated '$eps_fn'.\n";
}

if ($args{GIF}) {
    foreach my $eps_fn (`ls *$K*.eps`)
    {
        chomp $eps_fn;
        my $gif_fn = ($eps_fn =~ /^(.+).eps$/)[0] . ".gif";
        print "     converting '$eps_fn' to '$gif_fn'.\n";
        system "convert -density 120 $eps_fn $gif_fn";
    }
}

print "     done.\n\n";



sub file_write_or_die
{
    my ($fn, @l) = @_;
    open(FILE, ">$fn") or die "\n**** can't open '$fn' for writing.\n\n";
    print FILE @l;
    close(FILE);
}

sub limits 
{
    my @fns = @_;

    my $x_max_all = 0;
    my $y_max_all = 0;
    my $x2_max_all = 0;
    my $y2_max_all = 0;

    foreach my $fn (@fns) {
        open(FILE, "<$fn") or die "\n**** can't open '$fn' for reading.\n\n";
        my @lines = <FILE>;
        close(FILE);
        my @x = ();
        my @y = ();
        my @y2 = ();
        foreach my $ln (@lines) {
            if ($ln =~ /^\s*(\d+)\s+(\d+)\s+\S+\s+\S+\s+\S+\s+(\d+)\s/) {
                push @x, $1;
                push @y, $2;
                push @y2, $3;
            }
        }
        my $n = @y;
        my $i_max = 10;
        my $i2_max = 10;
        for (my $i = 11; $i < $n; $i++) {
            $i_max = $i   if ($y[$i]  > $y[$i_max]);
            $i2_max = $i  if ($y2[$i] > $y2[$i2_max]);
        }

        $x_max_all = $x[$i_max] if ($x[$i_max] > $x_max_all || $x_max_all == 0);
        $y_max_all = $y[$i_max] if ($y[$i_max] > $y_max_all || $y_max_all == 0);
        $x2_max_all = $x[$i2_max]  if ($x[$i2_max]  > $x2_max_all || $x2_max_all == 0);
        $y2_max_all = $y2[$i2_max] if ($y2[$i2_max] > $y2_max_all || $y2_max_all == 0);
    }
    return ($x_max_all, $y_max_all, $x2_max_all, $y2_max_all);
}




__DATA__
#!/usr/bin/env gnuplot

scale = 1.0


set encoding iso_8859_1

#set terminal postscript {landscape | portrait | eps | default}
#                        {enhanced | noenhanced}
#                        {color | monochrome} {solid | dashed}
#                        {<duplexing>}
#                        {"<fontname>"} {<fontsize>}
#
#set terminal postscript eps enhanced monochrome dashed 'Helvetica' 22*scale
#black = 1; grey = 0
set terminal postscript eps noenhanced color solid 'Helvetica' 22 * scale
black = 7; grey = 9



set border 31 lt black lw 4*scale
#set zeroaxis lt 7 lw 1
set grid  lt grey lw scale

size_h    = 0.80 * scale 
margin_l  = 0.25 * scale
margin_ih = 0.20 * scale
margin_r  = 0.10 * scale

size_v    = 1.00 * scale
margin_t  = 0.20 * scale
margin_iv = 0.20 * scale
margin_b  = 0.20 * scale


nh = 1
nv = 1


set lmargin  0; set rmargin  0
set tmargin  0; set bmargin  0


set tics scale -1.0 * scale


# ----------------------------------
# number of distinct kmers logx logy
# ----------------------------------

set output '__HEAD__num_distinct___K__-mers.__K__-mer_freq.log.log.eps'
set size margin_l + nh*size_h +(nh-1)*margin_ih + margin_r, \
         margin_b + nv*size_v +(nv-1)*margin_iv + margin_t 
set origin 0,0
set multiplot

  ih = 0; iv = 0;

  set log x
  set xlabel '__K__-mer frequency (copy number)' font 'Helvetica'
  set format x
  set xtics
  
  set log y
  set ylabel 'number of distinct __K__-mers' font 'Helvetica'
  set format y
  set ytics

  set grid xtics ytics lt grey
  unset key
  set size size_h, size_v
  set origin margin_l + ih*(size_h + margin_ih), margin_b+ iv*(size_v + margin_iv)
  plot [__MIN_FREQ__:__MAX_FREQ__] [] \
         __PLOT_LINES_1__

unset multiplot


# ----------------------------------
# number of distinct kmers logx liny
# ----------------------------------

set output '__HEAD__num_distinct___K__-mers.__K__-mer_freq.log.lin.eps'
set size margin_l + nh*size_h +(nh-1)*margin_ih + margin_r, \
         margin_b + nv*size_v +(nv-1)*margin_iv + margin_t 
set origin 0,0
set multiplot

  ih = 0; iv = 0;

  set log x
  set xlabel '__K__-mer frequency (copy number)' font 'Helvetica'
  set format x
  set xtics
  
  unset log y
  set ylabel 'number of distinct __K__-mers' font 'Helvetica'
  set format y
  set ytics

  set grid xtics ytics lt grey
  unset key
  set size size_h, size_v
  set origin margin_l + ih*(size_h + margin_ih), margin_b+ iv*(size_v + margin_iv)
  plot [1:1000] [] \
         __PLOT_LINES_1__

unset multiplot


# ----------------------------------
# number of distinct kmers linx liny
# ----------------------------------

set output '__HEAD__num_distinct___K__-mers.__K__-mer_freq.lin.lin.eps'
set size margin_l + nh*size_h +(nh-1)*margin_ih + margin_r, \
         margin_b + nv*size_v +(nv-1)*margin_iv + margin_t 
set origin 0,0
set multiplot

  ih = 0; iv = 0;

  unset log x
  set xlabel '__K__-mer frequency (copy number)' font 'Helvetica'
  set format x
  set xtics
  
  unset log y
  set ylabel 'number of distinct __K__-mers' font 'Helvetica'
  set format y
  set ytics

  set grid xtics ytics lt grey
  unset key
  set size size_h, size_v
  set origin margin_l + ih*(size_h + margin_ih), margin_b+ iv*(size_v + margin_iv)
  plot [0:4.0 * __X_MAX__.0] [0:1.2 *  __Y_MAX__.0] \
         __PLOT_LINES_1__

unset multiplot


# ------------------------------------
# fraction of distinct kmers logx logy
# ------------------------------------

set output '__HEAD__frac_distinct___K__-mers.__K__-mer_freq.log.log.eps'
set size margin_l + nh*size_h +(nh-1)*margin_ih + margin_r, \
         margin_b + nv*size_v +(nv-1)*margin_iv + margin_t 
set origin 0,0
set multiplot

  ih = 0; iv = 0;

  set log x
  set xlabel '__K__-mer frequency (copy number)' font 'Helvetica'
  set format x
  set xtics
  
  set log y
  set ylabel 'fraction of distinct __K__-mers' font 'Helvetica'
  set format y
  set ytics
  
  set grid xtics ytics lt grey
  unset key
  set size size_h, size_v
  set origin margin_l + ih*(size_h + margin_ih), margin_b+ iv*(size_v + margin_iv)
  plot [__MIN_FREQ__:__MAX_FREQ__] [:1] \
    __PLOT_LINES_2__

unset multiplot


# ------------------------------------
# fraction of distinct kmers logx liny
# ------------------------------------

set output '__HEAD__frac_distinct___K__-mers.__K__-mer_freq.log.lin.eps'
set size margin_l + nh*size_h +(nh-1)*margin_ih + margin_r, \
         margin_b + nv*size_v +(nv-1)*margin_iv + margin_t 
set origin 0,0
set multiplot

  ih = 0; iv = 0;

  set log x
  set xlabel '__K__-mer frequency (copy number)' font 'Helvetica'
  set format x
  set xtics
  
  unset log y
  set ylabel 'fraction of distinct __K__-mers' font 'Helvetica'
  set format y
  set ytics
  
  set grid xtics ytics lt grey
  unset key
  set size size_h, size_v
  set origin margin_l + ih*(size_h + margin_ih), margin_b+ iv*(size_v + margin_iv)
  plot [1:1000] [:1] \
    __PLOT_LINES_2__

unset multiplot


# -------------------------
# number of kmers logx logy
# -------------------------

set output '__HEAD__num___K__-mers.__K__-mer_freq.log.log.eps'
set size margin_l + nh*size_h +(nh-1)*margin_ih + margin_r, \
         margin_b + nv*size_v +(nv-1)*margin_iv + margin_t 
set origin 0,0
set multiplot

  ih = 0; iv = 0;

  set log x
  set xlabel '__K__-mer frequency (copy number)' font 'Helvetica'
  set format x
  set xtics
  
  set log y
  set ylabel 'number of __K__-mers' font 'Helvetica'
  set format y
  set ytics
  
  set grid xtics ytics lt grey
  unset key
  set size size_h, size_v
  set origin margin_l + ih*(size_h + margin_ih), margin_b+ iv*(size_v + margin_iv)
  plot [__MIN_FREQ__:__MAX_FREQ__] [:] \
    __PLOT_LINES_3__

unset multiplot


# -------------------------
# number of kmers logx liny
# -------------------------

set output '__HEAD__num___K__-mers.__K__-mer_freq.log.lin.eps'
set size margin_l + nh*size_h +(nh-1)*margin_ih + margin_r, \
         margin_b + nv*size_v +(nv-1)*margin_iv + margin_t 
set origin 0,0
set multiplot

  ih = 0; iv = 0;

  set log x
  set xlabel '__K__-mer frequency (copy number)' font 'Helvetica'
  set format x
  set xtics
  
  unset log y
  set ylabel 'number of __K__-mers' font 'Helvetica'
  set format y
  set ytics
  
  set grid xtics ytics lt grey
  unset key
  set size size_h, size_v
  set origin margin_l + ih*(size_h + margin_ih), margin_b+ iv*(size_v + margin_iv)
  plot [1:1000] [] \
    __PLOT_LINES_3__

unset multiplot


# -------------------------
# number of kmers logx liny
# -------------------------

set output '__HEAD__num___K__-mers.__K__-mer_freq.lin.lin.eps'
set size margin_l + nh*size_h +(nh-1)*margin_ih + margin_r, \
         margin_b + nv*size_v +(nv-1)*margin_iv + margin_t 
set origin 0,0
set multiplot

  ih = 0; iv = 0;

  unset log x
  set xlabel '__K__-mer frequency (copy number)' font 'Helvetica'
  set format x
  set xtics
  
  unset log y
  set ylabel 'number of __K__-mers' font 'Helvetica'
  set format y
  set ytics
  
  set grid xtics ytics lt grey
  unset key
  set size size_h, size_v
  set origin margin_l + ih*(size_h + margin_ih), margin_b+ iv*(size_v + margin_iv)
  plot [0:4.0 * __X2_MAX__.0] [0:1.2 * __Y2_MAX__.0] \
    __PLOT_LINES_3__

unset multiplot


# ---------------------------
# fraction of kmers logx logy
# ---------------------------

set output '__HEAD__frac___K__-mers.__K__-mer_freq.log.log.eps'
set size margin_l + nh*size_h +(nh-1)*margin_ih + margin_r, \
         margin_b + nv*size_v +(nv-1)*margin_iv + margin_t 
set origin 0,0
set multiplot

  ih = 0; iv = 0;

  set log x
  set xlabel '__K__-mer frequency (copy number)' font 'Helvetica'
  set format x
  set xtics
  
  set log y
  set ylabel 'fraction of __K__-mers' font 'Helvetica'
  set format y
  set ytics
  
  set grid xtics ytics lt grey
  unset key
  set size size_h, size_v
  set origin margin_l + ih*(size_h + margin_ih), margin_b+ iv*(size_v + margin_iv)
  plot [__MIN_FREQ__:__MAX_FREQ__] [:1] \
    __PLOT_LINES_4__

unset multiplot

# ---------------------------
# fraction of kmers logx liny
# ---------------------------

set output '__HEAD__frac___K__-mers.__K__-mer_freq.log.lin.eps'
set size margin_l + nh*size_h +(nh-1)*margin_ih + margin_r, \
         margin_b + nv*size_v +(nv-1)*margin_iv + margin_t 
set origin 0,0
set multiplot

  ih = 0; iv = 0;

  set log x
  set xlabel '__K__-mer frequency (copy number)' font 'Helvetica'
  set format x
  set xtics
  
  unset log y
  set ylabel 'fraction of __K__-mers' font 'Helvetica'
  set format y
  set ytics
  
  set grid xtics ytics lt grey
  unset key
  set size size_h, size_v
  set origin margin_l + ih*(size_h + margin_ih), margin_b+ iv*(size_v + margin_iv)
  plot [1:1000] [:1] \
    __PLOT_LINES_4__

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
