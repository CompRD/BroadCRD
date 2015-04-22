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
    (HEADS_IN => { value => undef,
		   help  => "Looks for files '<HEAD_IN>.prob' to plot."},
     HEAD_OUT => { value => "distribution",
		   help  => "Creates '<HEAD_OUT>.*.eps' files." },
     FRAC     => { value => 1.00,
		   help  => "Center quantile fraction to plot."},
     X_MIN    => { value => "",
		   help  => "Lower x limit."},
     X_MAX    => { value => "",
		   help  => "Upper x limit."},
     P_MAX    => { value => "",
		   help  => "Upper p limit."},
     TITLE    => { value => "",
		   help  => "Plot title."},
     GIF      => { value => 0,
		   help  => "Whether to convert .eps to .gif (better for powerpoint)."} );


my @fns = map $_ .= ".prob", array_from_ref_or_value($args{HEADS_IN});

map { die "can't find '$_'.\n" if (! -e "$_" || -d "$_"); } @fns;




# ---- Determining minimum and maximum

my $frac_lo = 0.5 * (1.0 - $args{FRAC});
my $frac_hi = 1.0 - $frac_lo;
my $min_all;
my $max_all;
foreach my $fn (@fns)
{
    my $min;
    my $max;
    open FILE, "<$fn";
    while (<FILE>) {
	if (/^\s*(\S+)\s+(\S+)\s+(\S+)\s*$/) {
	    $min = $1 if ($3 <= $frac_lo || ! defined $min);
	    $max = $1 if ($3 >= $frac_hi && ! defined $max);
	}
    }
    close FILE;

    $min_all = $min if (! defined $min_all || $min < $min_all);
    $max_all = $max if (! defined $max_all || $max > $max_all);
}

$min_all = $args{X_MIN} if ($args{X_MIN} ne "");
$max_all = $args{X_MAX} if ($args{X_MAX} ne "");





my $gp_fn = $args{HEAD_OUT} . ".gp";

my @lts = map "lc rgbcolor ". gnuplot_rgb($_), (1..@fns);


# columns:
# 1: x
# 2: p(x)
# 3: F(x)

my $pl1 = "";  # p(x)
for (my $i = 0; $i != @fns; $i++) 
{
    $pl1 .= ", " if ($pl1 ne "");
    #$pl1 .= "'$fns[$i]' u 1:2 not w l lt $lts[$i] lw 4*scale";
    $pl1 .= "'$fns[$i]' u 1:2 t '$i' w l $lts[$i] lw 4*scale";
}

my $pl2 = "";  # F(x)
for (my $i = 0; $i != @fns; $i++) 
{
    $pl2 .= ", " if ($pl2 ne "");
    #$pl2 .= "'$fns[$i]' u 1:3 not w l lt $lts[$i] lw 4*scale";
    $pl2 .= "'$fns[$i]' u 1:3 t '$i' w l $lts[$i] lw 4*scale";
}




my @gp = <DATA>; # get the gnuplot template from the __DATA__ section at the end 

foreach my $line (@gp)
{
    $line =~ s/__TITLE__/$args{TITLE}/g;
    $line =~ s/__HEAD_OUT__/$args{HEAD_OUT}/g;
    $line =~ s/__X_MIN__/$min_all/;
    $line =~ s/__X_MAX__/$max_all/;
    $line =~ s/__P_MAX__/$args{P_MAX}/;
    $line =~ s/__PLOT_LINES_1__/$pl1/;
    $line =~ s/__PLOT_LINES_2__/$pl2/;
}

file_write_or_die($gp_fn, @gp);

print "     generating plots.\n";

#system "cd $dn ; chmod +x $gp_fn ; $gp_fn";
system "chmod --quiet +x $gp_fn ; ./$gp_fn";

foreach my $eps_fn (`ls $args{HEAD_OUT}.*.eps`)
{
    chomp $eps_fn;
    print "     generated '$eps_fn'.\n";
}

if ($args{GIF}) {
    foreach my $eps_fn (`ls $args{HEAD_OUT}.*.eps`)
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


sub gnuplot_rgb
{
    my ($i) = @_;
    my @hex = ("0".."9", "A".."F");
    my $rgb = "'#" . ($hex[(19 * $i) % 15] .
		      $hex[(13 * $i) % 16] .
		      $hex[(5 + 21 * $i) % 15] .
		      $hex[(43 * $i) % 16] .
		      $hex[(10 + 23 * $i) % 15] .
		      $hex[(25 * $i) % 16]) . "'";
    #map $rgb .= $hex[int(rand(16))] . $hex[int(rand(16))] , (1..3);
    return $rgb;
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
margin_r  = 0.20 * scale

size_v    = 1.00 * scale
margin_t  = 0.20 * scale
margin_iv = 0.20 * scale
margin_b  = 0.20 * scale


nh = 1
nv = 1


set lmargin  0; set rmargin  0
set tmargin  0; set bmargin  0


set tics scale 1.0 * scale


# ----------------------------------
# p(x)
# ----------------------------------

set output '__HEAD_OUT__.p_x.eps'
set size margin_l + nh*size_h +(nh-1)*margin_ih + margin_r, \
         margin_b + nv*size_v +(nv-1)*margin_iv + margin_t 
set origin 0,0
set title '__TITLE__'
set multiplot

  ih = 0; iv = 0;

  unset log x
  set xlabel 'x' font 'Helvetica'
  set format x
  set xtics
  
  unset log y
  set ylabel 'p(x)' font 'Helvetica'
  set format y
  set ytics

  set grid xtics ytics lt grey
  set key left center at graph 1.03, 0.50 vertical Left reverse
  set size size_h, size_v
  set origin margin_l + ih*(size_h + margin_ih), margin_b+ iv*(size_v + margin_iv)
  plot [__X_MIN__:__X_MAX__] [0:__P_MAX__] \
         __PLOT_LINES_1__

unset multiplot


# ----------------------------------
# F(x)
# ----------------------------------

set output '__HEAD_OUT__.F_x.eps'
set size margin_l + nh*size_h +(nh-1)*margin_ih + margin_r, \
         margin_b + nv*size_v +(nv-1)*margin_iv + margin_t 
set origin 0,0
set title '__TITLE__'
set multiplot

  ih = 0; iv = 0;

  unset log x
  set xlabel 'x' font 'Helvetica'
  set format x
  set xtics
  
  unset log y
  set ylabel 'F(x)' font 'Helvetica'
  set format y
  set ytics

  set grid xtics ytics lt grey
  set key left center at graph 1.03, 0.50 vertical Left reverse
  set size size_h, size_v
  set origin margin_l + ih*(size_h + margin_ih), margin_b+ iv*(size_v + margin_iv)
  plot [__X_MIN__:__X_MAX__] [] \
         __PLOT_LINES_2__

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
