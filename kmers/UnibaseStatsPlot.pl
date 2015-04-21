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

# Local libraries, from elsewhere in the BroadCRD repository.
use lib "$FindBin::Bin/../";
use lib "$FindBin::Bin/";

use PerlRunTime;    # run_or_die
use PerlUtils;      # ISO_date, abort
use ArachneArgs;    # Command-line ARG=VALUE parser


my %args = getCommandArguments
    (STATS      => { value => undef,
                     help  => "List of the files to plot."},
     HEAD_OUT   => { value => "ub_stats",
                     help  => "Prefix for output files."},
     ORGANISM   => { value => "",
                     help  => "Organism name to show up at the top of each plot."},
     GIF        => { value => 0,
                     help  => "Whether to convert .eps to .gif (better for powerpoint)."} );

my $head_out = $args{HEAD_OUT};

my @fns = array_from_ref_or_value($args{STATS}) or die "no STATS found.\n";

print "     plotting unibase stats from:\n";
foreach my $fn (@fns)
{
    print "       $fn\n";
    die "can't find '$fn'.\n" if (! -e $fn || -d $fn);
}


my @titles = map /\/?([^\/]+)$/, @fns;

my $gp_fn = "$head_out.gp";

my @lts = (1,3,4,2,5,8,6,9);

# columns in kmer spectrum files
# 1: kmer_frequency
# 2: num_distinct_kmers
# 3: frac_distinct_kmers = (2)/num_distinct_kmers_total
# 4: cummulative(2)
# 5: cummulative(3)      = (4)/num_distinct_kmers_total
# 6: num_kmers           = (1)*(2)
# 7: frac_kmers          = (6)/num_kmers_total
# 8: cummulative(6)
# 9: cummulative(7)      = (8)/num_kmers_total


my $pl1 = "";  # unibase length  x  cumulative length
for (my $i = 0; $i != @fns; $i++) 
{
    $pl1 .= ", " if ($pl1 ne "");
    $pl1 .= "'$fns[$i]' u 4:5 t '$titles[$i]' w l lt $lts[$i] lw 4*scale";
}


my $pl2 = "";  # unibase length  x  cumulative number
for (my $i = 0; $i != @fns; $i++) 
{
    $pl2 .= ", " if ($pl2 ne "");
    $pl2 .= "'$fns[$i]' u 4:(column(1)+1) t '$titles[$i]' w l lt $lts[$i] lw 4*scale";
}


my @gp = <DATA>; # get the gnuplot template from the __DATA__ section at the end 


foreach my $line (@gp)
{
    $line =~ s/__ORG__/$args{ORGANISM}/g;
    $line =~ s/__HEAD__/$head_out/g;
    $line =~ s/__PLOT_LINES_1__/$pl1/g;
    $line =~ s/__PLOT_LINES_2__/$pl2/g;
}

file_write_or_die($gp_fn, @gp);

print "     generating plots.\n";

#system "cd $dn ; chmod +x $gp_fn ; $gp_fn";
system "chmod --quiet +x $gp_fn ; ./$gp_fn";

foreach my $eps_fn (`ls $head_out*.eps`)
{
    chomp $eps_fn;
    print "     generated '$eps_fn'.\n";
}

if ($args{GIF}) {
    foreach my $eps_fn (`ls $head_out*.eps`)
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

size_h    = 0.75 * scale 
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

#set key at graph 0.5, 1.05 center bottom Left reverse nobox 

set label '__ORG__' at graph 0.5, 1.13 center font 'Helvetica,24'

# ------------------------------------------
# unibase length  x  cumulative unibase length
# ------------------------------------------

set output '__HEAD__.length.log.lin.eps'
set size margin_l + nh*size_h +(nh-1)*margin_ih + margin_r, \
         margin_b + nv*size_v +(nv-1)*margin_iv + margin_t 
set origin 0,0
set multiplot

  ih = 0; iv = 0

  set title 'unibase lengths  x  cumulative length' font 'Helvetica,24'

  set log x
  set xlabel 'unibase lengths' font 'Helvetica' 
  set format x
  set xtics
  
  unset log y
  set ylabel 'cumulative unibase length' font 'Helvetica'
  set format y
  set ytics

  set key top right nobox 
  set grid xtics ytics lt grey
  set size size_h, size_v
  set origin margin_l + ih*(size_h + margin_ih), margin_b + iv*(size_v + margin_iv)
  plot [] [] \
         __PLOT_LINES_1__

unset multiplot

# ------------------------------------------
# unibase length  x  cumulative unibase number
# ------------------------------------------

set output '__HEAD__.number.log.log.eps'
set size margin_l + nh*size_h +(nh-1)*margin_ih + margin_r, \
         margin_b + nv*size_v +(nv-1)*margin_iv + margin_t 
set origin 0,0
set multiplot

  ih = 0; iv = 0

  set title 'unibase lengths  x  cumulative number' font 'Helvetica,24'

  set log x
  set xlabel 'unibase lengths' font 'Helvetica' 
  set format x
  set xtics
  
  set log y
  set ylabel 'cumulative unibase number' font 'Helvetica'
  set format y
  set ytics

  set key top right nobox 
  set grid xtics ytics lt grey
  set size size_h, size_v
  set origin margin_l + ih*(size_h + margin_ih), margin_b + iv*(size_v + margin_iv)
  plot [] [] \
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
