#!/usr/bin/perl -w
#############################################################################
#                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     #
#       This software and its documentation are copyright (2009) by the     #
#   Broad Institute.  All rights are reserved.  This software is supplied   #
#   without any warranty or guaranteed support whatsoever. The Broad        #
#   Institute is not responsible for its use, misuse, or functionality.     #
#############################################################################
# 
# This utility plots scaffolds centered around a specific contig, using gnuplot.
# It expects data to be present on file 'data.txt'.
# Usage:
#  
#  links-plot.pl <contig_id>
#  gv <contig_id>.ps

use strict;

die unless @ARGV;

my ($s1) = @ARGV;

my $dat_fn = "$s1.dat";

if (! -e $dat_fn) {
    system "grep '#' data.txt > $dat_fn";
    system "grep \"^$s1 \" data.txt >> $dat_fn";
}

open FILE, "<$dat_fn";
my @lines = <FILE>;
close FILE;

print shift @lines;    
foreach my $line (@lines) {
    printf("%6d %7d %6d %7d %7d %4d %4d %6.2f %6.2f %6.2f %6.2f %6.2f\n", 
           split " ", $line);
}




my $gp_fn = "$s1.gp";
my $ps_fn = "$s1.ps";
    
my @gp = <DATA>; # get the gnuplot template from the __DATA__ section at the end 
    

my $n = (split " ", `wc -l $dat_fn`)[0] - 1;

foreach my $line (@gp)
{
    $line =~ s/__N__/$n/;
    $line =~ s/__DATA_FN__/$dat_fn/;
    $line =~ s/__PS_FN__/$ps_fn/;
}
    
file_write_or_die($gp_fn, @gp);
    
print "     generating plot\n";
    
#system "cd $dn ; chmod +x $gp_fn ; $gp_fn";
system "chmod --quiet +x $gp_fn ; $gp_fn";
    
#print "     showing plot\n";
#exec "gv $ps_fn";




sub file_write_or_die
{
    my ($fn, @l) = @_;
    open(FILE, ">$fn") or die "\n**** can't open '$fn' for writing.\n\n";
    print FILE @l;
    close(FILE);
}



__DATA__
#!/usr/bin/env gnuplot

scale = 0.43

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
#set grid  lt grey lw scale

size_h    = 1.60 * scale 
margin_l  = 0.15 * scale
margin_ih = 0.10 * scale
margin_r  = 0.15 * scale

size_v    = 1.00 * scale
margin_t  = 0.20 * scale
margin_iv = 0.10 * scale
margin_b  = 0.20 * scale


nh = 1
nv = 1

set size margin_l + nh*size_h +(nh-1)*margin_ih + margin_r, \
         margin_b + nv*size_v +(nv-1)*margin_iv + margin_t 

set lmargin  0; set rmargin  0
set tmargin  0; set bmargin  0

#set nokey


set tics scale -1.0 * scale

set output '__PS_FN__'


n = __N__

min(a,b) = a < b ? a : b
max(a,b) = a > b ? a : b


sign(x) = (x >= 0) ? 1 :  -1;

signs(x) = (x >= 0) ? "+" : "-";




xs1min(size1)  = -0.5 * abs(size1)
xs1max(size1)  =  0.5 * abs(size1)
xshalf(size1) =  0.5 * abs(size1)

csep(size1, size2, sep) = (0.5*abs(size1) + sep + 0.5*abs(size2))  


xs2med(size1, size2, sep) = sign(size1) * csep(size1, size2, sep)

xs2min(size1, size2, sep) = xs2med(size1, size2, sep) - 0.5*abs(size2)
xs2max(size1, size2, sep) = xs2med(size1, size2, sep) + 0.5*abs(size2)



xl1med(size1, i10, i11) = xs1min(size1) + 0.5 * (i10 + i11)
c2end(size2, i20, i21) = (size2 > 0) ? i21 : -size2  - i20;


dx1lib(sep, sizelib, sizeread) = - sign(sep) * (sizelib + 2 * sizeread)
dx2lib(sep, sizelib, sizeread) =   sign(sep) * (sizelib + 2 * sizeread)


xlibmed(size, pct0, pct1) = (0.005 * (pct1 + pct0) - 0.5) * size
xlibhalf(size, pct0, pct1) = abs(0.005 * (pct1 - pct0) * size)


set bars small


set multiplot

  ih = 0; iv = 0;

  set xlabel 'seps' font 'Helvetica'
  set format x
  set xtics
  
  set ylabel '' font 'Helvetica'
  set format y ''
  unset ytics 
  #set mytics 1
  #set my2tics 1
  
  set grid xtics lt grey lw 0.5
  #set key left center at graph 1.03,1.10 vertical Left reverse height -40
  set key off
  set size size_h, size_v
  set origin margin_l + ih*(size_h + margin_ih), margin_b+ iv*(size_v + margin_iv)
   
  plot [] [-0.2:n+0.2] \
   '__DATA_FN__' u (0):(0.5*n):(xshalf($2)):(0.5*n - 0.1)        w boxxy not lw 2, \
   '' u (xs2med($2, $4, $5)):($0 + 0.5):(0.5*abs($4)):(0.4) w boxxy not lt 3,\
   '' u (xs2min($2, $4, $5)):($0 + 0.5):($6    ):(0.35)  w boxxy not lt 4,\
   '' u (xs2max($2, $4, $5)):($0 + 0.5):($6    ):(0.35)  w boxxy not lt 4,\
   '' u (0):(n):(sprintf("+ %d", $1)) w labels not, \
   '' u (xs2med($2, $4, $5)):($0 + 1.0):(sprintf("%s %d/%d/%.2f", signs($2*$4), $3, $7, $8))        w labels not, \
   '' u (xs2med($2, $4, $5) + xlibmed($4 * sign($2), $11, $12)):($0 +0.5):(xlibhalf($4, $11, $12)):(0.2) w boxxy not lt 5 lw 2,\
   '' u (xlibmed(abs($2), $9, $10)):($0 +0.5):(xlibhalf($2, $9, $10)):(0.2) w boxxy not lt 2 lw 2


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



