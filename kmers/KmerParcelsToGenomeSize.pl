#!/usr/bin/perl -w 
###########################################################################
#       This software and its documentation are copyright (2009) by the   #
#   Broad Institute/Massachusetts Institute of Technology.  All rights    #
#   are reserved.  This software is supplied without any warranty or      #
#   guaranteed support whatsoever. Neither the Broad Institute nor MIT    #
#   can be responsible for its use, misuse, or functionality.             #
###########################################################################
#
# KmerParcelsToGenomeSize
#
# Take kmer parcels statistics for a set of reads in "kmer_frequencies.count.dat" 
# and estimate the genome size, coverage, and other statistics. 
#
# 2010-06    Filipe Ribeiro       ribeiro@broadinstitute.org
#



use strict;
use FindBin;

# ---- Local libraries, from elsewhere in the BroadCRD repository.
use lib "$FindBin::Bin/../";
use lib "$FindBin::Bin/";
use ArachneArgs; # Command-line ARG=VALUE parser



my %args = getCommandArguments(PARCELS_DIR => undef,
                               READ_LEN => 101);
  


my $dir = $args{PARCELS_DIR};

my $Nbpr = $args{READ_LEN};  # num bases per read



die "**** Can't find '$dir'\n" unless (-d $dir);


# find value of K based on directory name

my ($K) = $dir =~ /\.(\d+)merParcels\/?$/;
print "K = $K\n";



# parse kmer frequencies data file 

my @kf = ();   # kmer frequencies
my @ndk = ();  # number of DISTINCT kmers for each kmer freq
my @cndk = (); # cummulative number of DISTINCT kmers for each kmer freq
my @nk = ();   # number of kmers for each kmer freq
my @cnk = ();  # cummulative number of kmers for each kmer freq

my $fn = "kmer_frequencies.count.dat";
open(FILE, "<$dir/$fn") or die "**** Can't open '$dir/$fn'\n";
while (<FILE>)
{
    # 1:kf  2:ndk 3:ndk/Tndk  4:tndk 5:tndk/Tndk  6:nk 7:nk/Tnk  8:cnk 9:cnk/Tnk 
    if (/^\s*(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
        push @kf, $1;
        push @ndk, $2;
        push @cndk, $4;
        push @nk, $6;
        push @cnk, $8;
    }
}
close FILE;







# total num kmers

my $nk_tot = $cnk[-1];  
printf "Total number of read kmers = %s\n", number_string($nk_tot);


my $ndk_tot = $cndk[-1];  
printf "Total number of distinct kmers = %s\n", number_string($ndk_tot);



# number of reads
my $Nr = $nk_tot / ($Nbpr - $K + 1);
printf "Total number of reads = %s\n", number_string($Nr);




# find first kmer frequency minimum

my $ikf = 1;
while ($ikf < @kf && $nk[$ikf] < $nk[$ikf - 1]) {
    $ikf++;
}
my $ikf_min = $ikf-1;

print "kfmin = $kf[$ikf_min]\n";




# find kmer frequency mode after minimum 

$ikf = $ikf_min;  
my $ikf_mode = $ikf;
while ($ikf < @kf) {
    $ikf_mode = $ikf if ($nk[$ikf] > $nk[$ikf_mode]);
    $ikf++;
}

my $kcoverage = $kf[$ikf_mode];

print "kfmax = $kf[$ikf_mode]\n";




# find real minimum between 1 and maximum

$ikf = $ikf_min;
while ($ikf < $ikf_mode) {
    $ikf_min = $ikf if ($nk[$ikf] < $nk[$ikf_min]);
    $ikf++;
}

print "kfmin = $kf[$ikf_min]\n";



# define limits around main peak 

my $ikf_lo = $ikf_min;
my $ikf_hi = $ikf_mode + 1*($ikf_mode - $ikf_min);

print "kf[ikf_lo = $ikf_lo] = $kf[$ikf_lo]\n";
print "kf[ikf_mode = $ikf_mode] = $kf[$ikf_mode]\n";
print "kf[ikf_hi = $ikf_hi] = $kf[$ikf_hi]\n";


# define maximum kf above which we neglect data 
# assumes 
#          ndk[f >= 2]  ~=  4 * ndk[2] / f^2
#
#          ndk[f_max] = 1   =>  f_max = sqrt(4 * ndk[2])
#
#          because  f  is actually  f * coverage  we need to multiply by the coverage freq. 
# 

my $kf_max = sqrt(4 * $ndk[2 * $ikf_mode] * $kf[$ikf_mode]) * $kf[$ikf_mode];
my $ikf_max = $ikf_hi;
while ($ikf_max < @kf - 1 && $kf[$ikf_max] < $kf_max) { $ikf_max++; };

print "kf[ikf_max = $ikf_max] = $kf[$ikf_max]\n";
print "\n";

# number of bad and good kmers

my $nk_bad_lo    = $cnk[$ikf_lo] - 0.5 * $nk[$ikf_lo];
my $nk_bad_hi    = $nk_tot - ($cnk[$ikf_max] - 0.5 * $nk[$ikf_max]);
my $nk_good      = $nk_tot - $nk_bad_lo - $nk_bad_hi;
my $nk_good_uniq =  (($cnk[$ikf_hi] - 0.5 * $nk[$ikf_hi]) - 
                     ($cnk[$ikf_lo] - 0.5 * $nk[$ikf_lo]));
my $nk_good_rep  = $nk_good - $nk_good_uniq;

printf("nk_bad_lo    = %s (%.1f%%)\n", number_string($nk_bad_lo), 100 * $nk_bad_lo / $nk_tot);
printf("nk_bad_hi    = %s (%.1f%%)\n", number_string($nk_bad_hi), 100 * $nk_bad_hi / $nk_tot);
printf("nk_good_uniq = %s (%.1f%%)\n", number_string($nk_good_uniq), 100 * $nk_good_uniq / $nk_tot);
printf("nk_good_rep  = %s (%.1f%%)\n", number_string($nk_good_rep), 100 * $nk_good_rep / $nk_tot);
printf("nk_good      = %s (%.1f%%)\n", number_string($nk_good), 100 * $nk_good / $nk_tot);
print "\n";



# number of bad and good distinct kmers 

my $ndk_bad_lo    = $cndk[$ikf_lo] - 0.5 * $ndk[$ikf_lo];
my $ndk_bad_hi    = $ndk_tot - ($cndk[$ikf_max] - 0.5 * $ndk[$ikf_max]);
my $ndk_good      = $ndk_tot - $ndk_bad_lo - $ndk_bad_hi;
my $ndk_good_uniq = (($cndk[$ikf_hi] - 0.5 * $ndk[$ikf_hi]) - 
                     ($cndk[$ikf_lo] - 0.5 * $ndk[$ikf_lo]));
my $ndk_good_rep  = $ndk_good - $ndk_good_uniq;

printf("ndk_bad_lo    = %s (%.1f%%)\n", number_string($ndk_bad_lo), 100 * $ndk_bad_lo / $ndk_tot);
printf("ndk_bad_hi    = %s (%.1f%%)\n", number_string($ndk_bad_hi), 100 * $ndk_bad_hi / $ndk_tot);
printf("ndk_good_uniq = %s (%.1f%%)\n", number_string($ndk_good_uniq), 100 * $ndk_good_uniq / $ndk_tot);
printf("ndk_good_rep  = %s (%.1f%%)\n", number_string($ndk_good_rep), 100 * $ndk_good_rep / $ndk_tot);
printf("ndk_good      = %s (%.1f%%)\n", number_string($ndk_good), 100 * $ndk_good / $ndk_tot);
print "\n";






# average kmer frequencies

my $kf_good_ave = $nk_good / $ndk_good;
#printf("Average kmer frequency for good distinct kmers = %.1f\n", $kf_good_ave);

my $kf_ave = $nk_tot / $ndk_good;
#printf("Average kmer frequency for all distinct kmers = %.1f\n", $kf_ave);

my $kf_uniq_ave = ($cnk[$ikf_hi] - $cnk[$ikf_lo]) / ($cndk[$ikf_hi] - $cndk[$ikf_lo]);  
printf("Average kmer frequency for all distinct and unique kmers = %.1f\n", $kf_uniq_ave);



# the genome size

my $G = 0;
my $coverage = 0;

if (0) {   # original estimate, kf_mode is a bit noisy
    $G = $nk_good / $kf[$ikf_mode];
    $coverage = $Nbpr * $Nr / $G;
    
}
else {    # second estimate, fluctuates about 1/7 of above estimate
    
    #$coverage = $kf_ave * $Nbpr / ($Nbpr - $K + 1); 
    #$G = $Nbpr * $Nr / $coverage;


    my $G0 = $nk_tot / $kf_ave;    # actually, G0 = ndk_good


    my $G_uniq = $ndk_good_uniq;
    printf "G_uniq = %s bases\n", number_string($G_uniq);

    $G = $G_uniq + $nk_good_rep / $kf_uniq_ave;

    #$G = $ndk_good;

    $coverage = $Nbpr * $Nr / $G;

}


printf "Genome size = %s bases\n", number_string($G);
printf "Coverage = %.1fx\n", $coverage;

#printf "Genome fraction that is repetitive     (at K = $K) = %.1f%%\n", 100 * $nk_good_rep / $nk_good;
#printf "Genome fraction that is not repetitive (at K = $K) = %.1f%%\n", 100 * $nk_good_uniq / $nk_good;









sub number_string 
{
    my $x = shift;
    return $x if ($x < 1e3);

    return sprintf("%.2f K", $x * 1e-3) if ($x < 1e4);
    return sprintf("%.1f K", $x * 1e-3) if ($x < 1e5);
    return sprintf("%.0f K", $x * 1e-3) if ($x < 1e6);

    return sprintf("%.2f M", $x * 1e-6) if ($x < 1e7);
    return sprintf("%.1f M", $x * 1e-6) if ($x < 1e8);
    return sprintf("%.0f M", $x * 1e-6) if ($x < 1e9);

    return sprintf("%.2f G", $x * 1e-9) if ($x < 1e10);
    return sprintf("%.1f G", $x * 1e-9) if ($x < 1e11);
    return sprintf("%.0f G", $x * 1e-9) if ($x < 1e12);
}





sub kf_sd 
{
    # kmer frequency standard deviation (a very rough estimate)
    
    my $kf_sd;
    if (1) {
        my $dikf = int(0.5 + sqrt($ikf_mode));
        my $ikf_l = $ikf_mode - $dikf;
        my $ikf_r = $ikf_mode + $dikf;
        
        # in a gaussian, sigma^2 = -f(mu) / f"(mu)
        
        $kf_sd = - $nk[$ikf_mode] * 2 * $dikf / ($nk[$ikf_l] + $nk[$ikf_r] - 2 * $nk[$ikf_mode]);
    }
    
    #printf("kmer frequency stddev for good distinct kmers = %.1f\n", $kf_sd);
}


