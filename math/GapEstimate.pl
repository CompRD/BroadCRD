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
    (HEAD_PROB  => { value => undef,
                     help  => "Looks for jump distribution in '<HEAD_PROB>.prob'."},
     S1         => { value => undef,
                     help  => "Size of contig 1." },
     S2         => { value => undef,
                     help  => "Size of contig 2." },
     G          => { value => undef,
                     help  => "Size of gap." },
     COV        => { value => 30,
                     help  => "Coverage." },
     CYCLES     => { value => 1,
                     help  => "Number of cycles." });



my $prob = prob_read($args{HEAD_PROB} . ".prob");

#my $n0 = 4000;
#my $p_0 = 0.01 / $n0;
#my $prob = {x0 => 0,
#	    p => [map { ($_ == 1000) ? 0.90 : ($_ == 3000) ? 0.09 : $p_0; } (0..$n0-1)],
#	    F => [map { $p_0 * ($_ + 1) + (($_ < 1000) ? 0.0 : ($_ < 3000) ? 0.90 : 0.99); } (0..$n0-1)]};

#prob_write($prob, "test.prob");

for (my $cycle = 0; $cycle != $args{CYCLES}; $cycle++) {

    my $S1 = $args{S1};
    my $S2 = $args{S2};
    my $G = $args{G};
    my $cov = $args{COV};
    my $len = 101;

    print "Simulating links\n" unless $cycle;
    my ($l12, $l1, $l2) = simulate_links($S1, $S2, $G, $cov, $len, $prob);


    print "Computing gap probability 0\n" unless $cycle;
    my $x0 = $prob->{x0};
    my $n = scalar @{$prob->{p}};


    my $p_G0 = prob_unif($prob);
    
    
    my $n12 = 0;
    my $n21 = 0;
    my $n11 = 0;
    my $n22 = 0;
    foreach my $l (@$l12) {
        my $logp0s = $p_G0->{logp};
        my $p0s = $p_G0->{p};
        my ($d1, $d2) = @$l;
        if ($d1 >= 0 && $d2 >= 0) {
            #map $logp0s->[$_] += prob_log(($x0 + $_) + $d1 + $d2, $prob), (0 .. $n - 1);
            #prob_log_normalize($p_G0);
            map $p0s->[$_] *= prob(($x0 + $_) + $d1 + $d2, $prob), (0 .. $n - 1);
            prob_normalize($p_G0);
            $n12 ++;
        }
        elsif ($d1 < 0 && $d2 < 0) {
            #map $logp0s->[$_] += prob_log(($x0 + $_) - $d1 - $d2, $prob), (0 .. $n - 1);
            #prob_log_normalize($p_G0);
            map $p0s->[$_] += prob(($x0 + $_) - $d1 - $d2, $prob), (0 .. $n - 1);
            prob_normalize($p_G0);
            $n21 ++;
        }
        elsif ($d1 >= 0 && $d2 < 0) {
            $n11++;
        }
        elsif ($d1 < 0 && $d2 >= 0) {
            $n22++;
        }
    }
    
    unless ($cycle) {
        printf "n10 = %d\n", scalar @$l1;
        print  "n11 = $n11\n";
        print  "n12 = $n12\n";
        printf "n20 = %d\n", scalar @$l2;
        print  "n21 = $n21\n";
        print  "n22 = $n22\n";
    }
    prob_write($p_G0, "G0.prob") unless $cycle;
    



    print "Computing gap probability 1\n" unless $cycle;
    my $p_G1 = prob_unif($prob);
    
    foreach my $d1 (@$l1) {
        my $logp1s = $p_G1->{logp};
        my $p1s = $p_G1->{p};
        if ($d1 >= 0) {
            my $p11 = prob_in($d1 - $S1, $d1 - $len, $prob);
            #map $logp1s->[$_] += log(1.0 - ($p11 + prob_in($d1 + ($x0 + $_), 
            #                                               $d1 + ($x0 + $_) + $S2 - $len - 1, 
            #                                               $prob))), (0 .. $n-1);
            map $p1s->[$_] *= 1.0 - ($p11 + prob_in($d1 + ($x0 + $_), 
                                                    $d1 + ($x0 + $_) + $S2 - $len - 1, 
                                                    $prob)), (0 .. $n-1);
        }
        else {
            my $p22 = prob_in($d1, $d1 + $S2 - $len - 1, $prob);
            #map $logp1s->[$_] += log(1.0 - ($p22 + prob_in($d1 - $S1 - ($x0 + $_), 
            #                                               $d1 - $len - ($x0 + $_), 
            #                                               $prob))), (0 .. $n-1);
            map $p1s->[$_] *= 1.0 - ($p22 + prob_in($d1 - $S1 - ($x0 + $_), 
                                                    $d1 - $len - ($x0 + $_), 
                                                    $prob)), (0 .. $n-1);
        }    
        #prob_log_normalize($p_G1);
        prob_normalize($p_G1);
    }        
    
    
    
    print "Computing gap probability 2\n" unless $cycle;
    my $p_G2 = prob_unif($prob);
    
    foreach my $d2 (@$l2) {
        my $logp2s = $p_G2->{logp};
        my $p2s = $p_G2->{p};
        if ($d2 >= 0) {
            my $p22 = prob_in($d2 - $S2, $d2 - $len, $prob);
            #map $logp2s->[$_] += log(1.0 - ($p22 + prob_in($d2 + ($x0 + $_), 
            #                                               $d2 + ($x0 + $_) + $S1 - $len - 1, 
            #                                               $prob))), (0 .. $n-1);
            map $p2s->[$_] *= 1.0 - ($p22 + prob_in($d2 + ($x0 + $_), 
                                                    $d2 + ($x0 + $_) + $S1 - $len - 1, 
                                                    $prob)), (0 .. $n-1);
        }
        else {
            my $p11 = prob_in($d2, $d2 + $S1 - $len - 1, $prob);
            #map $logp2s->[$_] += log(1.0 - ($p11 + prob_in($d2 - $S2 - ($x0 + $_), 
            #                                               $d2 - $len - ($x0 + $_), 
            #                                               $prob))), (0 .. $n-1);
            map $p2s->[$_] *= 1.0 - ($p11 + prob_in($d2 - $S2 - ($x0 + $_), 
                                                    $d2 - $len - ($x0 + $_), 
                                                    $prob)), (0 .. $n-1);
        }    
        #prob_log_normalize($p_G2);
        prob_normalize($p_G2);
    }
    print "\n" unless $cycle;
    
    prob_write($p_G1, "G1.prob");
    prob_write($p_G2, "G2.prob");

    
    my $p_G01 = prob_unif($prob);
    my $p01s = $p_G01->{p};
    map $p01s->[$_] = $p_G0->{p}[$_] * $p_G1->{p}[$_], (0 .. $n-1);
    prob_write($p_G01, "G01.prob");
    prob_print($p_G01, "p_G01");


    my $p_G02 = prob_unif($prob);
    my $p02s = $p_G02->{p};
    map $p02s->[$_] = $p_G0->{p}[$_] * $p_G2->{p}[$_], (0 .. $n-1);
    prob_write($p_G02, "G02.prob");
    prob_print($p_G02, "p_G02");
    



    my $p_Gw = prob_unif($prob);
    my $pws = $p_Gw->{p};
    map $pws->[$_] = $p_G0->{p}[$_] * $p_G1->{p}[$_] * $p_G2->{p}[$_], (0 .. $n-1);
    prob_write($p_Gw, "Gw.prob");
    prob_print($p_Gw, "p_Gw");
    


    my $p_G = prob_unif($prob);
    my $ps = $p_G->{p};
    map $ps->[$_] = sqrt($p_G01->{p}[$_] * $p_G02->{p}[$_]), (0 .. $n-1);
    prob_write($p_G, "G.prob");
    prob_print($p_G, "p_G");




    write_links($l12, $l1, $l2, $S1, $S2, $G, "links.pos");
}








sub prob_print
{
    my ($p, $label) = @_;

    prob_normalize($p);
    my $median = quantile(0.5, $p);

    printf("$label: median = %d    1sig = (%d, %d)   2sig = (%d, %d)   3sig = (%d, %d)    percentile = %f\n", 
	   $median, 
	   quantile(    0.158, $p) - $median, 
	   quantile(1 - 0.158, $p) - $median,
	   quantile(    0.022, $p) - $median, 
	   quantile(1 - 0.022, $p) - $median,
	   quantile(    0.001, $p) - $median, 
	   quantile(1 - 0.001, $p) - $median,
	   prob_le($args{G}, $p));

    print "\n";
}
    



sub prob
{
    my ($x, $prob) = @_;
    my $x0 = $prob->{x0};
    die "x0 = $x0\n" if (!defined $x);
    die "x = $x\n" if (!defined $x0);
    my $i = $x - $x0;

    return 0.0 if ($i < 0);
    my $p = $prob->{p};
    return 0.0 if ($i >= @$p);
    return $p->[$i];
}

sub prob_log
{
    my ($x, $prob) = @_;
    my $x0 = $prob->{x0};
    die "x0 = $x0\n" if (!defined $x);
    die "x = $x\n" if (!defined $x0);
    my $i = $x - $x0;

    return -1e1000 if ($i < 0);
    my $logp = $prob->{logp};
    return -1e1000 if ($i >= @$logp);
    return $logp->[$i];
}





sub prob_le
{
    my ($x, $prob) = @_;
    my $x0 = $prob->{x0};
    my $i = $x - $x0;
    return 0.0 if ($i < 0);
    my $F = $prob->{F};
    return 1.0 if ($i >= @$F);
    return $F->[$i];
}


sub prob_in
{
    my ($a, $b, $prob) = @_;
    die  unless ($b >= $a);
    my $p = prob_le($b, $prob) - prob_le($a - 1, $prob);
    die unless ($p >= 0);
    return $p;    
}



sub quantile
{
    my ($y, $prob) = @_;
    my $x0 = $prob->{x0};
    my $F = $prob->{F};
    my $n = scalar @$F;
    my $i = 0;
    while ($i < $n && $F->[$i] < $y) { $i++; }
    die if ($i == $n);

    return $x0 + $i;

    # ---- ths is off
    return $x0 if ($i == 0);

    my $dF = $F->[$i] - $F->[$i - 1];
    my $dFy = $F->[$i] - $y;

    return $x0 + ($i - 1) if ($dFy / $dF > 0.5); 
    return $x0 + $i;
}



sub rand_from_distrib
{
    my ($prob) = @_;
    return quantile(rand(), $prob);
}



sub d1 
{
    my ($S1, $S2, $G, $len, $end1) = @_;

    if ($end1 >= $len - 1 && $end1 < $S1) {
        return $S1 - $end1 - 1;
    } 
    elsif ($end1 >= $S1 + $G + $len - 1 && $end1 < $S1 + $G + $S2) {
        return $S1 + $G - $end1 - 1;
    } 
    return undef;
}



sub d2
{
    my ($S1, $S2, $G, $len, $end2) = @_;
    
    if ($end2 >= $S1 + $G && $end2 <= $S1 + $G + $S2 - $len) {
        return $end2 - $S1 - $G;
    } 
    elsif ($end2 >= 0 && $end2 <= $S1 - $len) {
        return $end2 - $S1;
    }
    return undef;
}
        



sub simulate_links
{
    my ($S1, $S2, $G, $cov, $len, $prob) = @_;

    my $Stot = $S1 + $S2 + $G;
    my $n = $cov * $Stot / $len / 2.0;

    my @l12 = ();
    my @l1 = ();
    my @l2 = ();

    printf "%d links\n", $n;

    # read 1
    for (my $i = 0; $i < $n/2; $i++) {
        my $end1 = int(0.5 + $len + rand($Stot));
	my $rnd = rand_from_distrib($prob);
        my $end2 = $end1 + $rnd;

	my $d1 = d1($S1, $S2, $G, $len, $end1);
	my $d2 = d2($S1, $S2, $G, $len, $end2);

	if (defined $d1 && defined $d2) {
		push @l12, [ $d1, $d2 ];
	}
	elsif (defined $d1) {
	    push @l1, $d1;
	}
	elsif (defined $d2) {
	    push @l2, $d2;
	}
    }
    # read 2
    for (my $i = 0; $i < $n/2; $i++) {
        my $end2 = int(0.5 + $len + rand($Stot));
	my $rnd = rand_from_distrib($prob);
        my $end1 = $end2 - $rnd;

	my $d1 = d1($S1, $S2, $G, $len, $end1);
	my $d2 = d2($S1, $S2, $G, $len, $end2);

	if (defined $d1 && defined $d2) {
		push @l12, [ $d1, $d2 ];
	}
	elsif (defined $d1) {
	    push @l1, $d1;
	}
	elsif (defined $d2) {
	    push @l2, $d2;
	}
    }

    return (\@l12, \@l1, \@l2);
}


sub simulate_links_old
{
    my ($S1, $S2, $G, $cov, $len, $prob) = @_;

    my $Stot = $S1 + $S2 + $G;
    my $n = $cov * $Stot / $len / 2.0;

    my @l12 = ();
    my @l1 = ();
    my @l2 = ();

    printf "%d links\n", $n;

    for (my $i = 0; $i < $n; $i++) {
	#print " $i";
	
        my $good_link = 1;
        
        my $end1 = int(0.5 + $len + rand($Stot));

	my $rnd = rand_from_distrib($prob);
	#print "rnd $rnd\n";
        my $end2 = $end1 + $rnd;

        if ($end2 > $Stot || $end2 < 0) {
            $good_link = 0;
            $end2 %= $Stot;
        }

	my $d1 = d1($S1, $S2, $G, $len, $end1);
	my $d2 = d2($S1, $S2, $G, $len, $end2);

	if (defined $d1 && defined $d2) {
	    if ($good_link) {
		push @l12, [ $d1, $d2 ];
	    }
	    else {
		push @l1, $d1;
		push @l2, $d2;
	    }
	}
	elsif (defined $d1) {
	    push @l1, $d1;
	}
	elsif (defined $d2) {
	    push @l2, $d2;
	}
    }
    #print "\n";

    return (\@l12, \@l1, \@l2);
}







sub prob_unif
{
    my ($prob) = @_;
    my $n = scalar @{$prob->{p}};
    my $p0 = 1.0 / $n;
    my $logp0 = log($p0);
    return { x0 => $prob->{x0}, 
	     p => [map $p0, (1..$n)],
             logp => [map $logp0, (1..$n)],
	     F => [map $_ / $n, (1..$n)] };
}





sub prob_read
{
    my ($fn) = @_;

    my $x0 = undef;
    my @p = ();
    my @logp = ();
    my @F = ();
    
    open FILE, "<$fn" or die "can't find '$fn'.\n";
    while (<FILE>) {
	if (/^\s*([^\#\s]\S*)\s+(\S+)\s+(\S+)\s*/) {
	    $x0 = $1 unless defined $x0;
            my $p = 1.0e-100 + $2;
	    push @p, $p;
            push @logp, log($p);
	    push @F, $3;
	}
    }
    close FILE;
    
    return { x0 => $x0, p => \@p, logp => \@logp, F => \@F };
}




sub prob_normalize
{
    my ($prob) = @_;

    my $p = $prob->{p};
    my $F = $prob->{F};
    $F->[0] = $p->[0];
    map $F->[$_] = $F->[$_ - 1] + $p->[$_], (1..@$p-1);
    my $sum = $F->[-1];
    map $_ /= $sum, @$p;
    map $_ /= $sum, @$F;
}


sub prob_log_normalize
{
    my ($prob) = @_;

    my $logp = $prob->{logp};
    my $p = $prob->{p};
    die unless defined $logp->[0];
    map $p->[$_] = exp($logp->[$_]), (0..@$p-1);

    my $sum = 0.0;
    map $sum += $_, @$p;
    map $_ /= $sum, @$p;

    my $logsum = log($sum);
    map $_ -= $logsum, @$logp;


}


sub prob_write
{
    my ($prob, $fn) = @_;

    prob_normalize($prob);
    open FILE, ">$fn" or die "can't open '$fn'.\n";
    my $F = 0.0;
    my $ps = $prob->{p};
    my $x0 = $prob->{x0};
    for (my $i = 0; $i != @$ps; $i ++) {
	my $p = $ps->[$i];
	$F += $p;
	printf FILE " %6d %12.8f %12.8f\n", $x0 + $i, $p, $F;
    }
    close FILE;
}




sub write_links
{
    my ($l12, $l1, $l2, $S1, $S2, $G, $fn) = @_;

    open FILE, ">$fn" or die "can't open '$fn'.\n";

    foreach my $l (@$l12) {
	my ($d1, $d2) = @$l;
	my $end1 = - $d1 + (($d1 >= 0) ? 0 : $G);
	my $end2 = + $d2 + (($d2 >= 0) ? $G : 0);
	printf FILE "%6d %6d\n", $end1, $end2;
    }
    print FILE "\n\n";

    foreach my $d1 (@$l1) {
	my $end1 = - $d1 + (($d1 >= 0) ? 0 : $G);
	printf FILE "%6d\n", $end1;
    }


    print FILE "\n\n";
    foreach my $d2 (@$l2) {
	my $end2 = + $d2 + (($d2 >= 0) ? $G : 0);
	printf FILE "%6d\n", $end2;
    }
    close FILE;
}








