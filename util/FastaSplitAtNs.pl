#!/usr/bin/perl -w


use strict;
use FindBin;

# ---- Local libraries, from elsewhere in the BroadCRD repository.
use lib "$FindBin::Bin/../";
use lib "$FindBin::Bin/";

use PerlRunTime;    # run_or_die
use PerlUtils;      # ISO_date, abort
use ArachneArgs;    # Command-line ARG=VALUE parser


sub Tag { return ISO_date() . " (PAPI): "; }

my %args = getCommandArguments
    (NUM_N       => { value => 25,
                      help  => "The minimum number of Ns." },
     FASTA_IN    => { value => undef,
                      help  => "Fasta input file." },
     FASTA_OUT   => { value => undef,
                      help  => "Fasta output file." });






my @lens = ();
my @seqs = ();
my @names = ();
open(FILE, "<$args{FASTA_IN}");
my $new_seq = 0;
while (<FILE>) {
    if (/^>(.+)$/) {
        push @names, $1;
        push @seqs, "";
        push @lens, 0;
    }
    elsif (/^(\w)+$/) {
        my $line = $_;
        chomp $line;
        $seqs[-1] .= $line;
        $lens[-1] += length $line;
    }
}
close(FILE);


print "N50_orig = ", N50(\@lens), "\n";


my $n = @seqs;

my $n_Ns = $args{NUM_N};
my $Ns = "N" x $n_Ns;

my @names_new = ();
my @seqs_new = ();
my @lens_new = ();

for (my $i = 0; $i != $n; $i++) {
    my $seq = $seqs[$i];
    my $name = $names[$i];
    #print "seq= $seq\n";
    my $j = 0;
    while ((my $ind = index($seq, $Ns, 0)) != -1) {
        push @names_new, "$name.$j";
        push @seqs_new, substr($seq, 0, $ind);
        push @lens_new, $ind;
        
        while (substr($seq, $ind, 1) eq "N") { $ind++; }

        $seq = substr($seq, $ind, length($seq) - $ind);
        $j++;
    }
    push @names_new, "$name.$j";
    push @seqs_new, $seq;
    push @lens_new, length $seq;
}



print "N50_new = ", N50(\@lens_new), "\n";



my $n_new = @names_new;
open(OUT, ">$args{FASTA_OUT}");
for (my $i = 0; $i != $n_new; $i++) {
    print OUT ">$names_new[$i]\n";
    my $seq = $seqs_new[$i];
    my $nb = int((length($seq) + 79) / 80);
    for (my $ib = 0; $ib < $nb; $ib++) {
        my $str = substr($seq, $ib * 80, 80);
        #if ($ib == $nb - 1) { print "str = $str\n"; }
        print OUT "$str\n";
    }
}










sub N50
{
    my ($lens0) = @_;
    my @lenss = sort { $a <=> $b } @$lens0;
    my $l_tot = 0;
    map $l_tot += $_, @lenss;

    my $l_50 = 0;
    my $i = 0;
    while ($l_50 * 2 < $l_tot) {
        $l_50 += $lenss[$i++];
    }
    return $lenss[$i-1];
}
