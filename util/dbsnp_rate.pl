#!/usr/bin/perl

$|=1;
use strict;

my ($dbsnp_file, $allele_file, $threshold_start, $threshold_step, $threshold_end, $output) = @ARGV;

my $temp_dbsnp_file    = "$output.temp_truth";
my $temp_allele_file   = "$output.temp_alleles";
my $temp_joined_file_1 = "$output.temp_joined_1";
my $temp_joined_file_2 = "$output.temp_joined_2";

open(TRUTH, $dbsnp_file) or die;

my $command;

$command = "cat $dbsnp_file | perl -pe 's/\\\// /g;' | awk '{ if (\$5 < \$6) { \$7 = \$5; \$5 = \$6; \$6 = \$7; } print \"chr\"\$1\":\"\$2-1,\$4,\$5\$6}' | perl -pe 's/chrX/23/; s/chrY/24/; s/chrM/chr0/;' | sort -u -k 1 | perl -pe 'if (\$_ =~ m/\\-/) { s/A/T/g; s/C/G/g; s/G/C/g; s/T/A/g; s/\\+//g;}' > $temp_dbsnp_file";
print "\n$command\n";
system($command);
if ($? != 0) { die("failed: $command"); }

$command = "cat $allele_file |  sort -k 1 > $temp_allele_file";
print "\n$command\n";
system($command);
if ($? != 0) { die("failed: $command"); }

$command = "join -j 1 $temp_dbsnp_file $temp_allele_file > $temp_joined_file_1";
print "\n$command\n";
system($command);
if ($? != 0) { die("failed: $command"); }

open(J_IN, $temp_joined_file_1) or die;
open(J_OUT, ">$temp_joined_file_2") or die;
sub alpha_sort
{
    my ($s, $strand) = @_;
    my @s = split("", $s);
    my @sorted_s = sort { $a <=> $b } @s;
    my $ans = join("", @sorted_s);
    return $ans;
}
while(<J_IN>)
{
    $_ =~ s/(^\s+)|(\s+$)//g;
    my @tokens     = split(/\s+/, $_);
    my $strand     = $tokens[1];
    my $genotype_1 = alpha_sort($tokens[2], $strand);
    my $genotype_2 = alpha_sort($tokens[4], $strand);
    $tokens[1] = "+"; 
    $tokens[2] = $genotype_1;
    $tokens[4] = $genotype_1;
    print J_OUT join(" ", @tokens) . "\n";
}
close(J_IN);
close(J_OUT);

open(OUT, ">$output") or die;
print OUT "threshold number_of_SNPs number_of_dbSNPs percent_of_SNPs_in_dbSNP\n";
for (my $threshold = $threshold_start; $threshold <= $threshold_end; $threshold += $threshold_step)
{
    my $number_of_SNPs   = `cat $temp_allele_file    | awk '{if (\$5 >= $threshold)                   { print; }}' | wc -l;`; $number_of_SNPs =~ s/(^\s+)|(\s+)$//g;
    my $number_of_dbSNPs = `cat $temp_joined_file_2  | awk '{if ((\$7 >= $threshold) && (\$3 == \$5)) { print; }}' | wc -l;`; $number_of_dbSNPs =~ s/(^\s+)|(\s+)$//g;

    print OUT sprintf("%0.01f %d %d %0.01f\n", $threshold, $number_of_SNPs, $number_of_dbSNPs, 100.0 * $number_of_dbSNPs / $number_of_SNPs);
}
close(OUT);

