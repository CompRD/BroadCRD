# PerlMath.pm
#
# A Perl module to provide Arachne-style math functions, including:
# -- N50
#
# To use this module, follow the directions as in PerlRunTime.
#
# Josh Burton
# August 2008 

package PerlMath;
use strict;

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);
require Exporter;
@ISA = qw(Exporter AutoLoader);
@EXPORT = qw(N50);






# Returns the N50 value of an array of positive integers.
# Returns -1 if the array is empty or if any element is not a positive integer.
sub N50 {
    
    # Verify that the data in the incoming array consists of positive integers
    return -1 unless @_;
    map { return -1 if /\D/ } @_; # i.e., return -1 if any element contains any characters other than [0-9]
    
    # Sort the array (in reverse - to optimize for speed)
    my @array = sort {$b <=> $a} @_;
    return -1 if $array[-1] == 0;
    
    # Find sum of array
    my $half_sum = 0;
    map {$half_sum += $_ } @array;
    $half_sum /= 2;
    
    # Find the N50 point in the array
    my $half = 0;
    foreach my $i (0..scalar @array - 1) {
	return -1 if ($array[$i] < 0);
	$half += $array[$i];
	return $array[$i] if ($half > $half_sum);
	return 0.5 * ($array[$i] + $array[$i+1]) if ($half == $half_sum);
    }
}





1;
