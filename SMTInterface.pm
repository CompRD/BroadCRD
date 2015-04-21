package SMTInterface;

use strict;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);

require Exporter;

@ISA = qw(Exporter AutoLoader);

@EXPORT = qw( getMetrics getAllMetrics printMetrics );

sub getMetrics {
	my ($metricfile) = @_;
	my %metrics;

	open(METRICS, "$metricfile") || die("Could not open '$metricfile'\n");
	while (my $metric = <METRICS>) {
		if ($metric =~ /(.+)=(.+)/) {
			my ($key, $value) = ($1, $2);

			if ($value =~ /^{({.*})}/) {
				my $multiarray = $1;
				$multiarray =~ s/},{/;/g;
				$multiarray =~ s/[{}]//g;
				my @arrays = split(';', $multiarray);

				my @values;
				for (my $i = 0; $i <= $#arrays; $i++) {
					my @elements = split(',', $arrays[$i]);
					$values[$i] = \@elements;
				}

				$metrics{$key} = \@values;
			} elsif ($value =~ /^{(.*)}$/) {
				my @values = split(',', $1);
				$metrics{$key} = \@values;
			} else {
				$metrics{$key} = $value;
			}
		}
	};
	close(METRICS);

	return %metrics;
}

sub getAllMetrics {
	my ($dir) = @_;
	my ($date, $fc) = $dir =~ /\/(\d{6})_(\w+)/;

	my %allmetrics;
	for (my $i = 1; $i <= 8; $i++) {
		my $metricfile = "$dir/$fc.$i.metrics";
		if (-e $metricfile) {
			my %metric = &getMetrics($metricfile);
			$allmetrics{$i} = \%metric;
		}
	}

	return %allmetrics;
}

sub printMetrics {
	my ($metricsref) = @_;
	my %metrics = %$metricsref;

	foreach my $key (sort { $a cmp $b } keys(%metrics)) {
		if (ref($metrics{$key}) eq 'ARRAY') {
			if (ref(${$metrics{$key}}[0]) eq 'ARRAY') {
				print "$key =>\n";
				foreach my $arrayref (@{$metrics{$key}}) {
					print "\t" . join(',', @$arrayref) . "\n";
				}
			} else {
				print "$key => " . join(',', @{$metrics{$key}}) . "\n";
			}
		} else {
			print "$key => $metrics{$key}\n";
		}
	}
}

1;
