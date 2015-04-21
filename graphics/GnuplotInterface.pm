package GnuplotInterface;

use File::Temp qw/ tempfile tempdir /;

sub new {
	my ($class, %args) = @_;

	my %self = (
		'title'  => '',
		'xlabel' => '',
		'ylabel' => '',
		'xrange' => '',
		'yrange' => '',
		'xscale' => '',
		'yscale' => '',
		'xtics'  => '',
		'ytics'  => '',
		'using'  => '',
		'output' => "gpi_" . time() . ".ps",
		%args,
	);

	bless(\%self, $class);

	open($self{'gp'}, "|gnuplot >/dev/null 2>/dev/null");

	(\%self)->setTitle($self{'title'});
	(\%self)->setXLabel($self{'xlabel'});
	(\%self)->setYLabel($self{'ylabel'});
	(\%self)->setXRange($self{'xrange'});
	(\%self)->setYRange($self{'yrange'});
	(\%self)->setXScale($self{'xscale'});
	(\%self)->setYScale($self{'yscale'});
	(\%self)->setXTics($self{'xtics'});
	(\%self)->setYTics($self{'ytics'});

	return \%self;
}

sub execute {
	my ($self, $cmd) = @_;
	print { ($self->{'gp'}) } "$cmd;\n";
}

sub setXLabel {
	my ($self, $xlabel) = @_;
	$self->execute("set xlabel '$xlabel'");
}

sub setYLabel {
	my ($self, $ylabel) = @_;
	$self->execute("set ylabel '$ylabel'");
}

sub setTitle {
	my ($self, $title) = @_;
	$self->execute("set title '$title'");
}

sub setOutput {
	my ($self, $output) = @_;

	my @path = split('/', $output);
	$path[-1] =~ s/\.(\w+)$//;
	$output = join("/", @path);
	$self->{'output'} = $output;

	$self->execute("set term postscript color solid");
	$self->execute("set output '$output.ps'");
}

sub setXRange {
	my ($self, $range) = @_;

	if ($range =~ /^?(\d+):?(\d+)$/) {
		$self->execute("set xrange [$range]");
	}
}

sub setYRange {
	my ($self, $range) = @_;

	if ($range =~ /^?(\d+):?(\d+)$/) {
		$self->execute("set yrange [$range]");
	}
}

sub setXScale {
	my ($self, $scale) = @_;

	if ($scale =~ /log/) {
		$self->execute("set logscale x");
	} else {
		$self->execute("set nologscale x");
	}
}

sub setYScale {
	my ($self, $scale) = @_;

	if ($scale =~ /log/) {
		$self->execute("set logscale y");
	} else {
		$self->execute("set nologscale y");
	}
}

sub setXTics {
	my ($self, $xtics) = @_;

	if (ref($xtics) eq 'ARRAY') {
		$self->execute("set xtics (" . join(",", @{$xtics}) . ")");
	}
}

sub setYTics {
	my ($self, $ytics) = @_;

	if (ref($ytics) eq 'ARRAY') {
		$self->execute("set ytics (" . join(",", @{$ytics}) . ")");
	}
}

sub plot {
	my ($self, @plotitems) = @_;
	my @plots;

	foreach my $plotitem (@plotitems) {
		my %plothash = %$plotitem;
		my %fullplothash = (
		                     'using' => $self->{'using'},
		                     'title' => '',
		                     'with'  => 'lines lw 2',
		                     %plothash,
		                   );

		if (exists($fullplothash{'file'}) && -e $fullplothash{'file'} && -s $fullplothash{'file'} > 0) {
			push(@plots, "'$fullplothash{'file'}'" .
					(($fullplothash{'using'}) eq '' ? "" : " using $fullplothash{'using'} ") .
					" title '$fullplothash{'title'}'" .
					(($fullplothash{'with'})  eq '' ? "" : " with $fullplothash{'with'} ")
			    );
		} elsif (exists($fullplothash{'data'}) && ref($fullplothash{'data'}) eq 'ARRAY') {
			chomp(my @lines = @{$fullplothash{'data'}});

			(my $fh, $fullplothash{'file'}) = tempfile("gnuplot_data_XXXX", "UNLINK" => 1);
			open(FILE, ">$fullplothash{'file'}");
			map { print FILE "$_\n"; } @lines;
			close(FILE);

			push(@plots, "'$fullplothash{'file'}'" .
					(($fullplothash{'using'}) eq '' ? "" : " using $fullplothash{'using'} ") .
					" title '$fullplothash{'title'}'" .
					(($fullplothash{'with'})  eq '' ? "" : " with $fullplothash{'with'} ")
			    );
		}
	}

	if ($#plots >= 0) {
		$self->setOutput($self->{'output'});
		my $output = $self->{'output'};
        unlink($output);

		$self->execute("plot " . join(", ", @plots));
		system("gs -q -sPAPERSIZE=letter -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dAutoRotatePages=/All -sOutputFile=$output.pdf $output.ps >/dev/null 2>&1");

		unlink("$output.ps");
	}
}

1;
