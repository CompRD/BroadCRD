package DistributedMake;

use strict;
use File::Temp qw/ tempfile tempdir /;
use File::Basename;

sub ParseHostsString {
    my ($hoststring) = @_;

    if ($hoststring !~ /\s+\+\s+/) {
        return undef;
    }

    my @hostobjs = split(/\s+\+\s+/, $hoststring);

    my @hosts;
    foreach my $hostobj (@hostobjs) {
        my ($multiplier, $server) = $hostobj =~ /(\d+)\*(\w+)/;
        for (my $i = 0; $i < $multiplier; $i++) {
            push(@hosts, $server);
        }
    }

    return \@hosts;
}

sub new {
    my ($class, %args) = @_;

    my %self = (
        'dryRun'         => 1,
        'numJobs'        => undef,
        'keepGoing'      => 0,
        'alwaysMake'     => 0,
        'debugging'      => 0,
        'ignoreErrors'   => 0,
        'printDirectory' => 0,
        'unlink'         => 1,
        'hosts'          => "",
        'target'         => 'all',

        %args,

        'targets'        => [],
        'hostindex'      => 0,
    );

    $self{'makefile'}  = new File::Temp(TEMPLATE => "/broad/shptmp/DistributedMake_XXXXXX", SUFFIX => ".makefile", UNLINK => $self{'unlink'}),
    $self{'hostarray'} = &ParseHostsString($self{'hosts'});

	bless \%self, $class;

	return \%self;
}

sub AddRule {
    my ($self, $targetsref, $dependenciesref, $cmdsref) = @_;
    my @targets      = (ref($targetsref)      eq 'ARRAY') ? @$targetsref      : ( $targetsref );
    my @dependencies = (ref($dependenciesref) eq 'ARRAY') ? @$dependenciesref : ( $dependenciesref );
    my @cmds         = (ref($cmdsref)         eq 'ARRAY') ? @$cmdsref         : ( $cmdsref );

    my $cmdprefix = "";
    if (defined($self->{'hostarray'})) {
        $cmdprefix = "ssh ${$self->{'hostarray'}}[$self->{'hostindex'}] ";

        $self->{'hostindex'}++;
        if ($self->{'hostindex'} == scalar(@{$self->{'hostarray'}}) - 1) {
            $self->{'hostindex'} = 0;
        }
    }

    my $rootdir = dirname($targets[0]);
    my $mkdircmd = (!-e $rootdir) ? "\n\ttest \"!-d $rootdir\" && mkdir -p $rootdir\n" : "";

    print { $self->{'makefile'} } "$targets[0]: " . join(" ", @dependencies) . "$mkdircmd\n\t$cmdprefix" . join("\n\t$cmdprefix", @cmds) . "\n";

    push(@{$self->{'targets'}}, $targets[0]);
}

sub Execute {
    my ($self, %overrides) = @_;

    print { $self->{'makefile'} } "all: " . join(" ", @{$self->{'targets'}}) . "\n";

	my %makeargs = (
        'dryRun'         => $self->{'dryRun'},
        'numJobs'        => $self->{'numJobs'},
        'keepGoing'      => $self->{'keepGoing'},
        'alwaysMake'     => $self->{'alwaysMake'},
        'debugging'      => $self->{'debugging'},
        'ignoreErrors'   => $self->{'ignoreErrors'},
        'printDirectory' => $self->{'printDirectory'},
        'target'         => $self->{'target'},
		%overrides,
	);

    my $numjobs = $makeargs{'numJobs'};
    if (!defined($numjobs)) {
        if (defined($self->{'hostarray'}) && scalar($self->{'hostarray'}) > 0) {
            $numjobs = scalar(@{$self->{'hostarray'}});
        } else {
            $numjobs = 1;
        }
    }

    my $makecmd = "make" .
                    ($makeargs{'dryRun'}         ? " -n" : "") .
                    ($makeargs{'keepGoing'}      ? " -k" : "") .
                    ($makeargs{'alwaysMake'}     ? " -B" : "") .
                    ($makeargs{'ignoreErrors'}   ? " -i" : "") .
                    ($makeargs{'printDirectory'} ? " -w" : "") .
                    ($makeargs{'debugging'} == 1 ? " -d" : "") .
                    ($makeargs{'debugging'} =~ /[abvijm]+/ ? " --debug=$makeargs{'debugging'}" : "") .
                    " -j $numjobs" .
                    " -f " . $self->{'makefile'}->filename .
                    " $makeargs{'target'}";

    print "$makecmd\n";
    system($makecmd);
    print "$makecmd\n";
}

1;
