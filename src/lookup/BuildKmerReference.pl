#!/util/bin/perl -w

use strict;

use lib "$ENV{'ARACHNE_DIR'}";
use lib "$ENV{'ARACHNE_DIR'}/LSF";

use ArachneArgs;
use LSFInterface;
use MetricsInterface;

sub doOrDie {
	my ($lsf, @cmds) = @_;

	my @ijids;

	foreach my $cmd (@cmds) {
		print "Running $cmd...\n";
		push(@ijids, $lsf->submit($cmd));
	}

	my @failures = $lsf->wait();

	foreach my $ijid (@ijids) {
		my $cmd = $lsf->getCommand($ijid);
		my $log = $lsf->getLog($ijid);
		my $elapsedtime = $lsf->getEndTime($ijid) - $lsf->getStartTime($ijid);
		my $status = $lsf->getCleanlyExited($ijid);
		
		if ($status == 1) {
			print "Success: $cmd ($elapsedtime seconds) $log\n";
		} else {
			print "Failed: $cmd ($elapsedtime seconds) $log\n";
		}
	}

	if (@failures) {
		print "Some commands failed.  Abort.\n";
		exit(3);
	}

	$lsf->reset();
}

my $starttime = time();

my %args = &getCommandArguments("HEAD"      => { 'value' => undef,       'help' => 'Prefix for reference.fastb file' },
                                "SZ"        => { 'value' => undef,       'help' => 'Bundle size' },
                                "CLEANUP"   => { 'value' => 1,           'help' => 'If true, remove temporary files when finished' },
                                "QUEUE"     => { 'value' => 'localhost', 'help' => 'Queue to which jobs should be submitted' });

my $lsf = new LSFInterface("outputDir" => "./lsf", "checkpointFile" => "cpfile", "queue" => $args{'QUEUE'}, "verbose" => 0);
my @failures;


&doOrDie($lsf, "BuildKmerSortBundles HEAD=$args{'HEAD'} SZ=$args{'SZ'}");

my @sortcmds;
my @bundles = sort { $a cmp $b } glob("$args{'HEAD'}.kmerr.*");
foreach my $bundle (@bundles) {
	my ($bundleid) = $bundle =~ /kmerr\.(\d)/;
	push(@sortcmds, "SortKmerBundle HEAD=$args{'HEAD'} BUNDLE=$bundleid");
}
&doOrDie($lsf, @sortcmds);

my ($startBundleIndex) = $bundles[0] =~ /kmerr\.(\d)/;
my ($endBundleIndex) = $bundles[-1] =~ /kmerr\.(\d)/;

&doOrDie($lsf, "MergeSortedKmerBundles HEAD=$args{'HEAD'} START_BUNDLE=$startBundleIndex END_BUNDLE=$endBundleIndex");
&doOrDie($lsf, "RemoveKmerDuplicates HEAD=$args{'HEAD'}");

if ($args{'CLEANUP'}) {
	map { unlink } glob("$args{'HEAD'}.kmer*");
}

print "Total: " . (time() - $starttime) . " seconds\n";
