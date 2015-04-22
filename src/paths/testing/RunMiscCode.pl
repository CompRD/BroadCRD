#! /usr/bin/perl -w

## run some modules in assembly directories 

##first argument should be a full name of a file containing a list of assembly locations

#example:

#PRE=/wga/scr1/ALLPATHS DATA=N.brichardi/fullbrichardi/fullbrichardi_v1.1 RUN=iainm.05Nov2010 SUBDIR=test
#PRE=/wga/scr1/ALLPATHS DATA=O.garnettii/fullbaby/fullbaby_v5.1 RUN=iainm.19jun SUBDIR=test


$locListFile = $ARGV[0];


my %h_locs;
my $ia = 0;
open( FHIN, $locListFile ) ||
    die "did not find file=$locListFile, stopped";
while(<FHIN>){
    if ( /^\#|^\s*$/ ){ next; }
    $ia++;
    ( $pre, $data, $run, $res ) = 
	($_=~/PRE=(\S+)\s+DATA=(\S+)\s+RUN=(\S+)\s+SUBDIR=(\S+)/ );
    $species = $data; $species =~s/^\///; $species =~s/\/.*//;

    $h_locs{$ia}{"species"}   = $species;
    $h_locs{$ia}{"pre"}       = $pre;
    $h_locs{$ia}{"data"}      = $data;
    $h_locs{$ia}{"run"}       = $run;
    $h_locs{$ia}{"res"}       = $res;
    $dataDir = $pre."/".$data;
    $runDir  = $dataDir."/".$run;
    $subDir  = $runDir."/ASSEMBLIES/".$res;
    $h_locs{$ia}{"dataDir"}   = $dataDir;
    $h_locs{$ia}{"runDir"}    = $runDir;
    $h_locs{$ia}{"subDir"}    = $subDir;
}
close FHIN;


foreach $ia (sort {$a<=>$b} keys %h_locs ){
    $species = $h_locs{$ia}{"species"};
    print "SPECIES=$species\n";

    # Compute graph unibase overlap adjacency graph if it does not exist
    my $graphFile = $h_locs{$ia}{"runDir"}."/all_reads.unibases.ovrlp_adjgraph.k96";    
    if ( ! -e $graphFile ){
	$cmd1 = "BuildUnipathAdjGraph PRE=".$h_locs{$ia}{"pre"}." DATA=".$h_locs{$ia}{"data"}." RUN=".$h_locs{$ia}{"run"}." K=96 BUILD_UNIPATH_GRAPH=False READS=all_reads";
	print $cmd1."\n";
	system($cmd1) == 0 || 
	    die "failed on command=$cmd1, stopped";
    }
    
    if ( ! -e $graphFile ){
	print "$graphFile not found\n";
	next;
    }

    $cmd2 = "ComputeGraphBubbliness GRAPH_IN= ".$h_locs{$ia}{"runDir"}."/all_reads.unibases.ovlp_adjgraph.k96";
    print $cmd2."\n";
    system($cmd2) == 0 || 
	die "failed on command=$cmd2, stopped";
}

