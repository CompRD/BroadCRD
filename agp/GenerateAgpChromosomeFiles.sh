#! /usr/bin/tcsh -f

###
# GenerateAgpChromosomeFiles
#
# Run CreateChromosomeFasta and CreateChromosomeQual on the given
# project (a "ForDistribution" subdir in DATA must exist). This script
# can be run on older projects, or on projects for which the
# ForDistribution code did not generate the agp.chromosome.fasta/qual
# files.
#
# DRAFT: the base name of the output dir of the module ForDistribution
##
set PRE    = "$ARACHNE_PRE"
set DATA   = "projects/Armadillo"
set RUN    = "run/work"
set SUB    = "assisted_2.18"
set DRAFT  = "Draft_v1"

##
# Derived path names.
set pdr      = "PRE=$PRE DATA=$DATA RUN=$RUN"
set fulldata = "$PRE/$DATA"
set fullrun  = "$PRE/$DATA/$RUN"
set fullsub  = "$fullrun/$SUB"

##
# Generate assembly.agp
if ( ! -e $fullsub/assembly.agp ) then
    SupersToAgp $pdr \
	SUBDIR=$SUB \
	OUTDIR=$DRAFT

    cp $fulldata/$DRAFT/assembly.agp $fullsub/assembly.agp
endif

##
# Create agp.chromosome.fasta
if ( ! -e $fulldata/$DRAFT/Draft_v1.agp.chromosome.fasta.gz ) then
    CreateChromosomeFasta $pdr \
	SUBDIR=$SUB \
	OUTDIR=$SUB \
	ASSEMBLY_NAME=$DRAFT \
	AGP_FILE=assembly.agp

    set sub_agp = $fullsub/Draft_v1.agp.chromosome.fasta.gz
    set draft_agp = $fulldata/$DRAFT/Draft_v1.agp.chromosome.fasta.gz
    mv $sub_agp $draft_agp
endif

##
# Create agp.chromosome.qual
if ( ! -e $fulldata/$DRAFT/Draft_v1.agp.chromosome.qual.gz ) then
    CreateChromosomeQual $pdr \
	SUBDIR=$SUB \
	OUTDIR=$SUB \
	ASSEMBLY_NAME=$DRAFT \
	AGP_FILE=assembly.agp

    set sub_agp = $fullsub/Draft_v1.agp.chromosome.qual.gz
    set draft_agp = $fulldata/$DRAFT/Draft_v1.agp.chromosome.qual.gz
    mv $sub_agp $draft_agp
endif
