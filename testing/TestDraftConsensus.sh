#! /usr/bin/tcsh -f


##
# TestDraftConsensus
#
# Given a SUBDIR in input it will generate consensus by using two methods:
#  1) MergeReadsConsensus, and
#  2) DraftConsensus/FixConsensus;
# and will compare the two by calling MarkBy and CmpSeq.
##


# Set (with example)
set DATA       = $1   # projects/Dog
set RUN        = $2   # run/work
set SUBDIR     = $3   # redog/wuffi19.b.closed.patched.patched
set OUTDIR     = $4   # redog/CompareConsensusTest
set CONTIG_ID  = $5   # 455

set pre_dir    = $ARACHNE_PRE/$DATA/$RUN
set sub_dir    = $pre_dir/$SUBDIR
set out_dir    = $pre_dir/$OUTDIR


# Clean up old runs.
rm -rf $sub_dir/Tilings
rm -rf $out_dir
mkdir -p $out_dir


# Generate consensus with MergeReads (MergeReads Consensus).
MergeReadsConsensus \
    DATA=$DATA \
    RUN=$RUN \
    SUBDIR=$SUBDIR \
    OUTDIR=$OUTDIR/MergeReadsConsensus \
    CONTIG_IDS=$CONTIG_ID


# Generate a draft consensus (Draft Consensus, I).
DraftConsensus \
    DATA=$DATA \
    RUN=$RUN \
    SUBDIR=$SUBDIR \
    OUTDIR=$OUTDIR/DraftConsensus \
    CONTIG_IDS=$CONTIG_ID


# Refine consensus (Draft Consensus, II).
FixConsensus \
    DATA=$DATA \
    RUN=$RUN \
    SUBDIR=$OUTDIR/DraftConsensus \
    TILINGS=$OUTDIR/DraftConsensus/draft_tiles \
    EXPERIMENTAL=True


# Generate locs indices (needed by MarkBy).
AddLocationIndices \
    DATA=$DATA \
    RUN=$RUN \
    SUBDIR=$OUTDIR/DraftConsensus \

AddLocationIndices \
    DATA=$DATA \
    RUN=$RUN \
    SUBDIR=$OUTDIR/MergeReadsConsensus \


# MarkBy the two runs.
MarkBy \
    DATA=$DATA \
    RUN=$RUN \
    SUB1=$OUTDIR/MergeReadsConsensus \
    SUB2=$OUTDIR/DraftConsensus \
    S1=0


# CmpSeq the two.
set COUNT1 = `head -1 $out_dir/MergeReadsConsensus/correspondence.log`
set COUNT2 = `head -1 $out_dir/DraftConsensus/correspondence.log`

CmpSeq \
    PRE=$out_dir \
    HEAD1=MergeReadsConsensus/mergedcontigs \
    HEAD2=DraftConsensus/mergedcontigs \
    START1=0 \
    START2=0 \
    COUNT1=$COUNT1 \
    COUNT2=$COUNT2 \
    K=24 \
    MAXCLIQ=1000 \
    MIN_OVERLAP=50 \
    BLOCK_SIZE=500 \
    ACTUAL_PASSES=100 \
    MAKEALIGNS_MAXERRS=200 \
    SW_GAP_METHOD=True \
    SW_GAP_MAX_GAP=10000 \
    SW_GAP_MAX_OFFSET_DIFF=10000 \
    SW_GAP_MIN_OVERLAP_FRAC=0.6 \
    MAKEALIGNS_END_STRETCH=40 \
    ABBREVIATE_ALIGNMENTS=True
