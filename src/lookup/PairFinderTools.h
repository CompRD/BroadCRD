/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef _PAIR_FINDER_TOOLS_H
#define _PAIR_FINDER_TOOLS_H

#include "CoreTools.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"

typedef enum { INNER, OUTER, LEFT, RIGHT, UNKNOWN } pair_direction;
typedef enum { INTRA, INTER } pair_placement;

/// Alignment pair <la1>, <la2> (assumed to be paired end reads from the same fragment)
/// is "reasonable" if the alignments are on the same contig and form an INNER pair ( --1--> <--2-- or --2--> <--1--).
bool IsReasonablePair(const look_align &la1, const look_align &la2) {
    return (la1.target_id == la2.target_id && la1.rc1 != la2.rc1 && 
	    ( (la1.rc1 == 0 && la1.StartOnTarget() < la2.StartOnTarget() ) || 
	      (la2.rc1 == 0 && la2.StartOnTarget() < la1.StartOnTarget() )
            )
	   );
}

/// Returns INTRA if the two aligns are on the same contig, or INTER otherwise.
inline pair_placement PairPlacement(const look_align &la1, const look_align &la2) {
    return (la1.target_id == la2.target_id) ? INTRA : INTER;
}

/// Returns the relative direction of the pair of aligns (assumed to be aligns for the paired end reads 
/// from the same fragment):
///   *  --1--> <--2-- or  --2--> <--1-- : INTER
///   *  <--1-- --2--> or <--2--  --1--> : OUTER
///   *  --1--> --2--> or  --2--> --1--> : RIGHT
///   *  <--1-- <--2-- or <--2--  <--1-- : LEFT
inline pair_direction PairDirection(const look_align &la1, const look_align &la2) {
    if (la1.rc1 == 0 && la2.rc1 == 1 && la1.StartOnTarget() < la2.StartOnTarget()) { return INNER; } //  --1--> <--2--
    if (la2.rc1 == 0 && la1.rc1 == 1 && la2.StartOnTarget() < la1.StartOnTarget()) { return INNER; } //  --2--> <--1--

    if (la1.rc1 == 1 && la2.rc1 == 0 && la1.StartOnTarget() < la2.StartOnTarget()) { return OUTER; } // <--1--   --2-->
    if (la2.rc1 == 1 && la1.rc1 == 0 && la2.StartOnTarget() < la1.StartOnTarget()) { return OUTER; } // <--2--   --1-->

    if (la1.rc1 == 0 && la2.rc1 == 0) { return RIGHT; } //  --1--> --2--> or --2--> --1-->
    if (la1.rc1 == 1 && la2.rc1 == 1) { return LEFT; }  // <--1-- <--2-- or <--2-- <--1--

    return UNKNOWN;
}

/// Assuming that the two alignments are for paired end reads from the same fragment, computes
/// the total length of the fragment (including reads, *not* the insert size!). In other words,
/// the returned fragment length is the distance on th etarget contig from the start of the 
/// leftmost alignment to the end of the rightmost alignment (regardless of the alignment orientations).
/// NOTE: this method does not check that the alignments are on the same contig!!
unsigned int FragmentSize(const look_align &la1, const look_align &la2) {
    return abs((la1.StartOnTarget() < la2.StartOnTarget()) ? (la1.StartOnTarget() - la2.EndOnTarget()) : (la2.StartOnTarget() - la1.EndOnTarget()));
}

vec<float> LoadFragmentDistribution(String distributionfile) {

    // Get maximum probability.

    float maxprob = 0.0;
    {    ifstream fdstream(distributionfile.c_str());
         String fdline;
         while(getline(fdstream, fdline)) {
             if (!fdline.Contains("#")) {
                 float prob = fdline.After(" ").Double();
                 maxprob = Max( maxprob, prob );    }    }    }

    ifstream fdstream(distributionfile.c_str());
    String fdline;
    vec<float> fd(1000);
    unsigned int maxgapsize = 0;
    const float minprob_mult = 0.001;

    while(getline(fdstream, fdline)) {
        if (!fdline.Contains("#")) {
            unsigned int gapsize = fdline.Before(" ").Int();
            float prob = fdline.After(" ").Double();

            // Remove outliers.

            if ( prob < minprob_mult * maxprob ) continue;

            if (gapsize >= fd.size()) {
                fd.resize(static_cast<int>(1.3*static_cast<float>(gapsize)));
            }

            fd[gapsize] = prob;

            if (gapsize > maxgapsize) {
                maxgapsize = gapsize;
            }
        }
    }

    float fdsum = Sum(fd);
    for ( int i = 0; i < fd.isize( ); i++ )
         fd[i] /= fdsum;

    fdstream.close();

    fd.resize(maxgapsize);

    return fd;
}


#endif
