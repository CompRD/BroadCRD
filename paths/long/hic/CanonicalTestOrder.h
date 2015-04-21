///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Aug 8, 2014 - <crdhelp@broadinstitute.org>
//

#ifndef CANONICALTESTORDER_H_
#define CANONICALTESTORDER_H_


// this is just to stuff some ugly, but related things into one place
//
// when testing pairs of lines A and B, we establish a canonical order
// for the TESTS, which is AB, BA, AB', B'A
//


class CanonicalTestOrder {
public:

// order below MUST MATCH testLinePair -- this
// needs to be a class
static String pairIndexToConfigString( size_t idx )
{
    if (idx==0) return "AB";
    else if ( idx == 1) return "BA";
    else if (idx == 2) return "AB'";
    else if ( idx==3) return "B'A";
    else FatalErr("weird index");
}

static pair<TaggedLineId,TaggedLineId> pairIndexToConfig( size_t idx, TaggedLineId const& A, TaggedLineId const& B)
{
    auto Bpr = B;
    Bpr.flipTag();
    if ( idx == 0 ) return make_pair(A,B);
    else if ( idx == 1 ) return make_pair(B,A);
    else if ( idx == 2 ) return make_pair(A,Bpr);
    else if ( idx == 3 ) return make_pair(Bpr,A);
    else FatalErr("weird index");
}

static pair<int, bool> pairIndexToPairConfig( size_t idx )
{
    if ( idx == 0 ) return make_pair(1, false); // AB
    else if ( idx == 1 ) return make_pair(-1, false); //BA
    else if ( idx == 2 ) return make_pair(1, true); // AB'
    else if ( idx == 3 ) return make_pair(-1, true); // B'A
    else FatalErr("weird index");
}

};




#endif /* CANONICALTESTORDER_H_ */
