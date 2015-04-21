///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * FindHomopolymers.h
 *
 *  Created on: Mar 27, 2013
 *      Author: tsharpe
 */

#ifndef FINDHOMOPOLYMERS_H_
#define FINDHOMOPOLYMERS_H_

#include "IteratorRange.h"
#include "Vec.h"
#include "dna/Bases.h"
#include "system/SysConf.h"
#include "feudal/BinaryStream.h"
#include "feudal/MasterVec.h"
#include "feudal/SerfVec.h"
#include "paths/long/SupportedHyperBasevector.h"
#include <algorithm>
#include <iostream>

struct HomopolymerRun
{
    HomopolymerRun()=default;

    HomopolymerRun( vec<int> const& edgeIds, unsigned initialOffset,
                      unsigned finalOffset, unsigned totalLen,
                      bool isCyclic, unsigned char base )
    : mEdgeIds(edgeIds), mInitialOffset(initialOffset),
      mFinalOffset(finalOffset), mTotalLen(totalLen), mIsCyclic(isCyclic),
      mBase(base)
    {}

    HomopolymerRun( int edgeId, unsigned initialOffset, unsigned finalOffset,
                      unsigned totalLen, bool isCyclic, unsigned char base )
    : mInitialOffset(initialOffset),
      mFinalOffset(finalOffset), mTotalLen(totalLen),
      mIsCyclic(isCyclic), mBase(base)
    { mEdgeIds.reserve(1); mEdgeIds.push_back(edgeId); }

    bool isSolo() const { return mEdgeIds.size()==1; }

    friend bool overlaps( HomopolymerRun const& hr1, HomopolymerRun const& hr2 )
    { vec<int> const& e1 = hr1.mEdgeIds;
      vec<int> const& e2 = hr2.mEdgeIds;
      return (e1.front() == e2.front() &&
                  hr1.mInitialOffset == hr2.mInitialOffset) ||
             (e1.back() == e2.back() &&
                  hr1.mFinalOffset == hr2.mFinalOffset) ||
             (e1.size() > 2 && e2.size() > 2 &&
                  std::find_first_of(e1.begin()+1,e1.end()-1,
                                     e2.begin()+1,e2.end()-1) != e1.end()-1); }

    friend bool operator<( HomopolymerRun const& hr1,HomopolymerRun const& hr2 )
    { if ( hr1.mEdgeIds < hr2.mEdgeIds ) return true;
      if ( hr2.mEdgeIds < hr1.mEdgeIds ) return false;
      if ( hr1.mInitialOffset < hr2.mInitialOffset ) return true;
      return false; }

    friend ostream& operator<<( ostream& os, HomopolymerRun const& hr )
    { os << 'p' << Base::val2Char(hr.mBase) << '[' << hr.mTotalLen;
      if ( hr.mIsCyclic ) os << '+';
      os << "] at " << hr.mEdgeIds.front() << '[' << hr.mInitialOffset << ']';
      if ( hr.mEdgeIds.size() > 1 )
        os << ',' << rangePrinter(hr.mEdgeIds.begin()+1,hr.mEdgeIds.end(),",");
      return os << '[' << hr.mFinalOffset << ']'; }

    vec<int> mEdgeIds; // the hbv edge indices for the path
    unsigned mInitialOffset; // the number of bases to skip at the beginning of the first path element
    unsigned mFinalOffset; // the number of bases to skip at the end of the last path element
    unsigned mTotalLen; // the total length of the homopolymer run
    bool mIsCyclic; // if a cycle was detected
    unsigned char mBase; // this is a homopolymer run for which base code
};
TRIVIALLY_SERIALIZABLE(HomopolymerRun);
typedef SerfVec<HomopolymerRun> VecHomopolymerRun;
typedef MasterVec<VecHomopolymerRun> VecHomopolymerRunVec;

class HomopolymerFinder
{
public:
    HomopolymerFinder( unsigned minRunLen = 20 )
    : mMinRunLen(minRunLen), mpHBV(0), mpOut(0)
    {}

    void find( VecHomopolymerRun* pOut, HyperBasevector const& hbv )
    { mpHBV = &hbv;
      mpOut = pOut;
      int nVs = hbv.N();
      for ( int idx = 0; idx != nVs; ++idx )
        explore(idx); }

    static void cluster( VecHomopolymerRunVec* pOut,
                            VecHomopolymerRun const& runs,
                            unsigned minClusterSize = 2 );

    /// checks to see that all paths from v1 contain a single base.  checking
    /// stops when v2 is encountered.
    static bool isHomopolymerCell( HyperBasevector const& hbv, int v1, int v2 );

private:
    void explore( int vvv );
    bool extend( int vvv, unsigned char base, unsigned off, unsigned curLen,
                    bool cyclic, vec<int>& path );

    // if any of the predecessor edges ends with the same base that this one
    // begins with, then it would be redundant to emit a match for the first run
    // of bases, even if that run is long enough -- we'll discover those runs
    // during extension.  so bump the starting point for the exploration past
    // the first run of identical bases if we're in that situation.
    bvec::const_iterator findStartPos( int vvv, bvec const& edge )
    { auto itr = edge.begin();
      unsigned char firstBase = *itr;
      int nVs = mpHBV->To(vvv).size();
      for ( int iii = 0; iii != nVs; ++iii )
        if ( mpHBV->EdgeObjectByIndexTo(vvv,iii).back() == firstBase )
          return std::find_if(itr,edge.end(),getNEPred(firstBase));
      return itr; }

    typedef std::binder1st<std::not_equal_to<unsigned char>> Pred_t;

    static Pred_t getNEPred( unsigned char base )
    { return std::bind1st(std::not_equal_to<unsigned char>(),base); }

    unsigned mMinRunLen;
    HyperBasevector const* mpHBV;
    VecHomopolymerRun* mpOut;
};

#endif /* FINDHOMOPOLYMERS_H_ */
