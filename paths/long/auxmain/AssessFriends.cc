///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file AssessFriends.cc
 * \author tsharpe
 * \date Apr 20, 2012
 *
 * \brief
 */
#include "MainTools.h"
#include "Intvector.h"
#include "paths/long/Friends.h"
#include <algorithm>
#include <vector>

namespace
{

struct GenomicLocus
{
    size_t refID;
    unsigned offset;
    unsigned length;

    friend bool operator<( GenomicLocus const& loc1,
                                GenomicLocus const& loc2 )
    { if ( loc1.refID < loc2.refID ) return true;
      if ( loc1.refID == loc2.refID && loc1.offset < loc2.offset )
          return true;
      return false; }

    friend int order( GenomicLocus const& loc1, GenomicLocus const& loc2 )
    { if ( loc1.refID < loc2.refID ) return -1;
      if ( loc1.refID > loc2.refID ) return 1;
      if ( loc1.offset + loc1.length <= loc2.offset ) return -1;
      if ( loc1.offset >= loc2.offset + loc2.length ) return 1;
      return 0; }

    friend unsigned overlap( GenomicLocus const& loc1,
                                    GenomicLocus const& loc2 )
    { return std::min(loc1.offset+loc1.length,loc2.offset+loc2.length)
                - std::max(loc1.offset,loc2.offset); }
};
typedef std::vector<GenomicLocus> VecGenomicLocus;

struct ReadLocus
{
    GenomicLocus locus;
    size_t readID;

    friend bool operator<( ReadLocus const& loc1, ReadLocus const& loc2 )
    { return loc1.locus < loc2.locus; }

    friend int order( ReadLocus const& loc1, ReadLocus const& loc2 )
    { return order(loc1.locus,loc2.locus); }

    friend unsigned overlap( ReadLocus const& loc1, ReadLocus const& loc2 )
    { return overlap(loc1.locus,loc2.locus); }
};
typedef std::vector<ReadLocus> VecReadLocus;

void computeTrueFriends( String const& locsFile, unsigned minOverlap,
                                VecULongVec* pVV )
{
    VecReadLocus locs;
    ReadLocus loc;
    loc.readID = 0;
    BinaryIteratingReader<VecGenomicLocus> rdr(locsFile);
    locs.reserve(rdr.remaining());
    while ( rdr.next(&loc.locus) )
    {
        locs.push_back(loc);
        loc.readID += 1;
    }
    std::sort(locs.begin(),locs.end());
    pVV->resize(locs.size());
    ULongVec friends;
    typedef VecReadLocus::iterator Itr;
    Itr start(locs.begin());
    for ( Itr itr(start), end(locs.end()); itr != end; ++itr )
    {
        while ( order(*start,*itr) < 0 )
            ++start;
        for ( Itr itr2(start); itr2 != end; ++itr2 )
        {
            int ord = order(*itr2,*itr);
            if ( ord > 0 )
                break;
            if ( !ord && itr2 != itr && overlap(*itr2,*itr) >= minOverlap )
                friends.push_back(itr2->readID);
        }
        std::sort(friends.begin(),friends.end());
        pVV[0][itr->readID] = friends;
        friends.clear();
    }
}

void assessAccuracy( VecULongVec const& trueFriends,
                        IAndOsVec& testFriends,
                        bool verbose )
{
    size_t nPlusTot = 0;
    size_t nMinusTot = 0;
    size_t nEqTot = 0;
    size_t nnn = trueFriends.size();
    for ( size_t idx = 0; idx != nnn; ++idx )
    {
        size_t nPlus = 0;
        size_t nMinus = 0;
        size_t nEq = 0;
        ULongVec const& trueV = trueFriends[idx];
        IAndOs& testV = testFriends[idx];
	std::sort(testV.begin(),testV.end());
        typedef ULongVec::const_iterator Itr;
        if ( verbose )
            std::cout << idx << ':';
        Itr itr1(trueV.begin());
        Itr end1(trueV.end());
        typedef IAndOs::const_iterator FItr;
        FItr itr2(testV.begin());
        FItr end2(testV.end());
        while ( itr1 != end1 || itr2 != end2 )
        {
            if ( itr1 == end1 )
            {
                do
                {
                    if ( verbose )
                        std::cout << " +" << itr2->getId();
                    nPlus += 1;
                }
                while ( ++itr2 != end2 );
            }
            else if ( itr2 == end2 )
            {
                do
                {
                    if ( verbose )
                        std::cout << " -" << *itr1;
                    nMinus += 1;
                }
                while ( ++itr1 != end1 );
            }
            else
            {
                if ( *itr1 < itr2->getId() )
                {
                    if ( verbose )
                        std::cout << " -" << *itr1;
                    nMinus += 1;
                    ++itr1;
                }
                else if ( *itr1 > itr2->getId() )
                {
                    if ( verbose )
                        std::cout << " +" << itr2->getId();
                    nPlus += 1;
                    ++itr2;
                }
                else
                {
                    nEq += 1;
                    ++itr1;
                    ++itr2;
                }
            }
        }

        if ( verbose )
            std::cout << "\n+" << nPlus << " -" << nMinus << " =" << nEq
                        << std::endl;

        nPlusTot += nPlus;
        nMinusTot += nMinus;
        nEqTot += nEq;
    }
    std::cout << "There were " << nPlusTot << " false friends." << std::endl;
    std::cout << "There were " << nMinusTot << " absent friends." << std::endl;
    std::cout << "There were " << nEqTot << " fine friends." << std::endl;
}

}
TRIVIALLY_SERIALIZABLE(GenomicLocus);

int main( int argc, char** argv )
{
    String const EMPTY;
    RunTime();
    BeginCommandArguments;
    CommandArgument_String_Doc(TRUTH, "Genomic loci of reads.");
    CommandArgument_String_Doc(TEST, "VecULongVec of friends.");
    CommandArgument_UnsignedInt_OrDefault_Doc(MIN_OVERLAP,50u,
            "Minimum interval overlap for friendship.");
    CommandArgument_Bool_OrDefault_Doc(VERBOSE,False,
            "Print per-read details about absent and false friends.")
    EndCommandArguments;

    VecULongVec trueFriends;
    computeTrueFriends(TRUTH,MIN_OVERLAP,&trueFriends);
    IAndOsVec testFriends;
    testFriends.ReadAll(TEST);
    assessAccuracy(trueFriends,testFriends,VERBOSE);
}

