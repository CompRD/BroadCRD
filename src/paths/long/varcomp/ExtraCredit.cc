///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * ExtraCredit.cc
 *
 *  Created on: Aug 27, 2013
 *      Author: tsharpe
 *
 *  A program to scan a file having the variants.all format (you can just cat
 *   together the variants.all files from each of the fosmids into a single big
 *   file, if you like).
 *  It spits out the regions where it is plausible that some caller is being
 *   robbed of credit by making emitting a set of calls different from, but
 *   equivalent to, our canonical set of fosmid events.
 *  It doesn't attempt to figure out anything about whether the calls are
 *   equivalent or not -- it just spits out the regions where this might be so.
 */
#include "MainTools.h"
#include "Basevector.h"
#include "TokenizeString.h"
#include "Vec.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include "paths/long/fosmid/Fosmids.h"

class Call
{
public:
    Call() = default;
    Call( size_t fosmidId, String const& context, size_t fosmidOffset,
            String const& chromosome, size_t chromoOffset,
            String const& refSeq, String const& altSeq, String const& info,
            int recognizedBy, String const& line )
    : mFosmidId(fosmidId), mContext(context), mFosmidOffset(fosmidOffset),
      mChromosome(chromosome), mChromoOffset(chromoOffset),
      mRefSeq(refSeq), mAltSeq(altSeq), mInfo(info),
      mRecognizedBy(recognizedBy), mLine(line)
    {}

    size_t getFosmidId() const { return mFosmidId; }
    size_t getChromoOffset() const { return mChromoOffset; }
    bool isFosmid() const { return mRecognizedBy & 1; }
    bool isIndel() const { return mRefSeq.size() != mAltSeq.size(); }

    // is fosmid indel, not universally recognized by all callers
    bool isFINUR() const
    { return isFosmid() && isIndel() && mRecognizedBy != 31; }

    // test for no common callers
    friend bool isDisjoint( Call const& c1, Call const& c2 )
    { return !(c1.mRecognizedBy & c2.mRecognizedBy); }

    friend bool operator<( Call const& c1, Call const& c2 )
    { if ( c1.mFosmidId < c2.mFosmidId ) return true;
      if ( c1.mFosmidId > c2.mFosmidId ) return false;
      if ( c1.mChromoOffset < c2.mChromoOffset ) return true;
      if ( c1.mChromoOffset > c2.mChromoOffset ) return false;
      if ( c1.isFosmid() && !c2.isFosmid() ) return true;
      return false; }

    friend ostream& operator<<( ostream& os, Call const& call )
    { return os << call.mLine; }

private:
    size_t mFosmidId;
    String mContext;
    size_t mFosmidOffset;
    String mChromosome;
    size_t mChromoOffset;
    String mRefSeq;
    String mAltSeq;
    String mInfo;
    int mRecognizedBy;
    String mLine;
};

int main( int argc, char** argv )
{
    RunTime();
    BeginCommandArguments;
//    CommandArgument_String_OrDefault(REFS, "refs.fastb");
//    CommandArgument_String_OrDefault(FOSMIDS, "fosmids.fastb");
    CommandArgument_String_OrDefault(CALLS, "variants.all");
    CommandArgument_UnsignedInt_OrDefault(REGION,100u);
    EndCommandArguments;

//    vecbvec refs(REFS);
//    vecbvec fosmids(FOSMIDS);
    std::ifstream calls(CALLS.c_str());
    String line;
    vec<String> fields;
    vec<Call> callVec;
    callVec.reserve(10000);
    char const* callers[]{"Fosmid","DISCOVAR","GATK-100","GATK-250","CORTEX"};
    while ( getline(calls,line) )
    {
        Tokenize(line,fields);
        ForceAssertEq(fields.size(),7ul);
        size_t fosmidId = fields[0].substr(1).Int();
        String const& context = fields[1];
        size_t fosmidOffset = fields[2].Int();
        String chromosome = fields[3].Before(":");
        size_t chromoOffset = fields[3].After(":").Int();
        String const& refSeq = fields[4];
        String const& altSeq = fields[5];
        String const& info = fields[6];
        int recognizedBy = 0;
        int mask = 1;
        for ( char const* caller : callers )
        {
            if ( info.Contains(caller) )
                recognizedBy |= mask;
            mask <<= 1;
        }
        callVec.emplace_back(fosmidId,context,fosmidOffset,chromosome,
                              chromoOffset,refSeq,altSeq,info,recognizedBy,
                              line);
    }
    std::sort(callVec.begin(),callVec.end());

    int g, rstart, rstop;
    String gid, loc;

    auto itr = callVec.begin();
    auto end = callVec.end();
    vec<Call const*> finurs;
    while ( itr != end )
    {
        auto itr2 = itr;
        size_t prevOffset = itr->getChromoOffset();
        if ( itr->isFINUR() )
            finurs.push_back(&*itr);
        while ( ++itr2 != end )
        {
            if ( itr2->getFosmidId() != itr->getFosmidId() )
                break;
            if ( itr2->getChromoOffset()-prevOffset > REGION )
                break;
            prevOffset = itr2->getChromoOffset();
            if ( itr2->isFINUR() )
                finurs.push_back(&*itr2);
        }

        bool goodGroup = false;
        for ( Call const* finur : finurs )
        {
            for ( auto itr1=itr; !goodGroup && itr1 != itr2; ++itr1 )
            {
                if ( isDisjoint(*finur,*itr1) )
                    goodGroup = true;
            }
            if ( goodGroup )
                break;
        }
        finurs.clear();

        if ( !goodGroup )
            itr = itr2;
        else
        {
            while ( itr != itr2 )
            {
                GetRegionInfo( ToString(itr->getFosmidId()), g, rstart, rstop, gid, loc );
                std::cout << *itr << " ("  << (itr->getChromoOffset()-rstart) << ")\n";
                ++itr;
            }
            std::cout << '\n' << std::endl;
        }
    }
}
