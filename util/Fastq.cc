///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file Fastq.cc
 * \author blau
 * \date Mar 5, 2014
 *
 * \brief Converts a Fastq file to a fastb and qualb.
 */
/*#include <cstdio>

#include "MainTools.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "String.h"
#include "system/ProcBuf.h"
#include <fstream>
#include "util/Fastq.h"
*/

#include <istream>
#include "system/System.h"
#include "String.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "PairsManager.h"
#include "util/Fastq.h"

namespace fastq{


// factored out from FastqToFastbQualb.cc at r48873
bool NextFastqEntry( std::istream& is, FastqEntry &entry,
                            vec<size_t> * p_n_amb, const char phred_offset )
{
    bool result = false;
    std::string seqheader;
    std::string seq;
    std::string qualheader;
    std::string qual;
    if ( std::getline(is,seqheader) && seqheader.size() )
    {
        std::getline(is,seq);
        std::getline(is,qualheader);
        std::getline(is,qual);

        // ---- count number of ambiguous bases
        typedef std::string::iterator Itr;
        for ( Itr itr(seq.begin()), end(seq.end()); itr != end; ++itr )
            (*p_n_amb)[GeneralizedBase::fromChar(*itr).getAmbiguityCount()-1]++;

        entry.sequence = bvec(seq.begin(),seq.end(),GenCharToRandomBaseMapper());

        entry.qualities.clear().reserve(qual.size());
        Swizzler swiz(phred_offset);
        for ( Itr itr(qual.begin()), end(qual.end()); itr != end; ++itr )
            entry.qualities.push_back(swiz(*itr));
        entry.name = seqheader;
        entry.name = entry.name.SafeAfter( "@" );
        entry.phred_min = swiz.getMin();
        entry.phred_max = swiz.getMax();
        result = true;
    }
    return result;
}

void Fastq2FastbQualbPair_Interleave( const String& fastq1, const String& fastq2, const char q_offset
                                    , const String& outhead){
    IncrementalWriter<basevector> bases_out((outhead + ".fastb").c_str());
    IncrementalWriter<qualvector> quals_out((outhead + ".qualb").c_str());
    uint64_t nPairs=0;
    FastQReader r1(fastq1,q_offset),r2(fastq2,q_offset);
    for(;r1.State()&&r2.State(); r1.Step(),r2.Step()){
        if(r1.Get().sequence.size() != r1.Get().qualities.size()){ FatalErr("#base != #qual for read 1 of pair "+ToString(nPairs)+"."); }
        if(r2.Get().sequence.size() != r2.Get().qualities.size()){ FatalErr("#base != #qual for read 2 of pair "+ToString(nPairs)+"."); }
        bases_out.add(r1.Get().sequence);
        bases_out.add(r2.Get().sequence);

        quals_out.add(r1.Get().qualities);
        quals_out.add(r2.Get().qualities);
        ++nPairs;
    }
    if(r1.State() != r2.State()){ FatalErr("The two fastq files have different number of entries."); }

    PairsManager pairs(nPairs * 2);
    int lib_id = pairs.addLibrary(0, 0, fastq1+"_"+fastq2);
    for (size_t i = 0; i < nPairs; i++) pairs.addPairToLib( i * 2, i * 2 + 1, lib_id);
    pairs.Write(outhead + ".pairs");
}
void Fastq2FastbQualbPair( const String& fastq, const char q_offset , const String& outhead){
    IncrementalWriter<basevector> bases_out((outhead + ".fastb").c_str());
    IncrementalWriter<qualvector> quals_out((outhead + ".qualb").c_str());
    uint64_t nReads=0;
    for(FastQReader rr(fastq,q_offset) ; rr.State() ; rr.Step()){

        if(rr.Get().sequence.size() != rr.Get().qualities.size()){ FatalErr("#base != #qual for read "+ToString(nReads)+"."); }
        bases_out.add(rr.Get().sequence);
        quals_out.add(rr.Get().qualities);
        ++nReads;
    }
    if(nReads%2!=0) FatalErr("odd number of reads.");

    PairsManager pairs(nReads);
    int lib_id = pairs.addLibrary(0, 0, fastq);
    for (size_t i = 0; i < nReads; i+=2) pairs.addPairToLib( i * 2, i * 2 + 1, lib_id);
    pairs.Write(outhead + ".pairs");
}



size_t ReadFastq( String const& filename, vecbasevector& bases_out,
        vecqualvector& quals_out, char const phred_offset,
        vec<String>* names_outp) {

    bases_out.clear();
    quals_out.clear();

    size_t count = 0;
    for(FastQReader rr(filename,phred_offset) ; rr.State() ; rr.Step()){
        if(rr.Get().sequence.size() != rr.Get().qualities.size()){
            FatalErr("#base != #qual for read "+
                    ToString(count)+" from "+filename); }
        bases_out.push_back(rr.Get().sequence);
        quals_out.push_back(rr.Get().qualities);
        if ( names_outp ) names_outp->push_back(rr.Get().name);
        count++;
    }

    return count;
}


}
