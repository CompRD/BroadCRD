///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * Dedupper.cc
 *
 *  Created on: Nov 15, 2013
 *      Author: neilw
 */
#include "Basevector.h"
#include "IteratorRange.h"
#include "PairsManager.h"
#include "Qualvector.h"
#include "dna/Bases.h"
#include "feudal/BinaryStream.h"
#include "feudal/VirtualMasterVec.h"
#include "feudal/HashSet.h"
#include "kmers/ReadPatherDefs.h"
#include "paths/HyperBasevector.h"
#include "paths/UnibaseUtils.h"
#include "system/SortInPlace.h"
#include "system/SpinLockedData.h"
#include "system/WorklistN.h"
#include <algorithm>
#include <atomic>
#include <fstream>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>
#include "MainTools.h"

namespace
{

unsigned const K_2 = 19;
unsigned const K = K_2 * 2;
typedef KMer<K> DoppelKmer;
typedef typename DoppelKmer::Hasher Hasher;
typedef std::equal_to<DoppelKmer> Comparator;
typedef HashSet<DoppelKmer,Hasher,Comparator> KmerSet;

class GoodLenFinder
{
public:
    GoodLenFinder( vecqualvector const& vqv, unsigned minQual,
                        std::vector<unsigned>* pGoodLens )
    : mVqv(vqv), mMinQual(minQual), mGoodLens(*pGoodLens) {}

    void operator()( size_t readId )
    {
      qualvector const& qv = mVqv[readId];
      unsigned minQual = mMinQual;
      mGoodLens[readId] =
        std::find_if(qv.begin(),qv.end(),
              [minQual](unsigned char qqq){return qqq<minQual;})-qv.begin(); }

private:
    vecqualvector const& mVqv;
    unsigned mMinQual;
    std::vector<unsigned>& mGoodLens;
};



} // end of anonymous namespace



void Dedupper( vecbasevector& bases, vecqualvector& quals, PairsManager& pairs,
        vec<unsigned>& mapq, const bool dedup,
        const bool trim_bases, const unsigned trim_qual,
        unsigned min_mapq, const bool keep_zero )
{
    size_t nvecs = bases.size();
    size_t ndups = 0;
    size_t badmapq = 0;

    ForceAssertEq( nvecs, quals.size() );
    ForceAssertEq( nvecs, pairs.nReads() );
    if (min_mapq) ForceAssertEq( nvecs, mapq.size() );

    // establish a trim size for the reads
    cout << Date() << ": processing quals." << std::endl;
    std::vector<unsigned> goodLens(nvecs);
    if ( trim_bases ) {
        size_t nThreads = (7*getConfiguredNumThreads()+7)/8;
        parallelForBatch(0ul,nvecs,100000,
                            GoodLenFinder(quals,trim_qual,&goodLens),
                            nThreads);
    } else {
        for (size_t i = 0; i < bases.size(); ++i )
            goodLens[i] = bases[i].size();
    }

    // find duplicates
    vec<Bool> deleteme( nvecs, false );
    float dupRate = 0.01;
    size_t dictSize = (1.0-dupRate)*bases.size()*K;
    KmerSet map( dictSize );

    // find low mapq
    if ( min_mapq ) {
        for ( size_t i = 0; i < bases.size(); ++i ) {
            if ( mapq[i] < min_mapq && ( mapq[i] > 0 || !keep_zero ) ) {
                if ( pairs.isPaired( i ) ) {
                    size_t id_other = pairs.getPartnerID( i );
                    if ( !deleteme[id_other] ) {
                        deleteme[id_other] = true;
                        badmapq++;
                    }
                }
                if ( ! deleteme[i] ) {
                    deleteme[i] = true;
                    badmapq++;
                }
            }
        }
    }

    // for each pair
    for ( size_t i = 0; dedup && i < pairs.nPairs(); ++i ) {

        // build kmer first
        size_t id1 = pairs.ID1(i);
        size_t id2 = pairs.ID2(i);

        if ( deleteme[id1] || deleteme[id2] )
            continue;

        if ( goodLens[id1] < K || goodLens[id2] < K )
            continue;

        DoppelKmer this_kmer;

        const basevector& b1 = bases[id1];
        const basevector& b2 = bases[id2];

        // the ordering here is purposeful such that the RC of the
        // fragment from which the kmers come would yield the same
        // canonical DoppelKmer
        for ( size_t i = 0; i < K_2; ++i )
            this_kmer.set(i, b1[i]);
        for ( size_t i = 0; i < K_2; ++i )
            this_kmer.set(K_2+i, b2[K_2-1-i]);

        this_kmer.canonicalize();

        // try to add to the dictionary
        if ( !map.add( this_kmer ) ) {
            deleteme[id1]=deleteme[id2]=true;
            ndups+=2;
        }
    }

    cout << Date() << ": found " << ndups << " duplicates of " << nvecs << endl;

    // trim bases
    if ( trim_bases ) {
        cout << Date() << ": trimming bases" << endl;
        for ( size_t i = 0; i < bases.size(); ++i ) {
            bases[i].resize( goodLens[i] );
        }
    }
    cout << Date() << ": average read length " << bases.SizeSum()/bases.size() << endl;

    // Delete pairs, quals, and bases
    cout << Date() << ": pairs before duplicate/lowq removal " << pairs.nPairs() << endl;
    cout << Date() << ": reads in pairs manager before duplicate/lowq removal " << pairs.nReads() << endl;
    pairs.removeReads( deleteme );
    cout << Date() << ": pairs after duplicate/lowq removal " << pairs.nPairs() << endl;
    cout << Date() << ": reads in pairs manager after duplicate/lowq removal " << pairs.nReads() << endl;
    bases.EraseIf(deleteme);
    quals.EraseIf(deleteme);
    cout << Date() << ": reads in bases after duplicate removal (will be the same) " << bases.size() << endl;
    if ( min_mapq ) {
        EraseIf(mapq,deleteme);
        cout << Date() << ": removed because of bad mapping quality " << badmapq << endl;
    }

    ForceAssertEq(bases.size(), quals.size());
    ForceAssertEq(bases.size(), pairs.nReads() );
}


int main( int argc, char** argv )
{
    String empty;
    RunTime();
    BeginCommandArguments;
    CommandArgument_String_Doc(IN_HEAD, "looks for IN_HEAD.{fastb,qualb,pairs}, "
          "in ParseStringSet format");
    CommandArgument_Bool_OrDefault_Doc(TRIM_BASES,False,"trim read after first below min quality base");
    CommandArgument_Bool_OrDefault_Doc(DEDUP,True,"whether to deduplicate");
    CommandArgument_UnsignedInt_OrDefault_Doc(TRIM_QUAL,20,"minimum acceptible quality score for trimming");
    CommandArgument_String_Doc(OUT_HEAD, "writes a single OUT_HEAD.{fastb.qualb,pairs}");
    CommandArgument_UnsignedInt_OrDefault(NUM_THREADS,0u);
    CommandArgument_UnsignedInt_OrDefault_Doc(MIN_MAPQ,0u, "exclude reads where either of the pair has lower mapq than specified (via .mqpq file)");
    CommandArgument_Bool_OrDefault_Doc(MIN_MAPQ_KEEP_ZERO, False, "whether to keep mapping quality zero reads, regardless of threshold");
    CommandArgument_LongLong_OrDefault(MAX_MEMORY_GB,0);
    EndCommandArguments;

    NUM_THREADS = configNumThreads(NUM_THREADS);
    SetMaxMemory(MAX_MEMORY_GB<<30);

    vecbasevector bases;
    vecqualvector quals;
    PairsManager pairs;
    vec<unsigned> mapq;

    cout << Date( ) << ": loading reads" << endl;
    vec<String> heads;
    ParseStringSet( IN_HEAD, heads );
    cout << Date( ) << ": found " << heads.size( ) << " heads" << endl;
    for ( int i = 0; i < heads.isize( ); i++ ) {
         bases.ReadAll( heads[i] + ".fastb", True );
         quals.ReadAll( heads[i] + ".qualb", True );
         pairs.Append( PairsManager( heads[i] + ".pairs" ) );
         if ( MIN_MAPQ ) {
             Ifstream( mapin, heads[i] + ".mapq" );
             while ( mapin.good() ) {
                 unsigned int tmp;
                 if ( mapin >> tmp ) {
                     mapq.push_back(tmp);
                 }
             }
             ForceAssertEq(mapq.size(), bases.size());
         }
    }

    cout << Date( ) << ": total reads = " << ToStringAddCommas( bases.size( ) )
         << endl;

    Dedupper( bases, quals, pairs, mapq, DEDUP, TRIM_BASES, TRIM_QUAL, MIN_MAPQ, MIN_MAPQ_KEEP_ZERO );

    cout << Date() << ": writing output " << OUT_HEAD << ".{fastb,qualb,pairs}" << endl;
    bases.WriteAll( OUT_HEAD + ".fastb" );
    quals.WriteAll( OUT_HEAD + ".qualb" );
    pairs.Write( OUT_HEAD + ".pairs" );

    if ( MIN_MAPQ ) {
        Ofstream( mapoutfile, OUT_HEAD + ".mapq" );
        ostream_iterator<unsigned> mapout( mapoutfile, "\n" );
        std::copy( mapq.begin(), mapq.end(), mapout );
    }

    cout << Date() << ": Done." << std::endl;
}
