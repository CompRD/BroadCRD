///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * Test.cc
 *
 *  Created on: Aug 15, 2013
 *      Author: tsharpe
 */

#include "MainTools.h"
#include "MapReduceEngine.h"
#include "Basevector.h"
#include "feudal/VirtualMasterVec.h"
#include "kmers/KMer.h"
#include <atomic>
#include <vector>

unsigned const K = 60u;
typedef KMer<K> Kmer;
typedef std::vector<std::atomic<size_t>> Spectrum;

class KSpecImpl
{
public:
    KSpecImpl( VirtualMasterVec<bvec> const& vbv, Spectrum* pOut )
    : mVBV(vbv), mpOut(pOut)
    { mBins.resize(mpOut->size(),0ul); }

    ~KSpecImpl()
    { auto oItr = mpOut->begin();
      for ( auto itr=mBins.begin(),end=mBins.end(); itr != end; ++itr,++oItr )
        oItr->fetch_add(*itr); }

    template <class OItr>
    void map( size_t idx, OItr oItr )
    { bvec const& bv = mVBV[idx];
      Kmer::kmerize(bv.begin(),bv.end(),oItr); }

    void reduce( Kmer const* beg, Kmer const* end )
    { size_t count = end-beg;
      if ( --count >= mBins.size() ) count = mBins.size()-1;
      mBins[count] += 1; }

    Kmer* overflow( Kmer* beg, Kmer* end )
    { return beg+std::min(size_t(end-beg),mBins.size()); }

private:
    VirtualMasterVec<bvec> mVBV;
    Spectrum* mpOut;
    std::vector<size_t> mBins;
};

int main( int argc, char** argv )
{
    RunTime();
    BeginCommandArguments;
    CommandArgument_String_Doc(READS, "The input reads.");
    CommandArgument_UnsignedInt_OrDefault_Doc(NBINS,1000u,"Number of spectrum bins.");
    CommandArgument_Int_OrDefault_Doc(NUM_THREADS,0,"Max number of threads.");
    CommandArgument_UnsignedInt_OrDefault(MAX_MEMORY_GB,0u);
    CommandArgument_Double_OrDefault(MEAN_USAGE,.9);
    EndCommandArguments;
    SetMaxMemory( (1ul*MAX_MEMORY_GB) << 30 );
    NUM_THREADS = configNumThreads(NUM_THREADS);

    typedef MapReduceEngine<KSpecImpl,Kmer,Kmer::Hasher> MRE;
    Spectrum spectrum(NBINS);
    VirtualMasterVec<bvec> vbv(READS);
    KSpecImpl impl(vbv,&spectrum);
    MRE mre(impl);
    if ( !mre.run(vbv.getKmerCount(K),0ul,vbv.size(),
                    MRE::VERBOSITY::NOISY,MEAN_USAGE) )
        FatalErr("Ran out of buffer space in MRE.");
    size_t bin = 0;
    for ( auto itr=spectrum.begin(),end=spectrum.end(); itr != end; ++itr )
    {
        ++bin;
        size_t count = *itr;
        if ( count ) std::cout << bin << ' ' << count << '\n';
    }
}
