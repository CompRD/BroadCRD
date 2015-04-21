///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file ParallelSortTest2.cc
 * \author tsharpe
 * \date Dec 3, 2012
 *
 * \brief
 */

#include "MainTools.h"
#include "Basevector.h"
#include "ParallelVecUtilities.h"
#include "system/SortInPlace.h"

int main( int argc, char *argv[] )
{
     RunTime();
     BeginCommandArguments;
     CommandArgument_UnsignedInt_OrDefault_Doc(G,1000000,"Genome size.");
     CommandArgument_UnsignedInt_OrDefault_Doc(K,31,"Kmer size.");
     CommandArgument_UnsignedInt_OrDefault_Doc(COV,200,"Coverage.");
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS,0,"Number of threads.");
     EndCommandArguments;
     NUM_THREADS = configNumThreads(NUM_THREADS);

     bvec genome;
     genome.reserve(G);
     for ( size_t count = 0; count != G; ++count )
         genome.push_back(random()%4u);

     size_t picks = 1UL*COV*G/K;
     vecbvec kmers;
     kmers.reserve(picks);
     size_t modulus = G-K+1;
     bvec tmp;
     for ( size_t count = 0; count != picks; ++count )
     {
         size_t idx = random()%modulus;
         bvec::const_iterator itr(genome.cbegin(idx));
         tmp.assign(itr,itr+K);
         kmers.push_back(tmp);
     }
     vecbvec kmers2(kmers);

     double clock = WallClockTime();
     sortInPlaceParallel(kmers.begin(),kmers.end(),NUM_THREADS);
     std::cout << "In place: " << TimeSince(clock) << std::endl;

     clock = WallClockTime();
     __gnu_parallel::sort(kmers2.begin(),kmers2.end());
     std::cout << "STL parallel: " << TimeSince(clock) << std::endl;

     ForceAssert(kmers==kmers2);
}
