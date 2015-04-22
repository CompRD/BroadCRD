///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * DumpHBV.cc
 *
 *  Created on: Dec 18, 2013
 *      Author: tsharpe
 */
#include "MainTools.h"
#include "feudal/BinaryStream.h"
#include "paths/HyperBasevector.h"
#include "system/SortInPlace.h"

int main( int argc, char** argv )
{
    RunTime();
    BeginCommandArguments;
    CommandArgument_String(HBV);
    EndCommandArguments;

    HyperBasevector hbv;
    BinaryReader::readFile(HBV,&hbv);
    unsigned BIGK = hbv.K();
#if 0
    vecbvec edges(hbv.Edges().begin(),hbv.Edges().end());
    edges.Sort();
    for ( bvec const& edge : edges )
        std::cout << edge << std::endl;
#endif
    vec<bvec> const& edges = hbv.Edges();
    size_t nKmers = 0;
    for ( bvec const& edge : edges )
        nKmers += edge.size()-BIGK+1;
    vecbvec kmers;
    kmers.reserve(nKmers);
    for ( bvec const& edge : edges )
    {
        auto beg = edge.begin();
        auto end = edge.end()-BIGK+1;
        while ( beg != end )
        {
            kmers.push_back(bvec(beg,beg+BIGK));
            ++beg;
        }
    }
    kmers.Sort();
    for ( bvec const& kmer : kmers )
        std::cout << kmer << std::endl;
}
