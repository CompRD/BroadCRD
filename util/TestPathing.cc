///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file TestPathing.cc
 * \author tsharpe
 * \date Aug 9, 2012
 *
 * \brief
 */
#include "MainTools.h"
#include "kmers/ReadPatherDefs.h"
#include "paths/HyperBasevector.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/Unipath.h"
#include <iostream>

int const K = 28;

void dumpCompoundPaths( String const& READS, String const& EDGES )
{
    PathCollection<K> pc(PathCollection<K>::getInfoFilename(EDGES),
                            UnipathGraph<K>::getInfoFilename(READS));
    size_t nEdges = pc.getNReads();
    for ( size_t idx = 0; idx < nEdges; ++idx )
    {
        EdgeList const& el = pc.getEdgeList(idx);
        if ( el.size() != 1 )
        {
            std::cout << idx << ": " << el << std::endl;
        }
    }
}

void pathEdges( String const& READS, String const& EDGES, int NUM_THREADS )
{
    UnipathGraph<K> ug(UnipathGraph<K>::getInfoFilename(READS));
    VirtualMasterVec<bvec> edgeVMV(EDGES);
    PathBuilder<K> pb(PathCollection<K>::getInfoFilename(EDGES),ug);
    pb.processReads(edgeVMV,ug.getDict(),NUM_THREADS,100);
    pb.writeEvidence(NUM_THREADS);
}

void writeOldEdges( String const& READS, String const& EDGES, int NUM_THREADS )
{
    vecbvec reads(READS);
    vecKmerPath paths, paths_rc;
    vec<tagged_rpint> pathsdb;
    ReadsToPathsCoreY(reads, K, paths, paths_rc, pathsdb, "", NUM_THREADS);
    vecKmerPath unipaths;
    vec<tagged_rpint> unipathsdb;
    Unipath(paths, paths_rc, pathsdb, unipaths, unipathsdb);
    digraph A;
    BuildUnipathAdjacencyGraph(paths,paths_rc,pathsdb,unipaths,unipathsdb,A);
    HyperKmerPath h;
    BuildUnipathAdjacencyHyperKmerPath(K, A, unipaths, h);
    KmerBaseBroker kbb(K, paths, paths_rc, pathsdb, reads);
    HyperBasevector hb( h, kbb );
    size_t nEdges = hb.EdgeObjectCount();
    IncrementalWriter<bvec> writer(EDGES,nEdges);
    for ( size_t idx = 0; idx < nEdges; ++idx )
        writer.add(hb.EdgeObject(idx));
    writer.close();
}

int main( int argc, char** argv )
{
    RunTime();
    BeginCommandArguments;
    CommandArgument_String_Doc(READS, "The input reads.");
    CommandArgument_String_Doc(EDGES, "The old graph's edges.");
    CommandArgument_Int_OrDefault_Doc(NUM_THREADS,0,"Max number of threads.")
    CommandArgument_UnsignedInt_OrDefault_Doc(COVERAGE,50u,"Coverage estimate.");
    EndCommandArguments;

    NUM_THREADS = configNumThreads(NUM_THREADS);

    dumpCompoundPaths(READS,EDGES);
}
