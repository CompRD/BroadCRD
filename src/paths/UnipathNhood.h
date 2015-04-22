/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef UNIPATH_NHOOD_H
#define UNIPATH_NHOOD_H

// File: UnipathNhood.h
//
// This file defines a toolkit for building a neighborhood (abbreviated "nhood")
// of unipaths around a given seed unipath, and identifying the reads that came from
// this neighborhood.

#include "CoreTools.h"
#include "Equiv.h"
#include "Basevector.h"
#include "Intvector.h"
#include "ReadLocation.h"
#include "ReadPairing.h"
#include "SemanticTypes.h"
#include "graph/Digraph.h"
#include "paths/KmerPath.h"
#include "paths/PdfEntry.h"
#include "paths/simulation/Placement.h"
#include "paths/UnipathNhood.h"
#include "paths/UnipathNhoodCommon.h"

// FuncDecl: BuildUnipathLinkGraph
//
// Build the unipath graph, in which vertices are 
// normal unipaths and edges come from read pairs.
//
// Instantiated for sepdev and fsepdev in the .cc file.
template<class T>
void BuildUnipathLinkGraph( 

     // inputs:

     const int K,                              // as in Kmer
     const vec<read_pairing>& pairs,           // read pairs
     const vec<read_location_short>& ulocs,    // locations of reads on unipaths
     const VecULongVec& ulocs_indexr,         // index to it by reads
     const vec<Bool>& normal,                  // is a given unipath normal
     const vec<nbases_t>& ulen,                // unipath lengths
     const vecKmerPath& paths,                 // the reads
     const vec<unipath_id_t>& to_rc,           // map unipath to its rc
     const int min_edge_multiplicity,          // else ignore edge
     const VecPdfEntryVec& cp,                 // unipath copy number

     // output:

     digraphE< Tsepdev<T> >& G,                      // the graph

     // optional:

     bool verbose = false
          );

/**
   FuncDecl: FindUnipathNhood

   Given a seed unipath, find all nearby <normal unipaths> to form a
   <neighborhood> around this seed.
*/
void FindUnipathNhood( 

     // Inputs:

     const int v,                             // seed unipath
     const digraphE<sepdev>& G,               // graph of all normal unipaths
     const vec<int>& ulen,                    // length of each unipath
     const VecPdfEntryVec& cp,             // copy number pdf for unipaths
     const vec<int>& predicted_copyno,        // predicted copy number for unipaths
     const vec<Bool>& branch,                 // which unipaths are at branches
     const vecKmerPath& paths,                // the read paths
     const vec<int>& path_lens,               // the read path lengths
     const vec<read_location_short>& ulocs,   // locations of reads on unipaths
     const vec<int>& uindex,                  // index to ulocs
     const vec<read_pairing>& pairs,          // all the read pairs
     const vec<int>& pairs_index,             // index to read pairs by read ids
     const Bool FILTER_NHOOD,                 // processing option: filter?
     const int MAX_COPY_NUMBER_OTHER,         // screens unipaths entering nhood
     const int NHOOD_RADIUS,                  // how far we go away from the seed
     const int MAX_TO_CUT,                    // cutpoints longer than this are kept
     const int MAX_DEV,                       // how stretchy sep from seed can get
     const int MAX_PATHS_IN_NHOOD,            // how big the nhood can get
     const Bool BUILD_NHOOD_FW_ONLY,          // forward only?

     // Output:

     vec<ustart>& processed                   // positions of unipaths in nhood

          );



/// FuncDecl: FillInTransitiveEdges
///
/// For each vertex of predicted-copy-number one in the graph, join it
/// to other vertices whose separation will be within radius of that
/// vertex.  But this naive version won't find connections which
/// require passing though multiple vertices further than distance
/// radius away.
///
/// An edge we construct will replace an existing edge with the same
/// endpoints if the new one has smaller deviation.
///
/// Do not introduce new edges having deviation > max_dev.
///
/// Do not replace an edge unless its deviation decreases by percent_improvement.
template<class T>   // defined for int and double
void FillInTransitiveEdges( digraphE< Tsepdev<T> >& graph, 
			    const int radius, 
                            const double max_dev,
                            const double percent_improvement,
			    const vec<int>& predicted_copyno,
			    const vec<nbases_t>& unipath_len, 
                            const int min_unipath = 0,
			    int verbosity = 0,
			    VecPlacementVec* locs_p = NULL );

/**
   FuncDecl: PopulateNhoodWithReads

   Find the reads which go in a particular neighborhood.   Also return their orientations and predicted positions.

   Input parameters:

     int v                                  -  seed for nhood: unipath id of the seed for this nhood
     const vec<ustart>& processed           -  this defines the nhood: for each unipath in the nhood,
                                               its id and its position relative to the seed 'v'
     int K                                  -  size of kmers in terms of which all paths are defined
     const vec<int> ulen,                   -  unipath lengths, for each unipath id.
     const vec<read_pairing>& pairs         -  read pairs: what reads go into each pair
     const vec<pair_id_t>& pairs_index            -  index to read pairs: for each read, the index in 'pairs'
                                               of the read_pairing structure describing its read and its partner
     const vecKmerPath& paths               -  kmer paths for reads (the <read paths>)
     const vec<read_location_short>& ulocs  -  locations of reads on unipaths.  each element describes the location
                                               of one read on one unipath; locations of reads on a given unipath
					       are contiguous in the array.  'uindex' maps each unipath to the block
					       of read alignments to that unipath in 'ulocs'.
     const vec<int>& uindex,                -  index to ulocs: locations of reads aligned to unipath w are in
                                               ulocs[] locations [ uindex[w], uindex[w+1] ).
     const int NHOOD_RADIUS_INTERNAL        -  how far from origin we should go
     const int MAX_DEV                      -  how sloppy read locations can get
     const Bool REACH_FW_ONLY               -  only go forward?

   Output parameters:

     vec< pair<read_id_t,orient_t> >& use   -  reads in nhood, with orientations
     vec<pos_rel_to_seed_t>& usestart                     -  predicted start of read relative to the seed 'v'
     vec<pair_id_t>& pairs_to_use                 -  the constituent pairs
   
*/
void PopulateNhoodWithReads(

     // Inputs:

     int v,                                 // seed for nhood
     const vec<ustart>& processed,          // this defines the nhood
     int K,                                 // as in Kmer
     const vec<int> ulen,                   // unipath lengths
     const vec<read_pairing>& pairs,        // read pairs
     const vec<pair_id_t>& pairs_index,           // index to read pairs
     const vecKmerPath& paths,              // kmer paths for reads
     const vec<read_location_short>& ulocs, // locations of reads on unipaths
     const vec<int>& uindex,                // index to ulocs
     const int NHOOD_RADIUS_INTERNAL,       // how far from origin we should go
     const int MAX_DEV,                     // how sloppy read locations can get
     const Bool REACH_FW_ONLY,              // only go forward?

     // Outputs:

     vec< pair<read_id_t,orient_t> >& use,  // reads in nhood, with orientations
     vec<pos_rel_to_seed_t>& usestart,      // predicted start of read rel. v
     vec<pair_id_t>& pairs_to_use                 // the constituent pairs

          );

// Find the short-insert read pairs such that each of its reads could be
// completely contained in a contig created from the reads in a populated nhood.

void GetShortInsertReads(

     // Inputs:

     vec< pair<read_id_t,orient_t> >& use,            // reads in nhood, with orientations
     const vec<tagged_rpint>& pathsdb,      // paths database for reads
     const vecKmerPath& paths,              // kmer paths for reads
     const vecKmerPath& paths_rc,           // kmer paths for rc of reads

     const vec<nbases_t>& PATH_KS,
     const vec< vec<tagged_rpint> >& extra_paths_db,
     const vec< vecKmerPath >& extra_paths,
     const vec< vecKmerPath >& extra_paths_rc,			   
     
     const vec<read_id_t>& partner,               // map read to partner
     const vec<Bool>& is_short_pair_read,   // is read an end of a short-insert pair
     

     // Output:

     vec< pair<read_id_t,orient_t> >& P               // the short insert reads, orientations

          );




Bool Linked( int x, int y, const vec<read_location_short>& ulocs,
     const vec<int>& uindex, const vec<read_pairing>& pairs, 
     const vec<int>& pairs_index );

// CalcLinkProbability: given normal unipaths x and y, with hypothetical separation 
// s +/- d (x --> y), estimate the probability that there is at least one link from 
// x to y with separation sep +/- dev such that |sep-s#| <= 2*dev, where s# is an
// instantiation of the normal random variable corresponding to s +/- d.

// Calculate the answer.
// Actually this still uses a little bit of Monte Carlo, but only
// asks for one normally-distributed random number per iteration.

double CalcLinkProbability( int x, int y, int s, int d, const vec<int>& ulen,
     const vec<read_location_short>& ulocs, const vec<int>& uindex,
     const VecPdfEntryVec& cp, const vec<read_pairing>& pairs,
     const vec<int>& pairs_index, const vecKmerPath& paths,
     const vec<int>& path_lens );

#endif
