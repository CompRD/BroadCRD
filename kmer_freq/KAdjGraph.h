/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef __INCLUDE_kmer_freq_KAdjGraph
#define __INCLUDE_kmer_freq_KAdjGraph



#include "CoreTools.h"
#include "math/Functions.h"
#include "kmers/KmerShape.h"
#include "kmers/KmerRecord.h"
#include "Bitvector.h"
#include "Intvector.h"
#include "CommonSemanticTypes.h"
#include "graph/Digraph.h"
#include "paths/KmerPathInterval.h"
#include "paths/KmerPath.h"
#include "feudal/BinaryStream.h"

// Semantic type: kseq_id_t
// Identifies a kmer.
// These are different from <kmer numbers> assigned
// by ReadsToPaths(); for that reason we are not
// calling this kmer_id_t, even though each value
// of type kseq_id_t identifies a kmer.
// Is the index of the kmer in KAdjGraph::kseqs.
typedef int kseq_id_t;
typedef IntVec KSeqIdVec;
typedef VecIntVec VecKSeqIdVec;

/**
    Struct: UnipathInfo

    Information about one unipath of a KAdjGraph.
 */
struct UnipathInfo {

  // Field: first_kseq;
  // The kseq_id of the first kmer of the unipath.
  kseq_id_t first_kseq;
  
  // Field: last_kseq;
  // The kseq_id of the last kmer of the unipath.
  kseq_id_t last_kseq;

  // Field: to
  // The number of kmers in the unipath.
  nkmers_t len;

  UnipathInfo() { }
  UnipathInfo( kseq_id_t _first_kseq, kseq_id_t _last_kseq, nkmers_t _len ):
    first_kseq( _first_kseq ), last_kseq( _last_kseq ), len( _len ) { }
  
};  // struct UnipathInfo

TRIVIALLY_SERIALIZABLE(UnipathInfo);

/**
   Class: KAdjGraph

   A graph of kmer adjacencies.  Kmers and their adjacencies
   are represented in sequence space (as opposed to by kmer numbers).
 */
template <nbases_t K>
class KAdjGraph
{
private:

  // Local type: kseq_set_id_t
  // Identifies a set of kseq_ids: if -1, the empty set;
  // if >= 0, a singleton set containing this kseq_id;
  // if < -1, the index in predsuccMult of the
  // kseq_id_t vector containing the (up to four)
  // predecessors or successors in the set.
  typedef int kseq_set_id_t;
  
public:

  KAdjGraph() { }

  // Constructor: KAdjGraph
  // Creates a KAdjGraph from a given set of kadjs.
  // *NOTE* that this destroys the _kadjs parameter!
  // (It gets swapped into KAdjGraph, and replaced
  // with an empty vecbasevector).  It's done this way
  // to avoid unneeded copying of large arrays.
  KAdjGraph( vec< kmer< K+1 > > & kadjs );

  /// MethodDecl: GetRc
  /// Given a kseq_id of a kmer, find the kseq_id for its rc.
  kseq_id_t GetRc( kseq_id_t k ) const;

  void ComputeUnibases( vecbasevector& unibases );

  /// MethodDecl: unipathAdjGraph
  /// Construct a <digraph> whose nodes are unipath ids (index of unipath in unipathInfos),
  /// and there is an edge (u1,u2) if there is an edge of the form
  /// (last-kmer-of(u1), first-kmer-of(u2)) in this KAdjGraph.
  void ComputeUnipathAdjGraph( digraph& unipathAdjGraph ) const;

  /// Test whether a given kmer has exactly one
  /// predecessor.
  Bool HasSolePred( kseq_id_t k ) const {
    return predsucc[ k ].first >= 0;
  }

  /// Test whether a given kmer has exactly one
  /// successor.
  Bool HasSoleSucc( kseq_id_t k ) const {
    return predsucc[ k ].second >= 0;
  }

  kseq_id_t GetSolePred( kseq_id_t k ) const {
    ForceAssert( HasSolePred( k ) );
    return predsucc[ k ].first;
  }

  kseq_id_t GetSoleSucc( kseq_id_t k ) const {
    ForceAssert( HasSoleSucc( k ) );
    return predsucc[ k ].second;
  }

  // Test whether the given kmer has the given successor
  // in the kmer graph
  Bool HasSucc( kseq_id_t k, kseq_id_t k_succ ) const {
    return predsucc[ k ].second == k_succ ||
      predsucc[ k ].second < -1  &&
      ShortVecContains( predsuccMult[ -predsucc[ k ].second - 2 ], k_succ );
  }
  
  // Test whether the given kmer has the given predecessor
  // in the kmer graph
  Bool HasPred( kseq_id_t k, kseq_id_t k_pred ) const {
    return predsucc[ k ].first == k_pred ||
      predsucc[ k ].first < -1  &&
      ShortVecContains( predsuccMult[ -predsucc[ k ].first - 2 ], k_pred );
  }

  void RemoveOrphans(int maxOrphanLength = 10);

  void ShaveGraph(int minLength = 100, int maxShaveLength = 20,
		  bool singleHairsOnly = False);

  void RemoveUnipathKAdjs(int unipath_index);

  void DisconnectUnipath(int unipath_index);

  void CheckKAdjGraphSanity() const;
  
private:  
  void GetPredsOrSuccs( kseq_set_id_t z, vec< kseq_id_t >& predsOrSuccs ) const {
    predsOrSuccs.clear();
    if ( z >= 0 )
      predsOrSuccs.push_back( z );
    else if ( z < -1 ) {
      KSeqIdVec const& predsOrSuccsVec = predsuccMult[ -z -2 ];
      predsOrSuccs.insert( predsOrSuccs.end(), predsOrSuccsVec.begin(), predsOrSuccsVec.end() );
    }
  }
  
public:

  void GetSuccs( kseq_id_t k, vec< kseq_id_t >& succs ) const {
    GetPredsOrSuccs( predsucc[ k ].second, succs );
  }
  
  void GetPreds( kseq_id_t k, vec< kseq_id_t >& preds ) const {
    GetPredsOrSuccs( predsucc[ k ].first, preds );
  }


  void WriteKmerGraphAround( vec< kseq_id_t > seeds, int radius, filename_t dotFileName );

  void GetKAdjs( vec< kmer< K+1 > > & kadjs );

  void GetKAdjs( vec< kmer_with_count< K+1 > > & kadjs );

  void writeBinary( BinaryWriter& writer ) const
  { writer.write(kseqs);
    writer.write(kseq2rc);
    writer.write(predsucc);
    writer.write(unipathInfos);
    writer.write(kseq2upath); }

  void readBinary( BinaryReader& reader )
  { reader.read(&kseqs);
    reader.read(&kseq2rc);
    reader.read(&predsucc);
    reader.read(&unipathInfos);
    reader.read(&kseq2upath); }

  static size_t externalSizeof() { return 0; }

private:
  
  /// Field: kseqs
  /// The sequences of all kmers.  Both fw and rc versions of each
  /// kmer are included; palindromic kmers are included only once.
  /// The type kseq_id_t denotes a kmer's index in this array.
  vec< kmer< K > > kseqs;

  /// Field: kseq2rc
  /// Maps each kseq to the kseq_id of its rc
  vec< kseq_id_t > kseq2rc;

  /// Field: predsucc
  /// For each kmer in kseq, the index in kseq of its predecessor and
  /// its successor in the unipath graph.  If no successor (predecessor),
  /// the corresponding value is -1.  If more than one successor (predecessor),
  /// then the value is -(i+2) where i is the index in 'predsuccMult' of the list of
  /// predecessors or successors.  This representation is chosen because
  /// most kmers lie in the middle of unipaths, and thus have exactly one predecessor
  /// and one successor.
  vec< pair< kseq_set_id_t, kseq_set_id_t > > predsucc;

  /// Field: predsuccMult
  /// For kmers in kseq which have multiple predecessors or successors ( up to four ),
  /// the list of these predecessors/successors, as indices into kseqs.
  /// See also documentation of predsucc above.
  VecKSeqIdVec predsuccMult;

  /// Field: unipathsInfos
  /// For each unipath, the kseq_ids of its first and last kmers,
  /// and its length.
  vec< UnipathInfo > unipathInfos;

  /// Field: kseq2upath
  /// For each kseq, the id of the unipath ( in unipathInfos )
  /// to which this kseq belongs.
  vec< unipath_id_t > kseq2upath;

  /// MethodDecl: AddKAdj
  /// Record the fact that the adjacency (k1,k2) is genomic.
  void AddKAdj( kseq_id_t k1, kseq_id_t k2 );

  /// MethodDecl: RemoveKAdj
  /// Record the fact that the adjacency (k1,k2) is non-genomic.
  void RemoveKAdj( kseq_id_t k1, kseq_id_t k2 );

  void ComputeUnipathInfos();
  
  /// MethodDecl: ShortVecContains
  /// Test whether a short vector ( of length up to 4 ) contains the given
  /// value.
  static Bool ShortVecContains( KSeqIdVec const& v, kseq_id_t k );
  static Bool ShortVecDoesNotContain( KSeqIdVec const& v, kseq_id_t k );
  
};  // class KAdjGraph

template <nbases_t K>
struct Serializability< KAdjGraph<K> >
{ typedef SelfSerializable type; };

#endif
// #ifndef __INCLUDE_kmer_freq_KAdjGraph
