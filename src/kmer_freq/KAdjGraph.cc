////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include <map>
#include <queue>

#include "CoreTools.h"
#include "Set.h"
#include "kmers/SortKmers.h"
#include "kmer_freq/KAdjGraph.h"
#include "paths/UnibaseUtils.h"

/**
   Constructor: KAdjGraph constructor

   Construct a KAdjGraph from the list of (K+1)-mers.
 */
template <nbases_t K>
KAdjGraph< K >::KAdjGraph( vec< kmer< K+1 > > & kadjs ) {

  cout << "Starting Kmer Adjacency size = " << kadjs.size() << endl;

  typedef int kadj_id_t;
  

  // Extract the kseqs that appear in the kadjs.
  // Each kadj has two kseqs: the first and the second.
  // For each kadj, record the ids of its two kseqs.
  
  vec< pair< kseq_id_t, kseq_id_t > > kadj2kseqIds( kadjs.size(), make_pair( -1, -1 ) );

  typedef typename KmerShapeDflt< K >::type mykshape_t;
  typedef kmer_record< K, 2 > mykrec_t;

  kseqs.reserve( kadjs.size() );
  kseq2rc.reserve( kadjs.size() );

  // turn kadjs into a vecbasevector.
  // proper solution: SortKmers() should be templatized on the iterator that 
  // iterates over the records, and for each returns the dataasints().
  // So, we need to define a Basevec concept, modeled by both
  // Basevector and kmer records, but implemented more efficiently
  // by the fixed-size kmer records.

  vecbasevector kadjs_basevecs( kadjs.size() );
  for ( kadj_id_t i = 0; i < kadjs.isize(); i++ )
    kadjs[ i ].GetBasevector( kadjs_basevecs[ i ] );

  FOR_KMER_OCCS_BEG( mykshape_t, mykrec_t, kadjs_basevecs, DONT_USE_STABLE_SORT, PALIND_FW_ONLY, from, to ) {
    basevector b, b_rc;
    from->GetBasevector( b );
    b_rc.ReverseComplement( b );

    kseq_id_t kseq_id_fw = kseqs.isize(), kseq_id_rc = kseqs.isize();
    kseqs.push_back( kmer< K >( b ) );
    if ( b != b_rc ) {
      kseq_id_rc = kseqs.isize();
      kseqs.push_back( kmer< K >( b_rc ) );
      kseq2rc.push_back( kseqs.isize()-1, kseqs.isize()-2 );
    } else
      kseq2rc.push_back( kseqs.isize()-1 );

    for ( typename vec< mykrec_t >::const_iterator kadjOcc = from; kadjOcc != to; kadjOcc++ ) {
      kadj_id_t kadjId = kadjOcc->GetId();
      int posInKadj = abs( kadjOcc->GetPos() ) - 1;

      ForceAssert( kadjOcc->EqualKmers( *from ) );
      ForceAssert( posInKadj == 0  ||  posInKadj == 1 );
      ForceAssertLe( 0, kadjId );
      ForceAssertLt( kadjId, kadjs.isize() );

      pair< kseq_id_t, kseq_id_t >& this_kadj2kseqIds = kadj2kseqIds[ kadjId ];
      ( posInKadj == 0 ? this_kadj2kseqIds.first : this_kadj2kseqIds.second ) =
	kadjOcc->IsFw() ? kseq_id_fw : kseq_id_rc;
    }
  }  // FOR_KMER_OCCS_BEG()
  FOR_KMER_OCCS_END( );

  Destroy( kadjs_basevecs );

  // Cross-reference the kseqs and kadjs:
  // for each kseq, store the predecessors and successors.

  predsucc.resize( kseqs.size(), make_pair( -1, -1 ) );

  cout << "\nCross-referencing the unique " << kseqs.size() << " seqs..." << endl;
  for ( kadj_id_t kadjId = 0; kadjId < kadjs.isize(); kadjId++ ) {
    const pair< kseq_id_t, kseq_id_t >& kadj = kadj2kseqIds[ kadjId ];
    if ( kadj.first != -1  &&  kadj.second != -1 ) {
      kseq_id_t k1 = kadj.first, k2 = kadj.second;
      AddKAdj( k1, k2 );
      AddKAdj( GetRc( k2 ), GetRc( k1 ) );
    }  // if this kadj links together two kseqs
  }  // for each kadj

  CheckKAdjGraphSanity();

  cout << "Computing unipath infos..." << endl;
  ComputeUnipathInfos();

}  // KAdjGraph constructor

template < nbases_t K >
void KAdjGraph< K >::CheckKAdjGraphSanity() const {

  cout << "Sanity Checking " << kseqs.isize() << " kseqs..." << endl;
  for ( kseq_id_t k = 0; k < kseqs.isize(); k++ ) {
    if ( !( k % 500000 ) )
      PRINT( k );
    {
      vec< kseq_id_t > succs;
      GetSuccs( k , succs );
      for ( int i = 0; i < succs.isize(); i++ ) {
	kseq_id_t ksucc = succs[ i ];
	if ( !HasSucc( k, ksucc ) ) {
	  PRINT5( k, i, succs.isize(), k, ksucc );
	}
	ForceAssert( HasSucc( k, ksucc ) );
	ForceAssert( HasPred( ksucc, k ) );
	ForceAssert( HasSucc( GetRc( ksucc ), GetRc( k ) ) );
	ForceAssert( HasPred( GetRc( k ), GetRc( ksucc ) ) );
      }
    }

    {
      vec< kseq_id_t > preds;

      GetPreds( k , preds );
      for ( int i = 0; i < preds.isize(); i++ ) {
	kseq_id_t kpred = preds[ i ];
	ForceAssert( HasPred( k, kpred ) );
	ForceAssert( HasSucc( kpred, k ) );
	ForceAssert( HasPred( GetRc( kpred ), GetRc( k ) ) );
	ForceAssert( HasSucc( GetRc( k ), GetRc( kpred ) ) );
      }
    }

    basevector bv, bv_rc;
    kseqs[ k ].GetBasevector( bv );
    kseqs[ GetRc( k ) ].GetBasevector( bv_rc );
    bv_rc.ReverseComplement();
    ForceAssert( bv_rc == bv );
  }
  
}

/// Method: GetRc
/// Given a kseq_id of a kmer, find the kseq_id for its rc.
template < nbases_t K >
kseq_id_t KAdjGraph< K >::GetRc( kseq_id_t k ) const {
  return kseq2rc[ k ];
}


template < nbases_t K >
Bool KAdjGraph< K >::ShortVecDoesNotContain( KSeqIdVec const& v, kseq_id_t k ) {
  return v.empty() ? True : 
    k != v.front() &&
    ( v.size() < 2  ||
      ( k != v[1]  &&
	(v.size() < 3  ||
	 ( k != v[2]  &&
	   ( v.size() < 4  ||
	     k != v[3] ) ) ) ) );
}

template < nbases_t K >
Bool KAdjGraph< K >::ShortVecContains( KSeqIdVec const& v, kseq_id_t k ) {
  return !ShortVecDoesNotContain( v, k );
}

/// Method: AddKAdj
/// Record the fact that the adjacency (k1,k2) is genomic.
template < nbases_t K >
void KAdjGraph< K >::AddKAdj( kseq_id_t k1, kseq_id_t k2 ) {
  {
    int& k1succs = predsucc[ k1 ].second;
    
    if ( k1succs == -1 ) {
      // we have not previously recorded a successor to k1,
      // so now we have exactly one successor.
      // that is represented by simply recording the successor's
      // kseq_id.
      k1succs = k2;
    } else {
      if ( k1succs == k2 ) {
	ForceAssert( HasPred( k2, k1 ) );
	return;
      } else {
	if ( k1succs >= 0 ) {
	  // previously, k1 had one successor.
	  // now we are adding a second, so we have to
	  // create an array to hold the list of k1's successors,
	  // and record in k1succs the index of that array in predsuccMult.
	  kseq_id_t old_k1succs = k1succs;
	  k1succs = - ( predsuccMult.size() + 2 );
	  predsuccMult.push_back_reserve( KSeqIdVec( 1, old_k1succs ) );
	}
	KSeqIdVec& k1succsVec = predsuccMult[ - k1succs - 2 ];
	if ( ShortVecContains( k1succsVec, k2 ) ) {
	  return;
	}
	k1succsVec.push_back( k2 );
      }
    }
  }

  {
    int& k2preds = predsucc[ k2 ].first;
    if ( k2preds == -1 ) {
      k2preds = k1;
    } else if ( k2preds == k1 ) {
      ForceAssert( HasSucc( k1, k2 ) );
      return;
    } else {
      if ( k2preds >= 0 ) {
	kseq_id_t old_k2preds = k2preds;
	k2preds = - ( predsuccMult.size() + 2 );
	predsuccMult.push_back_reserve( KSeqIdVec( 1, old_k2preds ) );
      }
      KSeqIdVec& k2predsVec = predsuccMult[ - k2preds - 2 ];
      if ( ShortVecContains( k2predsVec, k1 ) ) {
	return;
      }
      k2predsVec.push_back( k1 );
    }
  }
}  // KAdjGraph::AddKAdj()

/// Method: RemoveKAdj
/// Record the fact that the adjacency (k1,k2) is non-genomic.
template < nbases_t K >
void KAdjGraph< K >::RemoveKAdj( kseq_id_t k1, kseq_id_t k2 ) {

  ForceAssert( HasSucc( k1, k2 ) );
  ForceAssert( HasPred( k2, k1 ) );

  for ( int i = 0; i < 2; i++ ) {
    int& k1succs = ( i==0 ? predsucc[ k1 ].second : predsucc[ k2 ].first );
    kseq_id_t kOther = ( i==0 ? k2 : k1 );
    if ( k1succs == kOther )
      k1succs = -1;
    else {
      ForceAssert( k1succs < -1 );
      KSeqIdVec& k1succsVec = predsuccMult[ -k1succs - 2 ];
      if ( k1succsVec.size() == 2 ) {
	k1succs = k1succsVec[ k1succsVec[0] == kOther ? 1 : 0 ];
	k1succsVec.clear();
      } else {
	remove( k1succsVec.begin(), k1succsVec.end(), kOther );
	k1succsVec.resize(k1succsVec.size() - 1);
      }
    }
  }

  ForceAssert( !HasSucc( k1, k2 ) );
  ForceAssert( !HasPred( k2, k1 ) );
}


template < nbases_t K >
void KAdjGraph< K >::WriteKmerGraphAround( vec< kseq_id_t > seeds, int radius, filename_t dotFileName ) {
  vec< kseq_id_t > seeds_rc;
  For_( kseq_id_t, s, seeds )
    seeds_rc.push_back( GetRc( *s ) );
  seeds.append( seeds_rc );
  UniqueSort( seeds );
  PRINT( seeds );

  set< kseq_id_t > relevantNodes;
  
  typedef int distFromSeed_t;

  vec< kseq_id_t > active;
  map< kseq_id_t, distFromSeed_t > active2dist;
  
  For_( kseq_id_t, s, seeds ) {
    active.push_back( *s );
    active2dist.insert( make_pair( *s, 0 ) );
  }

  while( active.nonempty() ) {
    kseq_id_t k = active.back();
    active.pop_back();
    if ( !STLContains( relevantNodes, k ) ) {
      relevantNodes.insert( k );
      distFromSeed_t distFromSeed = active2dist[ k ];
      if ( distFromSeed < radius ) {
	for ( int ps = 0; ps < 2; ps++ ) {
	  
	  vec< kseq_id_t > succs;
	  if ( ps == 0 )
	    GetSuccs( k, succs );
	  else
	    GetPreds( k, succs );
	  PRINT3( ps, k, succs );
	  For_( kseq_id_t, ks, succs ) {
	    if ( !STLContains( relevantNodes, *ks ) ) {
	      map< kseq_id_t, distFromSeed_t >::iterator it = active2dist.find( *ks );
	      if ( it != active2dist.end() )
		it->second = min( it->second, distFromSeed+1 );
	      else {
		active2dist.insert( make_pair( *ks, distFromSeed+1 ) );
	      }
	      active.push_back( *ks );
	      
	    }
	  }
	}
      }
    }
  }

  Ofstream( df, dotFileName );

  vec< kseq_id_t > relNodes;
  relNodes.insert( relNodes.end(), relevantNodes.begin(), relevantNodes.end() );
  UniqueSort( relNodes );

  df << "digraph kseqs {\n";

  For_( kseq_id_t, r, relNodes ) {
    df << *r << " [ label=\"" << *r << "\\n" << kseqs[ *r ].ToString() << "\", shape=egg ];" << endl;
    vec< kseq_id_t > succs;
    GetSuccs( *r, succs );
    For_( kseq_id_t, s, succs )
      if ( BinMember( relNodes, *s ) )
	df << *r << " -> " << *s << endl;
  }
  
  df << "}\n";

  
}


template < nbases_t K >
void KAdjGraph< K >::ComputeUnipathInfos( ) {

  unipathInfos.clear();
  
  PRINT( kseqs.isize() );

  kseq2upath.clear();
  kseq2upath.resize( kseqs.size(), -1 );
  
  for ( kseq_id_t k = 0; k < kseqs.isize(); k++ ) {

    if ( kseq2upath[ k ] != -1 ) {
      continue;
    }

    unipath_id_t unipathId = unipathInfos.isize();
    kseq2upath[ k ] = unipathId;

    // Look forward and backward from this kmer,
    // moving in each direction along the kmer graph
    // until we hit a branch.

    nkmers_t unipathLen = 1;
    
    kseq_id_t kBef = k;
    kseq_id_t kBefPred = -1;
    while ( HasSolePred( kBef ) &&
	    kseq2upath[ kBefPred = GetSolePred( kBef ) ] == -1 &&
	    HasSoleSucc( kBefPred ) &&
	    GetSoleSucc( kBefPred ) == kBef ) {
      kseq2upath[ kBef = kBefPred ] = unipathId;
      unipathLen++;
    }

    kseq_id_t kAft = k;
    kseq_id_t kAftSucc = -1;
    while ( HasSoleSucc( kAft ) &&
	    kseq2upath[ kAftSucc = GetSoleSucc( kAft ) ] == -1 &&
	    HasSolePred( kAftSucc ) &&
	    GetSolePred( kAftSucc ) == kAft ) {
      kseq2upath[  kAft = kAftSucc ] = unipathId;
      unipathLen++;
    }

    unipathInfos.push( kBef, kAft, unipathLen );
  }

  // For circular unipaths, make sure that the start and end
  // on both copies is the same.

  for ( unipath_id_t u = 0; u < unipathInfos.isize(); u++ ) {
    const UnipathInfo& ui = unipathInfos[ u ];
    if ( HasSoleSucc( ui.last_kseq )  &&  GetSoleSucc( ui.last_kseq ) == ui.first_kseq &&
	 HasSolePred( ui.first_kseq )  &&  GetSolePred( ui.first_kseq ) == ui.last_kseq ) {
      unipath_id_t u_rc = kseq2upath[ GetRc( ui.first_kseq ) ];
      UnipathInfo& ui_rc = unipathInfos[ u_rc ];
      
      ForceAssert( HasSoleSucc( ui_rc.last_kseq ) );
      ForceAssert( GetSoleSucc( ui_rc.last_kseq ) == ui_rc.first_kseq );
      ForceAssert( HasSolePred( ui_rc.first_kseq ) );
      ForceAssert( GetSolePred( ui_rc.first_kseq ) == ui_rc.last_kseq );
      if ( u_rc > u ) {
	ui_rc.first_kseq = GetRc( ui.last_kseq );
	ui_rc.last_kseq = GetRc( ui.first_kseq );
      }
    }
  }

  vec< nkmers_t > unipathLens( unipathInfos.isize() );
  for ( unipath_id_t u = 0; u < unipathInfos.isize(); u++ )
    unipathLens[ u ] = unipathInfos[ u ].len;
  cout << "Unipath N50 = " << N50( unipathLens ) << endl;

  vec< int > succCounts( 5, 0 );
  vec< int > predCounts( 5, 0 );
  
  typedef pair< kseq_id_t, kseq_id_t > ps_t;
  For_( ps_t, ps, predsucc )
    for ( int i = 0; i < 2; i++ ) {
      kseq_id_t kseq_id = i==0 ? ps->first : ps->second; 
      ( i == 0 ? predCounts : succCounts )[ kseq_id == -1 ? 0 : ( kseq_id >= 0 ? 1 : predsuccMult[ -kseq_id - 2 ].size() ) ]++;
  }
  // PRINT2( succCounts, predCounts );

  {
    // check the sanity of unipath infos

    // check that every kmer became part of a unipath
    for ( kseq_id_t k = 0; k < kseqs.isize(); k++ ) {
      unipath_id_t u = kseq2upath[ k ];
      ForceAssertLe( 0, u );
      ForceAssertLt( u, unipathInfos.isize() );
    }

    cout << "Sanity Checking " << unipathInfos.isize() << " unipath infos..." << endl;
    for ( unipath_id_t u = 0; u < unipathInfos.isize(); u++ ) {
      if ( !( u % 10000 ) )
	PRINT( u );
      const UnipathInfo& ui = unipathInfos[ u ];
      ForceAssertEq( kseq2upath[ ui.first_kseq ], u );
      ForceAssertEq( kseq2upath[ ui.last_kseq ], u );

      unipath_id_t u_rc = kseq2upath[ GetRc( unipathInfos[ u ].first_kseq ) ];
      const UnipathInfo& ui_rc = unipathInfos[ u_rc ];

      ForceAssertEq( ui.len, ui_rc.len );
      ForceAssertEq( GetRc( ui.last_kseq ), ui_rc.first_kseq );
      ForceAssertEq( GetRc( ui.first_kseq ), ui_rc.last_kseq );

      nkmers_t kmersSeen = 0;
      
      kseq_id_t k = ui.first_kseq, k_rc = ui_rc.last_kseq;
      for ( ; k != ui.last_kseq; k = GetSoleSucc( k ), k_rc = GetSolePred( k_rc ) ) {
	ForceAssertEq( kseq2upath[ k ], u );
	ForceAssertEq( kseq2upath[ k_rc ], u_rc );
	
	ForceAssertEq( GetRc( k ), k_rc );
	ForceAssertEq( GetRc( k_rc ), k );
	kmersSeen++;
      }
      ForceAssertEq( k, ui.last_kseq );
      ForceAssertEq( k_rc, ui_rc.first_kseq );
      ForceAssertEq( GetRc( k ), k_rc );
      ForceAssertEq( GetRc( k_rc ), k );
      kmersSeen++;
      ForceAssertEq( kmersSeen, ui.len );
    }
  }  // check sanity of unipath infos
  
  
}  // KAdjGraph::ComputeUnipathInfos()


/// Method: unipathAdjGraph
/// Construct a <digraph> whose nodes are unipath ids (index of unipath in unipathInfos),
/// and there is an edge (u1,u2) if there is an edge of the form
/// (last-kmer-of(u1), first-kmer-of(u2)) in this KAdjGraph.
template < nbases_t K >
void KAdjGraph< K >::ComputeUnipathAdjGraph( digraph& unipathAdjGraph ) const {

  int nuni = unipathInfos.isize();
  unipathAdjGraph.Clear( );
  unipathAdjGraph.ToMutable( ).resize(nuni);
  unipathAdjGraph.FromMutable( ).resize(nuni);

  vec< kseq_id_t > succs;
  for ( unipath_id_t u = 0; u < nuni; u++ ) {
    const UnipathInfo& uinfo = unipathInfos[ u ];
    GetSuccs( uinfo.last_kseq, succs );
    For_( kseq_id_t, succ, succs )
      unipathAdjGraph.AddEdge( u, kseq2upath[ *succ ] );
  }
  
}  // KAdjGraph< K >::ComputeUnipathAdjGraph( )


template < nbases_t K >
void KAdjGraph< K >::ComputeUnibases( vecbasevector& unibases ) {
  unibases.resize( unipathInfos.size() );
  basevector btmp;
  for ( unipath_id_t u = 0; u < unipathInfos.isize(); u++ ) {
    const UnipathInfo& uinfo = unipathInfos[ u ];
    basevector& unibase = unibases[ u ];
    
    kseq_id_t k = uinfo.first_kseq;
    kseqs[ k ].GetBasevector( unibase );
    unibase.resize( uinfo.len + K - 1);
    for ( int i = 1; i < uinfo.len; i++ ) {
      k = GetSoleSucc( k );
      kseqs[ k ].GetBasevector( btmp );
      ForceAssertLe( K + i - 1, unibase.isize() );
      unibase.Set( K + i - 1, btmp[ K-1 ] );
    }
  }
}

template < nbases_t K >
void KAdjGraph< K >::RemoveOrphans(int maxOrphanLength) {

  vec<int> uToRemove;
  For_( UnipathInfo, u, unipathInfos ) {
    int ulen = u->len;
    if (ulen > maxOrphanLength)  // Unipath is too long
      continue;

    // Find neighbouring unipaths
    vec<kseq_id_t> neighbors;
    GetSuccs(u->last_kseq, neighbors);
    if (neighbors.empty())
      continue;
    GetPreds(u->first_kseq, neighbors);
    if (neighbors.empty())
      continue;
    uToRemove.push_back(kseq2upath[u->last_kseq]);
  }

  cout << "Removing " << uToRemove.isize() << " orphan unipaths" << endl;
  for (int i = 0; i < uToRemove.isize(); i++) {
    RemoveUnipathKAdjs(uToRemove[i]);
  }

  CheckKAdjGraphSanity();
  ComputeUnipathInfos();
}


template < nbases_t K >
void KAdjGraph< K >::ShaveGraph(int minLength, int maxShaveLength, bool singleHairsOnly) {

  ForceAssertLe(maxShaveLength, minLength);

  vec<int> uToRemove;
  // iterate over all unipaths
  For_( UnipathInfo, u, unipathInfos ) {
    int ulen = u->len;
    if (ulen < minLength)  // Unipath is too short
      continue;

    enum {forward = 0, backward = 1};
    for (int dir = 0; dir < 2; dir++ ) {
      vec<kseq_id_t> neighbors;
      if (dir == forward)
	GetSuccs(u->last_kseq, neighbors);
      else 
      	GetPreds(u->first_kseq, neighbors);

      vec<int> longUnipaths;
      vec<int> shortUnipaths;
      for (int i = 0; i < neighbors.isize(); i++) {
	int neighborIndex = kseq2upath[neighbors[i]];
	UnipathInfo& neighborUni = unipathInfos[neighborIndex];
	if (neighborUni.len >= minLength)
	  longUnipaths.push_back(neighborIndex);
	else
	  shortUnipaths.push_back(neighborIndex);
      }
      
      // Make sure we have at least 1 short edge to shave
      if (shortUnipaths.empty() || longUnipaths.empty())
	continue;
      
      for (int i = 0; i < shortUnipaths.isize(); i++) {
	set<int> seen;
	queue<int> unseen;
	int unipathIndex = shortUnipaths[i];
	unseen.push(unipathIndex);
	int hairLength = 0;
	bool loop = false;
	while(!unseen.empty() && hairLength < maxShaveLength) {
	  int currentIndex = unseen.front();
	  unseen.pop();
	  seen.insert(currentIndex);

	  UnipathInfo& currentUni = unipathInfos[currentIndex];
	  
	  hairLength += currentUni.len;

	  vec<kseq_id_t> neighbors;
	  if (dir == forward)
	    GetSuccs(currentUni.last_kseq, neighbors);
	  else 
	    GetPreds(currentUni.first_kseq, neighbors);
	  
	  for (int j = 0; j < neighbors.isize(); j++) {
	    int next = neighbors[j];
	    if (seen.find(next) != seen.end()) {
	      unseen.push(next);
	    } else
	      loop = True;
	  }
	  if (singleHairsOnly && !neighbors.empty())
	    break;


	}
	if (unseen.empty() && hairLength < maxShaveLength && !loop)
	  uToRemove.insert(uToRemove.end(), seen.begin(), seen.end());
      }
    }	      
  }
  UniqueSort(uToRemove);

  cout << "Removing " << uToRemove.isize() << " shaved unipaths" << endl;

  // Disconnect and remove shaved unipaths
  for (int i = 0; i < uToRemove.isize(); i++) {
    DisconnectUnipath(uToRemove[i]);
    RemoveUnipathKAdjs(uToRemove[i]);
  }

  for (int i = 0; i < uToRemove.isize(); i++)
    cout << "Removed Unipath " << uToRemove[i] << ", size "
	 << unipathInfos[uToRemove[i]].len << " kmers\n";
  
  //  CheckKAdjGraphSanity();
  ComputeUnipathInfos();
}


template < nbases_t K >
void KAdjGraph< K >::RemoveUnipathKAdjs(int unipath_index) {
  const UnipathInfo& unipath = unipathInfos[unipath_index];
  kseq_id_t current = unipath.first_kseq;
  for (int i = 0; i < unipath.len - 1; i++) {
    kseq_id_t next = GetSoleSucc(current);
    RemoveKAdj(current, next);
    current = next;
  }
}

template < nbases_t K >
void KAdjGraph< K >::DisconnectUnipath(int unipath_index) {
  const UnipathInfo& unipath = unipathInfos[unipath_index];

  // cout << "\nDisconnecting\n";
  kseq_id_t first = unipath.first_kseq;
  kseq_id_t last = unipath.last_kseq;
  // PRINT2(unipath_index, kseq2upath[kseq2rc[first]]);
  //  PRINT2(unipathInfos[unipath_index].len, unipathInfos[kseq2upath[kseq2rc[first]]].len);
  //  PRINT2(first, kseq2rc[first]);
  //  PRINT2(last, kseq2rc[last]);

  vec<kseq_id_t> preds;
  GetPreds(first, preds);
  //  cout << "Start [" << preds.isize() << "]";
  for (int i = 0; i < preds.isize(); i++) {
    kseq_id_t pred = preds[i];
    //    cout << "{" << pred << "," << first << "}\n";
    RemoveKAdj(pred, first);
  }  
  vec<kseq_id_t> succs;
  GetSuccs(last, succs);
  //  cout << "End [" << succs.isize() << "]";
  for (int i = 0; i < succs.isize(); i++) {
    kseq_id_t succ = succs[i];
    //    cout << "{" << last << "," << succ << "}\n";
    RemoveKAdj(last, succ);
  } 
  //  cout << "\n";
}

template <nbases_t K>
void KAdjGraph< K >::GetKAdjs( vec< kmer< K+1 > > & kadjs ) {

  int adjCount = 0;
  vec<kseq_id_t> succs;
  basevector b1(K + 1);
  basevector b2(K);
  for (int i = 0; i < kseqs.isize(); i++) {
    GetSuccs(i, succs);
    for (int j = 0; j < succs.isize(); j++) {
      kseqs[i].GetBasevector(b1);
      kseqs[succs[j]].GetBasevector(b2);
      b1.resize(K + 1);
      b1.Set(K, b2[K - 1]);
      kadjs.push_back(b1);
    }
  }
}

template <nbases_t K>
void KAdjGraph< K >::GetKAdjs( vec< kmer_with_count< K+1 > > & kadjs ) {

  int adjCount = 0;
  vec<kseq_id_t> succs;
  basevector b1(K + 1);
  basevector b2(K);
  for (int i = 0; i < kseqs.isize(); i++) {
    GetSuccs(i, succs);
    for (int j = 0; j < succs.isize(); j++) {
      kseqs[i].GetBasevector(b1);
      kseqs[succs[j]].GetBasevector(b2);
      b1.resize(K + 1);
      b1.Set(K, b2[K - 1]);
      b1.Canonicalize();
      kadjs.push(b1,1);
    }
  }
  UniqueSort(kadjs);
}

#define INSTANTIATE(K, dummy) template class KAdjGraph<K>

FOR_ALL_K_WITH_K_PLUS_1(INSTANTIATE,unused);

