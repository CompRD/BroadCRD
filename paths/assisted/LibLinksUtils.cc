/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "PairsManager.h"
#include "lookup/LookAlign.h"
#include "graph/Digraph.h"
#include "paths/Alignlet.h"
#include "paths/assisted/CProx.h"
#include "paths/assisted/LibLinksUtils.h"
#include "paths/assisted/Proxgraph.h"

// #define MASSIVE_LOG

/**
 * FillLibLinks
 *
 * Notice that each pair will contribute two links, as in this
 * example. Say reads r1 and r2 align fw on contigs v1 and v2. Then,
 * there are two possible ways ti link v1 and v2, as follows:
 *
 *            r1          r2
 *            ->..........<-
 *      ---------->   <----------
 *          v1            v2
 *  
 *            r2          r1
 *            ->..........<-
 *      ---------->   <----------
 *          v2            v1
 *
 * Ie, the pair implies a link between v1[+] and v2[-], and one
 * between v2[+] and v1[-].
 */
void FillLibLinks( const int MAX_CN,
		   const PairsManager &pairs,
		   const vec<bool> &isrep,
		   const vec<int> &cns,
		   const vec<alignlet> &aligns,
		   const vec<int> &idx,
		   vec<CLibLinks> &lib_links )
{
  lib_links.clear( );
  
  // Resize links (links[ii] = all links for library ii in pairs).
  lib_links.resize( pairs.nLibraries( ) );
  for (size_t ii=0; ii<pairs.nLibraries( ); ii++) {
    int gsep = pairs.getLibrarySep( ii );
    int gdev = pairs.getLibrarySD( ii );
    String name = pairs.getLibraryName( ii );
    lib_links[ii].SetLib( gsep, gdev, name );
  }
  
  // Count events.
  vec<size_t> n_events( pairs.nLibraries( ), 0 );
  for (size_t pid=0; pid<pairs.nPairs( ); pid++) {
    int r1 = pairs.ID1( pid );
    int r2 = pairs.ID2( pid );

    int idx1 = idx[r1];
    int idx2 = idx[r2];
    if ( idx1 < 0 || idx2 < 0 ) continue;
    if ( aligns[idx1].TargetId( ) == aligns[idx2].TargetId( ) ) continue;
    if ( isrep[ aligns[idx1].TargetId( ) ] ) continue;
    if ( isrep[ aligns[idx2].TargetId( ) ] ) continue;
    if ( cns[ aligns[idx1].TargetId( ) ] > MAX_CN ) continue;
    if ( cns[ aligns[idx2].TargetId( ) ] > MAX_CN ) continue;
    
    n_events[pairs.libraryID( pid )] += 2;
  }

  // Reserve memory.
  for (size_t ii=0; ii<pairs.nLibraries( ); ii++)
    lib_links[ii].Reserve( n_events[ii] );
  
  // Generate links.
  for (size_t pid=0; pid<pairs.nPairs( ); pid++) {
    int r1 = pairs.ID1( pid );
    int r2 = pairs.ID2( pid );

    int idx1 = idx[r1];
    int idx2 = idx[r2];
    if ( idx1 < 0 || idx2 < 0 ) continue;
    if ( aligns[idx1].TargetId( ) == aligns[idx2].TargetId( ) ) continue;
    if ( isrep[ aligns[idx1].TargetId( ) ] ) continue;
    if ( isrep[ aligns[idx2].TargetId( ) ] ) continue;
    if ( cns[ aligns[idx1].TargetId( ) ] > MAX_CN ) continue;
    if ( cns[ aligns[idx2].TargetId( ) ] > MAX_CN ) continue;
    
    int lib_id = pairs.libraryID( pid );
    int t1 = aligns[idx1].TargetId( );
    int t2 = aligns[idx2].TargetId( );
    int tlen1 = aligns[idx1].TargetLength( );
    int tlen2 = aligns[idx2].TargetLength( );
    int p1 = aligns[idx1].pos2( );
    int P1 = aligns[idx1].Pos2( );
    int p2 = aligns[idx2].pos2( );
    int P2 = aligns[idx2].Pos2( );
    bool fw1 = aligns[idx1].Fw1( );
    bool fw2 = aligns[idx2].Fw1( );

    int sep = 0;
    pair<int,int> set1;
    pair<int,int> set2;
    if ( fw1 && fw2 ) {
      sep = pairs.sep( pid ) - ( tlen1 - P1 ) - ( tlen2 - P2 );
      set1 = make_pair( 2*t1,  1 + 2*t2 );
      set2 = make_pair( 2*t2,  1 + 2*t1 );
    }
    else if ( fw1 && !fw2 ) {
      sep = pairs.sep( pid ) - ( tlen1 - P1 ) - p2;
      set1 = make_pair(     2*t1,      2*t2 );
      set2 = make_pair( 1 + 2*t2,  1 + 2*t1 );
    }
    else if ( ( !fw1 ) && fw2 ) {
      sep = pairs.sep( pid ) - ( tlen2 - P2 ) - p1;
      set1 = make_pair( 1 + 2*t1,  1 + 2*t2 );
      set2 = make_pair(     2*t2,      2*t1 );
    }
    else  {
      sep = pairs.sep( pid ) - ( tlen2 - P2 ) - p1;
      set1 = make_pair( 1 + 2*t1,  2*t2 );
      set2 = make_pair( 1 + 2*t2,  2*t1 );
    }
    
    lib_links[lib_id].AddLink( CRawLink( set1.first, set1.second, sep ) );
    lib_links[lib_id].AddLink( CRawLink( set2.first, set2.second, sep ) );

#ifdef MASSIVE_LOG
    cout << set1.first / 2 << ( set1.first % 2 == 0 ? "[+]" : "[-]" )
	 << " => "
	 << set1.second / 2 << ( set1.second % 2 == 0 ? "[+]" : "[-]" )
	 << "   p" << pid << " r(" << r1 << ", " << r2 << ")"
	 << "   " << ( fw1 ? "->" : "<-" )
	 << "[" << p1 << ", " << P1 << ")_" << tlen1
	 << "   " << ( fw2 ? "->" : "<-" )
	 << "[" << p2 << ", " << P2 << ")_" << tlen2
	 << "\n";

    cout << set2.first / 2 << ( set2.first % 2 == 0 ? "[+]" : "[-]" )
	 << " => "
	 << set2.second / 2 << ( set2.second % 2 == 0 ? "[+]" : "[-]" )
	 << "   p" << pid << " r(" << r2 << ", " << r1 << ")"
	 << "   " << ( fw2 ? "->" : "<-" )
	 << "[" << p2 << ", " << P2 << ")_" << tlen2
	 << "   " << ( fw1 ? "->" : "<-" )
	 << "[" << p1 << ", " << P1 << ")_" << tlen2
	 << "\n";

    cout << "\n";
#endif
    
  }
  
  // Sort raw links.
  for (size_t lib_id=0; lib_id<pairs.nLibraries( ); lib_id++) {
    lib_links[lib_id].SortLinks( );
    lib_links[lib_id].BuildFirst( cns.size( ) );
  }
  
}

/**
 * FillRefLinks
 *
 * Notice that each pair will contribute two links, as in FillLibLinks
 * above.
 */
void FillRefLinks( const int MIN_REF_GAP,
		   const int MAX_REF_GAP,
		   const int MAX_CN,
		   const vec<look_align> &hits,
		   const vec<int> &cns,
		   CLibLinks &ref_links )
{
  ref_links.Clear( );
  ref_links.SetLib( 0, 0, "ref_based" );

  vec<int> idx( cns.size( ), -1 );
  for (int ii=0; ii<hits.isize( ); ii++) {
    if ( cns[ hits[ii].query_id ] > MAX_CN ) continue;
    if ( idx[ hits[ii].query_id ] == -1 ) idx[ hits[ii].query_id ] = ii;
    else idx[ hits[ii].query_id ] = -2;
  }

  vec<bool> select( hits.size( ), false );
  for(size_t ii=0; ii<idx.size( ); ii++)
    if ( idx[ii] > -1 )
      select[ idx[ii] ] = true;
  
  for (size_t hid1=0; hid1<hits.size( ); hid1++) {
    if ( ! select[hid1] ) continue;

    int q1 = hits[hid1].query_id;
    int t1 = hits[hid1].target_id;
    int fw1 = hits[hid1].IsQueryFW( );
    int P1 = hits[hid1].Pos2( );

    for (size_t hid2=hid1+1; hid2<hits.size( ); hid2++) {
      if ( ! select[hid2] ) continue;

      int q2 = hits[hid2].query_id;
      if ( q1 == q2 ) continue;

      int t2 = hits[hid2].target_id;
      int fw2 = hits[hid2].IsQueryFW( );
      int p2 = hits[hid2].pos2( );
      if ( t1 != t2 ) continue;

      int sep = p2 - P1;
      if ( sep < MIN_REF_GAP ) continue;
      if ( sep > MAX_REF_GAP ) break; 
      
      int id1 = 2 * q1 + ( fw1 ? 0 : 1 );
      int id2 = 2 * q2 + ( fw2 ? 0 : 1 );
      int id1rc = 2 * q1 + ( fw1 ? 1 : 0 );
      int id2rc = 2 * q2 + ( fw2 ? 1 : 0 );
      ref_links.AddLink( CRawLink( id1, id2, sep ) );
      ref_links.AddLink( CRawLink( id2rc, id1rc, sep ) );
    }
  }
  
  ref_links.SortLinks( );
  ref_links.BuildFirst( cns.size( ) );

}

/**
 * LibLinksToDigraphE
 */
void LibLinksToDigraphE( const int n_contigs,
			 const vec<CLibLinks> *jlinks,
			 const vec<CLibLinks> *Jlinks,
			 const CLibLinks *rlinks,
			 const CProxInfo *info,
			 proxgraph &graphE )
{

  // The core descriptors for a digraphE.
  vec< vec<int> > from( 2 * n_contigs );
  vec< vec<int> > to( 2 * n_contigs );
  vec< vec<int> > fromEO( 2 * n_contigs );
  vec< vec<int> > toEO( 2 * n_contigs );
  vec<CProx> edges;
  
  // Count libraries.
  int n_jlibs = jlinks ? jlinks->isize( ) : 0;
  int n_Jlibs = Jlinks ? Jlinks->isize( ) : 0;

  // Find all mates, tag generate maps to CProxInfo tags.
  vec<int> jtags;
  vec<int> Jtags;
  vec< pair<int,int> > mates;
  {
    vec< vec< pair<int,int> > > jmates;
    vec< vec< pair<int,int> > > Jmates;
    vec< pair<int,int> > rmates;
    if ( jlinks ) {
      jtags.resize( jlinks->size( ) );
      jmates.resize( jlinks->size( ) );
      for (int ii=0; ii<jlinks->isize( ); ii++) {
	jtags[ii] = info->TagId( 0, ii );
	(*jlinks)[ii].AllMates( jmates[ii] );
      }
    }
    if ( Jlinks ) {
      Jtags.resize( Jlinks->size( ) );
      Jmates.resize( Jlinks->size( ) );
      for (int ii=0; ii<Jlinks->isize( ); ii++) {
	Jtags[ii] = info->TagId( 1, ii );
	(*Jlinks)[ii].AllMates( Jmates[ii] );
      }
    }
    if ( rlinks )
      rlinks->AllMates( rmates );

    size_t ntot = 0;
    for (size_t ii=0; ii<jmates.size( ); ii++)
      ntot += jmates[ii].size( );
    for (size_t ii=0; ii<Jmates.size( ); ii++)
      ntot += Jmates[ii].size( );
    ntot += rmates.size( );
    mates.reserve( ntot );
    
    for (size_t ii=0; ii<jmates.size( ); ii++)
      copy( jmates[ii].begin( ), jmates[ii].end( ), back_inserter( mates ) );
    for (size_t ii=0; ii<Jmates.size( ); ii++)
      copy( Jmates[ii].begin( ), Jmates[ii].end( ), back_inserter( mates ) );
    copy( rmates.begin( ), rmates.end( ), back_inserter( mates ) );
    sort( mates.begin( ), mates.end( ) );
    mates.erase( unique( mates.begin( ), mates.end( ) ), mates.end( ) );
  }
  
  // Reserve memory.
  edges.reserve( mates.size( ) );
  
  // Loop over all mates.
  for (size_t mate_id=0; mate_id<mates.size( ); mate_id++) {
    int v = mates[mate_id].first;
    int w = mates[mate_id].second;
    
    CProx prox( info );
    for (int ii=0; ii<jtags.isize( ); ii++)
      prox.AddLinkGaps( jtags[ii], (*jlinks)[ii].AllSeps( v, w ) );
    for (int ii=0; ii<Jtags.isize( ); ii++)
      prox.AddLinkGaps( Jtags[ii], (*Jlinks)[ii].AllSeps( v, w ) );
    if ( rlinks ) {
      vec<int> seps = rlinks->AllSeps( v, w );
      if ( seps.size( ) > 0 ) prox.AddRefGap( seps[0] );
    }
    prox.ComputeLinkGap( );
    
    // Decide if the two edges are joined (by links and/or ref).
    int ref_gap = prox.RefGap( );
    pair<int,int> link_gap = prox.LinkGap( );
    if ( ref_gap == INT_MAX && link_gap.first == INT_MAX ) continue;

    // Join edges.
    from[v].push_back( w );
    fromEO[v].push_back( edges.size( ) );
    to[w].push_back( v );
    toEO[w].push_back( edges.size( ) );
    edges.push_back( prox );
  }

  // Build digraphE (make sure proxgraph is symmetric).
  graphE.Initialize( from, to, edges, toEO, fromEO );
  ForceAssert( IsSymmetric( graphE ) );
  
}

/**
 * LibLinksToDigraphE
 */
void AddBundlesToDigraphE(const digraphE<CLinkBundle> &bundles,
			  const CProxInfo *info,
			  proxgraph &graphE )
{
  vec<int> to_left;
  vec<int> to_right;
  bundles.ToLeft( to_left );
  bundles.ToRight( to_right );
  int minLinks = info->MinLinks();
  int minGap = info->MinGap();
  int maxGap = info->MaxGap();
    
  for (int edge_id=0; edge_id < bundles.EdgeObjectCount( ); ++edge_id) {
    CLinkBundle bundle = bundles.EdgeObject(edge_id);

    // ignore edges which don't meet our standards
    if (minLinks > 0 && bundle.Weight() < minLinks) continue;
    if (bundle.Sep() < minGap || bundle.Sep() > maxGap) continue;

    int v1 = to_left[edge_id];
    int v2 = to_right[edge_id];
    vec<int> edges = graphE.EdgesBetween(v1, v2);
    if (edges.size() > 0) {
      int e = edges[0];
      CProx prox = graphE.EdgeObject(e);
      prox.AddBundle(bundle);
      graphE.SetEdgeObject(e, prox);
    } else {
      CProx prox(info);
      prox.AddBundle(bundle);
      graphE.AddEdge(v1, v2, prox);
    }
  }
}

    
