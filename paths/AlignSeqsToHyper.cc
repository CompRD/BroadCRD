///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "paths/AlignSeqsToHyper.h"

#include "Alignment.h"
#include "Basevector.h"
#include "Charvector.h"
#include "Intvector.h"
#include "CoreTools.h"
#include "PrintAlignment.h"
#include "ReadLocationLG.h"
#include "graph/Digraph.h"
#include "lookup/LookAlign.h"
#include "lookup/LookupTableBuilder.h"
#include "lookup/PerfectLookup.h"
#include "pairwise_aligners/PerfectAlignerLG.h"
#include "paths/HyperBasevector.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/SeqOnHyper.h"


// Local class: simple_align
//
// In a simple_align, the offset is the NEGATIVE of the usual offset.  Thus if
// sequence 1 starts at position p on sequence 2, then the offset is p.

class simple_align {

     public:

     simple_align( ) { }
     simple_align( int id1, int id2, int offset, Bool rc1 )
          : id1(id1), id2(id2), offset(offset), rc1(rc1) { }

     int id1, id2;
     int offset;
     Bool rc1;

};



typedef multimap<unsigned int,size_t> EdgeLenToIdMap;

int FillInCachedAligns( const HyperBasevector& hb, const unsigned int minEdgeLenToCache,
                        const vecbasevector& cachedEdges, 
                        const vec<partial_perfect_align>& cachedPartialAligns, 
                        vec<partial_perfect_align>& newPartialAligns, 
                        vec<Bool>& edgeWasCached ) 
{
  EdgeLenToIdMap cachedEdgeLenMap;
  for ( size_t e = 0; e < cachedEdges.size(); ++e ) {
    if ( cachedEdges[e].size() < minEdgeLenToCache ) continue;
    cachedEdgeLenMap.insert( make_pair(cachedEdges[e].size(), e) );
  }
  
  vec< vec<int> > edgeToAlignsMap( cachedEdges.size() );
  for ( unsigned int i = 0; i < cachedPartialAligns.size(); ++i ) {
    int edge_id = cachedPartialAligns[i].edge_id;
    if ( cachedEdges[edge_id].size() < minEdgeLenToCache ) continue;
    edgeToAlignsMap[edge_id].push_back(i);
  }
  
  int numCachedEdges = 0;
  for ( int e = 0; e < hb.EdgeObjectCount(); ++e ) {
    if ( hb.EdgeObject(e).size() < minEdgeLenToCache ) continue;
    pair<EdgeLenToIdMap::iterator,EdgeLenToIdMap::iterator> rng_matchingSize;
    rng_matchingSize = cachedEdgeLenMap.equal_range( hb.EdgeObject(e).size() );
    while ( rng_matchingSize.first != rng_matchingSize.second ) {
      int cachedEdgeIdx = rng_matchingSize.first->second;
      if ( cachedEdges[cachedEdgeIdx] == hb.EdgeObject(e) ) {
        edgeWasCached[e] = True;
        ++numCachedEdges;
        vec<int>& alignIdxs = edgeToAlignsMap[cachedEdgeIdx];
        for ( unsigned int i = 0; i < alignIdxs.size(); ++i ) {
          newPartialAligns.push_back( cachedPartialAligns[alignIdxs[i]] );
          newPartialAligns.back().edge_id = e;
        }
        break;
      }
      ++rng_matchingSize.first;
    }
  }
  return numCachedEdges;
}

// Finds all prefect alignments to the HyperKmerPath - PerfectLookup Method (Fastest)
// Note that hb must be HyperBasevector(h)
void AlignSeqsToHyper( const HyperKmerPath& h, const HyperBasevector& hb,
		       const vecbasevector& seqs, const String& sub_dir,
		       vec<CompressedSeqOnHyper>& csaligns, Bool verbose,
                       vecbasevector* p_cachedEdges, vec<partial_perfect_align>* p_cachedAligns )
{
  if ( verbose )
    cout << Date( ) << ": Aligning sequences to HyperKmerPath " << endl;

  vec<partial_perfect_align> aligns;

  vec<Bool> edgeWasCached( hb.EdgeObjectCount(), False );
  const unsigned int minEdgeLenToCache = h.K()+1;

  int numCachedEdges = 0;
  if ( p_cachedEdges && p_cachedAligns ) {
    numCachedEdges = FillInCachedAligns( hb, minEdgeLenToCache,
                                         *p_cachedEdges, *p_cachedAligns, 
                                         aligns, edgeWasCached );
    if ( verbose )
      cout << Date( ) << ": Restored " << aligns.size() << " aligns to "
	   << numCachedEdges << " edges from cache" << endl;
  }

  if ( hb.EdgeObjectCount() > numCachedEdges ) {
    
    if ( verbose ) {
      cout << Date( ) << ": Creating LookupTable for "
	   << hb.EdgeObjectCount() - numCachedEdges << " edges" << endl;
      PrintMemUsage( );
    }
    
    // Build vecbasevector of edges to align to
    vecbasevector edgesToAlign;
    vec<int> edgeIdMap;
    for ( int i = 0; i < hb.EdgeObjectCount( ); i++ ) {
      if ( ! edgeWasCached[i] ) {
	edgesToAlign.push_back_reserve( hb.EdgeObject(i) );
	edgeIdMap.push_back(i);
      }
    }
    
    // Build Lookup Table
    String table_name = sub_dir + "/" + "edges.lookup";
    const unsigned int K = 12;
    LookupTableBuilder(edgesToAlign, table_name, False, K);
    
    Destroy( edgesToAlign );
    
    // Perform Alignment
    const int chunk_size = 10000000;    // 10,000,000 reads at a time
    size_t chunks = (seqs.size() + chunk_size - 1) / chunk_size;
    
    if ( verbose ) {
      cout << Date( ) << ": Aligning in " << chunks << " chunk(s)" << endl;
      PrintMemUsage( );
    }
    
    for (size_t i = 0; i < chunks; ++i) {
      vec<look_align> lookAligns;
      using std::min;
      PerfectLookup(K,  seqs, table_name, lookAligns, FW_OR_RC,
		    i * chunk_size, min( (i+1) * chunk_size, seqs.size() ) - 1,
		    False, h.K() - 1);
      
      if ( verbose ) {
	cout << Date( ) << ": Chunk " << i+1 << " : Found "
	     << lookAligns.size() << " potential aligns " << endl;
	PrintMemUsage( );
      }
      
      for ( unsigned int j = 0; j < lookAligns.size(); ++j ) {
	look_align& la = lookAligns[j];
	aligns.push_back( partial_perfect_align() );
	int read_id = la.query_id;
	aligns.back().read_id = ( la.rc1 ? -read_id-1 : read_id );
	aligns.back().edge_id = edgeIdMap[la.target_id];
	aligns.back().read_length = la.query_length;
	aligns.back().pos1 = la.pos1();
	aligns.back().pos2 = la.pos2();
	aligns.back().length = la.a.Lengths(0);
      }
    }
    
    // Remove lookup table.
    Remove( table_name );
    
  }
  
  // Try to extend partial alignments, compress all full alignments.
  if ( verbose ) {
    cout << Date( ) << ": Found " << aligns.size() << " partial aligns " << endl;
    PrintMemUsage( );
  }
  
  vecbasevector allEdges;
  for ( int i = 0; i < hb.EdgeObjectCount( ); i++ )
    allEdges.push_back_reserve( hb.EdgeObject(i) );
  
  csaligns.clear( );
  ExtendAndCompressAlignments( h, hb, seqs, allEdges, aligns, csaligns );
  if ( verbose ) {
    cout << Date( ) << ": Extended and compressed " << csaligns.size() << " aligns " << endl;
    PrintMemUsage( );
  }
  
  // Cache edges aligned to and the resulting aligns.
  if ( p_cachedEdges ) {
    p_cachedEdges->swap( allEdges );
  }
  Destroy( allEdges );
  
  if ( p_cachedAligns )
    p_cachedAligns->swap( aligns );
  Destroy( aligns );
  
  if ( verbose ) {
    cout << Date( ) << ": Found " << csaligns.size() << " compressed Aligns" << endl;
    PrintMemUsage( );
  }
}



// Finds alignments to the HyperBasevector - PerfectAlignerLG Method (Slow)
// This version returns the aligns in the SeqOnHyper form
void AlignSeqsToHyper( const HyperBasevector& hb, const vecbasevector& seqs,
     vec<SeqOnHyper>& saligns )
{    vec<alignment_plus> aligns;
     int K = hb.K( );
     vecbasevector seqs_edges = seqs;
     for ( int i = 0; i < hb.EdgeObjectCount( ); i++ )
          seqs_edges.push_back_reserve( hb.EdgeObject(i) );
     PerfectAlignerLG pal( hb.K( ), PerfectAlignerLG::findProperOnly );
     pal.Align( seqs_edges, aligns, 1, seqs.size( ) );
     AlignSeqsToHyper( hb, seqs, aligns, saligns );    }


// Finds alignments to the HyperBasevector - PerfectAlignerLG Method (Slow)
// This version returns the aligns in the CompressedSeqOnHyper form
void AlignSeqsToHyper( const HyperBasevector& hb, const vecbasevector& seqs,
     vec<CompressedSeqOnHyper>& csaligns )
{    vec<alignment_plus> aligns;
     int K = hb.K( );
     vecbasevector seqs_edges = seqs;
     for ( int i = 0; i < hb.EdgeObjectCount( ); i++ )
          seqs_edges.push_back_reserve( hb.EdgeObject(i) );
     PerfectAlignerLG pal( hb.K( ), PerfectAlignerLG::findProperOnly );
     pal.Align( seqs_edges, aligns, 1, seqs.size( ) );
     AlignSeqsToHyper( hb, seqs, aligns, csaligns );    }


// Finds alignments to the HyperBasevector - Alternative Method (Faster)
void AlignSeqsToHyper( const HyperKmerPath& h, const HyperBasevector& hb,
     const vecbasevector& seqs, const vecbasevector& unibases, 
     const vecKmerPath& unipaths, const vec<ReadLocationLG>& ulocs, 
     const VecULongVec& ulocs_indexr, const String& sub_dir, const String& tmprun,
     vec<alignment_plus>& aligns )
{
     // Define some variables.

     String brun_dir = sub_dir + "/" + tmprun;
     int N = seqs.size( ), K = h.K( );
     String KS = ToString(K);

     // Align unipaths to the edges.

     int nuni = unipaths.size( );
     vecbasevector unipaths_edges;
     for ( int i = 0; i < nuni; i++ )
          unipaths_edges.push_back_reserve( unibases[i] );
     for ( int i = 0; i < hb.EdgeObjectCount( ); i++ )
          unipaths_edges.push_back_reserve( hb.EdgeObject(i) );
     PerfectAlignerLG pal( K, PerfectAlignerLG::findProperOnly );
     vec<alignment_plus> uealigns;
     pal.Align( unipaths_edges, uealigns, 1, unipaths.size( ) );
     vec< vec<int> > uealigns_index(nuni);
     for ( int i = 0; i < uealigns.isize( ); i++ )
     {    const alignment_plus& ap = uealigns[i];
          uealigns_index[ ap.Id1( ) ].push_back(i);    }

     // Find the reads that are placed on a single unipath and infer their
     // locations on the edges.

     vec<Bool> translated( N, False );
     vec<simple_align> ALIGNS;
     vec<simple_align> simpaligns;
     for ( int id = 0; id < N; id++ )
     {    if ( ulocs_indexr[id].size( ) != 2 ) continue;
          Bool placed = False, broken = False;
          simpaligns.clear( );
          for ( ULongVec::size_type j = 0; j < ulocs_indexr[id].size( ); j++ )
          {    const ReadLocationLG& rl = ulocs[ ulocs_indexr[id][j] ];
               int u = rl.Contig( );
               for ( int x = 0; x < uealigns_index[u].isize( ); x++ )
               {    const alignment_plus& ap = uealigns[ uealigns_index[u][x] ];
                    if ( ap.Rc2( ) ) continue;
                    int e = ap.Id2( );
                    if ( IntervalOverlap( rl.Start( ), 
                         rl.Start( ) + seqs[id].isize( ), 
                         ap.a.pos1( ), ap.a.Pos1( ) ) == 0 )
                    {    continue;    }
                    if ( !( ap.a.pos1( ) <= rl.Start( ) ) ) broken = True;
                    if ( !( rl.Start( ) + seqs[id].isize( ) <= ap.a.Pos1( ) ) )
                         broken = True;
                    if (broken) continue;
                    if ( !broken ) placed = True;
                    int read_start_on_edge 
                         = rl.Start( ) + ap.a.pos2( ) - ap.a.pos1( );
                    simpaligns.push_back( simple_align( id, e, 
                         read_start_on_edge, rl.Rc( ) ) );    }    }
          if ( placed && !broken ) 
          {    translated[id] = True;
               ALIGNS.append(simpaligns);    }    }
     Destroy(simpaligns);
     
     // Form a vecbasevector, containing one entry per isomorphism class of 
     // residual reads (where a read is regarded as isomorphic to its reverse
     // complement), and parallel VecIntVec and vecBoolvector structures that
     // record the reads and orientations corresponding to the isomorphism classes.

     int Nr = N - Sum(translated);
     vec<basevector> reads;
     vec<int> ids;
     vec<Bool> rc;
     reads.reserve(Nr), ids.reserve(Nr), rc.reserve(Nr);
     for ( int i = 0; i < N; i++ )
     {    if ( translated[i] ) continue;
          ids.push_back(i);
          basevector b = seqs[i];
          b.ReverseComplement( );
          if ( b < seqs[i] )
          {    reads.push_back(b);
               rc.push_back(True);    }
          else
          {    reads.push_back( seqs[i] );
               rc.push_back(False);    }     }
     SortSync( reads, ids, rc );
     vecbasevector ereads;
     VecIntVec eids;
     vecBoolvector erc;
     if ( true )
     {
         IntVec cids;
         Boolvector crc;
         for ( int i = 0; i < Nr; i++ )
         {    int j;
              for ( j = i + 1; j < Nr; j++ )
                   if ( reads[j] != reads[i] ) break;
              cids.clear( ), crc.clear( );
              for ( int k = i; k < j; k++ )
              {    cids.push_back( ids[k] );
                   crc.push_back( rc[k] );    }
              ereads.push_back_reserve( reads[i] );
              eids.push_back_reserve(cids);
              erc.push_back_reserve(crc);
              i = j - 1;    }
     }

     // Generate KmerPaths for the equivalence classes of residual reads and the
     // edges.


     vecbasevector ereads_plus_edges = ereads;
     for ( int i = 0; i < hb.EdgeObjectCount( ); i++ ) {
       ereads_plus_edges.push_back_reserve( hb.EdgeObject(i) );
     }

     vecKmerPath paths;
     ReadsToPathsCoreY( ereads_plus_edges, K, paths );

     // vecbasevector ereads_plus_edges = ereads;
     // for ( int i = 0; i < hb.EdgeObjectCount( ); i++ )
     //      ereads_plus_edges.push_back_reserve( hb.EdgeObject(i) );
     // ereads_plus_edges.WriteAll( brun_dir + "/reads.fastb" );
     // ReadsToPathsCore( sub_dir, ".", tmprun, "", K, False, True, "",
     //      False, True, "/dev/null" );

     // Find the alignments between the remaining reads and the edges.

     // vecKmerPath paths( brun_dir + "/reads.paths.k" + KS );
     vecKmerPath edges;
     for ( int i = 0; i < h.EdgeObjectCount( ); i++ )
          edges.push_back_reserve( paths[ ereads.size( ) + i ] );
     vec<big_tagged_rpint> edgedb;
     CreateDatabase( edges, edgedb );
     for ( size_t i = 0; i < ereads.size( ); i++ )
     {    for ( int pass = 1; pass <= 2; pass++ )
          {    Bool rc = ( pass == 2 );
               KmerPath p = paths[i];
               if (rc) p.Reverse( );
               vec< pair<int,longlong> > e_offsets;   // (read to edge)
               for ( int j = 0; j < p.NSegments( ); j++ )
               {    vec<longlong> places;
                    Contains( edgedb, p.Segment(j), places );    
                    for ( int m = 0; m < places.isize( ); m++ )
                    {    const big_tagged_rpint& t = edgedb[ places[m] ];
                         int e = t.PathId( );
                         if ( !ProperOverlapExt( p, edges[e], j, t.PathPos( ) ) )
                              continue;
                         longlong o = t.Start( ) - p.Start(j);
                         for ( int x = 0; x < j; x++ )
                              o += p.Length(x);
                         for ( int x = 0; x < t.PathPos( ); x++ )
                              o -= edges[e].Length(x);
                         e_offsets.push_back( make_pair( e, o ) );    }    }    
               UniqueSort(e_offsets);    
               for ( int j = 0; j < e_offsets.isize( ); j++ )
               {    int e = e_offsets[j].first, o = -e_offsets[j].second;
                    for ( IntVec::size_type m = 0; m < eids[i].size( ); m++ )
                    {    int id = eids[i][m];
                         Bool isrc = erc[i][m] ^ rc;
                         ALIGNS.push_back( 
                              simple_align( id, e, o, isrc ) );    }    }    }    }
     
     // Generate alignment_plus objects.  This is stupid.

     aligns.clear( );
     aligns.reserve( ALIGNS.size( ) );
     for ( int i = 0; i < ALIGNS.isize( ); i++ )
     {    const simple_align& s = ALIGNS[i];
          int len1 = seqs[s.id1].size( ), len2 = hb.EdgeLength(s.id2);
          int pos1, pos2;
          int offset = -s.offset;
          if ( offset < 0 )
          {    pos1 = 0;
               pos2 = -offset;    }
          else
          {    pos1 = offset;
               pos2 = 0;    }
          int overlap = IntervalOverlap( 0, len1, offset, offset + len2 );
	  align a;
	  avector<int> g(1, 0), l(1);
          l(0) = overlap;
          a.Set( pos1, pos2, g, l );
          if (s.rc1) a.ReverseThis( len1, len2 );
	  alignment al;
          al.Set( a, 0 );
          alignment_plus ap( s.id1, s.id2, len1, len2, s.rc1, al, 0.0 );
          aligns.push_back(ap);    }    }


// Finds alignments to the HyperBasevector - Extend alignment_plus Method
// Extends previously found alignments along the HyperBasevector.
void AlignSeqsToHyper( const HyperBasevector& hb, const vecbasevector& seqs,
     const vec<alignment_plus>& aligns, vec<SeqOnHyper>* p_saligns,
     vec<CompressedSeqOnHyper>* p_csaligns)
{    vec<int> to_right_vertex( hb.EdgeObjectCount( ) );
     for ( int i = 0; i < hb.N( ); i++ )
     {    for ( int j = 0; j < hb.To(i).isize( ); j++ )
          {    int e = hb.EdgeObjectIndexByIndexTo( i, j );
               to_right_vertex[e] = i;    }    }
     vec<int> to_left_vertex( hb.EdgeObjectCount( ) );
     for ( int i = 0; i < hb.N( ); i++ )
     {    for ( int j = 0; j < hb.From(i).isize( ); j++ )
          {    int e = hb.EdgeObjectIndexByIndexFrom( i, j );
               to_left_vertex[e] = i;    }    }
     int K = hb.K( );
     vec< vec<int> > aligns_index( seqs.size( ) );
     for ( int i = 0; i < aligns.isize( ); i++ )
          aligns_index[ aligns[i].Id1( ) ].push_back(i);
     vec<int> L( hb.EdgeObjectCount( ) );
     for ( int i = 0; i < hb.EdgeObjectCount( ); i++ )
          L[i] = hb.EdgeLength(i);
     if( p_saligns ) {
       p_saligns->clear( );
       p_saligns->reserve( aligns.size( ) );
     }
     if( p_csaligns ) {
       p_csaligns->clear( );
       p_csaligns->reserve( aligns.size( ) );
     }
     for ( size_t id1 = 0; id1 < seqs.size( ); id1++ )
     {    int len1 = seqs[id1].size( );
          for ( int pass = 1; pass <= 2; pass++ )
          {    
               // Find all the alignments of sequence id1, forward if on the first
               // pass, otherwise reverse.

               vec<alignment_plus> places;
               for ( int j = 0; j < aligns_index[id1].isize( ); j++ )
               {    const alignment_plus& ap = aligns[ aligns_index[id1][j] ];
                    if ( pass == 1 && ap.Rc2( ) ) continue;
                    if ( pass == 2 && !ap.Rc2( ) ) continue;

                    // Handle the case where the read goes end-to-end.

                    if ( ap.a.pos1( ) == 0 && ap.a.Pos1( ) == len1 )
                    {    int id2 = ap.Id2( );
                         SeqOnHyper s( False, -1, vec<PartialSeqOnHyper>(1) );
                         PartialSeqOnHyper& p = s.Part(0);
                         p.Setpos1(0);
                         int pos2;
                         if ( !ap.Rc2( ) ) pos2 = ap.a.pos2( );
                         else pos2 = L[id2] - ap.a.Pos2( );
                         p.Setpos2(pos2);
                         p.SetId2(id2);
                         p.SetLen(len1);
                         s.SetId1(id1);
                         s.SetRc1( ap.Rc2( ) );
			 
			 if( p_saligns )
			   p_saligns->push_back(s);
			 if( p_csaligns ) {
			   p_csaligns->push_back( CompressedSeqOnHyper() );
			   p_csaligns->back().Compress( s, hb, &seqs[id1],
							to_right_vertex );
			 }

                         continue;    }

                    // Handle the other cases.

                    if ( pass == 2 ) 
		    {    alignment_plus apx = ap;
                         apx.a.ReverseThis( len1, L[ apx.Id2( ) ] );
                         places.push_back(apx);    }
                    else places.push_back(ap);    }

               // Define a graph G whose vertices are the alignments and for which
               // there is an edge if one alignment follows another.

               if ( places.empty( ) ) continue;
	       vec< vec<int> > from, to;
               from.resize( places.size( ) ), to.resize( places.size( ) );
	       vec<int> start0;
               for ( int j = 0; j < places.isize( ); j++ )
                    if ( places[j].a.pos2( ) == 0 ) start0.push_back(j);
               for ( int j1 = 0; j1 < places.isize( ); j1++ )
               {    const alignment_plus &ap1 = places[j1];
                    if ( ap1.a.Pos2( ) != L[ ap1.Id2( ) ] ) continue;
                    int r = to_right_vertex[ ap1.Id2( ) ];
                    for ( int i2 = 0; i2 < start0.isize( ); i2++ )
                    {    int j2 = start0[i2];
                         if ( j2 == j1 ) continue;
                         const alignment_plus &ap2 = places[j2];
                         if ( r != to_left_vertex[ ap2.Id2( ) ] ) continue;

                         if ( ap1.a.Pos1( ) - ap2.a.pos1( ) != K - 1 ) continue;
                         from[j1].push_back(j2), to[j2].push_back(j1);    }    }
               digraph G( from, to );
	       vec<int> goodsources, goodsinks;
               for ( int j = 0; j < G.N( ); j++ )
               {    if ( places[j].a.pos1( ) == 0 ) goodsources.push_back(j);
                    if ( places[j].a.Pos1( ) == len1 ) goodsinks.push_back(j);    }

               // Go through all the alignments that are sources of the graph
               // and start at the beginning of the read.
          
               for ( int j1 = 0; j1 < goodsources.isize( ); j1++ )
               {    for ( int j2 = 0; j2 < goodsinks.isize( ); j2++ )
                    {    vec< vec<int> > Gpaths;
                         G.AllPaths( goodsources[j1], goodsinks[j2], Gpaths );
                         for ( int m = 0; m < Gpaths.isize( ); m++ )
                         {    const vec<int>& p = Gpaths[m];
                              vec<PartialSeqOnHyper> v;
                              for ( int i = 0; i < p.isize( ); i++ )
                              {    const alignment_plus& ap = places[ p[i] ];
                                   v.push_back( 
                                        PartialSeqOnHyper( ap.Id2( ), ap.a.pos1( ),
                                        ap.a.pos2( ), ap.a.extent1( ) ) );    }
			      SeqOnHyper s( pass == 2, id1, v );
			      if( p_saligns )
				p_saligns->push_back(s);
			      if( p_csaligns ) {
				p_csaligns->push_back( CompressedSeqOnHyper() );
				p_csaligns->back().Compress( s, hb, &seqs[id1],
							     to_right_vertex );
			      }
			 }
		    }
	       }
	  }
     }
}


// Wrapper around AlignSeqToHyper (Extend alignment_plus Method)
// This version returns the aligns in the SeqOnHyper from only
void AlignSeqsToHyper( const HyperBasevector& hb, const vecbasevector& seqs,
     const vec<alignment_plus>& aligns, vec<SeqOnHyper>& saligns ) { 
  AlignSeqsToHyper( hb, seqs, aligns, &saligns, NULL );
}


// Wrapper around AlignSeqToHyper (Extend alignment_plus Method)
// This version returns the aligns in the CompressedSeqOnHyper from only
void AlignSeqsToHyper( const HyperBasevector& hb, const vecbasevector& seqs,
     const vec<alignment_plus>& aligns, vec<CompressedSeqOnHyper>& csaligns ) { 
  AlignSeqsToHyper( hb, seqs, aligns, NULL, &csaligns );
}


// Compresses alignments, extending any partial alignments along the hyperkmerpath.

void ExtendAndCompressAlignments ( const HyperKmerPath& h, const HyperBasevector& hb,
		     const vecbasevector& reads, const vecbasevector& edges,
		     const vec<partial_perfect_align>& aligns,
		     vec<CompressedSeqOnHyper>& csaligns )
{
  const int K = h.K();
  const int K_1 = K - 1;    // Use K-1 often, lets make it clear when we do.

  const int max_depth = 100;  // Max no. of edges to extend over (excluding initial edge)

  // Obtain edge lists
  vec<int> to_right_vertex( h.EdgeObjectCount( ) );
  h.ToRight(to_right_vertex);
  vec<int> to_left_vertex( h.EdgeObjectCount( ) );
  h.ToLeft(to_left_vertex);

  unsigned int num_aligns = aligns.size();
  for (unsigned int i = 0; i < num_aligns; ++i) {
    const partial_perfect_align& ppa = aligns[i];

    // Initial Alignment Info
    int read_id = ppa.read_id;
    int edge_id = ppa.edge_id;
    int read_length = ppa.read_length;
    int la_pos1 = ppa.pos1;
    int la_pos2 = ppa.pos2;
    int la_length = ppa.length;
    bool rc1 = ( read_id < 0 );
    if ( rc1 ) read_id = -read_id-1;

    if (la_length == read_length) {
      // Read entirely contained on edge - very easy!
      csaligns.push_back(CompressedSeqOnHyper(rc1, read_id,
					      edge_id, 0, la_pos2, read_length,
					      &reads[read_id], to_right_vertex ));
    } else {
      // Partial alignment - try extending along HyperKmerPath

      // Obtain correctly orientated read
      const basevector& s = reads[read_id];
      basevector src;
      if ( rc1 ) src.ReverseComplement(s);
      const basevector& S = ( rc1 ? src : s );

      // What directions do we need to expand the alignment?
      Bool extend_left = la_pos1 != 0;
      Bool extend_right = (la_pos1 + la_length < read_length);

      // When extending in both directions, need both halves for CompressedSeqOnHyper
      vec<vec<PartialSeqOnHyper> > left_v;
      vec<vec<PartialSeqOnHyper> > right_v;

      // Get list of to and from vertices
      const vec<vec<int> >& from_edges = hb.FromEdgeObj();
      const vec<vec<int> >& to_edges = hb.ToEdgeObj();

      // Information we need to keep track of during extension search
      vec<int> current_index(max_depth + 1);
      vec<int> used_on_edge(max_depth + 1);
      vec<int> current_edge(max_depth + 1);
      
      current_index[0] = 0;
      used_on_edge[0] = la_length;
      current_edge[0] = edge_id;

      // Try extending left then extending right
      for (int direction = 0; direction < 2; direction++) {
	const bool left = (direction == 0);
	const bool right = !left;

	// Only try to extend in directions that we need to
	if (left && !extend_left)
	  continue;
	if (right && !extend_right)
	  continue;
	if (right && extend_left && left_v.empty())
	  continue;

	// Reset search information in preparation to start extension
	current_index[1] = 0;
	int current_depth = 1;
	int used = la_length;
	int remaining = (left ? la_pos1 : read_length - used - la_pos1 );

	// Obtain list of edges and vertices for the current extension direction
	const vec<int>& vertex_list = (left ? to_left_vertex : to_right_vertex);
	const vec<vec<int> >& edge_list = (left ? hb.ToEdgeObj() : hb.FromEdgeObj());
	
	// Obtain list of edges out in the current extension direction
	vec<int> from = edge_list[vertex_list[edge_id]];

	vec<PartialSeqOnHyper> v;
	v.push_back(PartialSeqOnHyper(edge_id, la_pos1, la_pos2, used ));

	// Search for valid extensions in current direction
	for (;;) {
	  int index = current_index[current_depth];
	  if (index < from.isize()) {
	    int edge = from[index];
	    current_edge[current_depth] = edge;
	    const basevector& e = edges[edge];
	    int edge_length = e.isize();
	    int on_edge = Min(remaining, e.isize() - K_1);

	    // See whether or not the bases on this edge match the bases on the
	    // read, up to the end of the edge (in the appropriate direction).
	    // If they do, we have a valid extension.
	    Bool mismatch = false;
	    if (left)  // Left
	      for ( int y = 1; y <= on_edge + K_1; y++ ) {
		if ( S[remaining - y + K_1 ] != e[edge_length - y] ) {
		  mismatch = True;
		  break;
		}
	      }
	    else       // Right
	      for ( int y = 0; y < on_edge + K_1; y++ ) {
		if ( S[la_pos1 + used + y - K_1 ] != e[y] ) {
		  mismatch = True;
		  break;
		}
	      }

	    if (mismatch) {
	      current_index[current_depth]++;
	      continue; 
	    } else {
	      int pos1, pos2, overlap;
	      if (left) {
		pos1 = remaining - on_edge;
		pos2 = edge_length - K_1 - on_edge;
		overlap = on_edge + K_1;
	      } else {
		pos1 = Max(0, used - K_1 + la_pos1);
		pos2 = Max(0, K_1 - used);
		overlap = (K_1 - pos2) + on_edge;
	      }
	      
	      // Record the alignment extension that we have found.
	      v.push_back(PartialSeqOnHyper(edge, pos1, pos2, overlap));
	      if (on_edge == remaining) {

		if (left)
		  v.ReverseMe();
		
		if (left && extend_right) {
		  left_v.push_back(v);
		} else if (right && extend_left) {
		  right_v.push_back(v);
		} else {
		  // We have found a valid and complete alignment of this read!
		  // Make a CompressedSeqOnHyper for this alignment.
		  SeqOnHyper soh(rc1, read_id, v);
		  // soh.PrintVisual(cout, reads, edges);
		  csaligns.push_back(CompressedSeqOnHyper(soh, hb, &reads[read_id],
							  to_right_vertex)); 
		}
		if (left)
		  v.ReverseMe();

		v.resize(current_depth);
		current_index[current_depth]++;
	      } else if (current_depth == max_depth) {
		// Gone as deep into the tree as we want to go
		break;
	      } else {
		used += on_edge;
		remaining -= on_edge;
		used_on_edge[current_depth] = on_edge;
		from = edge_list[vertex_list[edge]];
		current_depth++;
		current_index[current_depth] = 0;
	      }
	    }
	  } else {
	    if (current_depth == 1)
	      break;
	    else {
	      current_depth--;
	      current_index[current_depth]++;
	      used -= used_on_edge[current_depth];
	      remaining += used_on_edge[current_depth];
	      v.resize(current_depth);
	      from = edge_list[vertex_list[current_edge[current_depth - 1]]];
	    }
	  }
	}
      }

      // Join left and right extensions
      if (extend_left && extend_right && !right_v.empty()) {
	for (int left_index = 0; left_index < left_v.isize(); ++left_index) {
	  vec<PartialSeqOnHyper> v = left_v[left_index];
	  for (int right_index = 0; right_index < right_v.isize(); ++right_index) {
	    v.append(right_v[right_index], 1, right_v[right_index].size());
	    SeqOnHyper soh(rc1, read_id, v);
	    //	    soh.PrintVisual(cout, reads, edges);
	    csaligns.push_back(CompressedSeqOnHyper(soh, hb, &reads[read_id],
						    to_right_vertex)); 
	    v.resize(left_v[left_index].size());
	  }
	}
      }
    }
  }
}


// Builds an index into the alignments based on the alignment Id1 value
void BuildAlignmentIndex (const vec<CompressedSeqOnHyper>& csaligns, int nreads, 
			  vec<vec<int> >& csaligns_index, const bool verbose ) {

  if ( verbose ) {
    cout << Date( ) << ": building alignment index" << endl;
    PrintMemUsage( );
  }

  csaligns_index.clear_and_resize(nreads);
  int csaligns_size = csaligns.isize();
  for ( int i = 0; i < csaligns_size; i++ )
    csaligns_index[ csaligns[i].Id1( ) ].push_back(i);
}


// Useful only for debugging purposes:
int VerifyAlignSeqsToHyper( const HyperBasevector& hb, 
			    const vecbasevector& seqs ) {

  int errors = 0;

  // Copied from the three-argument version of AlignSeqsToHyper:

  vec<alignment_plus> aligns;
  int K = hb.K( );
  vecbasevector seqs_edges = seqs;
  for ( int i = 0; i < hb.EdgeObjectCount( ); i++ )
    seqs_edges.push_back_reserve( hb.EdgeObject(i) );
  PerfectAlignerLG pal( hb.K( ), PerfectAlignerLG::findProperOnly );
  pal.Align( seqs_edges, aligns, 1, seqs.size( ) );
  
  // Generate both regular and compressed SeqOnHypers:
  vec<SeqOnHyper> saligns;
  vec<CompressedSeqOnHyper> csaligns;

  AlignSeqsToHyper( hb, seqs, aligns, &saligns, &csaligns );

  PRINT3(aligns.size(), saligns.size(), csaligns.size());
  ForceAssertEq( saligns.size(), csaligns.size() );

  // Check whether they match:

  SeqOnHyper decompressed;

  for(int i=0; i < saligns.isize(); i++) {
    cout << "Align " << i << " of basevec " 
	 << saligns[i].Id1() << " " << seqs[saligns[i].Id1()].ToString()
	 << "\nwas: " << saligns[i] << endl;
    csaligns[i].DecompressInto( decompressed, hb );

    cout << "now: " << decompressed << endl;

    if( decompressed == saligns[i] )
      cout << "Good\n" << endl;
    else {
      cout << "ERROR!\n" << endl;
      errors++;
    }
  }

  return errors;

}


void CompareAligns( const HyperBasevector& hb,
		    const vecbasevector& reads, const vecbasevector& edges,
		    const vec<CompressedSeqOnHyper>& csaligns,
		    const vec<CompressedSeqOnHyper>& csaligns_old)
{
  unsigned int nreads = reads.size( );
  vec<vec<int> > new_index(nreads);
  vec<vec<int> > old_index(nreads);

  for (unsigned int i = 0; i < csaligns.size(); ++i)
    new_index[csaligns[i].Id1()].push_back(i);
  for (unsigned int i = 0; i < csaligns_old.size(); ++i)
    old_index[csaligns_old[i].Id1()].push_back(i);

  int new_align = 0;
  int old_align = 0;

  for (unsigned int i = 0; i < nreads; ++i ) {
    vec <int> new_list(new_index[i].size());
    vec <int> old_list(old_index[i].size());

    for (int j = 0; j < new_index[i].isize(); ++j) 
      new_list[j] = j;
    for (int k = 0; k < old_index[i].isize(); ++k)
      old_list[k] = k;

    for (int j = 0 ; j < new_index[i].isize(); ++j) {
      for (int k = 0; k < old_index[i].isize(); ++k) {

	if (old_list[k] != -1)
	  if (csaligns[new_index[i][j]] == csaligns_old[old_index[i][k]]) {
	    new_list[j] = -1;
	    old_list[k] = -1;
	    break;
	  }
      }
    }

    for (int j = 0; j < new_list.isize(); ++j) {
      if (new_list[j] != -1) {
 	cout << "NEW ALIGN of " << i  << endl;
 	SeqOnHyper prob_align;
 	csaligns[new_index[i][j]].DecompressInto(prob_align, hb);
 	prob_align.PrintStats(cout, reads, edges);
 	prob_align.PrintVisual(cout, reads, edges);
 	new_align++;
 	PRINT2(new_list.size(),old_list.size()  );
	cout << "All new alignments" << endl;
	for (int k = 0 ; k < new_index[i].isize(); ++k) {
	  SeqOnHyper prob_align;
	  csaligns[new_index[i][k]].DecompressInto(prob_align, hb);
	  prob_align.PrintStats(cout, reads, edges);
	  prob_align.PrintVisual(cout, reads, edges);
	}
	cout << "All Old Alignments" << endl;
	for (int k = 0 ; k < old_index[i].isize(); ++k) {
	  SeqOnHyper prob_align;
	  csaligns_old[old_index[i][k]].DecompressInto(prob_align, hb);
	  prob_align.PrintStats(cout, reads, edges);
	  prob_align.PrintVisual(cout, reads, edges);
	}
	cout << "END" << endl;
      }
    }

    for (int j = 0; j < old_index[i].isize(); ++j) {
      if (old_list[j] != -1) {
	cout << "OLD ALIGN of " << i  << endl;
	SeqOnHyper prob_align;
	csaligns_old[old_index[i][j]].DecompressInto(prob_align, hb);
	prob_align.PrintStats(cout, reads, edges);
	prob_align.PrintVisual(cout, reads, edges);
	old_align++;

	PRINT2(new_list.size(),old_list.size()  );
	cout << "All new alignments" << endl;
	for (int k = 0 ; k < new_index[i].isize(); ++k) {
	  SeqOnHyper prob_align;
	  csaligns[new_index[i][k]].DecompressInto(prob_align, hb);
	  prob_align.PrintStats(cout, reads, edges);
	  prob_align.PrintVisual(cout, reads, edges);
	}
	cout << "All Old Alignments" << endl;
	for (int k = 0 ; k < old_index[i].isize(); ++k) {
	  SeqOnHyper prob_align;
	  csaligns_old[old_index[i][k]].DecompressInto(prob_align, hb);
	  prob_align.PrintStats(cout, reads, edges);
	  prob_align.PrintVisual(cout, reads, edges);
	}
	cout << "END" << endl;
	
      }
    }
  }
  PRINT2(old_align, new_align);

}
