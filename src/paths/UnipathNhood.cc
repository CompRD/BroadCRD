/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "Intvector.h"
#include "ReadLocation.h"
#include "ReadPairing.h"
#include "TaskTimer.h"
#include "VecUtilities.h"
#include "graph/Digraph.h"
#include "math/Functions.h"
#include "paths/KmerPath.h"
#include "paths/UnipathNhood.h"
#include "paths/UnipathNhoodCommon.h"
#include "paths/simulation/Placement.h"
#include "random/NormalRandom.h"
#include "random/NormalDistribution.h"
#include "random/Random.h"


// Make it easy and uniform to hold on to the read partner data.

struct offset_interval {
  offset_interval(const int& len1, const int& len2, 
		  const int& dist2, const int& sd2)
    : len1(len1), len2(len2), dist2(dist2), sd2(sd2) {}

  int len1;  // length of interval 1, the other unipath (dist unspecified)
  int len2;  // length of interval 2, the read partner
  int dist2; // mean value of the offset to interval 2
  int sd2;   // standard deviation of dist2
};

// overload IntervalOverlapProbability to use offset_intervals.
// Just needs dist1 = actual separation between unipaths

inline float OverlapProb( const int& dist1, const offset_interval& oi) {
  return IntervalOverlapProbability( oi.len1, oi.len2,
				     oi.dist2-dist1, oi.sd2 );
}



/**
   FuncDecl: FindUnipathNhood

   Given a seed unipath, find all nearby <normal unipaths> to form a
   <neighborhood> around this seed.
*/
void FindUnipathNhood( 

     // Inputs:

     const unipath_id_t v,                    // seed vertex
     const digraphE<sepdev>& G,               // graph of all normal unipaths
     const vec<nbases_t>& ulen,               // length of each unipath
     const VecPdfEntryVec& cp,                // copy number pdf for unipaths
     const vec<int>& predicted_copyno,        // predicted copy number for unipaths
     const vec<Bool>& branch,                 // which unipaths are at branches
     const vecKmerPath& paths,                // the read paths
     const vec<int>& path_lens,               // the lengths of the read paths
     const vec<read_location_short>& ulocs,   // locations of reads on unipaths
     const vec<int>& uindex,                  // index to ulocs
     const vec<read_pairing>& pairs,          // all the read pairs
     const vec<int>& pairs_index,             // index to read pairs by reads
     const Bool FILTER_NHOOD,                 // processing option: filter?
     const int MAX_COPY_NUMBER_OTHER,         // screens unipaths entering nhood
     const int NHOOD_RADIUS,                  // how far we go away from the seed
     const int MAX_TO_CUT,                    // cutpoints longer than this are kept
     const int MAX_DEV,                       // how stretchy sep from seed can get
     const int MAX_PATHS_IN_NHOOD,            // how big the nhood can get
     const Bool BUILD_NHOOD_FW_ONLY,          // forward only?

     // Output:

     vec<ustart>& processed                   // positions of unipaths in nhood

          )

{    
  static vec<int> zero_dev;
  if ( zero_dev.empty( ) ) zero_dev.push_back(0);
  static vec<ustart> unprocessed;
  static vec<int> allowed;
  allowed.clear( );
  int upasses = ( FILTER_NHOOD ? 2 : 1 );
  for ( int upass = 1; upass <= upasses; upass++ ) {
    processed.clear( ), unprocessed.clear( );
    unprocessed.push_back( ustart( v, 0, zero_dev ) );
    while( unprocessed.nonempty( ) ) {
      pop_heap( unprocessed.begin(), unprocessed.end(), ustart::OrderByDescendingMeanDev() );
      int w = unprocessed.back( ).Uid( );
      int wstart = unprocessed.back( ).Start( );
      vec<int> wdev = unprocessed.back( ).Dev( );
      // If we've processed this unipath before, that ustart had a
      // MeanDev no worse than this one, so skip this one.
      Bool used = False;
      for ( int j = 0; j < processed.isize( ); j++ )
        if ( w == processed[j].Uid( ) ) {
          used = True;
          break;
        }
      if ( used ) {
        unprocessed.pop_back();
        continue;
      }

      processed.push_back( unprocessed.back( ) );
      unprocessed.pop_back();

      for ( int i = 0; i < G.From(w).isize( ); i++ ) {
        int x = G.From(w)[i];
        if ( predicted_copyno[x] > Max(MAX_COPY_NUMBER_OTHER, predicted_copyno[v]) ) 
             continue;    
        if ( predicted_copyno[v] == 1 && MAX_COPY_NUMBER_OTHER <= 1 
             && branch[x] && ulen[x] <= 100 )
        {    continue;    }
        if ( upass == 2 && !BinMember( allowed, x ) ) continue;
        int wxsep = G.EdgeObjectByIndexFrom( w, i ).Sep( );
        int xstart = wstart + ulen[w] + wxsep;
        if ( xstart > ulen[v] + NHOOD_RADIUS ) continue;

        static vec<int> dev;
        dev = wdev;
        dev.push_back( G.EdgeObjectByIndexFrom( w, i ).Dev( ) );
        ustart new_ustart( x, xstart, dev );
        if ( new_ustart.MeanDev( ) > MAX_DEV ) continue;

        Bool used = False;
        for ( int j = 0; j < processed.isize( ); j++ )
          if ( x == processed[j].Uid( ) ) {
            used = True;
            break;
          }
        if (used) continue;

        unprocessed.push_back( new_ustart );
        push_heap( unprocessed.begin(), unprocessed.end(), 
             ustart::OrderByDescendingMeanDev() );
        
        if ( processed.isize( ) + unprocessed.isize( ) >= MAX_PATHS_IN_NHOOD ) break;
      }
      if ( processed.isize( ) + unprocessed.isize( ) >= MAX_PATHS_IN_NHOOD ) break;

      if ( !BUILD_NHOOD_FW_ONLY ) {
        for ( int i = 0; i < G.To(w).isize( ); i++ ) {
          int x = G.To(w)[i];
          if ( predicted_copyno[x] > 
               Max( MAX_COPY_NUMBER_OTHER, predicted_copyno[v] ) )
          {    continue;    }
          if ( predicted_copyno[v] == 1 && MAX_COPY_NUMBER_OTHER <= 1 
               && branch[x] && ulen[x] <= 100 )
          {    continue;    }
          if ( upass == 2 && !BinMember( allowed, x ) ) continue;
          int xwsep = G.EdgeObjectByIndexTo( w, i ).Sep( );
          int xstart = wstart + - xwsep - ulen[x];
          if ( xstart + ulen[x] < -NHOOD_RADIUS ) continue;

          static vec<int> dev;
          dev = wdev;
          dev.push_back( G.EdgeObjectByIndexTo( w, i ).Dev( ) );
          ustart new_ustart( x, xstart, dev );
          if ( new_ustart.MeanDev( ) > MAX_DEV ) continue;

          Bool used = False;
          for ( int j = 0; j < processed.isize( ); j++ )
            if ( x == processed[j].Uid( ) ) {
              used = True;
              break;
            }
          if (used) continue;

          unprocessed.push_back( new_ustart );
          push_heap( unprocessed.begin(), unprocessed.end(), 
               ustart::OrderByDescendingMeanDev() );

          if ( processed.isize( ) + unprocessed.isize( ) >= MAX_PATHS_IN_NHOOD ) 
               break;
        }
      }
      if ( processed.isize( ) + unprocessed.isize( ) >= MAX_PATHS_IN_NHOOD ) break;
    }

    // If there are any unprocessed unipaths left, that are not already in 
    // processed, add them in.

    for ( int j = 0; j < unprocessed.isize( ); j++ )
    {    Bool used = False;
         for ( int z = 0; z < processed.isize( ); z++ )
              if ( processed[z].Uid( ) == unprocessed[j].Uid( ) ) used = True;
         if ( !used ) processed.push_back( unprocessed[j] );    }
    Sort(processed);
    
    // Note: it is likely that by more effectively using the available 
    // information, we could come up with better estimates for the positions
    // of the unipaths in the neighborhood.
    
    // Screen the neighborhood to remove unipaths which are incompatible
    // with the seed neighborhood (v).  This is an interim step.  We will
    // want to systematically address incompatibilities.  Also remove cut
    // points.
    
    if ( upass == 1 && FILTER_NHOOD ) {
      static vec<ustart> processedx;
      processedx.clear( );
      for ( int j = 0; j < processed.isize( ); j++ ) {
        int start = processed[j].Start( ), w = processed[j].Uid( );
        int dev = processed[j].MeanDev( );
        if ( start >= 0 ) {
          Bool linked = Linked( v, w, ulocs, uindex, pairs, pairs_index );
          double p = linked ? 0 
            : CalcLinkProbability( v, w, start - ulen[v], dev, ulen, ulocs, uindex, 
                               cp, pairs, pairs_index, paths, path_lens );
          if ( linked || p < 0.95 )
            processedx.push_back( processed[j] );    
        }
        else {
          Bool linked = Linked( w, v, ulocs, uindex, pairs, pairs_index );
          double p = linked ? 0 
            : CalcLinkProbability( w, v, -start - ulen[w], dev, ulen, ulocs, uindex, 
                               cp, pairs, pairs_index, paths, path_lens );
          if ( linked || p < 0.95 )
            processedx.push_back( processed[j] );
        }
      }
      processed = processedx;
      for ( int j = 0; j < processed.isize( ); j++ )
        allowed.push_back( processed[j].Uid( ) );
      UniqueSort(allowed);
      continue;
    }

    // For each graph cut point corresponding to a unipath of length
    // less than a middling insert length, remove the vertices that do
    // not lie in the connected component of the neighborhood seed v.
    
    if ( upass == 2 ) {
      // Create a graph that contains all the selected unipaths and
      // the edges between them.
      static vec<int> processedv;
      processedv.clear( );
      for ( int j = 0; j < processed.isize( ); j++ )
        processedv.push_back( int(processed[j].Uid( )) );
      digraphE<sepdev> N( G.Subgraph(processedv) );

      // Find the vertices (other than the seed) that, if removed,
      // would disconnect the graph.
      static vec<int> cuts;
      N.CutPoints(cuts);
      static vec<Bool> to_remove;
      to_remove.resize_and_set( N.N( ), False );
      int seed_vx = Position( processedv, int(v) );
      for ( int j = 0; j < cuts.isize( ); j++ ) {
        int c = cuts[j];
        if ( c == seed_vx ) continue;
        int cut_unipath = processedv[c];
        if ( ulen[ cut_unipath ] > MAX_TO_CUT ) continue;
        equiv_rel e( N.N( ) );
        for ( int vx = 0; vx < N.N( ); vx++ ) {
          if ( vx == c ) continue;
          for ( int j = 0; j < N.From(vx).isize( ); j++ ) {
            int w = N.From(vx)[j];
            if ( w == c ) continue;
            e.Join(vx, w);
          }
        }
        static vec<int> keep;
        e.Orbit( seed_vx, keep );
        Sort(keep);
        for ( int i = 0; i < N.N( ); i++ ) {
          if ( i != c && !BinMember( keep, i ) )
            to_remove[i] = True;
        }
      }
      EraseIf( processed, to_remove );    

//       Sort( cuts );
//       ofstream dot( "nhood.dot" );
//       dot << "digraph \"nhood_" << v << "\" {" << endl;
//       for ( int vx = 0; vx < N.N( ); vx++ ) {
//         dot << "  " << vx 
//             << " [label=\"" << processedv[vx] << "\\n" << ulen[ processedv[vx] ] << "\"";
//         if ( processedv[vx] == v )
//           dot << ", fillcolor=blue";
//         if ( to_remove[vx] )
//           dot << ", color=red";
//         if ( BinMember( cuts, vx ) )
//           dot << ", shape=diamond";
//         dot << "];\n";
//       }
//       for ( int vx = 0; vx < N.N( ); vx++ ) {
//         for ( int j = 0; j < N.From(vx).isize( ); j++ ) {
//           int w = N.From(vx)[j];
//           const sepdev& E = N.EdgeObjectByIndexFrom( vx, j );
//           dot << "  " << vx << " -> " << w
//               << " [label=\"" << E.Sep() << "~" << E.Dev() << "\"];\n";
//         }
//       }
//       dot << "}" << endl;
//       dot.close();

    }
  }
}

void BuildRawEdges(
     const int K,                              // as in Kmer
     const vec<read_pairing>& pairs,           // read pairs
     const vec<read_location_short>& ulocs,    // locations of reads on unipaths
     const VecULongVec& ulocs_indexr,          // index to it by reads
     const vec<Bool>& normal,                  // is a given unipath normal
     const vec<int>& ulen,                     // unipath lengths
     const vecKmerPath& paths,                 // the reads
     const vec<int>& to_rc,                    // map unipath to its rc
     const vec<Bool>& CN1,                     // if copy number one
     const vec< vec<Bool> >& dangerous,        // tagged bases on unipaths
     vec<edge>& raw_edges,                     // result
     bool verbose = false )
{    vec<edgeplus> plus;
     for ( int i = 0; i < pairs.isize( ); i++ )
     {    int id1 = pairs[i].id1, id2 = pairs[i].id2;
          static vec< pair<int,int> > uid_offset1, uid_offset2;
          uid_offset1.clear( ), uid_offset2.clear( );
          for ( ULongVec::size_type j = 0; j < ulocs_indexr[id1].size( ); j++ )
          {    const read_location_short& rl = ulocs[ ulocs_indexr[id1][j] ];
               int uid1 = rl.Contig( );
               if ( !normal[uid1] ) continue;
               if ( CN1[uid1] )
               {    Bool untrusted = False;
                    int start = rl.Start( );
                    int stop = start + paths[id1].KmerCount( ) + K - 1;
                    start = Max( 0, start );
                    stop = Min( ulen[uid1] + K - 1, stop );
                    for ( int x = start; x < stop; x++ )
                         if ( dangerous[uid1][x] ) untrusted = True;
                    if (untrusted) continue;    }
               if ( rl.Fw( ) )
                    uid_offset1.push_back( make_pair( uid1, -rl.Start( ) ) );    }
          for ( ULongVec::size_type j = 0; j < ulocs_indexr[id2].size( ); j++ )
          {    const read_location_short& rl = ulocs[ ulocs_indexr[id2][j] ];
               int uid2 = rl.Contig( );
               if ( !normal[uid2] ) continue;
               if ( CN1[uid2] )
               {    Bool untrusted = False;
                    int start = rl.Start( );
                    int stop = start + paths[id2].KmerCount( ) + K - 1;
                    start = Max( 0, start );
                    stop = Min( ulen[uid2] + K - 1, stop );
                    for ( int x = start; x < stop; x++ )
                         if ( dangerous[uid2][x] ) untrusted = True;
                    if (untrusted) continue;    }
               if ( rl.Rc( ) )
                    uid_offset2.push_back( make_pair( uid2, -rl.Start( ) ) );    }
          for ( int i1 = 0; i1 < uid_offset1.isize( ); i1++ )
          {    for ( int i2 = 0; i2 < uid_offset2.isize( ); i2++ )
               {    int uid1 = uid_offset1[i1].first, uid2 = uid_offset2[i2].first;
                    int off1 = uid_offset1[i1].second, off2 = uid_offset2[i2].second;
                    if ( uid1 == uid2 ) continue;

                    // Compute separation between uid1 and uid2, in kmers.

                    int sep = pairs[i].sep + K - 1
                         - ( off1 + ulen[uid1] - paths[id1].KmerCount( ) ) + off2;
                    int dev = pairs[i].sd;

                    // Don't allow an overlap of more than three deviations.

                    if ( sep - (K-1) < -3 * dev ) continue;

		    if (verbose)
		    {    cout << "Raw edges " << uid1 << "->" << uid2
			      << " and " << to_rc[uid2] << "->" << to_rc[uid1]
			      << ", sep=" << sep 
                              << ", dev=" << dev
                              << ", off1=" << off1 << ", off2=" << off2 << endl;    }

                    int roff1 = - ( ulen[uid1] + off1 - paths[id1].KmerCount( ) );
                    int roff2 = - ( ulen[uid2] + off2 - paths[id2].KmerCount( ) );
                    plus.push( uid1, uid2, sep, dev, off1, off2 );
                    plus.push( to_rc[uid2], to_rc[uid1], sep, dev, 
                         roff2, roff1 );    }    }    }

     // If two links have identical offsets on the unipaths (suggesting that they
     // may be duplicate molecules), treat them as a single link.

     UniqueSort(plus);
     for ( int i = 0; i < plus.isize( ); i++ )
          raw_edges.push_back( edge( plus[i] ) );    }

template<class T>
void BuildUnipathLinkGraph( 

     // inputs:

     const int K,                              // as in Kmer
     const vec<read_pairing>& pairs,           // read pairs
     const vec<read_location_short>& ulocs,    // locations of reads on unipaths
     const VecULongVec& ulocs_indexr,         // index to it by reads
     const vec<Bool>& normal,                  // is a given unipath normal
     const vec<int>& ulen,                     // unipath lengths
     const vecKmerPath& paths,                 // the reads
     const vec<int>& to_rc,                    // map unipath to its rc
     const int min_edge_multiplicity,          // else ignore edge
     const VecPdfEntryVec& cp,                 // unipath copy number

     // output:

     digraphE< Tsepdev<T> >& G,                // the graph

     // optional:

     bool verbose                              // default false

          )
{
     int nuni = ulen.size( );

     // Find regions of supposedly copy-number-one unipaths that appear to have
     // higher copy number.  These are tagged as 'dangerous'.  Here is the 
     // heuristic:
     //
     // 1. Consider every read placement on a copy-number-one unipath, for which
     // the partner has room to land, after stretching the link by 3 deviations.
     //
     // 2. If the partner lands on the unipath, increment a counter 'good' for
     // the bases under it.  Otherwise increment a counter 'bad'.
     //
     // 3. Any base on a copy-number-one unipath having a nonzero bad count
     // that is at least 20% of its good count is tagged as dangerous.

     vec<Bool> CN1( nuni, False );
     vec< vec<int> > good(nuni), bad(nuni);
     vec< vec<Bool> > dangerous(nuni);
     for ( int u = 0; u < nuni; u++ )
     {    int copyno = -1;
          double maxp = 0;
          for ( PdfEntryVec::size_type j = 0; j < cp[u].size( ); j++ )
          {    if ( cp[u][j].second > maxp )
               {    copyno = cp[u][j].first; 
                    maxp = cp[u][j].second;    }    }
          if ( copyno != 1 ) continue;
          CN1[u] = True;
          good[u].resize( ulen[u] + K - 1, 0 );
          bad[u].resize( ulen[u] + K - 1, 0 );
          dangerous[u].resize( ulen[u] + K - 1, False );    }
     for ( int i = 0; i < pairs.isize( ); i++ )
     {    const read_pairing& p = pairs[i];
          int id1 = p.id1, id2 = p.id2;
          static vec< pair<int,int> > uid_start1, uid_start2;
          uid_start1.clear( ), uid_start2.clear( );
          for ( ULongVec::size_type j = 0; j < ulocs_indexr[id1].size( ); j++ )
          {    const read_location_short& rl = ulocs[ ulocs_indexr[id1][j] ];
               if ( rl.Fw( ) )
               {    uid_start1.push_back( 
                         make_pair( rl.Contig( ), rl.Start( ) ) );    }    }
          for ( ULongVec::size_type j = 0; j < ulocs_indexr[id2].size( ); j++ )
          {    const read_location_short& rl = ulocs[ ulocs_indexr[id2][j] ];
               if ( rl.Rc( ) )
               {    uid_start2.push_back( 
                         make_pair( rl.Contig( ), rl.Start( ) ) );    }    }
          if ( !uid_start1.solo( ) || !uid_start2.solo( ) ) continue;
          int uid1 = uid_start1[0].first, uid2 = uid_start2[0].first;
          int start1 = uid_start1[0].second, start2 = uid_start2[0].second;
          if ( !CN1[uid1] ) continue;
          int len1 = paths[id1].KmerCount( ) + K - 1;
          int len2 = paths[id2].KmerCount( ) + K - 1;
          if ( start1 + len1 + p.sep + 3 * p.sd + len2 > ulen[uid1] + K - 1 )
               continue;
          int ruid1 = to_rc[uid1];
          for ( int j = 0; j < len1; j++ )
          {    if ( start1 + j >= 0 && start1 + j <= ulen[uid1] + K - 1 )
               {    if ( uid2 == uid1 ) 
                    {    ++good[uid1][ start1 + j ];
                         ++good[ruid1][ ulen[uid1] + K - 1 - (start1 + j) - 1 ];    }
                    else 
                    {    ++bad[uid1][ start1 + j ];
                         ++bad[ruid1][ ulen[uid1] + K - 1 - (start1 + j) - 1 ];    
                         }    }    }    }
     if (verbose)
     {    cout << "\nBases on 'copy-number-one' unipaths that may actually have "
          << "higher copy number:\n";    }
     for ( int u = 0; u < nuni; u++ )
     {    if ( !CN1[u] ) continue;
          for ( int j = 0; j < ulen[u] + K - 1; j++ )
          {    if ( bad[u][j] > 0 && 5 * bad[u][j] >= good[u][j] )
               {    dangerous[u][j] = True;
                    if (verbose)
                    {    cout << u << "." << j << ", good = " << good[u][j]
                              << ", bad = " << bad[u][j] << "\n";    }    }    }    }
     if (verbose) cout << "\n";

     // Find edges in graph (to be condensed).

     vec<edge> raw_edges;
     BuildRawEdges( K, pairs, ulocs, ulocs_indexr, normal, ulen, paths,
          to_rc, CN1, dangerous, raw_edges, verbose );

     // Condense edges.

     Sort(raw_edges);
     vec< Tedge<T> > condensed_edges;
     for ( int i = 0; i < raw_edges.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < raw_edges.isize( ); j++ )
          {    if ( raw_edges[j].uid1 != raw_edges[i].uid1 ) break;
               if ( raw_edges[j].uid2 != raw_edges[i].uid2 ) break;    }
          if ( j - i >= min_edge_multiplicity )
          {    int sep_sum = 0, nsep = 0, u;
               for ( u = i; u < j; u++ ) 
               {    if ( raw_edges[u].dev() > 2.0 * raw_edges[i].dev() ) 
                    {    if( nsep >= min_edge_multiplicity ) break;
		         // otherwise, try again, with the first larger i
		         // that will allow us get to a larger u.
		         sep_sum = nsep = 0;
		         for(i++; raw_edges[u].dev() > 2.0 * raw_edges[i].dev(); i++)
		              ;
		         u=i;    }
		     sep_sum += raw_edges[u].sep();
		     nsep++;    }
	       if ( nsep < min_edge_multiplicity ) 
               {    i = j - 1;
                    continue;    }

               double sep_ave = double(sep_sum) / double(nsep);
	       double dev = raw_edges[u-1].dev( );

	       // Do the edges (from i to u-1) really look like they all
	       // belong to the same genomic copies?  Hmm: what fraction 
	       // are within two dev's?

	       int good_edges=0;
	       for(int w=i; w<u; w++)
               {    if( abs(raw_edges[w].sep() - sep_ave) < 2.0 * dev )
		         good_edges++;    }
	       const double min_good_fraction = 0.75;
	       if( good_edges < min_good_fraction * (u-i) ) 
               {    // Abandon this edge unless, say, 75% of the raw edges are good.
		    if (verbose)
		    {    cout << "No edge " << raw_edges[i].uid1
			      << "->" << raw_edges[i].uid2 << ": "
			      << u-i-good_edges << " out of " << u-i
			      << " raw edges were outside of 2 dev's"
			      << endl;    }    }
	       else 
               {     // We have a good edge.
		     if (verbose)
		     {    cout << "Condensed graph edge " << raw_edges[i].uid1
			       << "->" << raw_edges[i].uid2 << " made from "
			       << nsep << " raw edges with average sep=" << sep_ave
			       << endl;    }

		      Tedge<T> this_edge( raw_edges[i].uid1, raw_edges[i].uid2,
                          sep_ave, dev );
		      condensed_edges.push_back( this_edge );    }    }
          i = j - 1;    }

     // Build graph.

     BuildGraphFromEdges( condensed_edges, nuni, G );    }

#define BULG_DEF(T) \
   template void BuildUnipathLinkGraph(const int, const vec<read_pairing>&, \
   const vec<read_location_short>&, const VecULongVec&, const vec<Bool>&,   \
   const vec<int>&, const vecKmerPath&, const vec<int>&, const int,         \
   const VecPdfEntryVec&, digraphE< Tsepdev<T> >&, bool )
BULG_DEF(int);
BULG_DEF(double);



template<class T>
void FillInTransitiveEdges( digraphE< Tsepdev<T> >& graph, 
			    const int radius, 
                            const double max_dev,
                            const double percent_improvement,
			    const vec<int>& predicted_copyno,
			    const vec<int>& unipath_len, 
                            const int min_unipath,
			    int verbosity,
			    VecPlacementVec* locs_p ) {

  double multiplier = 1.0 + percent_improvement / 100.0;

  // The Sep() of the edge v->w is the distance from the 
  // *right* end of v to the *left* end of w.  Ugh.
  // Let's transform the graph so that sep is relative to left 
  // endpoints throughout, and then transform back at the end.
  for(int v=0; v<graph.N(); v++)
    for(int i=0; i<graph.From(v).isize(); i++)
      graph.EdgeObjectByIndexFromMutable(v,i).AddToSep( unipath_len[v] );

  // Is the copy number of this vertex <= that of all adjacent vertices?
  vec<Bool> locally_minimal_copyno( graph.N(), true );
  for(int v=0; v<graph.N(); v++) {
    for(int i=0; i<graph.From(v).isize(); i++)
      locally_minimal_copyno[v] &=
	predicted_copyno[v] <= predicted_copyno[ graph.From(v)[i] ] ;
    for(int i=0; i<graph.To(v).isize(); i++)
      locally_minimal_copyno[v] &=
	predicted_copyno[v] <= predicted_copyno[ graph.To(v)[i] ] ;
  }


  vec<Bool> to_check( locally_minimal_copyno );

  if( verbosity ) {
    // DEBUGGING:
    cout << "Locally minimal copyno: ";
    for(int i=0; i<graph.N(); i++)
      if( locally_minimal_copyno[i] ) cout << i << "/" << predicted_copyno[i] << " ";
    cout << endl;
  }

//   vec<Bool> to_check( graph.N() );
//   for(int v=0; v<graph.N(); v++)
//     to_check[v] = (predicted_copyno[v]==1);


  // These hold info on vertices adjacent to the current vertex.
  // Storing in advance means we don't have to worry about things
  // changing under us as we modify the graph.
  vec<int> adj_vxs;
  vec< Tsepdev<T> > adj_seps;

  if (verbosity) cout << Date( ) << ": filling in transitive edges." << endl;

  while( to_check.CountValue(true) > 0 ) {

    if (verbosity) 
    {    int max_degree = 0;
         for ( int v = 0; v < graph.N( ); v++ )
         {    max_degree = Max( max_degree, 
                  graph.From(v).isize( ) + graph.To(v).isize( ) );    }
         cout << Date( ) << ": " << Sum(to_check) 
              << " vertices left to check, max degree = " << max_degree << endl;    }

    for(int vx=0; vx<graph.N(); vx++) {
      if( ! to_check[vx] ) continue;
      to_check[vx] = false;
      if ( unipath_len[vx] < min_unipath ) continue;

      // These vecs hold the info from both To(vx) and From(vx);
      // edges in To(vx) have their sense reversed (sep -> -sep).
      // Storing in advance means we don't have to worry about things
      // changing under us as we modify the graph.
      adj_vxs.clear();
      adj_seps.clear();

      for(int i=0; i<graph.To(vx).isize(); i++) {
        if ( unipath_len[ graph.To(vx)[i] ] < min_unipath ) continue;
	adj_vxs.push_back( graph.To(vx)[i] );
	adj_seps.push_back( graph.EdgeObjectByIndexTo(vx,i) );
	if( verbosity>2 ) {
	  cout << "* Edge from " << graph.To(vx)[i] << " to " << vx
	       << ", sep=" << graph.EdgeObjectByIndexTo(vx,i).Sep();
	  if( locs_p && IsBadEdge(graph.To(vx)[i], vx, 
				  graph.EdgeObjectByIndexTo(vx,i), *locs_p) )
	    cout << " (BAD EDGE!)";
	  cout << endl;
	}
	adj_seps.back().Flip();
      }
      for(int i=0; i<graph.From(vx).isize(); i++) {
        if ( unipath_len[ graph.From(vx)[i] ] < min_unipath ) continue;
	adj_vxs.push_back( graph.From(vx)[i] );
	adj_seps.push_back( graph.EdgeObjectByIndexFrom(vx,i) );
	if( verbosity>2 ) {
	  cout << "* Edge from " << vx << " to " << graph.From(vx)[i]
	       << ", sep=" << graph.EdgeObjectByIndexFrom(vx,i).Sep();
	  if( locs_p && IsBadEdge(vx, graph.From(vx)[i], 
				  graph.EdgeObjectByIndexFrom(vx,i), *locs_p) )
	    cout << " (BAD EDGE!)";
	  cout << endl;
	}
      }

      // Find vertices adjacent to vx and within radius of each other
      for(int i=0; i<adj_vxs.isize(); i++) {
	int vi = adj_vxs[i];
        if ( adj_seps[i].Dev( ) > max_dev ) continue;
	for(int j=0; j<adj_vxs.isize(); j++) {
	  int vj = adj_vxs[j];
          if ( adj_seps[j].Dev( ) > max_dev ) continue;

          // We are only interested in creating new edges that are
          // to and/or from a vertex of locally minimal copy number.
	  //
	  // No, wait, we can't do that right now: we might want to
	  // pick things with non-locally-minimal copyno.  Right?
	  // If there's a really long copyno=2 unipath which is
	  // linked to 1's, we may need to pick it anyway, I think.
	  //
	  // This will make things much slower, not to mention harder
	  // to think about.  Deep questions remain.

// 	  if( ! locally_minimal_copyno[vi] && ! locally_minimal_copyno[vj] )
// 	    continue;

	  double diff = adj_seps[j].Sep() - adj_seps[i].Sep();
	  if( vi != vj && diff <= radius && (diff>0 || (diff==0 && i>j)) ) {


	    Tsepdev<T> sd( diff,
                           sqrt( adj_seps[i].Dev()*adj_seps[i].Dev() +
                                 adj_seps[j].Dev()*adj_seps[j].Dev() ) );

            if ( sd.Dev( ) > max_dev ) continue;

	    //
	    // Is this better than any existing edge from vi to vj?
	    // If so, replace that edge with this one.
	    //
	    // For now, "better" = "smaller dev".
	    // If there may be multiple edges between vi and vj,
	    // eg reflecting different genomic positions, then
	    // this would need some more sophisticated heuristic
	    // to decide whether this edge represented the same
	    // area of the genome or not.  That seems hard.
	    //

	    int existing = BinPosition( graph.From(vi), vj );

	    // Finds any edge vi->vj, if one already exists.
	    bool add_edge = (existing == -1);

	    // For now, assume there is at most one edge vi->vj.
	    if( existing != -1 &&
		graph.EdgeObjectByIndexFrom(vi,existing).Dev() 
                     > multiplier * sd.Dev() ) {
	      // Delete the existing edge, and plan to replace it.
	      if( verbosity > 1 )
		cout << "Deleting edge " << vi << "->" << vj 
		     << " (sep=" << graph.EdgeObjectByIndexFrom(vi,existing).Sep()
		     << ", dev=" << graph.EdgeObjectByIndexFrom(vi,existing).Dev()
		     << ")" << endl;
	      graph.DeleteEdgeFrom(vi,existing);
	      add_edge = true;
	    }

	    if( add_edge ) {

	      if( verbosity > 1 ) {
		cout << "Adding edge " << vi << "->" << vj << " with sep="
		     << sd.Sep() << ", dev=" << sd.Dev() 
                     << ", via edge " << vx << endl;
		if( sd.Sep() < unipath_len[vi] )
		  cout << "  NOTE: unipath " << vi << " has length "
		       << unipath_len[vi] << ", so final sep will be "
		       << sd.Sep()-unipath_len[vi] 
		       << " (other way around would yield "
		       << -sd.Sep()-unipath_len[vj] << ")"
		       << endl;
		vec<int> offsets;
		if( locs_p && IsBadEdge(vi, vj, sd, *locs_p, &offsets) ) {
		  // Alert if this edge is wrong according to true locs.
		  const PlacementVec& iplaces = (*locs_p)[vi];
		  const PlacementVec& jplaces = (*locs_p)[vj];

		  cout << "BAD TRANSITIVE EDGE: "
		       << "This new edge is genomically incorrect."
		       << "\n   True placements of pivot vertex " << vx << ": ";
		  for(PlacementVec::size_type xx=0; xx<(*locs_p)[vx].size(); xx++)
		    cout << (*locs_p)[vx][xx] << " ";
		  cout << "\n   adj_seps[i].Sep()=" << adj_seps[i].Sep()
		       << ", adj_seps[j].Sep()=" << adj_seps[j].Sep();
		  cout << "\n   True placements of unipath vi=" << vi << ": ";
		  for(PlacementVec::size_type ii=0; ii<iplaces.size(); ii++)
		    cout << iplaces[ii] << " ";
		  cout << "\n   True placements of unipath vj=" << vj << ": ";
		  for(PlacementVec::size_type jj=0; jj<jplaces.size(); jj++)
		    cout << jplaces[jj] << " ";
		  if( offsets.empty() )
		    cout << "\n   In fact, vi and vj never appear in the same genome contig!";
		  else
		    cout << "\n   Closest true sep is off by " << offsets.front()
			 << " (" << ToString(offsets.front()/double(sd.Dev()))
			 << " dev's).";
		  cout << endl;
		}
	      } // end of verbosity>1 stuff

	      graph.AddEdge(vi,vj,sd);

	      // Update lmc.  Not entirely moral, since (I think) things now
	      // can depend on edge addition order, ie unipath numbering order.
	      locally_minimal_copyno[vi] &=
		predicted_copyno[vi] <= predicted_copyno[vj];
	      locally_minimal_copyno[vj] &=
		predicted_copyno[vj] <= predicted_copyno[vi];

	      to_check[vi] |= locally_minimal_copyno[vi];
	      to_check[vj] |= locally_minimal_copyno[vj];
	    }
	  }

	} // end loop over j
      } // end loop over i
    } // end loop over vx
  } // end while loop which runs until to_check is all false


  // Finally, transform the graph back, so that sep is relative to right
  // endpoints of From vertices, undoing what we did at the beginning.
  for(int v=0; v<graph.N(); v++)
    for(int i=0; i<graph.From(v).isize(); i++)
      graph.EdgeObjectByIndexFromMutable(v,i).AddToSep( -unipath_len[v] );

  if (verbosity) cout << Date( ) << ": done filling in transitive edges." << endl;

}

template
void FillInTransitiveEdges( digraphE<sepdev>&, const int, const double,
                            const double, const vec<int>&, const vec<int>&, 
                            const int, int, VecPlacementVec* );
template
void FillInTransitiveEdges( digraphE<fsepdev>&, const int, const double,
                            const double, const vec<int>&, const vec<int>&, 
                            const int, int, VecPlacementVec* );


// Debugging helper for FillInTransitiveEdges:
template<class SEPDEV> bool IsBadEdge( int vi, int vj, SEPDEV sd,
        VecPlacementVec& locs,
				       vec<int>* offsets_p = NULL ) {
  vec<int> offsets_local;
  vec<int>& offsets = (offsets_p ? *offsets_p : offsets_local);
  const PlacementVec& iplaces = locs[vi];
  const PlacementVec& jplaces = locs[vj];

  for(PlacementVec::size_type ii=0; ii<iplaces.size(); ii++)
    for(PlacementVec::size_type jj=0; jj<jplaces.size(); jj++)
      if( jplaces[jj].GenomeId() == iplaces[ii].GenomeId()
	  && jplaces[jj].Rc() == iplaces[ii].Rc() ) {
	// I think this is the right way to take rc-ness into account.
	// Note both swapping order of i/j and pos vs. Pos.
	if( jplaces[jj].Rc() )
	  offsets.push_back( int(abs( iplaces[ii].Pos() 
				      - jplaces[jj].Pos() 
				      - sd.Sep() )) );
	else
	  offsets.push_back( int(abs( jplaces[jj].pos() 
				      - iplaces[ii].pos() 
				      - sd.Sep() )) );
      }
  return( offsets.empty() || Min(offsets) > 5*sd.Dev() );
}

/**
   Function: PopulateNhoodWithReads

   Find the reads which go in a particular neighborhood.   Also return their orientations and predicted positions.

   Rough algorithm:

   We want to identify reads where at least one read of the pair is within <inner radius> of the <neighborhood seed>.
   For each unipath within the <outer radius> of the neighborhood seed, we take reads aligned to that unipath.
   From the estimated distance of the unipath to the nhood seed, and the position of the read on the unipath, we
   estimate how far the read is from the nhood seed; if within the inner radius, accept both the read and its partner.
   

   Input parameters:

     int v                                  -  seed for nhood: unipath id of the seed for this nhood
     const vec<ustart>& processed           -  this defines the nhood: for each unipath in the nhood,
                                               its id and its position relative to the seed 'v'
     int K                                  -  size of kmers in terms of which all paths are defined
     const vec<int> ulen,                   -  unipath lengths, for each unipath id.
     const vec<read_pairing>& pairs         -  read pairs: what reads go into each pair
     const vec<int>& pairs_index            -  index to read pairs: for each read, the index in 'pairs'
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

     vec< pair<read_id_t,orient_t> >& use        -  reads in nhood, with orientations
     vec<pos_rel_to_seed_t>& usestart            -  predicted start of read relative to the seed 'v'
     vec<pair_id_t>& pairs_to_use                -  the constituent pairs
   
*/
void PopulateNhoodWithReads(

     // Inputs:

     unipath_id_t v,                                 // seed for nhood
     const vec<ustart>& processed,          // this defines the nhood
     nbases_t K,                                 // as in Kmer
     const vec<nbases_t> ulen,                   // unipath lengths
     const vec<read_pairing>& pairs,        // read pairs
     const vec<pair_id_t>& pairs_index,           // index to read pairs
     const vecKmerPath& paths,              // kmer paths for reads
     const vec<read_location_short>& ulocs, // locations of reads on unipaths
     const vec<int>& uindex,                // index to ulocs
     const int NHOOD_RADIUS_INTERNAL,       // how far from origin we should go
     const int MAX_DEV,                     // how sloppy read locations can get
     const Bool REACH_FW_ONLY,              // only go forward?

     // Outputs:

     vec< pair<read_id_t,orient_t> >& use,            // reads in nhood, with orientations
     vec<pos_rel_to_seed_t>& usestart,                    // predicted start of read rel. v
     vec<pair_id_t>& pairs_to_use                 // the constituent pairs

          )

{    const double max_read_dist_devs = 2.0;
     use.clear( ), usestart.clear( ), pairs_to_use.clear( );
     for ( int j = 0; j < processed.isize( ); j++ )
     {
       // take one unipath from the nhood, take its start relative to the seed 'v',
       // and its unipath id.
          int start = processed[j].Start( );
	  unipath_id_t w = processed[j].Uid( );
          vec<int> dev = processed[j].Dev( );
	  // Look at all the reads aligned to this unipath-from-the-nhood.
          int indw1 = uindex[w], indw2 = uindex[w+1];
          for ( int t = indw1; t < indw2; t++ )
          {    read_id_t id1 = ulocs[t].ReadId( );
               Bool rc1 = ulocs[t].Rc( );
	       // Does this read fall roughly on this unipath?
               int start1 = start + ulocs[t].Start( );
               int stop1 = start1 + paths[id1].KmerCount( ) - 1;
               double dist1 = Distance( start1, stop1, 
                    -NHOOD_RADIUS_INTERNAL, NHOOD_RADIUS_INTERNAL + ulen[v] );
               double dev1 = processed[j].MeanDev( );
               Bool used_id1 = False;
               if ( dev1 == 0 || dist1/dev1 <= max_read_dist_devs )
               {    use.push_back( make_pair( id1, rc1 ? ORIENT_RC : ORIENT_FW ) );
                    usestart.push_back(start1);
                    used_id1 = True;    }
               if ( REACH_FW_ONLY && rc1 ) continue;
               pair_id_t pi = pairs_index[id1];
               if ( pi < 0 ) continue;
               read_id_t id2 = pairs[pi].Partner(id1);
               int sep12 = pairs[pi].sep, dev12 = pairs[pi].sd;
               double dev2 = double(dev12) * double(dev12);
               for ( int t = 0; t < dev.isize( ); t++ )
                    dev2 += double(dev[t]) * double(dev[t]);
               dev2 = sqrt(dev2);
               if ( dev2 > MAX_DEV ) continue;
               int start2;
               if ( !rc1 ) start2 = stop1 + sep12 + K - 1;
               else start2 = start1 - (sep12 + K - 1) - paths[id2].KmerCount( );
               int stop2 = start2 + paths[id2].KmerCount( ) - 1;
               double dist2 = Distance( start2, stop2, 
                    -NHOOD_RADIUS_INTERNAL, NHOOD_RADIUS_INTERNAL + ulen[v] );
               if ( dev2 == 0 || dist2/dev2 <= max_read_dist_devs )
               {    use.push_back( make_pair( id2, !rc1 ? ORIENT_RC : ORIENT_FW ) );    
                    usestart.push_back(start2);
                    if (used_id1) pairs_to_use.push_back(pi);    }    }    }
     SortSync(use, usestart), UniqueSort(pairs_to_use);
     vec<Bool> to_delete( use.isize( ), False );
     for ( int i = 0; i < use.isize( ) - 1; i++ )
          if ( use[i] == use[i+1] ) to_delete[i] = True;
     EraseIf( use, to_delete ), EraseIf( usestart, to_delete );    }

// Note: "use" now contains reads which extend substantially beyond
// NHOOD_RADIUS_INTERNAL, just because of cumulative pairing inaccuracies.
// It might be possible and desirable to eliminate some of these reads at
// this point.

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

          )

{
     // 1. Form a condensed list L of the KmerPathIntervals which appear in S
     // (where S = use).

     static vec<KmerPathInterval> L0, L;
     L0.clear( ), L.clear( );
     for ( int i = 0; i < use.isize( ); i++ )
     {    int id = use[i].first;
          if ( !use[i].second )
          {    for ( int j = 0; j < paths[id].NSegments( ); j++ )
                    L0.push_back( paths[id].Segment(j) );    }
          else
          {    for ( int j = 0; j < paths_rc[id].NSegments( ); j++ )
                    L0.push_back( paths_rc[id].Segment(j) );    }    }
     sort( L0.begin( ), L0.end( ), cmp_start );
     for ( int i = 0; i < L0.isize( ); i++ )
     {    longlong Istart = L0[i].Start( ), Istop = L0[i].Stop( );
          int j;
          for ( j = i + 1; j < L0.isize( ); j++ )
          {    if ( L0[j].Start( ) > Istop ) break;
               Istop = Max( Istop, L0[j].Stop( ) );    }
          L.push_back( KmerPathInterval( Istart, Istop ) );
          i = j - 1;    }

     // 2. Find all short-insert reads R that share kmers with L, and
     // whose partners do too.

     static vec< pair<read_id_t,Bool> > R0, R;
     R0.clear( ), R.clear( );
     for ( int i = 0; i < L.isize( ); i++ )
     {    static vec<longlong> con;
          Contains( pathsdb, L[i], con );
          for ( int j = 0; j < con.isize( ); j++ )
          {    const tagged_rpint& t = pathsdb[ con[j] ];
               read_id_t id = t.ReadId( );
               Bool rc = ( t.PathId( ) < 0 );
               if ( is_short_pair_read[id] ) 
                    R0.push_back( make_pair( id, rc ) );    }    }
     UniqueSort(R0);
     for ( int i = 0; i < R0.isize( ); i++ )
     {    Bool fw = !R0[i].second;
          if ( BinMember( R0, make_pair( partner[ R0[i].first ], fw ) ) )
               R.push_back( R0[i] );    }

     // 3. For each read in R, work from left to right, trying to cover it
     // by reads in S, yielding a new set Q.

     static vec< pair<int,Bool> > Q;
     Q.clear( );
     static vec<tagged_rpint> Spathsdb;
     Spathsdb.clear( );
     vec< vec<tagged_rpint> > extra_Spathsdb( PATH_KS.size() );
     
     for ( int i = 0; i < use.isize( ); i++ )
     {    if ( !use[i].second )
          {
	    paths[ use[i].first ].AppendToDatabase(Spathsdb, use[i].first);
	    for ( int j = 0; j < PATH_KS.isize(); j++ )
	      extra_paths[j][ use[i].first ].AppendToDatabase(extra_Spathsdb[j], use[i].first);
	  }
          else
          {    paths_rc[ use[i].first ]
                    .AppendToDatabase(Spathsdb, -use[i].first-1);

	    for ( int j = 0; j < PATH_KS.isize(); j++ )
	      extra_paths_rc[j][ use[i].first ]
		.AppendToDatabase(extra_Spathsdb[j], -use[i].first-1);	  
	  
	  }    }
     Prepare(Spathsdb);
     for ( int j = 0; j < PATH_KS.isize(); j++ )
       Prepare(extra_Spathsdb[j]);
     for ( int i = 0; i < R.isize( ); i++ )
     {    int id = R[i].first;
          Bool rc = R[i].second;
          const KmerPath& p = ( !rc ? paths[id] : paths_rc[id] );
	  Bool foundCover = False;
	  static int helped = 0, didntHelp = 0, tried = 0, shortButNotLong = 0;
	  for ( int j = PATH_KS.isize()-1; !foundCover && j >= 0; j-- )
	    if ( SubContig( ( !rc ? extra_paths[j][id] : extra_paths_rc[j][id] ),
			    extra_paths[j], extra_paths_rc[j], extra_Spathsdb[j] ) )
	      foundCover = True;
	  (foundCover ? helped : didntHelp)++;
	  Bool foundShortCover = False;
          if ( foundCover || ( foundShortCover = SubContig( p, paths, paths_rc, Spathsdb ) ) ) Q.push_back( R[i] );

	  if ( !foundCover && foundShortCover )
	    shortButNotLong++;

	  //if ( !(tried++ % 500000) )
	  //  cout << " helped=" << helped << " didntHelp=" << didntHelp << " shortButNotLong=" << shortButNotLong << endl;
     }

     // 4. Replace Q by the set P of reads in it whose partner are also in it.

     P.clear( );
     for ( int i = 0; i < Q.isize( ); i++ )
     {    Bool fw = !Q[i].second;
          if ( BinMember( Q, make_pair( partner[ Q[i].first ], fw ) ) )
               P.push_back( Q[i] );    }    }





Bool Linked( int x, int y, const vec<read_location_short>& ulocs,
     const vec<int>& uindex, const vec<read_pairing>& pairs, 
     const vec<int>& pairs_index )
{    vec<int> yreads;
     yreads.clear( );
     int indx1 = uindex[x], indx2 = uindex[x+1];
     int indy1 = uindex[y], indy2 = uindex[y+1];
     for ( int i = indy1; i < indy2; i++ )
          yreads.push_back( ulocs[i].ReadId( ) );
     for ( int j1 = indx1; j1 < indx2; j1++ )
     {    const read_location_short& rl1 = ulocs[j1];
          if ( !rl1.Fw( ) ) continue;
          int id1 = rl1.ReadId( );
          int pi = pairs_index[id1];
          if ( pi < 0 ) continue;
          int id2 = pairs[pi].Partner(id1);
          int j2 = BinPosition( yreads, id2 );
          if ( j2 < 0 ) continue;
          const read_location_short& rl2 = ulocs[ indy1 + j2 ];
          if ( !rl2.Rc( ) ) continue;
          return True;    }
     return False;    }



double CalcLinkProbability( int x, int y, int s, int d, const vec<int>& ulen,
     const vec<read_location_short>& ulocs, const vec<int>& uindex,
     const VecPdfEntryVec& cp, const vec<read_pairing>& pairs,
     const vec<int>& pairs_index, const vecKmerPath& paths,
     const vec<int>& path_lens ) {

  vec<offset_interval> x_partners;
  vec<offset_interval> y_partners;

  int indx1 = uindex[x], indx2 = uindex[x+1]; // interval of reads on unipath x
  int indy1 = uindex[y], indy2 = uindex[y+1]; // interval of reads on unipath y

  // Loop over the reads on unipath x, gather partner data:
  for ( int j1 = indx1; j1 < indx2; j1++ ) {
    const read_location_short& rl1 = ulocs[j1];
    if ( !rl1.Fw( ) ) continue;
    int id1 = rl1.ReadId( );
    int pi = pairs_index[id1];
    if ( pi < 0 ) continue;
    int sep = pairs[pi].sep, dev = pairs[pi].sd;
    int id2 = pairs[pi].Partner(id1);
    int id2_x_offset =  // right end of unipath x to left end of read id2
      -ulen[x] + rl1.Start( ) + path_lens[id1] + sep;

    x_partners.push_back( offset_interval(ulen[y], path_lens[id2], 
					  id2_x_offset, dev) );
  }

  // Loop over the reads on unipath y, gather partner data:
  for ( int j2 = indy1; j2 < indy2; j2++ ) {
    const read_location_short& rl2 = ulocs[j2];
    if ( !rl2.Rc( ) ) continue;
    int id2 = rl2.ReadId( );
    int pi = pairs_index[id2];
    if ( pi < 0 ) continue;
    int sep = pairs[pi].sep, dev = pairs[pi].sd;
    int id1 = pairs[pi].Partner(id2);
    int id1_y_offset =  // left=+; left end of y to right end of id1
      -rl2.Start( ) + sep;

    y_partners.push_back( offset_interval(ulen[x], path_lens[id1], 
					  id1_y_offset, dev) );
  }

  const PdfEntryVec& copy_pdf_x = cp[x];
  const PdfEntryVec& copy_pdf_y = cp[y];

  vec<float> xcopy_no_link_prob(copy_pdf_x.size());
  vec<float> ycopy_no_link_prob(copy_pdf_y.size());

  float sum_no_link_prob = 0.0;

  // NormalRandom random_unipath_sep(s,d);
  // tries /= 10;

  const int TRIES = 10;
  
  // Random numbers from normal distribution N(0,1):

  static double N[TRIES] = {0.0883064, 1.29813, -0.540512, -0.279401, 1.26195,
       0.938126, -0.918371, -1.83504, 0.347259, -0.418964};

  for( int t=0; t < TRIES; t++ ) {

    // For each i from 0 to copy_pdf_x.size(), xcopy_no_link_prob[i]
    // accumulates the product which is the probability that x has copy 
    // number copy_pdf_x[i].first and there are no links from x hitting y
    for(PdfEntryVec::size_type i=0; i<copy_pdf_x.size(); i++)
      xcopy_no_link_prob[i] = copy_pdf_x[i].second;
    for(PdfEntryVec::size_type j=0; j<copy_pdf_y.size(); j++)
      ycopy_no_link_prob[j] = copy_pdf_y[j].second;

    // Pick a unipath separation from the appropriate distribution.

//    int unipath_sep = int(round( s + d*FastNormal() ));
    // int unipath_sep = int(round(random_unipath_sep.value()));
    int unipath_sep = int( round( double(s) + ( N[t] * double(d) ) ) );

    float p;
    for( int xp=0; xp<x_partners.isize(); xp++ ) {
      p = OverlapProb( unipath_sep, x_partners[xp] );
      for( PdfEntryVec::size_type cx=0; cx<copy_pdf_x.size(); cx++ )
	xcopy_no_link_prob[cx] *= (1 - p/copy_pdf_x[cx].first);
    }
    for( int yp=0; yp<y_partners.isize(); yp++ ) {
      p = OverlapProb( unipath_sep, y_partners[yp] );
      for( PdfEntryVec::size_type cy=0; cy<copy_pdf_y.size(); cy++ )
	ycopy_no_link_prob[cy] *= (1 - p/copy_pdf_y[cy].first);
    }

    // Combine all the partial probabilities:
    float x_no_link_prob = accumulate( xcopy_no_link_prob.begin(),
				       xcopy_no_link_prob.end(), 0.0);
    float y_no_link_prob = accumulate( ycopy_no_link_prob.begin(),
				       ycopy_no_link_prob.end(), 0.0);
    float no_link_prob = x_no_link_prob * y_no_link_prob;
    
    sum_no_link_prob += no_link_prob;
  }

  float avg_no_link_prob = sum_no_link_prob / TRIES;

  return ( 1 - avg_no_link_prob );

}
