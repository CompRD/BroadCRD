 /////////////////////////////////////////////////////////////////////////////
 //                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
 //       This software and its documentation are copyright (2009) by the   //
 //   Broad Institute/Massachusetts Institute of Technology.  All rights    //
 //   are reserved.  This software is supplied without any warranty or      //
 //   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
 //   can be responsible for its use, misuse, or functionality.             //
 /////////////////////////////////////////////////////////////////////////////

#include "Superb.h"
#include "math/HoInterval.h"
#include "paths/CleanseScaffoldGraph.h"
#include "paths/reporting/CLinkBundle.h"

 /**
  * class CBadness
  *
  * Container for a descriptor of the types of "badness" for a center
  * super. Mainly used to appropriately sort the badness events and
  * decide which super to fix first.
  */
 class CBadness{

 public:

   CBadness( ) : id_ ( -1 ), len_ ( -1 ), wcenter_ ( 0 ), wothers_ ( 0 ) { }

   CBadness( int id, int len, int wcenter, int wothers ) :
     id_ ( id ), len_ ( len ), wcenter_ ( wcenter ), wothers_ ( wothers ) { }

   int Id( ) const { return id_; }
   int Len( ) const { return len_; }
   int WCenter( ) const { return wcenter_; }
   int WOthers( ) const { return wothers_; }

   void IncrementWCenter( ) { wcenter_ += 1; }
   void IncrementWOthers( ) { wothers_ += 1; }

   friend bool operator< ( const CBadness &left, const CBadness &right ) {
     if ( left.WCenter( ) != right.WCenter( ) )
       return ( left.WCenter( ) > right.WCenter( ) );   // note >
     if ( left.WOthers( ) != right.WOthers( ) )
       return ( left.WOthers( ) > right.WOthers( ) );   // note >
     if ( left.Len( ) != right.Len( ) )
       return ( left.Len( ) < right.Len( ) );
     return( left.Id( ) < right.Id( ) );  // for completeness
   }

   friend ostream &operator<< ( ostream &out, const CBadness &bad ) {
     out << "s_" << bad.Id( )
	 << "   " << bad.Len( ) << " bp,"
	 << "   " << bad.WCenter( ) << " overlaps_with_center,"
	 << "   " << bad.WOthers( ) << " overlaps_with_others";
     return out;
   }


 public:

   int id_;        // id of super
   int len_ ;      // length of super
   int wcenter_;   // number of overlaps with center
   int wothers_;   // number of overlaps between non-center supers

 };

/**
  * IsolateSuper
  *
  * Cut all links to and from this super. Notice each super is
  * associated to two vertices (2*super_id, 1 + 2*super_id).
  */
 void IsolateSuper( const int super_id,
		    const vec<superb> &supers,
		    digraphE<CLinkBundle> &graph,
		    ContigsManager &manager,
		    ostream &log,
		    const bool VERBOSE,
		    const bool LIGHT )
 {
   int v1 = 2 * super_id;
   int v2 = 1 + v1;
   if ( VERBOSE )
     log << "isolating s" << super_id
	 << " (" << supers[super_id].TrueLength( ) << " bp)\n";
   
   if ( ! LIGHT )
     for (int ii=0; ii<supers[super_id].Ntigs( ); ii++)
       manager.IsolateContig( supers[super_id].Tig( ii ) );
			      
   graph.DeleteEdgesAtVertex( v1 );
   graph.DeleteEdgesAtVertex( v2 );
   //graph.RemoveDeadEdgeObjects( );
 }

 /**
  * CutContig
  *
  * It assumes the given super is a singleton, and it cuts the contig
  * of the super at the two specified cut points. Data structures will
  * be updated.
  */
 void CutContig( int super_id,
		 int cutter_1,
		 int cutter_2,
		 vec<superb> &supers,
		 digraphE<CLinkBundle> &graph,
		 ContigsManager &manager,
		 ostream &log,
		 bool VERBOSE )
 {
   ForceAssertEq( supers[super_id].Ntigs( ), 1 );
   size_t contig_id = supers[super_id].Tig( 0 );
   int contig_len = supers[super_id].Len( 0 );

   // Split fasta and update aligns.
   vec<size_t> new_ids = manager.SplitContig( contig_id, cutter_1, cutter_2 );
   const vec<fastavector> &contigs = manager.Contigs( );

   // Update supers.
   ForceAssertEq( new_ids[0], contig_id );
   supers[super_id].SetLen( 0, contigs[contig_id].size( ) );

   vec<int> new_super_ids;
   for (int ii=1; ii<new_ids.isize( ); ii++) {
     int new_contig_id = new_ids[ii];
     int new_contig_len = contigs[new_contig_id].size( );
     superb new_super;
     new_super.PlaceFirstTig( new_contig_id, new_contig_len );
     supers.push_back( new_super );
     new_super_ids.push_back( supers.isize( ) - 1 );
   }

   // Update graph (two vertices per new super).
   int nadd = 2 * ( new_ids.isize( ) - 1 );
   graph.AddVertices( nadd );

   vec<int> from_ids = graph.From( 2 * super_id );
   for (int ii=0; ii<from_ids.isize( ); ii++) {
     int newv = 2 * new_super_ids.back( );
     int oldv = 2 * super_id;
     int w = from_ids[ii];
     vec<int> all_edge_ids = graph.EdgesBetween( oldv, w );
     for (int jj=0; jj<all_edge_ids.isize( ); jj++)
       graph.GiveEdgeNewFromVx( all_edge_ids[jj], oldv, newv );
   }

   vec<int> to_ids = graph.To( 1 + 2 * super_id );
   for (int ii=0;  ii<to_ids.isize( ); ii++) {
     int neww = 1 + 2 * new_super_ids.back( );
     int oldw = 1 + 2 * super_id;
     int v = to_ids[ii];
     vec<int> all_edge_ids = graph.EdgesBetween( v, oldw );
     for (int jj=0; jj<all_edge_ids.isize( ); jj++)
       graph.GiveEdgeNewToVx( all_edge_ids[jj], oldw, neww );
   }

   // Log event.
   if ( VERBOSE ) {
     log << "cutting s" << super_id
	 << "=c" << contig_id
	 << " (" << contig_len
	 << " bp) into: ";
     for (int ii=0; ii<new_ids.isize( ); ii++) {
       log << "  s" << ( ii == 0 ? super_id : new_super_ids[ii-1] )
	   << "=c" << new_ids[ii]
	   << " (" << contigs[ new_ids[ii] ].size( )
	   << " bp);" << ( ii == new_ids.isize( ) - 1 ? "\n" : "" );
     }
   }

 }	

 /**
  * BreakSuper
  *
  * Try to break a super at misassembly point(s). If these cannot be
  * reliably determined, just isolate the super from all the others.
  * It returns true on success (isolated - but not broken - contigs are
  * not considered a success).
  */
 bool BreakSuper( const int super_id,
		  const int MIN_LINKS,
		  vec<superb> &supers,
		  digraphE<CLinkBundle> &graph,
		  ContigsManager &manager,
		  ostream &log,
		  const vec<CSloc> &locs,
		  const bool VERBOSE,
		  const bool LIGHT )
 {
   // Arguments.
   const int SMALL_GAP = 1000;   // Used to merge close spread windows

   // Find scaffold locs (these are returned as a sorted vector).
   int slen = supers[super_id].TrueLength( );

   // Detect spread windows.
   vec<ho_interval> wins;
   {
     vec<ho_interval> all_wins;
     all_wins.reserve( locs.size( ) );
     for (size_t ii=0; ii<locs.size( ); ii++) {
       ho_interval new_win = locs[ii].Spread( );
       if ( new_win.Start( ) != 0 || new_win.Stop( ) != 0 ) 
	 all_wins.push_back( new_win );
     }
     sort( all_wins.begin( ), all_wins.end( ) );

     wins.reserve( all_wins.size( ) );
     if ( all_wins.size( ) > 0 ) wins.push_back( all_wins[0] );
     for (size_t ii=1; ii<all_wins.size( ); ii++) {
       if ( all_wins[ii].Start( ) - SMALL_GAP - wins.back( ).Stop( ) < 0 ) {
	 int new_stop = Max( wins.back( ).Stop( ), all_wins[ii].Stop( ) );
	 wins[wins.size( )-1].SetStop( new_stop );
	 continue;
       }
       wins.push_back( all_wins[ii] );
     }
   }  

   // Cannot find cutting points: isolate super.
   if ( wins.size( ) < 2 ) {
     IsolateSuper( super_id, supers, graph, manager, log, VERBOSE, LIGHT );
     return false;
   }

   // For now skip non-singleton supers.
   if ( supers[super_id].Ntigs( ) != 1 ) {
     if ( VERBOSE )
       log << "ERROR! super_" << super_id << ": this is not a singleton\n";
     IsolateSuper( super_id, supers, graph, manager, log, VERBOSE, LIGHT );
     return false;
   }

   // Unspecified problem with the cutter points.
   int cutter_1 = wins[0].Stop( );
   int cutter_2 = wins.back( ).Start( );
   if ( cutter_1 < 0 || cutter_1 >= slen || cutter_2 < 0 || cutter_2 >= slen ) {
     if ( VERBOSE )
       log << "ERROR! super_" << super_id << ": problem with cutter points\n";
     IsolateSuper( super_id, supers, graph, manager, log, VERBOSE, LIGHT );
     return false;
   }

   if ( LIGHT ) {
     IsolateSuper( super_id, supers, graph, manager, log, VERBOSE, LIGHT );
     return false;
   }   

   CutContig( super_id, cutter_1, cutter_2, supers, graph, manager,
	      log, VERBOSE );
	      
   return true;
   
 }

/**
 * CleanseScaffoldGraph
 */
int CleanseScaffoldGraph( digraphE<CLinkBundle> &graph,
			  vec<superb> &supers,
			  ContigsManager &manager,
			  ostream &out,
			  const bool VERBOSE,
			  const bool LIGHT,
			  const int MIN_LINKS,
			  const int MAX_OVERLAP )
{
  int ncuts = 0, nisolates = 0;
  vec< vec<CSloc> > locs_cache;
  out << Date( ) << ": starting CleanseScaffoldGraph" << endl;
  
  u_int bad_count;

  do {
    u_int n_scaffolds = supers.isize();
    if (locs_cache.size() != n_scaffolds)
      locs_cache.resize(n_scaffolds);

    // compute initial location neighborhoods centered on each super
    for (u_int s1=0; s1<n_scaffolds; s1++) {
      // Find scaffold locs with center s1.
      BuildScaffoldLocs( locs_cache[s1], s1, graph, supers, MIN_LINKS );
    }

    // Bad events.
    vec<CBadness> badness;

    for (u_int s1=0; s1<n_scaffolds; s1++) {
      vec<CSloc> &locs = locs_cache[s1];
      int len1 = supers[s1].TrueLength( );
      // Look for overlaps.
      for (int ii=0; ii<locs.isize( ); ii++) {
	for (int jj=ii+1; jj<locs.isize( ); jj++) {
	  if ( locs[ii].OverlapWith( locs[jj] ) < MAX_OVERLAP ) break;
	  if ( badness.size( ) < 1 || badness.back( ).Id( ) != (int)s1 )
	    badness.push_back( CBadness( s1, len1, 0, 0 ) );
	  if ( locs[ii].Begin( ) == 0 || locs[jj].Begin( ) == 0 )
	    badness[ badness.size( ) - 1 ].IncrementWCenter( );
	  else
	    badness[ badness.size( ) - 1 ].IncrementWOthers( );
	}
      }
    }
    sort( badness.begin( ), badness.end( ) );
    bad_count = badness.size();

    vec<bool> used(n_scaffolds, false);

    for (u_int b = 0; b < bad_count; ++b) {
         
      // Isolate or break bad edge.
      int sid = badness[b].Id( );
      if (used[sid]) continue;

      // mark the immediate neighborhood around this scaffold as used
      // during this iteration, so they aren't considered before recomputing
      used[sid] = true;
      {
	int v0 = sid*2;
	int v1 = v0 + 1;
	vec<int> links = graph.From(v0);
	links.append(graph.To(v0));
	links.append(graph.From(v1));
	links.append(graph.To(v1));
	for (u_int i = 0; i < links.size(); ++i)
	  used[links[i]/2] = true;
      }

      u_int nscaff_old = supers.size();

      // do the work (isolate or cut)
      if ( BreakSuper( sid, MIN_LINKS, supers, graph, manager,
		       out, locs_cache[sid], VERBOSE, LIGHT ) )
	++ncuts;
      else
	++nisolates;
     
      // if new scaffolds were created, mark them as used for this iteration
      u_int nscaff_new = supers.size();
      if (nscaff_new > nscaff_old) {
	used.resize(nscaff_new);
	for (u_int s = nscaff_old; s < nscaff_new; ++s)
	  used[s] = true;
      }

    }
    // Clean up graph
    graph.RemoveDeadEdgeObjects( );

    if (VERBOSE) PRINT2(ncuts, nisolates);
  } while (bad_count > 0);

  out << Date( ) << ": CleanseScaffoldGraph done ("
      << ncuts << " cuts)"
      << endl;

  // Return.
  return ncuts;
}

