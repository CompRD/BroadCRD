///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// LocalView.  Given some contigs, generate a view of the local data.
// Work in progress!
//
// Note that this can be very slow.  Nearly all the time is spent loading reads.
// This become very fast once these data are in cache.

#include "Alignment.h"
#include "Basevector.h"
#include "Equiv.h"
#include "FastIfstream.h"
#include "FeudalMimic.h"
#include "MainTools.h"
#include "ParseSet.h"
#include "Superb.h"
#include "graph/Digraph.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"
#include "paths/ReadLoc.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/Unipath.h"
#include "util/ReadTracker.h"

template<class T> ho_interval Bounds( const vec<T>& v, const T& x )
{    int low = lower_bound( v.begin( ), v.end( ), x ) - v.begin( );
     int high = upper_bound( v.begin( ), v.end( ), x ) - v.begin( );
     return ho_interval( low, high );    }

int main(int argc, char *argv[])
{
     RunTime( );
  
     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(SUBDIR, "test");
     CommandArgument_String_OrDefault(ASSEMBLY, "linear_scaffolds0.clean.patched");
     CommandArgument_String_OrDefault_Doc(TIGS, "", 
          "if unspecified, process all contigs;"
          " otherwise it is one of the following: \n"
          "(a) a list of contig ids (in ParseIntSet format) or \n"
          "(b) the letter s followed by a list of scaffolds or \n"
          "(c) s<scaffold id>.<list of indices of contigs in the scaffold");
     CommandArgument_String_OrDefault_Doc(READLOCS_PREFIX, "",
          "if specified, file extension is .READLOCS_PREFIX.readlocs "
          "instead of .readlocs");
     CommandArgument_Bool_OrDefault(PARTNERS, True);
     CommandArgument_Bool_OrDefault_Doc(ADD_CONTIGS, True,
          "add the contig sequence to the graph");
     CommandArgument_Int_OrDefault(K, 80);
     CommandArgument_Int_OrDefault(MIN_COMPONENT, 1000);
     CommandArgument_Bool_OrDefault(USE_CORRECTED_JUMPS, True);
     CommandArgument_Bool_OrDefault(REMOVE_HANGING_ENDS, True);
     CommandArgument_String_Doc(TMP, "name of a temporary file");

     // Output control.

     CommandArgument_String(DOT);
     CommandArgument_String(COMP_DOT);
     CommandArgument_Bool_OrDefault_Doc(DOT_LABEL_EDGES, True,
          "number graph edges");
     CommandArgument_String_Doc(EDGE_FILE, "generate fasta edge file of this name");

     EndCommandArguments;

     // Define directories.

     String data_dir = PRE + "/" + DATA, run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR; 

     // Parse TIGS.

     vec<int> tigs;
     {    if (TIGS == "") 
          {    int n_tigs = MastervecFileObjectCount( sub_dir + "/" + ASSEMBLY
                    + ".contigs.fastb" );
               for ( int j = 0; j < n_tigs; j++ )
                    tigs.push_back(j);    }
          else if (TIGS.Contains("s", 0)) 
          {    TIGS = TIGS.After("s");
               vec<superb> scaffolds;
               ReadSuperbs(sub_dir + "/" + ASSEMBLY + ".superb", scaffolds);
               if (TIGS.Contains(".")) 
               {    int scaffold = TIGS.Before(".").Int();
                    ForceAssertLt(scaffold, scaffolds.isize());
                    vec<int> spos;
                    ParseIntSet(TIGS.After("."), spos);
                    for (int j = 0; j < spos.isize(); j++)
                         tigs.push_back(scaffolds[scaffold].Tig(spos[j]));    }
               else 
               {    vec<int> s;
                    ParseIntSet(TIGS, s);
                    for (int i = 0; i < s.isize(); i++) 
                    {    int scaffold = s[i];
                         ForceAssertLt(scaffold, scaffolds.isize());
                         for (int j = 0; j < scaffolds[scaffold].Ntigs(); j++)
                              tigs.push_back(scaffolds[scaffold].Tig(j));
                              }    }    }
          else ParseIntSet(TIGS, tigs);    }

     // Print the contigs.

     String head = sub_dir + "/" + ASSEMBLY;
     vecbasevector fewtigs;
     fewtigs.Read( head + ".contigs.fastb", tigs );
     for ( int i = 0; i < tigs.isize( ); i++ )
          fewtigs[i].Print( cout, "contig_" + ToString( tigs[i] ) );

     // Figure out which reads we need.

     if ( READLOCS_PREFIX != "" ) head += "." + READLOCS_PREFIX;
     read_locs_on_disk locs_file( head, run_dir );
     vec< pair<int64_t,Bool> > frag_ids, jump_ids, long_ids;
     map<int64_t,int64_t> frag_to_partner, jump_to_partner, long_to_partner;
     for ( int it = 0; it < tigs.isize( ); it++ )
     {    vec<read_loc> locs;
          locs_file.LoadContig( tigs[it], locs );
          for ( int j = 0; j < locs.isize( ); j++ )
          {    const read_loc& rl = locs[j];
               int64_t id1 = rl.ReadId( ), id2 = rl.PartnerReadId( );
               if ( rl.Frag( ) )
               {    frag_ids.push( id1, rl.Fw( ) );
                    if (PARTNERS) frag_ids.push( id2, rl.Rc( ) );
                    frag_to_partner[id1] = id2;
                    frag_to_partner[id2] = id1;    }
               if ( rl.Jump( ) )
               {    jump_ids.push( id1, rl.Fw( ) );
                    if (PARTNERS) jump_ids.push( id2, rl.Rc( ) );
                    jump_to_partner[id1] = id2;
                    jump_to_partner[id2] = id1;    }
               if ( rl.LongJump( ) )
               {    long_ids.push( id1, rl.Fw( ) );
                    if (PARTNERS) long_ids.push( id2, rl.Rc( ) );
                    long_to_partner[id1] = id2;
                    long_to_partner[id2] = id1;    }    }    }
     UniqueSort(frag_ids), UniqueSort(jump_ids), UniqueSort(long_ids);
     vec<int64_t> frag_ids0, jump_ids0, long_ids0;
     for ( int i = 0; i < frag_ids.isize( ); i++ )
          frag_ids0.push_back( frag_ids[i].first );
     for ( int i = 0; i < jump_ids.isize( ); i++ )
          jump_ids0.push_back( jump_ids[i].first );
     for ( int i = 0; i < long_ids.isize( ); i++ )
          long_ids0.push_back( long_ids[i].first );

     // Map jump reads to corrected jump reads.

     vec<int64_t> ec_ids;
     vec<int> ec_ids_local;
     if (USE_CORRECTED_JUMPS)
     {    ReadTracker rt;
          int64_t njumps_ec 
               = MastervecFileObjectCount( run_dir + "/jump_reads_ec.fastb" );
          cout << Date( ) << ": loading read tracker" << endl;
          rt.Load( run_dir + "/jump_reads_ec" );
          cout << Date( ) << ": done" << endl;
          for ( int64_t ec_id = 0; ec_id < njumps_ec; ec_id++ )
          {    int64_t filt_id = rt.GetReadIndex(ec_id);
               int p = BinPosition( jump_ids0, filt_id );
               if ( p >= 0 ) 
               {    ec_ids.push_back(ec_id);
                    ec_ids_local.push_back(p);    }    }    }

     // Load the reads.

     vecbasevector frag_reads, jump_reads, long_reads;
     const Bool FRAG_EDIT = True;
     if ( frag_ids.size( ) > 0 )
     {    if ( !FRAG_EDIT )
               frag_reads.Read( run_dir + "/frag_reads_filt_cpd.fastb", frag_ids0 );
          else frag_reads.Read( run_dir + "/frag_reads_edit.fastb", frag_ids0 );    }
     for ( int i = 0; i < frag_ids.isize( ); i++ )
          if ( !frag_ids[i].second ) frag_reads[i].ReverseComplement( );
     if ( jump_ids.size( ) > 0 )
     {    jump_reads.Read( run_dir + "/jump_reads_filt_cpd.fastb", jump_ids0 );
          if (USE_CORRECTED_JUMPS)
          {    vecbasevector jump_reads_ec;
               cout << Date( ) << ": loading ec jumps" << endl;
               jump_reads_ec.Read( run_dir + "/jump_reads_ec.fastb", ec_ids );
               cout << Date( ) << ": done" << endl;
               for ( int i = 0; i < ec_ids.isize( ); i++ )
                    jump_reads[ ec_ids_local[i] ] = jump_reads_ec[i];    }    }
     for ( int i = 0; i < jump_ids.isize( ); i++ )
          if ( !jump_ids[i].second ) jump_reads[i].ReverseComplement( );
     if ( long_ids.size( ) > 0 )
          long_reads.Read( run_dir + "/long_jump_reads_filt.fastb", long_ids0 );
     for ( int i = 0; i < long_ids.isize( ); i++ )
          if ( !long_ids[i].second ) long_reads[i].ReverseComplement( );

     // Start to build graph.

     vecbasevector all(frag_reads);
     all.Append(jump_reads), all.Append(long_reads);
     if (ADD_CONTIGS) all.Append(fewtigs);

     // Generate unipath graph.

     vecKmerPath paths, paths_rc, unipaths;
     vec<tagged_rpint> pathsdb, unipathsdb;
     ReadsToPathsCoreY( all, K, paths );
     CreateDatabase( paths, paths_rc, pathsdb );
     Unipath( paths, paths_rc, pathsdb, unipaths, unipathsdb );
     digraph A;
     BuildUnipathAdjacencyGraph( paths, paths_rc, pathsdb, unipaths, unipathsdb, A );
     HyperKmerPath h;
     BuildUnipathAdjacencyHyperKmerPath( K, A, unipaths, h );

     // Clean up graph.

     if (REMOVE_HANGING_ENDS) RemoveHangingEnds( h, &KmerPath::KmerCount, 250, 5.0 );
     h.RemoveSmallComponents(MIN_COMPONENT);
     h.RemoveDeadEdgeObjects( );
     h.RemoveEdgelessVertices( );
     h.RemoveUnneededVertices( );
     h.RemoveDeadEdgeObjects( );

     // Dump edges.

     if ( EDGE_FILE != "" )
     {    KmerBaseBroker kbb( K, paths, paths_rc, pathsdb, all );
          Ofstream( out, EDGE_FILE );
          for ( int i = 0; i < h.EdgeObjectCount( ); i++ )
               kbb.Seq( h.EdgeObject(i) ).Print( out, "edge_" + ToString(i) );    }

     // Figure out which reads lie on which edges.

     vec< pair<int64_t,int> > frag_read_edge, jump_read_edge, long_read_edge;
     for ( int e = 0; e < h.EdgeObjectCount( ); e++ )
     {    const KmerPath& p = h.EdgeObject(e);
          for ( int i = 0; i < p.NSegments( ); i++ )
          {    vec<longlong> con;
               Contains( pathsdb, p.Segment(i), con );
               for ( int j = 0; j < con.isize( ); j++ )
               {    const tagged_rpint& t = pathsdb[ con[j] ];
                    int64_t id = t.ReadId( );
                    if ( id < (int64_t) frag_ids.size( ) )
                    {    frag_read_edge.push( frag_ids[id].first, e );
                         continue;    }
                    id -= frag_ids.size( );
                    if ( id < (int64_t) jump_ids.size( ) )
                    {    jump_read_edge.push( jump_ids[id].first, e );
                         continue;    }
                    id -= jump_ids.size( );
                    if ( id < (int64_t) long_ids.size( ) )
                         long_read_edge.push( long_ids[id].first, e );    }    }    }
     UniqueSort(frag_read_edge);
     UniqueSort(jump_read_edge); 
     UniqueSort(long_read_edge);
     vec<int64_t> frag1, jump1, long1;
     for ( int i = 0; i < frag_read_edge.isize( ); i++ )
          frag1.push_back( frag_read_edge[i].first );
     for ( int i = 0; i < jump_read_edge.isize( ); i++ )
          jump1.push_back( jump_read_edge[i].first );
     for ( int i = 0; i < long_read_edge.isize( ); i++ )
          long1.push_back( long_read_edge[i].first );

     // Show read locations.
                    
     vec< vec<int> > link_forward( h.EdgeObjectCount( ) );
     vec< vec<int> > link_backward( h.EdgeObjectCount( ) );
     for ( int e = 0; e < h.EdgeObjectCount( ); e++ )
     {    const KmerPath& p = h.EdgeObject(e);
          cout << "\nedge " << e << "\n";
          for ( int i = 0; i < p.NSegments( ); i++ )
          {    vec<longlong> con;
               Contains( pathsdb, p.Segment(i), con );
               for ( int j = 0; j < con.isize( ); j++ )
               {    const tagged_rpint& t = pathsdb[ con[j] ];
                    int64_t id = t.ReadId( );
                    if ( id < (int64_t) frag_ids.size( ) )
                    {    int64_t ID = frag_ids[id].first;
                         cout << "frag_" << ToString(ID) << "_" 
                              << ( frag_ids[id].second ? "fw" : "rc" ) << " (";
                         int64_t PID = frag_to_partner[ID];
                         ho_interval b = Bounds(frag1, ID), bp = Bounds(frag1, PID);
                         for ( int l = b.Start( ); l < b.Stop( ); l++ )
                              cout << frag_read_edge[l].second << " ";
                         cout << "/";
                         for ( int l = bp.Start( ); l < bp.Stop( ); l++ )
                         {    int e2 = frag_read_edge[l].second;
                              cout << " " << e2;
                              if ( frag_ids[id].second ) 
                                   link_forward[e].push_back(e2);
                              else link_backward[e].push_back(e2);    }
                         cout << ")\n";
                         continue;    }
                    id -= frag_ids.size( );
                    if ( id < (int64_t) jump_ids.size( ) )
                    {    int64_t ID = jump_ids[id].first;
                         cout << "jump_" << ToString(ID) << "_" 
                              << ( jump_ids[id].second ? "fw" : "rc" ) << " (";
                         int64_t PID = jump_to_partner[ID];
                         ho_interval b = Bounds(jump1, ID), bp = Bounds(jump1, PID);
                         for ( int l = b.Start( ); l < b.Stop( ); l++ )
                              cout << jump_read_edge[l].second << " ";
                         cout << "/";
                         for ( int l = bp.Start( ); l < bp.Stop( ); l++ )
                         {    int e2 = jump_read_edge[l].second;
                              cout << " " << e2;
                              if ( jump_ids[id].second ) 
                                   link_forward[e].push_back(e2);
                              else link_backward[e].push_back(e2);    }
                         cout << ")\n";
                         continue;    }
                    id -= jump_ids.size( );
                    if ( id < (int64_t) long_ids.size( ) )
                    {    int64_t ID = long_ids[id].first;
                         cout << "long_" << ToString(ID) << "_" 
                              << ( long_ids[id].second ? "fw" : "rc" ) << " (";
                         int64_t PID = long_to_partner[ID];
                         ho_interval b = Bounds(long1, ID), bp = Bounds(long1, PID);
                         for ( int l = b.Start( ); l < b.Stop( ); l++ )
                              cout << long_read_edge[l].second << " ";
                         cout << "/";
                         for ( int l = bp.Start( ); l < bp.Stop( ); l++ )
                         {    int e2 = long_read_edge[l].second;
                              cout << " " << e2;
                              if ( long_ids[id].second ) 
                                   link_forward[e].push_back(e2);
                              else link_backward[e].push_back(e2);    }
                         cout << ")\n";    }    }    }
          UniqueSort( link_forward[e] ), UniqueSort(link_backward[e] );
          cout << "edge " << e << " links forward to";
          for ( int j = 0; j < link_forward[e].isize( ); j++ )
               cout << " " << link_forward[e][j];
          cout << "\nedge " << e << " links backward to";
          for ( int j = 0; j < link_backward[e].isize( ); j++ )
               cout << " " << link_backward[e][j];
          cout << "\n";    }

     // Find links between components.

     equiv_rel e;
     h.ComponentRelation(e);
     vec<int> reps;
     e.OrbitRepsAlt(reps);
     vec<int> to_left, to_right;
     h.ToLeft(to_left), h.ToRight(to_right);
     vec< pair<int,int> > links;
     for ( int e1 = 0; e1 < h.EdgeObjectCount( ); e1++ )
     {    int v1 = to_left[e1];
          int c1 = BinPosition( reps, e.ClassId(v1) );
          for ( int j = 0; j < link_forward[e1].isize( ); j++ )
          {    int e2 = link_forward[e1][j];
               int v2 = to_left[e2];
               int c2 = BinPosition( reps, e.ClassId(v2) );
               if ( c2 != c1 ) links.push( c1, c2 );    }    }
     UniqueSort(links);
     cout << "\n";
     for ( int i = 0; i < links.isize( ); i++ )
     {    cout << "contig " << links[i].first << " --> contig "
               << links[i].second << "\n";    }

     // Build component overlap graph.

     vec< pair<int,int> > edge_overlaps;
     SystemSucceed( "CmpSeq PRE=. FILE1=" + EDGE_FILE + " FILE2=" 
          + EDGE_FILE + " OUT_IMPROPER_MATCHES=True MIN_OVERLAP=20 FW_ONLY=True "
          + "BINARY_ALIGNMENTS_FILE=" + TMP + " > /dev/null" );
     vec<alignment_plus> aligns;
     READX( TMP, aligns );
     vec< vec<int> > from( reps.size( ) ), to( reps.size( ) );
     for ( int i = 0; i < aligns.isize( ); i++ )
     {    const alignment_plus& ap = aligns[i];
          int id1 = ap.Id1( ), id2 = ap.Id2( );
          int v1 = to_right[id1], v2 = to_left[id2];
          int c1 = BinPosition( reps, e.ClassId(v1) );
          int c2 = BinPosition( reps, e.ClassId(v2) );
          if ( c1 == c2 ) continue;
          vec<int> o1, o2;
          e.Orbit( reps[c1], o1 ), e.Orbit( reps[c2], o2 );
          if ( !h.Sink(v1) || !h.Source(v2) ) continue;

          if ( o1.size( ) == 2 && ap.pos2( ) > 0 ) continue;
          if ( o2.size( ) == 2 && 
               ap.Pos1( ) != h.EdgeObject(id1).KmerCount( ) + K - 1 ) 
          {    continue;    }

          /*
          ... -----------------
                      ------------------- ...
          */

          if ( o1.size( ) ==2 && o2.size( ) == 2 )
          {    if ( ap.pos1( ) == 0 
                    && ap.Pos2( ) == h.EdgeObject(id2).KmerCount( ) + K - 1 )
               {    continue;    }    }
          from[c1].push_back(c2), to[c2].push_back(c1);    }
     for ( int i = 0; i < reps.isize( ); i++ )
     {    UniqueSort( from[i] );
          UniqueSort( to[i] );    }
     digraph G( from, to );
     Ofstream( yout, COMP_DOT );
     G.PrettyDOT( yout, False, True );

     // Generate dot file.

     Ofstream( out, DOT );
     const Bool DOT_LABEL_CONTIGS = True;
     const Bool DOT_LABEL_VERTICES = False;
     h.PrintSummaryDOT0w( out, DOT_LABEL_CONTIGS, DOT_LABEL_VERTICES,
          DOT_LABEL_EDGES );    }
