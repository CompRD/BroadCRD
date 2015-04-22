///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// This is an experimental version of Snorgle phase 2.  Needs to be kept in sync
// with it.
//
// A mess, could be simplified.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "MainTools.h"
#include "Qualvector.h"
#include "VecUtilities.h"
#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/FriendAligns.h"
#include "paths/long/LoadCorrectCore.h"
#include "paths/long/ReadPath.h"
#include "paths/long/ReadStack.h"
#include "paths/long/large/explore/SnorgleTools.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_Int_OrDefault_Doc(TRACE_PID, -1,
          "for phase 2, analyze just this PID; also follow in phase 1");
     CommandArgument_Int_OrDefault_Doc(TRACE_RID, -1,
          "sets TRACE_PID = TRACE_RID/2");
     CommandArgument_String_OrDefault(INSTANCE, "51400.newchem");
     CommandArgument_Bool_OrDefault(PRINT_PATHS, False);
     EndCommandArguments;

     if ( TRACE_RID >= 0 ) TRACE_PID = TRACE_RID/2;

     ForceAssert( TRACE_PID >= 0 );

     String work_dir = "/wga/scr4/jaffe/GapToy/" + INSTANCE;

     {    Bool FP_LOGGING = True;
          Bool PRINT_FRIENDS = True;
          Bool PRINT_STACKS = True;
          Bool PRINT_JSTACK = True;

          int64_t id1 = TRACE_PID * 2;
          int64_t id2 = id1 + 1;

          cout << Date( ) << ": reading friends" << endl;
          vec< quad<int64_t,int64_t,int,Bool> > xfriends;
          vec<int64_t> ss;
          vec<int> subset = {int(id1),int(id2)+1}; // should use int64_t!!
          // ss.ReadSubset( work_dir + "/xfriends.index", subset ); doesn't work!?
          vec<int> subset1 = {int(id1)}, subset2 = {int(id2)+1};
          ss.ReadSubset( work_dir + "/xfriends.index", subset1 );
          ss.ReadSubset( work_dir + "/xfriends.index", subset2, True );
          PRINT2( ss[0], ss[1] );
          xfriends.ReadRange( work_dir + "/xfriends", ss[0], ss[1] );

          cout << Date( ) << ": defining ids" << endl;
          vec<int64_t> all;
          all.push_back( id1, id2 );
          for ( int i = 0; i < (int) xfriends.size( ); i++ )
               all.push_back( xfriends[i].second );
          UniqueSort(all);
          PRINT( all.size( ) );

          cout << Date( ) << ": setting up virtuals" << endl;
          vecbasevector basesx;
          VecPQVec pqualsx;
          basesx.Read( work_dir + "/data/frag_reads_orig.fastb", all );
          pqualsx.Read( work_dir + "/data/frag_reads_orig.qualp", all );

          // Set up bins and closures.

          cout << Date( ) << ": setting up bins and closures" << endl;
          const int64_t bins = 1;
          vec<int64_t> starts( bins + 1 );
          for ( int j = 0; j <= bins; j++ )
               starts[j] = ( j * xfriends.jsize( ) ) / bins;
          for ( int j = bins - 1; j > 0; j-- )
          {    while( starts[j] > 0 && xfriends[ starts[j] ].first/2 
                    == xfriends[ starts[j]-1 ].first/2 )
               {    starts[j]--;    }    }
          vec<vec<basevector>> closures(bins);
          vec<vec<int64_t>> closure_ids(bins);

          // Go through the bins.

          cout << Date( ) << ": starting main loop" << endl;
          int64_t m = 0;
          {    int64_t start = starts[m], stop = starts[m+1];
               int64_t last_id = -1;
               readstack last_stack, stack, jstack;
               PairsManager pairs;
               vecbasevector bases;
               vecqualvector quals;
               Friends aligns;
               vec<Bool> suspect;
               basevector last_con, con, jcon;
               vec<basevector> p;
               vec<int> ids, off;
               for ( int64_t i = start; i < stop; i++ )
               {    // Reminder: xfriends[i] is data (id1,id2,offset,fw2).
                    int64_t ii = i;
                    int id = xfriends[ii].first;
                    if ( TRACE_PID >= 0 && id/2 != TRACE_PID ) continue;

                    // Restrict to those reads for which id1 = id.

                    int64_t j;
                    for ( j = i + 1; j < stop; j++ )
                         if ( xfriends[j].first != id ) break;

                    {    cout << "\n" << xfriends[i].first << " has " << j - i 
                              << " friends" << endl;    }
          
                    // if ( j - i <= 20 )
                    if (PRINT_FRIENDS)
                    {    for ( int64_t k = i; k < j; k++ )
                         {    cout << ( xfriends[k].fourth ? "+" : "-" )
                                   << xfriends[k].second 
                                   << "." << xfriends[k].third << endl;    }    }
     
                    // Load reads.

                    ids = {id};
                    for ( int64_t k = i; k < j; k++ )
                         ids.push_back( xfriends[k].second );
                    bases.clear( );
                    quals.clear( );

                    for ( int k = 0; k < ids.isize( ); k++ )
                    {    int64_t bid = BinPosition( all, ids[k] );
                         bases.push_back( basesx[bid] );
                         qvec q;
                         pqualsx[bid].unpack( &q );
                         quals.push_back(q);    }

                    // Cap quality scores.

                    vec<Bool> done( quals.size( ), False );
                    CapQualityScores( quals, done );

                    // Build read stack.

                    aligns.clear( );
                    for ( int64_t k = i; k < j; k++ )
                    {    aligns.push_back( Friend( k - i + 1, xfriends[k].third, 
                              !xfriends[k].fourth ) );    }
                    stack.Initialize( 0, aligns, 0, aligns.size( ), 
                         readstack::right_extended, bases, quals, pairs, False );
                    for ( int64_t k = i; k < j; k++ )
                    {    stack.SetId( k-i, ids[k-i] );
                         stack.SetPid( k-i, ids[k-i]/2 );    }
                    stack.Raise1(0);
                    stack.MotifDiff(1,suspect);
                    stack.Erase(suspect);

                    const int q_solid = 30;
                    stack.HighQualDiff( q_solid, 1, suspect );
                    PRINT2( int(suspect[0]), int(suspect[1]) );
                    stack.Erase(suspect);
     
                    if ( xfriends[i].first % 2 == 1 ) stack.Reverse( );
     
                    if (PRINT_STACKS) stack.Print(cout);
     
                    // Look for overlaps.
          
                    con = stack.Consensus1( );
                    if ( xfriends[i].first % 2 == 1 
                         && xfriends[i].first == last_id + 1 )
                    {    const int offset_verb = 3;
                         off = GetOffsets1( last_stack, stack, offset_verb, 0 );
                         cout << "off = " << printSeq(off) << endl;
                         const int max_offsets = 2;
                         if ( TRACE_PID >= 0 
                              && ( off.empty( ) || off.isize( ) > max_offsets ) )
                         {    cout << "\n";
                              con.Print( cout, "consensus" );
                              basevector conl = last_stack.Consensus1( );
                              conl.Print( cout, "last_consensus" );    
                              cout << "\n";    }

                         if ( off.isize( ) <= max_offsets )
                         {    for ( int o = 0; o < off.isize( ); o++ )
                              {    cout << "\nusing offset " << off[o] << endl;
                                   jstack = last_stack;
                                   jstack.Merge( stack, off[o] );

                                   // What does this do and does it assume
                                   // reads have same length?

                                   if ( -2*off[o] >= jstack.Cols( ) ) continue;

                                   // Trim so that stack does not extend beyond
                                   // the originating fragment.

                                   int left = 0, right = jstack.Cols( );
                                   if ( off[o] < 0 ) left = -off[o];
                                   int ext_right = last_stack.Cols( ) -
                                        ( off[o] + stack.Cols( ) );
                                   if ( ext_right > 0 ) right -= ext_right;
                                   jstack.Trim( left, right );

                                   // Ignore very short fragments.

                                   if ( jstack.Cols( ) < 100 ) continue;

                                   jstack.SortByPid( 
                                        jstack.Pid(0), 0, last_stack.Rows( ) );
                                   jstack.Unique( );
                                   // jstack.Raise1(0), jstack.Raise1(1);
                                   jstack.HighQualDiff( q_solid, 2, suspect );
                                   // if ( suspect[0] || suspect[1] ) continue;
                                   {    PRINT2( int(suspect[0]), int(suspect[1]) );
                                        PRINT( jstack.Rows( ) );    }
                                   jstack.Erase(suspect);
                                   PRINT( jstack.Rows( ) );
               
                                   // Clean stack.
          
                                   for ( int j = 0; j < jstack.Cols( ); j++ )
                                   {    vec<int> qsum( 4, 0 ); 
                                        vec<int> ids( 4, vec<int>::IDENTITY );
                                        for ( int l = 0; l < jstack.Rows( ); l++ )
                                        {    if ( jstack.Def( l, j ) )
                                             {    qsum[ jstack.Base( l, j ) ]
                                                       += jstack.Qual( l, j );    
                                                       }    }
                                        ReverseSortSync( qsum, ids );
                                        if ( qsum[0] >= 100 
                                             && qsum[0] >= 10 * qsum[1]
                                             && qsum[1] < 100 )
                                        {    for ( int l = 0; l < jstack.Rows( ); 
                                                  l++ )
                                             {    if ( jstack.Def( l, j ) 
                                                       && jstack.Base( l, j ) 
                                                            != ids[0] )
                                                  {    jstack.SetBase(l, j, ids[0]);
                                                       jstack.SetQual( l, j, 0 );    
                                                            }    }    }    }
     
                                   // Print stack.
     
                                   if (PRINT_JSTACK) 
                                   {    cout << "\njoint stack:\n";
                                        jstack.Print(cout);    }
                                   jcon = jstack.Consensus1( );
                                   jcon.Print( cout, "joint_consensus" );    
                                   FindPaths( jstack, p, FP_LOGGING );
                                   const int max_paths = 5;
                                   if ( p.isize( ) <= max_paths )
                                   {    for ( int c = 0; c < p.isize( ); c++ )
                                        {    closures[m].push_back( p[c] );    
                                             closure_ids[m].push_back(id);    }    }
                                   {    
                                        #pragma omp critical
                                        {    cout << "id = " << id << ", off[0] = "
                                                  << off[0] << ", found " 
                                                  << p.size( ) << " paths";
                                             if ( p.isize( ) > max_paths )
                                                  cout << " (too many)";
                                             cout << "\n";
                                             if ( p.size( ) <= 12 )
                                             {    for ( int j = 0; j < p.isize( ); 
                                                       j++ )
                                                  {    p[j].Print( cout, 
                                                            "p" + ToString(j+1) );
                                                       }    }    }    }    }    }

                         }

                    last_id = xfriends[i].first;
                    last_con = con;
                    last_stack = stack;
                    i = j - 1;    }    }

          // Done.

          cout << Date( ) << ": done" << endl;
          Scram(0);    }    }
