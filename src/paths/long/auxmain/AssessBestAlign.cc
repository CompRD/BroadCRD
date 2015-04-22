///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// AssessBestAlign: Take the output of CallRefTrace BEST_GLOBAL_OUT=..., make an
// alignment, parse it into "EVENTS", then assess the events using raw read data.

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "PackAlign.h"
#include "PrintAlignment.h"
#include "lookup/LookAlign.h"
#include "math/HoInterval.h"
#include "paths/long/AssessBestAlignCore.h"
#include "paths/long/fosmid/FosmidPool.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "util/SearchFastb2Core.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(BEST_GLOBAL, 
          "fastb file from CallRefTrace BEST_GLOBAL_OUT=...");
     CommandArgument_String_OrDefault_Doc(TMP1, "", 
          "tmp dir for first set of reads");
     CommandArgument_String_OrDefault_Doc(TMP2, "", 
          "tmp dir for second set of reads");
     CommandArgument_String_OrDefault_Doc(NAME1, "", "name for first set of reads");
     CommandArgument_String_OrDefault_Doc(NAME2, "", "name for second set of reads");
     CommandArgument_String_OrDefault_Doc(ANNOTATIONS, "", "file of annotations");
     CommandArgument_Int(X);
     CommandArgument_Bool_OrDefault_Doc(PERFS, False,
          "print perfect interval sizes and exit");
     EndCommandArguments;

     // Load sequences.

     vecbasevector best(BEST_GLOBAL);
     basevector b1 = best[0], b2 = best[1];

     // Horrible criminal kludge to get reference coordinates.  BAAAAAD!

     String chr;
     int ref_start;
     {    vec<String> regions;
          vec< vec< pair<String,String> > > junctions, breaks, edits;
          ParseFosmidPoolMetainfo( regions, junctions, breaks, edits );
          chr = regions[X].Before( ":" );
          int target;
          if ( chr == "X" ) target = 22;
          else target = chr.Int( ) - 1;
          vecbasevector G(1);
          G[0] = b2;
          G.WriteAll( "temptemp.fastb" );
          SystemSucceed( "QueryLookupTable K=12 MM=12 MC=0.15 "
               "L=/wga/scr4/bigrefs/human19/genome.lookup.lookup PARSEABLE=True "
               "SEQS=temptemp.fastb "
               "TARGETS_TO_PROCESS=" + ToString(target) + " > temptemp.aligns" );
          vec<look_align> aligns;
          LoadLookAligns( "temptemp.aligns", aligns );
          ref_start = aligns[0].pos2( );    }

     // Load annotations.

     vec< pair<ho_interval,String> > ann;
     if ( ANNOTATIONS != "" )
     {    fast_ifstream in(ANNOTATIONS);
          String line;
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               int start = line.Before( "-" ).Int( );
               int stop = line.Between( "-", "," ).Int( );
               String a = line.After( ", " );
               ann.push( ho_interval( start, stop ), a );    }    }

     // Generate alignment.

     alignment al;
     SmithWatAffineParallel( b1, b2, al );
     // SmithWatAffine( b1, b2, al );
     // cout << "\nBest alignment of assembly:\n";
     // PrintVisualAlignment( True, cout, b1, b2, al );    
     align a(al);
     if (PERFS)
     {    vec<ho_interval> perf2;
          a.PerfectIntervals2( b1, b2, perf2 );
          for ( int i = 0; i < perf2.isize( ); i++ )
               cout << perf2[i].Length( ) << endl;
          Scram(0);    }
     
     // Decompose alignment into events.

     vec< pair<int,int> > P1, P2;
     DecomposeAlign( a, b1, b2, P1, P2 );

     vecbasevector bb;
     for ( int i = 0; i < P1.isize( ); i++ )
     {    basevector B1( b1, P1[i].first, P1[i].second - P1[i].first );
          basevector B2( b2, P2[i].first, P2[i].second - P2[i].first );
          bb.push_back(B1), bb.push_back(B2);    }
     if ( TMP1 != "" && TMP2 != "" )
     {    bb.WriteAll( TMP1 + "/test.fastb" );
          bb.WriteAll( TMP2 + "/test.fastb" );    }
     const int K = 20;
     vec< vec< triple<int64_t,int64_t,int> > > aligns(2);
     if ( TMP1 != "" && TMP2 != "" )
     {    for ( int pass = 0; pass < 2; pass++ )
          {    String TMP = ( pass == 0 ? TMP1 : TMP2 );
               SearchFastb2( TMP + "/test.fastb", TMP + "/frag_reads_orig.fastb", 
                    K, &aligns[pass], NULL, -1, 0.90, False );    }    }
     int total_subs = 0, total_indels = 0;
     for ( int z = 0; z < P1.isize( ); z++ )
     {    basevector B1( b1, P1[z].first, P1[z].second - P1[z].first );
          basevector B2( b2, P2[z].first, P2[z].second - P2[z].first );
          align A = a.TrimmedTo1( P1[z].first, P1[z].second - P1[z].first );

          int p1 = 0, p2 = 0;
          for ( int j = 0; j < A.Nblocks( ); j++ ) 
          {    if ( A.Gaps(j) > 0 )  
               {    total_indels++;
                    p2 += A.Gaps(j);    }
               if ( A.Gaps(j) < 0 ) 
               {    total_indels++;
                    p1 -= A.Gaps(j);    }
               for ( int x = 0; x < A.Lengths(j); x++ ) 
               {    if ( B1[p1] != B2[p2] ) total_subs++;
                    ++p1;
                    ++p2;    

                    // This must be a workaround for a bug:

                    if ( p2 == B2.isize( ) ) break;    }    }

          vec<int> mgg = A.MutationsGap1Gap2( B1, b2 );
          int start = P2[z].first, stop = P2[z].second;
          cout << "\nMETA-EVENT " << z+1 << " (p2 = " << start << "-" << stop 
               << ", subs = " << mgg[0] << ", dels = "
               << mgg[1] << ", ins = " << mgg[2] << ")\n";
          cout << "chr" << chr << ":" << ref_start + start
               << "-" << ref_start + stop << endl;

          // Print alignment and remove spurious whitespace from it.

          PrintVisualAlignmentClean( False, cout, B1, b2, A );

          // Find evidence.

          if ( TMP1 != "" && TMP2 != "" )
          {    for ( int pass = 0; pass < 2; pass++ )
               {    vec< vec<int> > hits(2);
                    for ( int i = 0; i < aligns[pass].isize( ); i++ )
                    {    if ( aligns[pass][i].first == 2*z )
                              hits[0].push_back( aligns[pass][i].second );
                         if ( aligns[pass][i].first == 2*z + 1 )
                              hits[1].push_back( aligns[pass][i].second );    }
                    for ( int j = 0; j < 2; j++ )
                    {    UniqueSort( hits[j] );
                         cout << hits[j].size( ) << " hits of "
                              << ( j == 0 ? "assembly" : "reference" ) 
                              << " sequence " << "to " 
                              << ( pass == 0 ? NAME1 : NAME2 ) 
                              << " reads\n";    }    }    }

          // Present annotations.

          Bool first = False;
          for ( int i = 0; i < ann.isize( ); i++ )
          {    ho_interval win( start, stop );
               if ( Overlap( ann[i].first, win ) > 0 )
               {    if (first) cout << "\n";
                    first = True;
                    cout << ann[i].first << ", " << ann[i].second 
                         << endl;    }    }    }

     // Summarize results.

     cout << "\ntotal meta-events: " << P1.size( ) << endl;
     cout << "total indel events: " << total_indels << endl;
     cout << "total substitution events: " << total_subs << endl;    }
