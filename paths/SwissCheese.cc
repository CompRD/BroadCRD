///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// SwissCheese.  Make random holes in an assembly.  The maximum number N of holes is
// specified, as well as the HOLE_SIZES, and the minimum proximity MIN_PROX of
// two holes to each other or the edge of a contig.

#include "Fastavector.h"
#include "MainTools.h"
#include "ParseSet.h"
#include "Superb.h"
#include "efasta/EfastaTools.h"
#include "math/HoInterval.h"
#include "paths/AssemblyCleanupTools.h"
#include "random/Random.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(SUBDIR, "test");
     CommandArgument_String(SCAFFOLDS_IN);
     CommandArgument_String_OrDefault(SCAFFOLDS_OUT, SCAFFOLDS_IN + ".swiss_cheese");
     CommandArgument_Int_OrDefault(NHOLES, 1000);
     CommandArgument_String_OrDefault(HOLE_SIZES, "{5000,6000,7000,8000}");
     CommandArgument_Int_OrDefault(MIN_PROX, 20000);
     CommandArgument_Bool_OrDefault(WRITE, True);
     EndCommandArguments;

     // Define directories, etc.

     String data_dir = PRE + "/" + DATA; 
     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;

     // Parse arguments.

     vec<int> hole_sizes;
     ParseIntSet( HOLE_SIZES, hole_sizes );

     // Load assembly.

     cout << Date( ) << ": loading data" << endl;
     String supers_file = sub_dir + "/" + SCAFFOLDS_IN + ".superb";
     String efasta_file = sub_dir + "/" + SCAFFOLDS_IN + ".contigs.efasta";
     vec<superb> scaffolds;
     ReadSuperbs( supers_file, scaffolds );
     VecEFasta tigse;
     LoadEfastaIntoStrings( efasta_file, tigse );
     int ntigs = tigse.size( );

     // Set up for random location selection.

     vec<int64_t> starts( ntigs + 1 );
     starts[0] = 0;
     for ( int t = 0; t < ntigs; t++ )
          starts[t+1] = starts[t] + tigse[t].size( );
     int64_t N = starts.back( );

     // Define holes.

     cout << Date( ) << ": looking for hole locations" << endl;
     vec< vec<ho_interval> > holes(ntigs);
     int h;
     for ( h = 0; h < NHOLES; h++ )
     {    const int max_tries = 10000;
          Bool found = False;
          for ( int tr = 0; tr < max_tries; tr++ )
          {    int64_t start = big_random( ) % N;
               int t = -1;
               for ( t = 0; t < ntigs; t++ )
                    if ( starts[t+1] > start ) break;
               int pos = starts[t+1] - start;
               int hole_size = hole_sizes[ randomx( ) % hole_sizes.size( ) ];
               int hstart = pos, hstop = pos + hole_size;
               if ( hstart < MIN_PROX ) continue;
               if ( (int) tigse[t].size( ) - hstop < MIN_PROX ) continue;
               ho_interval h( hstart, hstop );
               Bool too_close = False;
               for ( int j = 0; j < holes[t].isize( ); j++ )
               {    if ( Distance( h, holes[t][j] ) < MIN_PROX )
                    {    too_close = True;
                         break;    }    }
               if (too_close) continue;
               int lbracks = 0, rbracks = 0;
               for ( int j = 0; j < hstart; j++ )
               {    if ( tigse[t][j] == '{' ) lbracks++;
                    if ( tigse[t][j] == '}' ) lbracks--;    }
               for ( int j = hstart; j < (int ) tigse[t].size( ); j++ )
               {    if ( tigse[t][j] == '{' ) rbracks++;
                    if ( tigse[t][j] == '}' ) rbracks--;    }
               if ( lbracks != 0 || rbracks != 0 ) continue;
               holes[t].push_back(h);
               found = True;
               break;    }
          if ( !found ) break;    }
     cout << Date( ) << ": Identified " << h << " hole locations." << endl;

     // Go through the scaffolds and edit them.  We set the deviation of the 
     // new gaps to 10% of the hole size.

     cout << Date( ) << ": modifying scaffolds" << endl;
     VecEFasta tigse2;
     for ( int s = 0; s < scaffolds.isize( ); s++ )
     {    superb& S = scaffolds[s];
          for ( int j = S.Ntigs( ) - 1; j >= 0; j-- )
          {    int t = S.Tig(j);
               Sort( holes[t] );
               const vec<ho_interval>& H = holes[t];
               VecEFasta T2( H.size( ) + 1 );
               if ( H.empty( ) ) T2[0] = tigse[t];
               else
               {    T2.front( ) = tigse[t].substr( 0, H[0].Start( ) );
                    for ( int l = 0; l < H.isize( ) - 1; l++ )
                    {    T2[l+1] = tigse[t].substr( H[l].Stop( ),
                              H[l+1].Start( ) - H[l].Stop( ) );    }
                    T2.back( ) = tigse[t].substr( H.back( ).Stop( ),
                         tigse[t].isize( ) - H.back( ).Stop( ) );    }
               superb X;
               X.SetNtigs( T2.size( ) );
               for ( int l = 0; l < X.Ntigs( ); l++ )
               {    X.SetTig( l, tigse2.size( ) + l );
                    X.SetLen( l, T2[l].Length1( ) );
                    if ( l < X.Ngaps( ) )
                    {    X.SetGap( l, H[l].Length( ) );
                         X.SetDev( l, H[l].Length( ) / 10 );    }     }
               S.ReplaceTigBySuper( j, X );
               tigse2.append(T2.begin(),T2.end());    }    }

     // Write assembly and edits.

     if (WRITE)
     {    cout << Date( ) << ": writing assembly and edits" << endl;
          Assembly A( scaffolds, tigse2 );
          A.WriteAll( sub_dir + "/" + SCAFFOLDS_OUT );    }    }
