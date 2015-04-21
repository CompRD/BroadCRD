///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Jumpomatic: call Jumpster.

#include "Basevector.h"
#include "MainTools.h"
#include "Qualvector.h"
#include "paths/GetNexts.h"
#include "paths/Jumpster.h"
#include "paths/KmerBaseBroker.h"
#include "paths/LongReadTools.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/UnibaseUtils.h"
#include "paths/Unipath.h"
#include "paths/long/CreateGenome.h"
#include <omp.h>
// MakeDepend: library OMP

Bool cmp213( const triple< basevector, basevector, vec<int> >& T1,
     const triple< basevector, basevector, vec<int> >& T2 )
{    if ( T1.second < T2.second ) return True;
     if ( T1.second > T2.second ) return False;
     if ( T1.first < T2.first ) return True;
     if ( T1.first > T2.first ) return False;
     return T1.third < T2.third;    }

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String(HEAD);
     CommandArgument_Int(K);
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0,
       "Number of threads to use (use all available processors if set to 0)");
     CommandArgument_Bool_OrDefault(WRITE, True);
     CommandArgument_Bool_OrDefault(VERBOSE, False);
     CommandArgument_Bool_OrDefault(EVALUATE, False);
     CommandArgument_Bool_OrDefault(PRINT_GRAPH, False);
     EndCommandArguments;

     // Thread control.
     
     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads( NUM_THREADS );
     
     // Define directories.
     
     String data_dir = PRE + "/" + DATA;
     String run_dir = PRE + "/" + DATA + "/" + RUN;

     // Load data.

     String KS = ToString(K);
     vecbasevector unibases( run_dir + "/" + HEAD + ".unibases.k" + KS );
     vecbasevector jumps( run_dir + "/jump_reads_filt.fastb" );
     vecqualvector quals( run_dir + "/jump_reads_filt.qualb" );

     // Normalize T runs.  We look for 10 or more Ts, flanked by at least 8 bases
     // on both ends, and we assume that for given flanks, at least 5 such instances
     // occur.  If we see flank1...flank2 and flank1...flank2' 
     // with flank1 != flank2' then we dump both, and so the same thing for the
     // the flank order swapped.  Sequences that look too much like homopolymers
     // on the ends (governed by max_diff) are excluded.  For fixed flanks if
     // we see too much length varation of the T runs (governed by max_ratio) we
     // exclude the case.  Otherwise we replace the T run by the median length.

     const int mint = 10;
     const int flank = 8;
     const int min_count = 5;
     const int min_diff = 2;
     const double max_ratio = 1.5;
     vec< triple<basevector,basevector,int> > truns;
     for ( size_t id = 0; id < jumps.size( ); id++ )
     {    basevector x = jumps[id];
          for ( int pass = 1; pass <= 2; pass++ )
          {    if ( pass == 2 ) x.ReverseComplement( );
               for ( int j = 0; j < x.isize( ); j++ )
               {    if ( as_base( x[j] ) != 'T' ) continue;
                    int k = j + 1;
                    while(1)
                    {    if ( k == x.isize( ) ) break;
                         if ( as_base( x[k] ) != 'T' ) break;
                         k++;    }
                    int nt = k - j;
                    if ( nt >= mint && j >= flank && x.isize( ) - k >= flank )
                    {    basevector left, right;
                         left.SetToSubOf( x, j - flank, flank );
                         right.SetToSubOf( x, k, flank );
                         truns.push( left, right, nt );    }
                    j = k - 1;    }    }    }
     Sort(truns);
     vec< triple< basevector, basevector, vec<int> > > T;
     for ( int j = 0; j < truns.isize( ); j++ )
     {    int k;
          for ( k = j + 1; k < truns.isize( ); k++ )
          {    if ( truns[k].first != truns[j].first ) break;
               if ( truns[k].second != truns[j].second ) break;    }
          if ( k - j >= min_count )
          {    vec<int> x;
               for ( int l = j; l < k; l++ )
                    x.push_back( truns[l].third );
               T.push( truns[j].first, truns[j].second, x );    }
          j = k - 1;    }
     vec<Bool> to_delete1( T.size( ), False ), to_delete2( T.size( ), False );
     for ( int j = 0; j < T.isize( ); j++ )
     {    if ( j > 0 && T[j].first == T[j-1].first ) to_delete1[j] = True;
          if ( j < T.isize( ) - 1 && T[j].first == T[j+1].first )
               to_delete1[j] = True;    }
     EraseIf( T, to_delete1 );
     sort( T.begin( ), T.end( ), cmp213 );
     for ( int j = 0; j < T.isize( ); j++ )
     {    if ( j > 0 && T[j].second == T[j-1].second ) to_delete2[j] = True;
          if ( j < T.isize( ) - 1 && T[j].second == T[j+1].second )
               to_delete2[j] = True;    }
     EraseIf( T, to_delete2 );
     vec< triple< basevector, basevector, int > > TX;
     for ( int j = 0; j < T.isize( ); j++ )
     {    vec<int> m1( 4, 0 ), m2( 4, 0 );
          for ( int l = 0; l < T[j].first.isize( ); l++ )
               m1[ T[j].first[l] ]++;
          for ( int l = 0; l < T[j].second.isize( ); l++ )
               m2[ T[j].second[l] ]++;
          Sort(m1), Sort(m2);
          int diff1 = m1[0] + m1[1] + m1[2], diff2 = m2[0] + m2[1] + m2[2];
          if ( T[j].third.isize( ) >= min_count 
               && diff1 >= min_diff && diff2 >= min_diff )
          {    vec<int> len, c;
               for ( int l = 0; l < T[j].third.isize( ); l++ )
               {    int m;
                    for ( m = l + 1; m < T[j].third.isize( ); m++ )
                         if ( T[j].third[m] != T[j].third[l] ) break;
                    len.push_back( T[j].third[l] );
                    c.push_back( m - l );
                    l = m - 1;    }
               if ( len.size( ) > 1 && len.back( ) <= max_ratio * len.front( ) )
               {    int m = Median( T[j].third );
                    if (VERBOSE)
                    {    cout << T[j].first << " (";
                         for ( int l = 0; l < len.isize( ); l++ )
                              cout << " " << len[l] << "[" << c[l] << "]";
                         cout << " --> " << m;
                         cout << " ) " << T[j].second << "\n";    }
                    TX.push( T[j].first, T[j].second, m );    }    }    }
     for ( size_t id = 0; id < jumps.size( ); id++ )
     {    basevector& x = jumps[id];
          qualvector& q = quals[id];
          int reverses = 0;
          rerun:
          if ( reverses % 2 == 1 ) 
          {    x.ReverseComplement( );
               q.ReverseMe( );
               reverses++;    }
          for ( int pass = 1; pass <= 2; pass++ )
          {    x.ReverseComplement( );
               q.ReverseMe( );
               reverses++;
               for ( int j = 0; j < x.isize( ); j++ )
               {    if ( as_base( x[j] ) != 'T' ) continue;
                    int k = j + 1;
                    while(1)
                    {    if ( k == x.isize( ) ) break;
                         if ( as_base( x[k] ) != 'T' ) break;
                         k++;    }
                    int nt = k - j;
                    if ( nt >= mint && j >= flank && x.isize( ) - k >= flank )
                    {    basevector left, right;
                         left.SetToSubOf( x, j - flank, flank );
                         right.SetToSubOf( x, k, flank );
                         for ( int r = 0; r < TX.isize( ); r++ )
                         {    if ( TX[r].first == left && TX[r].second == right
                                   && TX[r].third != nt )
                              {    basevector x1, x2, x3;
                                   qualvector q1, q2, q3;
                                   x1.SetToSubOf( x, 0, j );
                                   q1.SetToSubOf( q, 0, j );
                                   x2.resize( TX[r].third );
                                   q2.resize( TX[r].third );
                                   for ( int l = 0; l < x2.isize( ); l++ )
                                   {    x2.Set( l, as_char( 'T' ) );
                                        q2[l] = 0;    }
                                   x3.SetToSubOf( x, k, x.isize( ) - k );
                                   q3.SetToSubOf( q, k, (int) q.size( ) - k );
                                   x = Cat( x1, x2 );
                                   q = Cat( q1, q2 );
                                   x = Cat( x, x3 );
                                   q = Cat( q, q3 );
                                   goto rerun;    }    }    }    
                    j = k - 1;   }    }    }

     // Get nexts.

     vec< vec<int> > nexts;
     GetNexts( K, unibases, nexts );
     if (PRINT_GRAPH)
     {    vec<int> to_rc;
          UnibaseInvolution( unibases, to_rc );
          cout << "\nunibase graph:\n\n";
          for ( size_t u = 0; u < unibases.size( ); u++ )
          {    cout << u << "[l=" << unibases[u].isize( ) - K + 1
                    << ",rc=" << to_rc[u] << "] -->";
               for ( int j = 0; j < nexts[u].isize( ); j++ )
                    cout << " " << nexts[u][j];
               cout << "\n";    }
          cout << "\n";    }

     // Load and hash genome.

     vecbasevector genome, genome2;
     const int LG = 12;
     VecIntPairVec Glocs;
     int ng;
     if (EVALUATE)
     {    genome.ReadAll( data_dir + "/genome.fastb" );
          ng = genome.size( );
          genome2.resize( genome.size( ) ); 
          for ( size_t j = 0; j < genome.size( ); j++ )
               genome2[j] = Cat( genome[j], genome[j] );
          CreateGlocs(  genome2, LG, Glocs );   }

     // Bridge gaps.

     vec<basevector> bridges;
     Jumpster( K, unibases, nexts, jumps, quals, bridges, VERBOSE, 
          genome2, LG, Glocs );

     // Integrate bridges and write output.

     if (WRITE) 
     {    cout << Date( ) << ": building new unipaths" << endl;
          vecbasevector all(unibases);
          for ( int j = 0; j < bridges.isize( ); j++ )
          {    all.push_back_reserve( bridges[j] );
               bridges[j].ReverseComplement( );
               all.push_back_reserve( bridges[j] );    }
          for ( int id1 = 0; id1 < (int) unibases.size( ); id1++ ) 
          {    for (int j = 0; j < nexts[id1].isize( ); j++) 
               {    int id2 = nexts[id1][j];
	            basevector b = unibases[id1];
	            b.resize( b.size( ) + 1 );
	            b.Set( b.size( ) - 1, unibases[id2][K-1] );
	            all.push_back_reserve(b);    }    }
          vecKmerPath newpaths, newpathsrc, newunipaths;
          vec<tagged_rpint> newpathsdb, newunipathsdb;
          Mkdir777( run_dir + "/Jumpomatic.tmp" );
          ReadsToPathsCoreY( all, K, newpaths, newpathsrc, newpathsdb,
               run_dir + "/Jumpomatic.tmp", NUM_THREADS );
          Unipath( newpaths, newpathsrc, newpathsdb, newunipaths, newunipathsdb );
          KmerBaseBroker newkbb( K, newpaths, newpathsrc, newpathsdb, all );
          vecbasevector newunibases;
          for (size_t i = 0; i < newunipaths.size(); i++)
               newunibases.push_back_reserve( newkbb.Seq(newunipaths[i]) );
          cout << Date( ) << ": writing output files" << endl;
          String outhead = run_dir + "/" + HEAD + ".jumpomatic";
          newpaths.WriteAll( outhead + ".paths.k" + KS );
          newpathsrc.WriteAll( outhead + ".paths_rc.k" + KS );
          newunipaths.WriteAll( outhead + ".unipaths.k" + KS );
          BinaryWriter::writeFile( outhead + ".pathsdb.k" + KS, newpathsdb );
          BinaryWriter::writeFile( outhead + ".unipathsdb.k" + KS, newunipathsdb );
          newunibases.WriteAll( outhead + ".unibases.k" + KS );    }
     cout << Date( ) << ": done" << endl;    }
