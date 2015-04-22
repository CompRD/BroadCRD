///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MaizeDemo.  Given the maize genome reference sequence, suppose that we can
// close gaps of bounded size between unipaths; compute the resulting N50.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/KmerPath.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/LongReadsToPaths.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_Int(TOP); // 1000000000 = all
     EndCommandArguments;

     double all_clock = WallClockTime( );

     // Load maize genome reference sequence.

     cout << Date( ) << ": loading genome" << endl;
     vecbasevector genome( "/wga/scr4/bigrefs/maize/genome.fastb" );
     PRINT( genome.size( ) );

     // Shrink genome for rapid testing.

     PRINT(TOP);
     for ( int g = 0; g < (int) genome.size( ); g++ )
          if ( genome[g].isize( ) > TOP ) PRINT( genome[g].size( ) );
     for ( int g = 0; g < (int) genome.size( ); g++ )
          if ( genome[g].isize( ) > TOP ) genome[g].resize(TOP);

     // Build K=200 hyperbasevector from the genome.  We chop the genome up into
     // chunks because it makes the code a lot faster.

     double clock = WallClockTime( );
     const int K = 200;
     HyperBasevector hb;
     {    vecbasevector genomex;
          const int chunk = 10000;
          for ( int g = 0; g < (int) genome.size( ); g++ )
          {    const basevector& G = genome[g];
               for ( int p = 0; p < G.isize( ); p += chunk - K )
               {    basevector b( G, p, Min( chunk, G.isize( ) - p ) );
                    genomex.push_back(b);    }    }
          HyperKmerPath h;
          vecKmerPath paths, paths_rc;
          cout << Date() << ": calling LongReadsToPaths" << endl;
          unsigned const COVERAGE = 60u;
          LongReadsToPaths(genomex,K,COVERAGE,False,False,&hb,&h,&paths,&paths_rc);
          cout << Date( ) << ": done" << endl;
          cout << TimeSince(clock) << " used" << endl;    }

     // Map the hyperbasevector edges back to the genome.  Presumably a
     // redundant calculation.  We only consider edges of size >= 1000 that have
     // a unique place.

     const int min_len = 1000;
     vecbasevector all(genome), genome_rc(genome);
     for ( int g = 0; g < (int) genome.size( ); g++ )
          genome_rc[g].ReverseComplement( );
     all.Append(genome_rc);
     vec<int> tids;
     for ( int e = 0; e < hb.E( ); e++ )
     {    if ( hb.O(e).isize( ) >= min_len ) 
          {    basevector b( hb.O(e), 0, K );
               all.push_back(b);    
               tids.push_back(e);    }    }
     cout << Date( ) << ": making lookup" << endl;
     vec< triple<kmer<K>,int,int> > kmers_plus;
     MakeKmerLookup0( all, kmers_plus );
     cout << Date( ) << ": building places" << endl;
     vec< triple<int,int,int> > places;
     for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
     {    int64_t j;
          for ( j = i + 1; j < kmers_plus.jsize( ); j++ )
               if ( kmers_plus[j].first != kmers_plus[i].first ) break;
          if ( j - i == 2 && kmers_plus[i].second < (int) genome.size( )
               && kmers_plus[i+1].second >= 2 * (int) genome.size( ) )
          {    int gid = kmers_plus[i].second;
               int gpos = kmers_plus[i].third;
               int tid = tids[ kmers_plus[i+1].second - 2 * (int) genome.size( ) ];
               places.push( gid, gpos, tid );    }
          i = j - 1;    }
     cout << Date( ) << ": sorting" << endl;
     Sort(places);
     PRINT( places.size( ) );
     for ( int i = 0; i < places.isize( ); i++ )
     {    int gid = places[i].first, gpos = places[i].second;
          int tid = places[i].third;
          // cout << tid << "[l=" << hb.Bases(tid) << "] at "
          //      << gid << "." << gpos << endl;    
               }

     // Compute N50 edge size, and the N50s that result if gaps of various
     // sizes are closed.

     vec<int> len1;
     for ( int i = 0; i < places.isize( ); i++ )
     {    int tid = places[i].third;
          len1.push_back( hb.Bases(tid) );    }
     Sort(len1);
     cout << "N50.1 = " << N50(len1) << "\n" << endl;
     cout << "\n";
     vec<int> bounds = {2000,2500,3000,3500,4000,4500,5000};
     for ( int b : bounds )
     {    vec<int> len2;
          for ( int i = 0; i < places.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < places.isize( ); j++ )
               {    if ( places[j].first != places[i].first ) break;
                    int sep = places[j].second 
                         - places[j-1].second - hb.Bases( places[j-1].third );
                    if ( sep > b ) break;    }
               int len = places[j-1].second + hb.Bases( places[j-1].third )
                    - places[i].second;
               len2.push_back(len);
               i = j - 1;    }
          Sort(len2);
          cout << "N50(sep bound " << b << ") = " << N50(len2) << endl;    }
     cout << "\n";

     // Compute distribution of separation.

     vec<int> seps;
     for ( int i = 0; i < places.isize( ) - 1; i++ )
     {    if ( places[i].first == places[i+1].first )
          {    seps.push_back( places[i+1].second 
                    - places[i].second - hb.Bases( places[i].third ) );    }    }
     Sort(seps);
     // for ( int i = 0; i < seps.isize( ); i++ )
     //      PRINT2( i, seps[i] );
     int bin1 = 0, bin2 = 0, bin3 = 0;
     for ( int i = 0; i < seps.isize( ); i++ )
     {    if ( seps[i] <= 2000 ) bin1++;
          else if ( seps[i] <= 5000 ) bin2++;
          else bin3++;    }
     cout << PERCENT_RATIO( 3, bin1, seps.isize( ) )
          << " of gaps are <= 2 kb" << endl;
     cout << PERCENT_RATIO( 3, bin2, seps.isize( ) )
          << " of gaps are between 2 kb and 5 kb" << endl;
     cout << PERCENT_RATIO( 3, bin3, seps.isize( ) )
          << " of gaps are > 5 kb" << endl;

     // Write the gap sequences.

     cout << "\nGAP SEQUENCES\n";
     int count = 0;
     vec<int> N;
     #pragma omp parallel for
     for ( int i = 1; i < places.isize( ); i++ )
     {    int g = places[i].first;
          if ( g != places[i-1].first ) continue;
          int start = places[i-1].second + hb.Bases( places[i-1].third );
          int stop = places[i].second;
          int gap = stop - start;
          if ( gap <= 0 || gap > 5000 ) continue;
          basevector b( genome[g], start, gap );
          const int K = 100;
          HyperBasevector hb;
          vecbasevector genomex;
          genomex.push_back(b);
          HyperKmerPath h;
          vecKmerPath paths, paths_rc;
          unsigned const COVERAGE = 60u;
          LongReadsToPaths(genomex,K,COVERAGE,False,False,&hb,&h,&paths,&paths_rc);
          #pragma omp critical
          {    count++;
               cout << "\n";
               PRINT3( count, b.size( ), h.E( ) );
               hb.PrintSummaryDOT0w( cout, False, False, True );
               N.push_back( h.E( ) / 2 );    }    }
     Sort(N);
     cout << "\n";
     PRINT( N.size( ) );
     for ( int i = 0; i < N.isize( ); i++ )
     {    int j = N.NextDiff(i);
          cout << N[i] << " [" << j-i << "]\n";
          i = j - 1;    }

     // Done.

     cout << Date( ) << ": done" << endl;
     cout << "peak mem usage = " << PeakMemUsageGB( ) << " GB" << endl;
     cout << "total time used = " << TimeSince(all_clock) << endl;    }
