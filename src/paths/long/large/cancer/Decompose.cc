///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Try to decompose tumor-only edges.  For the moment filtering on long edges.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "FastIfstream.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "kmers/KmerRecord.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/GapToyTools.h"
#include "paths/long/large/Lines.h"
#include "paths/long/large/cancer/SomaticNormal.h"

int main( )
{    RunTime( );

     // Define hardcoded paths.

     // String work_dir = "/wga/scr4/jaffe/GapToy/50935.HCC1954+BL";
     String work_dir = "/wga/scr4/jaffe/GapToy/50921.HCC1143+BL";
     String DIR = work_dir + "/a.final";
     String fin_dir = DIR;
     Bool ALTREF = True;
     String suffix = ( ALTREF ? "_alt" : "" );

     // Load assembly.

     HyperBasevectorX hb;
     BinaryReader::readFile( DIR + "/a.hbx", &hb );
     vec<int> inv;
     BinaryReader::readFile( DIR + "/a.inv", &inv );
     vec<vec<int>> count;
     BinaryReader::readFile( DIR + "/a.countsb", &count );
     vec< vec< pair<int,int> > > hits;
     BinaryReader::readFile( fin_dir + "/a.aligns" + suffix, &hits );
     vec<vec<vec<vec<int>>>> lines;
     BinaryReader::readFile( DIR + "/a.lines", &lines );
     vec<int> tol;
     GetTol( hb, lines, tol );

     // Load genome record names.

     String gnf = DIR + "/../genome.names" + suffix;
     fast_ifstream in(gnf);
     String line;
     vec<String> genome_names;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          genome_names.push_back(line);    }

     // Find hits in bad regions.

     vec<Bool> bad( hb.E( ), False );
     FindBadEdges( hb, inv, lines, hits, bad, ALTREF );

     // Setup.
     
     const int L = 60;
     const int K = hb.K( );
     vec< triple<kmer<L>,int,int> > kmers_plus;
     vecbasevector edges( hb.E( ) );

     // Define edge trims.

     vec<int> ltrim( hb.E( ), 0 ), rtrim( hb.E( ), 0 );
     for ( int e = 0; e < hb.E( ); e++ )
     {    int v = hb.ToLeft(e), w = hb.ToRight(e);
          if ( hb.From(w).size( ) == 1 ) rtrim[e] = K-L;
          if ( hb.To(v).size( ) == 1 ) 
          {    if ( K-L + rtrim[e] <= hb.Bases(e) ) ltrim[e] = K-L;    }    }

     // Create lookup table.

     cout << Date( ) << ": making lookup table" << endl;
     vec<int64_t> starts;
     starts.push_back(0);
     for ( int e = 0; e < hb.E( ); e++ )
     {    const basevector& u = hb.EdgeObject(e);
          if ( count[1][e] == 0 ) 
          {    starts.push_back( starts.back( ) );
               continue;    }
          int nu = u.isize( ) - ltrim[e] - rtrim[e];
          starts.push_back( starts.back( ) + Max( 0, nu - L + 1 ) );    }
     kmers_plus.resize( starts.back( ) );
     #pragma omp parallel for
     for ( int e = 0; e < hb.E( ); e++ )
     {    if ( count[1][e] == 0 ) continue;
          const basevector& u = hb.EdgeObject(e);
          kmer<L> x;
          for ( int j = ltrim[e]; j <= u.isize( ) - L - rtrim[e]; j++ )
          {    int64_t r = starts[e] + j - ltrim[e];
               x.SetToSubOf( u, j ); 
               kmers_plus[r].first = x;
               kmers_plus[r].second = e; 
               kmers_plus[r].third = j;    }    }
     ParallelSort(kmers_plus);

     // Go through tumor-only edges.

     int ir = 0;
     for ( int e = 0; e < hb.E( ); e++ )
     {    
          // Temporary filter.

          if ( count[0][e] < 10 ) continue;
          if ( count[1][e] > 0 ) continue;
          if ( hb.Kmers(e) < 500 ) continue;
          if ( inv[e] < e ) continue;
          if ( bad[e] ) continue;

          // Find matches.

          int v = hb.ToLeft(e), w = hb.ToRight(e);
          vec< triple<int,int,int> > P;
          const basevector& E = hb.EdgeObject(e);
          Bool fail = False;
          for ( int j = 0; j <= E.isize( ) - L; j++ )
          {    kmer<L> x;
               x.SetToSubOf( E, j );
               int64_t low = LowerBound1( kmers_plus, x );
               int64_t high = UpperBound1( kmers_plus, x );

               vec<int64_t> candidates;
               Bool local = False;
               for ( int64_t m = low; m < high; m++ )
               {    candidates.push_back(m);
                    int f = kmers_plus[m].second;
                    if ( hb.ToLeft(f) == v || hb.ToRight(f) == v ) local = True;
                    if ( hb.ToLeft(f) == w || hb.ToRight(f) == w ) local = True;    }
               if (local)
               {    vec<Bool> to_delete( candidates.size( ), False );
                    for ( int64_t m = low; m < high; m++ )
                    {    Bool locm = False;
                         int f = kmers_plus[m].second;
                         if ( hb.ToLeft(f) == v || hb.ToRight(f) == v ) locm = True;
                         if ( hb.ToLeft(f) == w || hb.ToRight(f) == w ) locm = True;
                         if ( !locm ) to_delete[m-low] = True;    }
                    EraseIf( candidates, to_delete );     }

               if ( candidates.solo( ) )
               {    int64_t m = candidates[0];
                    P.push( j, kmers_plus[m].second, 
                         kmers_plus[m].third );    }    }

          // Collate matches.

          vec< quad<int,int,int,int> > Q; // (tstart, len, eid, estart)
          for ( int i = 0; i < P.isize( ); i++ )
          {    int l;
               for ( l = i + 1; l < P.isize( ); l++ )
               {    if ( P[l].second != P[i].second ) break;
                    int tshift = P[l].first - P[l-1].first;
                    int eshift = P[l].third - P[l-1].third;
                    if ( tshift != eshift || tshift > L ) break;    }
               Q.push( P[i].first, P[l-1].first - P[i].first + L,
                    P[i].second, P[i].third );
               i = l - 1;    }
          Bool cov = True;
          if ( Q.empty( ) || Q.front( ).first > 0 ) cov = False;
          else if ( Q.back( ).first + Q.back( ).second < E.isize( ) ) cov = False;
          if (cov)
          {    cout << "\n" << "[" << ++ir << "]\n" << e << " " << inv[e] << ":\n";
               for ( int i = 0; i < Q.isize( ); i++ )
               {    if ( i > 0 )
                    {    int start1 = Q[i-1].first, start2 = Q[i].first; 
                         int len1 = Q[i-1].second, len2 = Q[i].second;
                         int stop1 = Q[i-1].first + len1, stop2 = Q[i].first + len2;
                         int e1 = Q[i-1].third, e2 = Q[i].third;
                         int gap = start2 - stop1;
                         int estart1 = Q[i-1].fourth, estop1 = Q[i-1].fourth + len1;
                         int estart2 = Q[i].fourth, estop2 = Q[i].fourth + len2;

                         if ( gap == -L+1 && hb.ToRight(e1) == hb.ToLeft(e2)
                              && ( hb.Bases(e1) - estop1 ) + estart2 == K - L );
                         else if ( gap <= 0 ) cout << "(" << gap << ")\n";    
                         else
                         {    cout << "(";
                              for ( int z = stop1; z < start2; z++ )
                                   cout << as_base( E[z] );
                              cout << ")\n";    }    }
                    cout << Q[i].third << "." << Q[i].fourth << "-"
                         << Q[i].fourth + Q[i].second << " [" << Q[i].second << "]";
                    cout << PrintHits( Q[i].third, hits, hb, inv, genome_names );
                    cout << "\n";    }    }    }    }
