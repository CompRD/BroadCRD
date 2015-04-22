///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Trying to match ONT basecalling.
//
// Need to run this from a rebuilt machine.
//
// Current performance (fraction of bases in 10-mers matching reference):
// ONT   - 16.1%
// Broad - 13.6%

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "FastIfstream.h"
#include "FetchReads.h"
#include "lookup/LookAlign.h"
#include "math/HoInterval.h"
#include "MainTools.h"
#include "TokenizeString.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "math/Functions.h"

// /wga/scr4/vendor/on/lambda-2014-06-11
// ran convert2_gen

double Align( 
     const vec<pair<double,double>>& R, // read (events) -- (mean,dev)
     const vec<pair<double,double>>& A, // kmers (model) -- (mean,dev)
     ostringstream& out,                // for logging
     vec<int>& kmers,                   // kmer sequence in best path
     const vec<vec<int>>& next,         // next[j][l] (l=0..3) denotes next kmers
     vec<int>& p,                       // best path
     vec<double>& P                     // PARAMETERS
     )
{
     int nr = R.size( ), na = A.size( );
     int N = 2 + nr*na; // number of vertices
     PRINT3_TO( out, nr, na, N );

     // VERTEX                                     INDEX
     // 1. begin, end                              0, 1
     // 2. r_i x a_j                               2 + (i * na) + j
     //
     // EDGE                                       WEIGHT
     // 1. begin --> r_0 x a_j                     0
     // 2. r_nr-1 x a_j --> end                    delta( nr-1, j )
     // 3. r_i x a_j --> r_i+1 x a_j+              delta( i, j ),  j+ in next[j][l]
     // 4. r_i x a_j --> r_i x a_j+                1.5 + delta( i, j+ )
     // 5. r_i x a_j --> r_i+1 x a_j               3.0 + delta( i, j )

     // experiment 1
     // varied addend of edge type 4
     // 1.0   15.2%
     // 1.5   16.7% (best, using)
     // 2.0   16.1%
     // 2.5   13.6%

     // experiment 2
     // varied addend of edge type 5
     // 1.0   16.7%
     // 1.5   17.0% (best)
     // 2.0   16.9%
     // 2.5   16.8%
     // 3.0   16.7% (using)
     // 3.5   16.7%
     // 4.0   16.6%

     auto pos = [=]( int i, int j ){ return 2 + (i * na) + j; };
     auto r_ind = [=]( int pos ){ return (pos - 2) / na; };
     auto a_ind = [=]( int pos ){ return (pos - 2) % na; };

     auto delta = [&]( int i, int j )
     {    return ( P[0] * R[i].first + P[1] - (P[2]*i)/10000 - A[j].first ) 
               * ( P[0] * R[i].first + P[1] - (P[2]*i)/10000 - A[j].first )
               / ( R[i].second * R[i].second + A[j].second * A[j].second );    };

     auto deltaR = [&]( int i, int j )
     {    return ( R[i].first - R[j].first ) * ( R[i].first - R[j].first )
               / ( R[i].second * R[i].second + R[j].second * R[j].second );    };

     auto vert_name = [=]( int v )
     {    if ( v == 0 ) return String("begin");
          if ( v == 1 ) return String("end");
          int i = r_ind(v);
          int j = v - 2 - (i*na);
          return "r" + ToString(i) + ".a" + ToString(j);    };

     // Define G.

     auto f = [&]( const int v )
     {    vec< pair<int,double> > x;
          if ( v == 0 )                             // edge type 1
          {    for ( int j = 0; j < na; j++ )
                    x.push( pos(0, j), 0 );    }
          else if ( v >= 2 )
          {    int i = r_ind(v), j = a_ind(v);

               if ( i == nr-1 )                     // edge type 2
                    x.push( 1, delta( nr-1, j ) );

               if ( i < nr - 1 )                    // edge type 3
               {    for ( int l = 0; l < 4; l++ )
                         x.push( pos(i+1, next[j][l]), delta( i, j ) );    }

               for ( int l = 0; l < 4; l++ )        // edge type 4
                    x.push( pos(i, next[j][l]), P[3] + delta( i, next[j][l] ) );


               // OLD
               // if ( i < nr - 1 )                    // edge type 5
               //      x.push( pos(i+1, j), P[4] + P[5] * delta( i, j ) );    }

               if ( i < nr - 1 )                    // edge type 5
                    x.push( pos(i+1, j), P[4] + P[5] * deltaR( i, i+1 ) );    }


               /*
               for ( int l = 0; l < 4; l++ )        // edge type 4
                    x.push( pos(i, next[j][l]), 
                         1.4 + 0.1 * delta( i, j ) + 1.1 * delta( i, next[j][l] ) );
               if ( i < nr - 1 )                    // edge type 5
                    x.push( pos(i+1, j), 2.8 + 0.5 * delta( i, j ) );    }    
               */

          return x;    };

     digraphE_V1<double> G( N, f );

     // Find shortest path from begin to end.

     // vec<int> p;
     double d = G.ShortestPath( 0, 1, p );

     // Generate kmer sequence.

     kmers.resize( p.isize( ) - 2 );
     for ( int i = 1; i < p.isize( ) - 1; i++ )
          kmers[i-1] = a_ind( p[i] );

     // Print result.

     // simple printer
     for ( int j = 0; j < p.isize( ); j++ )
     {    if ( j > 0 ) out << "--> ";
          out << vert_name( p[j] ) << "\n";    }

     /*
     for ( int j = 0; j < p.isize( ); j++ )
     {    if ( j == 0 ) out << "begin\n";
          else
          {    out << "--> ";
               int k;
               for ( k = j + 1; k < p.isize( ); k++ )
                    if ( p[k] - p[k-1] != na + 1 ) break;
               if ( k-1 > j )
               {    out << vert_name( p[j] ) << " ... " 
                         << vert_name( p[k-1] ) << "\n";    }
               else out << vert_name( p[j] ) << "\n";
               j = k - 1;    }    }
     */

     out << "d = " << setprecision(5) << d << endl;
     out << "\n";
     return d;    }

double Delta( const double x, const double xd, const double y, const double yd )
{     return (x-y)*(x-y) / (xd*xd + yd*yd);    }

double Square( double x )
{     return x * x;    }

void Stuff( const String& rn )
{
     String rn_loc = "/wga/scr4/vendor/on/2014-07-22/ecoli_r7/convert/" + rn;

     // Get truth data.

     SystemSucceed( "QueryLookupTable K=12 MM=12 MC=0.02 "
          "SEQS=" + rn_loc + ".fast5.fasta "
          "L=/wga/scr4/macro/ecoli/scs.lookup FW_ONLY=True "
          "SMITH_WAT=True PARSEABLE=True BW_ADD=300 > woof.aligns" );
     vec<look_align> aligns;
     LoadLookAligns( "woof.aligns", aligns );
     PRINT( aligns.size( ) );
     if ( !aligns.solo( ) ) return;
     vecbasevector G( "/wga/scr4/macro/ecoli/scs.fastb" );
     vecbasevector R;
     FetchReads( R, 0, rn_loc + ".fast5.fasta" );
     int gid = aligns[0].target_id;
     String GS = G[gid].ToString( );
     align a( aligns[0].a );
     vec<ho_interval> perf1, perf2;
     a.PerfectIntervals1( R[0], G[gid], perf1 );
     a.PerfectIntervals2( R[0], G[gid], perf2 );
     int start = a.pos2( ), stop = a.Pos2( );
     PRINT2( start, stop );

     // Print perfect matches.

     Bool print_perfect_matches = False;
     if (print_perfect_matches)
     {    for ( int i = 0; i < perf1.isize( ); i++ )
          {    if ( perf1[i].Length( ) < 5 ) continue;
               cout << perf1[i].Start( ) << " " << perf2[i].Start( ) << " --> ";
               for ( int j = perf1[i].Start( ); j < perf1[i].Stop( ); j++ )
                    cout << as_base( R[0][j] );
               cout << "\n";    }    }

     // Score the alignment.

     int nperf1 = 0;
     for ( int i = 0; i < perf1.isize( ); i++ )
     {    if ( perf1[i].Length( ) < 10 ) continue;
          nperf1 += perf1[i].Length( );    }
     cout << PERCENT_RATIO( 3, nperf1, R[0].isize( ) )
          << " of bases in matching 10-mers" << endl;

     // Load data.

     cout << Date( ) << ": loading" << endl;
     vec< triple<String,double,double> > model; // (kmer,mean,dev)
     vec<String> kmers;
     fast_ifstream in1( rn_loc + ".fast5.model" );
     String line;
     vec<String> tokens;
     while(1)
     {    getline( in1, line );
          if ( in1.fail( ) ) break;
          Tokenize( line, ',', tokens );
          model.push( tokens[0], tokens[2].Double( ), tokens[3].Double( ) );
          kmers.push( tokens[0] );    }
     kmers.resize(1024), model.resize(1024);
     SortSync( kmers, model );

     vec<vec<int>> ghits(1024);
     for ( int j = start; j <= stop - 5; j++ )
     {    String kmer = GS.substr( j, 5 );
          int p = BinPosition( kmers, kmer );
          ghits[p].push_back( j - start );    }

     vec< triple<String,double,double> > events; // (kmer,mean,dev)
     vec<double> model_level;
     vec<int> move;


     // Use existing ONT base calls.

     fast_ifstream in2( rn_loc + ".fast5.events" );
     int lcount = 0;
     while(1)
     {    getline( in2, line );
          if ( in2.fail( ) ) break;
          lcount++;

          /*
          // FOCUSING ON TINY STRETCH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if ( lcount < 10 ) continue;
          if ( lcount > 25 ) break;
          */

          Tokenize( line, ',', tokens );
          events.push( tokens[4], tokens[0].Double( ), tokens[2].Double( ) );
          // events.push( tokens[8], tokens[0].Double( ), tokens[2].Double( ) );
          model_level.push( tokens[5].Double( ) );
          move.push_back( tokens[6].Int( ) );    }

     // Do ONT base calling, then extract events from the calls.

     /*
     SystemSucceed( "/wga/dev/local/bin/python "
          "/wga/scr4/vendor/on/Dragonet/Dragonet-1.1.2/bin/basecall_1d "
          "/wga/scr4/vendor/on/2014-07-22/ecoli_r7/reads/" + rn + ".fast5 "
          "hhh.fast5" );
     SystemSucceed( "/wga/dev/jaffe/BroadCRD/paths/long/magic/scripts/"
          "dump_hdf_tag.py call_ev hhh.fast5 > hhh.events" );
     // SystemSucceed( "/wga/dev/jaffe/BroadCRD/paths/long/magic/scripts/"
     //      "dump_hdf_tag.py basecall hhh.fast5 > hhh.fasta" );
     fast_ifstream in2( "hhh.events" );
     getline( in2, line );
     while(1)
     {    getline( in2, line );
          if ( in2.fail( ) ) break;
          Tokenize( line, ',', tokens );
          events.push( tokens[4], tokens[2].Double( ), tokens[3].Double( ) );
          // events.push( tokens[8], tokens[0].Double( ), tokens[2].Double( ) );
          model_level.push( tokens[5].Double( ) );
          move.push_back( tokens[6].Int( ) );    }
     */

     events.resize( events.isize( ) - 1 );

     // Substitute better model.

     /*
     fast_ifstream bin( 
          "/wga/scr4/vendor/on/lambda-2014-06-11/R6_vs_R7-template-row.csv" );
     getline( bin, line );
     while(1)
     {    getline( bin, line );
          if ( bin.fail( ) ) break;
          Tokenize( line, ',', tokens );
          String kmer = tokens[0];
          double level_mean = tokens[4].Double( ), level_stdv = tokens[5].Double( );
          int p = BinPosition( kmers, kmer );
          model[p].second = level_mean, model[p].third = level_stdv;    }
     */

     // Substitute inferred model.

     /*
     fast_ifstream woof( "/wga/scr4/vendor/on/lambda-2014-06-11/model.inferred" );
     while(1)
     {    getline( woof, line );
          if ( woof.fail( ) ) break;
          Tokenize( line, ' ', tokens );
          int p = BinPosition( kmers, tokens[0] );
          model[p].second = tokens[1].Double( );
          model[p].third = tokens[2].Double( );    }
     */

     // Define next.

     vec<vec<int>> next( model.size( ) );
     for ( int i = 0; i < model.isize( ); i++ )
     {    for ( int m = 0; m < 4; m++ )
          {    int ip = BinPosition( kmers, 
                    model[i].first.substr( 1, 4 ) + String( as_base(m) ) );
               next[i].push_back(ip);    }    }

     // Normalize by subtracting means.

     vec<double> modelx, eventsx;
     for ( int i = 0; i < model.isize( ); i++ )
          modelx.push_back( model[i].second );
     for ( int i = 0; i < events.isize( ); i++ )
          eventsx.push_back( events[i].second );
     double model_mean = Mean(modelx);
     double event_mean = Mean(eventsx);
     PRINT2( model_mean, event_mean );
     for ( int i = 0; i < model.isize( ); i++ )
          model[i].second -= model_mean;
     for ( int i = 0; i < events.isize( ); i++ )
          events[i].second -= event_mean;

     vec<double> P = { 1, 0, 0.0, 1.5, 1.0, 1.0 };
     double best_score = 0;
     for ( int apass = 0; apass < 10000; apass++ )
     {
          vec<double> Pnew = P;

          if ( apass > 0  )
          {    
               int j1 = random( ) % P.size( );
               PRINT(j1);
               int sg1 = ( 2 * ( random( ) % 2 ) ) - 1;
               double delta = 0.02 * ( 1 + ( random( ) % 6 ) );
               if ( j1 != 1 && Pnew[j1] < delta ) sg1 = 1;
               Pnew[j1] += sg1 * delta;    

               /*
               int j2 = random( ) % P.size( );
               PRINT(j2);
               int sg2 = ( 2 * ( random( ) % 2 ) ) - 1;
               delta = 0.01 * ( 1 + ( random( ) % 3 ) );
               if ( j2 != 1 && Pnew[j2] < delta ) sg2 = 1;
               Pnew[j2] += sg2 * delta;    
               */

               cout << "trying " << printSeq(Pnew) << endl;    }

     // Call Align.

     cout << "\n" << Date( ) << ": calling Align" << endl;
     {
     vec< pair<double,double> > R, A;
     for ( int i = 0; i < events.isize( ); i++ )
          R.push( events[i].second, events[i].third );
     for ( int i = 0; i < model.isize( ); i++ )
          A.push( model[i].second, model[i].third);
     ostringstream out;
     vec<int> kmer_ids;
     vec<int> path;
     Align( R, A, out, kmer_ids, next, path, Pnew );
     // cout << out.str( );
     // for ( int i = 0; i < kmer_ids.isize( ); i++ )
     //      cout << i << "\t" << kmers[ kmer_ids[i] ] << endl;














     // Trying to reproduce old output.  Note duplication with code below.

     /*
     cout << "\n";
     {
     int equals = 0;
     int xcount = 0;
     int movesum = 0;
     for ( int i = 1; i < path.isize( ) - 1; i++ )
     {    int e = ( path[i] - 2 ) / A.size( );
          movesum += move[e];
          String kmer = kmers[ kmer_ids[i-1] ];
          int p = BinPosition( kmers, kmer );
          int q = BinPosition( kmers, events[e].first );
          double obs = events[e].second;
          cout << movesum
               << "\t"
               << obs << "\t";

          {    int v = path[i];
               String vn;
               if ( v == 0 ) vn = "begin";
               else if ( v == 1 ) vn = "end";
               else
               {    int i = (v-2) / A.size( );
                    int j = v - 2 - (i * A.size( ));
                    vn = "r" + ToString(i) + "." + kmer;    }
               cout << vn;    }

          cout << "[" << ToString( 
                    // Square( Abs( obs - model[p].second ) )
                    Delta( events[e].second, events[e].third,
                         model_level[e], model[p].third )
                    , 3 ) 
               << "] "
               << '\t' << events[e].first 
               << "[" << ToString( 
                    // Square( Abs( obs - model[q].second ) )
                    Delta( events[e].second, events[e].third,
                         model_level[e], model[q].third )
                    , 3 ) 
               << "]";

          cout << "\t";
          if ( p == q ) 
          {    cout << "equal";
               equals++;    }

          vec<int> h;
          for ( int j = 0; j < ghits[q].isize( ); j++ )
          {    if ( Abs( ghits[q][j] - movesum ) > 10 + movesum/10 ) continue;
               h.push_back( ghits[q][j] );    }
          cout << "\t" << printSeq(h);
          cout << endl;
          if ( e == 0 || events[e].first != events[e-1].first ) xcount++;    }

     cout << "\nequals = " << equals << endl;
     }
     */










     // Convert to fasta.  Not quite right.

     {    Ofstream( pout, "qqq.fasta" );
          pout << ">\n";
          pout << kmers[ kmer_ids[0] ];
          for ( int i = 1; i < kmer_ids.isize( ); i++ )
               if ( kmer_ids[i] != kmer_ids[i-1] ) pout << kmers[ kmer_ids[i] ][4];
          pout << "\n";    }

     // Assess.

     SystemSucceed( "QueryLookupTable K=12 MM=12 MC=0.01 SEQS=qqq.fasta "
          "L=/wga/scr4/macro/ecoli/scs.lookup FW_ONLY=True "
          "SMITH_WAT=True PARSEABLE=True BW_ADD=300 > woof3.aligns" );
     vec<look_align> aligns3;
     LoadLookAligns( "woof3.aligns", aligns3 );
     PRINT( aligns3.size( ) );
     if ( !aligns3.solo( ) ) 
     {    cout << "no solo alignment" << endl;
          continue;    }
     vecbasevector G( "/wga/scr4/macro/ecoli/scs.fastb" );
     vecbasevector R3;
     FetchReads( R3, 0, "qqq.fasta" );
     int gid3 = aligns3[0].target_id;
     align a3( aligns3[0].a );
     {
     vec<ho_interval> perf1;
     a3.PerfectIntervals1( R3[0], G[gid3], perf1 );
     int nperf3 = 0;
     for ( int i = 0; i < perf1.isize( ); i++ )
     {    if ( perf1[i].Length( ) < 10 ) continue;
          nperf3 += perf1[i].Length( );    }
     cout << PERCENT_RATIO( 3, nperf3, R3[0].isize( ) )
          << " of bases in matching 10-mers" << endl;

               double score = double(nperf3) / R3[0].size( );
               PRINT(score);
               if ( score > best_score ) 
               {    best_score = score;
                    P = Pnew;    
                    cout << "score = " << score << ", P = " << printSeq(P) << endl;
                         }    }    }

     cout << "\n" << Date( ) << ": done" << endl;
     }
     Scram(0);

     // Define graph.

     cout << Date( ) << ": defining graph" << endl;
     vec<String> verts;
     verts.push_back( "begin" );
     int n = events.size( );
     for ( int i = 0; i < model.isize( ); i++ )
     for ( int j = 0; j < events.isize( ); j++ )
          verts.push( model[i].first + ".e" + ToString(j+1) );
     verts.push( "end" );
     Sort(verts);
     int N = verts.size( );
     PRINT(N);
     vec<vec<int>> from(N), to(N);
     vec<vec<int>> from_edge_obj(N), to_edge_obj(N);
     vec<double> edges;

     // Create vertex index.

     vec< vec<int> > vi( model.size( ) );
     for ( int i = 0; i < model.isize( ); i++ )
          vi[i].resize( events.size( ) );
     #pragma omp parallel for
     for ( int i = 0; i < verts.isize( ); i++ )
     {    if ( verts[i] == "begin" || verts[i] == "end" ) continue;
          String kmer = verts[i].Before( "." );
          int p = BinPosition( kmers, kmer );
          int e = verts[i].After( "e" ).Int( ) - 1;
          vi[p][e] = i;    }

     for ( int i = 0; i < model.isize( ); i++ )
     {    // begin --> kmer.e1
          from[ BinPosition( verts, String("begin") ) ].push_back( 
               BinPosition( verts, model[i].first + ".e1" ) );
          to[ BinPosition( verts, model[i].first + ".e1" ) ]
               .push_back( BinPosition( verts, String("begin") ) );
          from_edge_obj[ BinPosition( verts, String("begin") ) ].push_back(
               edges.size( ) );
          to_edge_obj[ BinPosition( verts, model[i].first + ".e1" ) ].push_back(
               edges.size( ) );
          edges.push_back(0);
          // kmer.en --> end
          from[ BinPosition( verts, model[i].first + ".e" + ToString(n) ) ].
               push_back( BinPosition( verts, String("end") ) );
          to[ BinPosition( verts, String("end") ) ]
               .push_back(
                    BinPosition( verts, model[i].first + ".e" + ToString(n) ) );
          from_edge_obj[ BinPosition( verts, model[i].first + ".e" + ToString(n) ) ]
               .push_back( edges.size( ) );
          to_edge_obj[ BinPosition( verts, String("end") ) ].push_back( 
               edges.size( ) );
          edges.push_back(
               Delta( events[n-1].second, events[n-1].third,
                    model[i].second
                    /* model_level[n-1] */, model[i].third ) );    }
     cout << Date( ) << ": woof" << endl;
     PRINT2( model.size( ), events.size( ) );
     for ( int i = 0; i < model.isize( ); i++ )
     for ( int j = 0; j < events.isize( ); j++ )
     {    int v1 = vi[i][j];
          int v2 = ( j < events.isize( ) - 1 ? vi[i][j+1] : -1 );
          // Type 3, manufactured event
          if ( j < events.isize( ) - 1 )
          {    from[v1].push_back(v2), to[v2].push_back(v1);
               from_edge_obj[v1].push_back( edges.size( ) );
               to_edge_obj[v2].push_back( edges.size( ) );
               // double add_penalty = 2.0;
               double add_penalty = 3.0; // perhaps marginally better
               edges.push_back( add_penalty +
                    Delta( events[j].second, events[j].third,
                         model[i].second
                         /* model_level[j] */, model[i].third ) );    }    }
     cout << Date( ) << ": snort" << endl;
     for ( int i = 0; i < model.isize( ); i++ )
     for ( int m = 0; m < 4; m++ )
     {    int ip = BinPosition( kmers, 
               model[i].first.substr( 1, 4 ) + String( as_base(m) ) );
          for ( int j = 0; j < events.isize( ); j++ )
          {    int v1 = vi[i][j];
               int v2 = ( j < events.isize( ) - 1 ? vi[i][j+1] : -1 );
               int v1p = vi[ip][j];
               int v2p = ( j < events.isize( ) - 1 ? vi[ip][j+1] : -1 );
               // Type 1
               if ( j < events.isize( ) - 1 )
               {    from[v1].push_back(v2p), to[v2p].push_back(v1);
                    from_edge_obj[v1].push_back( edges.size( ) );
                    to_edge_obj[v2p].push_back( edges.size( ) );
                    edges.push_back(
                         Delta( events[j].second, events[j].third,
                              model[i].second
                              /* model_level[j] */, model[i].third ) );    }

               // Type 2, don't see model as event
               from[v1].push_back(v1p), to[v1p].push_back(v1);
               from_edge_obj[v1].push_back( edges.size( ) );
               to_edge_obj[v1p].push_back( edges.size( ) );

               // DEFAULT
               double drop_penalty = 2.5;
               edges.push_back( drop_penalty + 
                    Delta( events[j].second, events[j].third,
                         model[ip].second
                         /* model_level[j] */, model[ip].third ) );
               //

               // RUNS FOREVER
               // double drop_penalty = 1.0;
               // edges.push_back( drop_penalty + 
               //      Delta( events[j].second, events[j].third,
               //           /* model[i].second */ model_level[j], model[i].third )
               //      - Delta( events[j].second, events[j].third,
               //           /* model[ip].second */ model_level[j], model[ip].third ) );
               //

                    }    }




     /*
     for ( int i = 0; i < model.isize( ); i++ )
     for ( int m1 = 0; m1 < 4; m1++ )
     for ( int m2 = 0; m2 < 4; m2++ )
     {    int ip = BinPosition( kmers, 
               model[i].first.substr( 2, 3 ) + String( as_base(m1) ) 
                    + String( as_base(m2) ) );
          for ( int j = 0; j < events.isize( ) - 1; j++ )
          {    int v1 = vi[i][j];
               int v2 = ( j < events.isize( ) - 1 ? vi[i][j+1] : -1 );
               int v1p = vi[ip][j];
               int v2p = ( j < events.isize( ) - 1 ? vi[ip][j+1] : -1 );
               // Type 2, don't see model as event
               from[v1].push_back(v2p), to[v2p].push_back(v1);
               from_edge_obj[v1].push_back( edges.size( ) );
               to_edge_obj[v2p].push_back( edges.size( ) );
               double drop_penalty = 2.5;
               edges.push_back( drop_penalty + 
                    Delta( events[j].second, events[j].third,
                         // model[i].second
                         model_level[j], model[i].third ) );    }    }
     */

     cout << Date( ) << ": wheeze" << endl;
     #pragma omp parallel for
     for ( int i = 0; i < N; i++ )
     {    SortSync( from[i], from_edge_obj[i] );
          SortSync( to[i], to_edge_obj[i] );    }

     digraphE<double> D( from, to, edges, to_edge_obj, from_edge_obj );

     cout << Date( ) << ": getting shortest path\n" << endl;
     vec<int> pathx;
     D.ShortestPath( BinPosition( verts, String("begin") ), 
          BinPosition( verts, String("end") ), pathx );

     int equals = 0;
     int xcount = 0;
     int movesum = 0;
     for ( int i = 1; i < pathx.isize( ) - 1; i++ )
     {    int e = verts[ pathx[i] ].After( "e" ).Int( ) - 1;
          movesum += move[e];
          String kmer = verts[ pathx[i] ].Before( "." );
          int p = BinPosition( kmers, kmer );
          int q = BinPosition( kmers, events[e].first );
          double obs = events[e].second;
          cout 
               // << xcount 

               /*
               << e
               << "\t"
               */

               << movesum
               << "\t"
               << obs << "\t" << verts[ pathx[i] ] 
               << "[" << ToString( 
                    // Square( Abs( obs - model[p].second ) )
                    Delta( events[e].second, events[e].third,
                         /* model[p].second */ model_level[e], model[p].third )
                    , 3 ) 
               << "] "
               << '\t' << events[e].first 
               << "[" << ToString( 
                    // Square( Abs( obs - model[q].second ) )
                    Delta( events[e].second, events[e].third,
                         /* model[q].second */ model_level[e], model[q].third )
                    , 3 ) 
               << "]";

          cout << "\t";
          if ( p == q ) 
          {    cout << "equal";
               equals++;    }

          vec<int> h;
          for ( int j = 0; j < ghits[q].isize( ); j++ )
          {    if ( Abs( ghits[q][j] - movesum ) > 10 + movesum/10 ) continue;
               h.push_back( ghits[q][j] );    }
          cout << "\t" << printSeq(h);
          cout << endl;
          if ( e == 0 || events[e].first != events[e-1].first ) xcount++;    }

     cout << "\nequals = " << equals << endl;

     // Generate fasta file.

     {    Ofstream( out, "ppp.fasta" );
          out << ">\n";
          int count = 0;
          for ( int i = 1; i < pathx.isize( ) - 1; i++ )
          {    if ( i == 1 || verts[ pathx[i] ].Before( "." ) 
                    != verts[ pathx[i-1] ].Before( "." ) )
               {    if ( count > 0 && count % 80 == 0 ) out << "\n";
                    count++;
                    out << verts[ pathx[i] ][4];    }    }
          out << endl;    }

     // Align our calls to reference.

     SystemSucceed( "QueryLookupTable K=12 MM=12 MC=0.01 "
          "SEQS=ppp.fasta "
          "L=/wga/scr4/macro/ecoli/scs.lookup FW_ONLY=True "
          "SMITH_WAT=True PARSEABLE=True BW_ADD=300 > woof2.aligns" );

     // Compare results.

     vec<look_align> aligns2;
     LoadLookAligns( "woof2.aligns", aligns2 );
     PRINT( aligns2.size( ) );
     if ( !aligns2.solo( ) ) return;
     // vecbasevector G( "/wga/scr4/macro/ecoli/scs.fastb" );
     vecbasevector R2;
     FetchReads( R2, 0, "ppp.fasta" );
     int gid2 = aligns2[0].target_id;
     align a2( aligns2[0].a );
     {
     vec<ho_interval> perf1;
     a2.PerfectIntervals1( R2[0], G[gid2], perf1 );
     int nperf2 = 0;
     for ( int i = 0; i < perf1.isize( ); i++ )
     {    if ( perf1[i].Length( ) < 10 ) continue;
          nperf2 += perf1[i].Length( );    }
     cout << PERCENT_RATIO( 3, nperf2, R2[0].isize( ) )
          << " of bases in matching 10-mers" << endl;
     }

     // Done.

     // Scram(0);    
          }

int main( )
{    RunTime( );

     vec<String> rns;
 
     /*
     vec<String> all = AllFiles( "/wga/scr4/vendor/on/2014-07-22/ecoli_r7/convert" );
     for ( int i = 0; i < all.isize( ); i++ )
     {    if ( all[i].Contains( ".fast5.fasta", -1 ) )
                    rns.push_back( all[i].RevBefore( ".fast5.fasta" ) );    }
     */

     // Note that alignments are required to be fw so "no aligns" does not mean
     // that the read is junk.

     // rns.push_back( "strands_500024_ch177_file6" );   // 16.1%; 16.7%
     // rns.push_back( "strands_500024_ch206_file5" );   // 22.6%; no alignment
     // perturbs to 18.9%: P = 1,-2.95,0.15,1.65,2.65,0.45

     // rns.push_back( "strands_500024_ch96_file3" );    // no aligns
     // rns.push_back( "strands_500024_ch8_file19" );    // 16.0%; no alignment
     // rns.push_back( "strands_500024_ch121_file18" );  // no aligns
     // rns.push_back( "strands_500024_ch121_file22" );  // 11.9%; no alignment
     // rns.push_back( "strands_500024_ch121_file24" );  // 11.3%; no alignment
     // rns.push_back( "strands_500024_ch136_file15" );  // no aligns
     // rns.push_back( "strands_500024_ch136_file2" );   // no aligns
     // rns.push_back( "strands_500024_ch142_file17" );  // no aligns
     // rns.push_back( "strands_500024_ch149_file16" );  // no aligns
     // rns.push_back( "strands_500024_ch199_file5" );   // 2 aligns
     // rns.push_back( "strands_500024_ch147_file7" );   // no aligns

     // CHANNEL 177 (with MC lowered from 0.2 to 0.1)
     // rns.push_back( "strands_500024_ch177_file0" );   // 2 aligns; short read
     // rns.push_back( "strands_500024_ch177_file16" );  // no aligns

     // rns.push_back( "strands_500024_ch177_file2" );   // 21.3%; 9.78%
     // tried to optimize but not finding solo alignment

     // rns.push_back( "strands_500024_ch177_file4" );  // 21.8%
     // perturbs to 17.9%: P = 1, 0.5,  0.75, 1.85, 2.2, 0.05

     // rns.push_back( "strands_500024_ch177_file12" );  // 20.4%; 14.1%
     // tried to optimize but not finding solo alignment

     // rns.push_back( "strands_500024_ch177_file20" );  // 19.0%; 9.56%
     // perturbs to 16.3%: P = 1, 0.65, 0.05, 1.5,  2.1, 0.45

     // rns.push_back( "strands_500024_ch177_file21" );     // 23.7%; 8.08%
     // perturbs to 20.8%: P = 1, 0.9,  0.0,  1.5,  2.3, 0.65

     // rns.push_back( "strands_500024_ch177_file24" );     // 18.3%
     // tried to optimize but not finding solo alignment

     // rns.push_back( "strands_500024_ch177_file7" );     // 26.1%
     // tried to optimize but not finding solo alignment (maybe terrible)

     // rns.push_back( "strands_500024_ch177_file6" );     // 16.1%
     // perturbs to 20.7%: 1,0.4,0.02,1.42,1.98,0.04

     // rns.push_back( "strands_500024_ch177_file12" );  // 20.4%; 16.6%

     rns.push_back( "strands_500024_ch121_file22" );

     for ( int i = 0; i < rns.isize( ); i++ )
     {    PRINT( rns[i] );
          Stuff( rns[i] );    }    }
