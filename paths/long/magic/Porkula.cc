///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Trying to match ONT basecalling.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "FastIfstream.h"
#include "FetchReads.h"
#include "math/HoInterval.h"
#include "MainTools.h"
#include "PrintAlignment.h"
#include "TokenizeString.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "pairwise_aligners/ClusterAligner.h"
#include "paths/LongReadTools.h"
#include "paths/long/CreateGenome.h"

double Align( 
     const vec<pair<double,double>>& R, // read (events) -- (mean,dev)
     const vec<pair<double,double>>& A, // kmers (model) -- (mean,dev)
     ostringstream& out,                // for logging
     vec<int>& kmers,                   // kmer sequence in best path
     const vec<vec<int>>& next,         // next[j][l] (l=0..3) denotes next kmers
     vec<int>& p                        // best path
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

     auto pos = [=]( int i, int j ){ return 2 + (i * na) + j; };
     auto r_ind = [=]( int pos ){ return (pos - 2) / na; };
     auto a_ind = [=]( int pos ){ return (pos - 2) % na; };

     auto delta = [&]( int i, int j )
     {    return ( R[i].first - A[j].first ) * ( R[i].first - A[j].first )
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
                    x.push( pos(i, next[j][l]), 1.5 + delta( i, next[j][l] ) );

               if ( i < nr - 1 )                    // edge type 5
                    x.push( pos(i+1, j), 1.0 + deltaR( i, i+1 ) );    }

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

class verb_control {
     public:
     Bool compare_base_calling;
};

double Evaluate( const basevector& q, const String& query, const int id, 
     const int nreads, const String& instance, int& start, int& stop, String& GS, 
     const vecbasevector& GG, const int K, const VecIntPairVec& Glocs, 
     ostringstream& sout, const Bool VISUAL, vec<ho_interval>& perfs )

{    // Align the query.

     vec<look_align> aligns;
     int BW_ADD = 300;
     const int MIN_CLUSTER = 20;
     const int MAX_OFFSET_DIFF = 1000;
     const int MISMATCH_PENALTY = 2;
     const int GAP_PENALTY = 3;
     Bool FW_ONLY = True;
     ClusterAligner( q, GG, K, Glocs, aligns, FW_ONLY, BW_ADD, MIN_CLUSTER, 
          MAX_OFFSET_DIFF, MISMATCH_PENALTY, GAP_PENALTY );
     if ( aligns.nonempty( ) && instance == /* "2" */ "1c" )
     {    sout << "LOC: " << id << " at " << aligns[0].target_id
               << "." << aligns[0].a.pos2( ) << "-" << aligns[0].a.Pos2( )
               << endl;    }
     if (VISUAL)
     {    for ( int i = 0; i < aligns.isize( ); i++ )
          {    sout << "\nalignment " << i+1 << endl;
               int g = aligns[i].target_id;
               const align& a = aligns[i].a;
               sout << "at " << g << "." << a.pos2( ) << "-" << a.Pos2( ) << endl;
               PrintVisualAlignment( True, sout, q, GG[g], a );    }    }

     vec<Bool> to_remove( aligns.size( ), False );
     for ( int i = 0; i < aligns.isize( ); i++ )
          if ( !aligns[i].FullLength( ) ) to_remove[i] = True;
     EraseIf( aligns, to_remove );

     PRINT_TO( sout, aligns.size( ) );
     if ( aligns.empty( ) )
     {    sout << instance << ": 0% of bases in matching 10-mers" << endl;
          return 0;    }

     // Find perfect matches.

     // vecbasevector R;
     int gid = aligns[0].target_id;
     GS = GG[gid].ToString( );
     align a( aligns[0].a );
     vec<ho_interval> perf1, perf2;
     a.PerfectIntervals1( q, GG[gid], perf1 );
     a.PerfectIntervals2( q, GG[gid], perf2 );
     start = a.pos2( ), stop = a.Pos2( );
     PRINT2_TO( sout, start, stop );

     // Print perfect matches.

     Bool print_perfect_matches = False;
     if (print_perfect_matches)
     {    for ( int i = 0; i < perf1.isize( ); i++ )
          {    if ( perf1[i].Length( ) < 5 ) continue;
               sout << perf1[i].Start( ) << " " << perf2[i].Start( ) << " --> ";
               for ( int j = perf1[i].Start( ); j < perf1[i].Stop( ); j++ )
                    sout << as_base( q[j] );
               sout << "\n";    }    }

     Bool PRINT_PERFECT_INTERVALS = True;
     if (PRINT_PERFECT_INTERVALS)
     {    
          #pragma omp critical
          {    for ( int i = 0; i < perf2.isize( ); i++ )
               {    if ( perf2[i].Length( ) < 6 ) continue;
                    int start = perf2[i].Start( ), stop = perf2[i].Stop( );
                    String idx = "+" + ToString(id);
                    if ( gid == 1 )
                    {    idx = "-" + ToString(id-nreads);
                         start = GG[0].isize( ) - perf2[i].Stop( );
                         stop = GG[0].isize( ) - perf2[i].Start( );    }
                    perfs.push( start, stop );
                    sout << "PERF: " << start << " " << stop << " " << stop - start
                         << " (read " << idx << ") ";
                    for ( int j = start; j < stop; j++ )
                         sout << as_base( GG[0][j] );
                    sout << endl;    }    }    }

     // Score the alignment.

     int nperf1 = 0;
     for ( int i = 0; i < perf1.isize( ); i++ )
     {    if ( perf1[i].Length( ) < 10 ) continue;
          nperf1 += perf1[i].Length( );    }
     sout << instance << ": " << PERCENT_RATIO( 3, nperf1, q.isize( ) )
          << " of bases in matching 10-mers" << endl;
     return double(nperf1) / q.size( );    }

void Stuff( const String& DATA_DIR, const String& rn, const int id, ostringstream& sout,
     const verb_control& vc, const vecbasevector& GG, const VecIntPairVec& Glocs,
     vec< pair<double,double> >& scores, const Bool VISUAL, const Bool PRINT_BASES,
     vecbasevector& calls, const int nreads, vec<ho_interval>& perfs,
     const Bool CONSENSUS )
{
     String rn_loc = DATA_DIR + "/" + rn;

     // Get truth data.

     int start, stop;
     String GS;

     const int K = 12;
     vecbasevector Q;
     FetchReads( Q, 0, rn_loc + ".fast5.fasta" );
     if ( Q.size( ) == 0 ) return;
     double score1 = Evaluate( Q[0], rn_loc + ".fast5.fasta", id, nreads,
          "1", start, stop, GS, GG, K, Glocs, sout, False, perfs );
     scores[id].first = score1;

     if (CONSENSUS)
     {    String con_file = DATA_DIR + "/con_" + rn + ".fast5.fasta";
          if ( IsRegularFile(con_file) )
          {    vecbasevector Q;
               FetchReads( Q, 0, con_file );
               if ( Q.size( ) > 0 )
               {    double score1 = Evaluate( Q[0], con_file, id, nreads,
                         "1c", start, stop, GS, GG, K, Glocs, sout, VISUAL, perfs );
                    scores[id].first = score1;  // NIW -- this was missing
                    calls[id] = Q[0]; // ******************************************
                    return; // ****************************************************
                         }    }
          else FatalErr("Consensus file " + con_file + " missing"); }

     // Load data.

     sout << Date( ) << ": loading" << endl;
     vec< triple<String,double,double> > model; // (kmer,mean,dev)
     vec<String> kmers;
     fast_ifstream in1( rn_loc + ".fast5.model" );
     String line;
     vec<String> tokens;
     sout << Date( ) << ": start while loop" << endl; // XXXXXXXXXXXXXXXXXXXXXXXXXXX
     while(1)
     {    getline( in1, line );
          if ( in1.fail( ) ) break;
          Tokenize( line, ',', tokens );
          if ( tokens.empty( ) ) break;
          model.push( tokens[0], tokens[2].Double( ), tokens[3].Double( ) );
          kmers.push( tokens[0] );    }
     sout << Date( ) << ": resizing" << endl; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     kmers.resize(1024), model.resize(1024);
     SortSync( kmers, model );

     sout << Date( ) << ": forming ghits" << endl; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     vec<vec<int>> ghits(1024);
     if ( GS.nonempty( ) )
     {    for ( int j = start; j <= stop - 5; j++ )
          {    String kmer = GS.substr( j, 5 );
               int p = BinPosition( kmers, kmer );
               ghits[p].push_back( j - start );    }    }

     vec< triple<String,double,double> > events; // (kmer,mean,dev)
     vec<double> model_level;
     vec<int> move;

     // Use existing ONT base calls.

     sout << Date( ) << ": reading events file" << endl; // XXXXXXXXXXXXXXXXXXXXXXXX
     fast_ifstream in2( rn_loc + ".fast5.events" );
     int lcount = 0;
     while(1)
     {    getline( in2, line );
          if ( in2.fail( ) ) break;
          lcount++;

          // FOCUSING ON TINY STRETCH!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          // if ( lcount > 4000 ) break;

          Tokenize( line, ',', tokens );
          if ( tokens.empty( ) ) break;
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

     // events.resize( events.isize( ) - 1 );

     // Define next.

     vec<vec<int>> next( model.size( ) );
     for ( int i = 0; i < model.isize( ); i++ )
     {    for ( int m = 0; m < 4; m++ )
          {    int ip = BinPosition( kmers, 
                    model[i].first.substr( 1, 4 ) + String( as_base(m) ) );
               next[i].push_back(ip);    }    }

     // Normalize by subtracting means.

     sout << Date( ) << ": normalizing" << endl; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     vec<double> modelx, eventsx;
     for ( int i = 0; i < model.isize( ); i++ )
          modelx.push_back( model[i].second );
     for ( int i = 0; i < events.isize( ); i++ )
          eventsx.push_back( events[i].second );
     double model_mean = Mean(modelx), event_mean = Mean(eventsx);
     PRINT2_TO( sout, model_mean, event_mean );
     for ( int i = 0; i < model.isize( ); i++ )
          model[i].second -= model_mean;
     for ( int i = 0; i < events.isize( ); i++ )
          events[i].second -= event_mean;

     // Plot.

     /*
     {    Ofstream( out, "data" );
          int window = 1000;
          for ( int i = 0; i <= eventsx.isize( ) - window; i++ )
          {    double sum = 0.0;
               for ( int j = 0; j < window; j++ )
                    sum += eventsx[i+j];
               sum /= window;
               out << i << " " << sum << endl;    }    }
     */

     // Call Align.

     sout << "\n" << Date( ) << ": calling Align" << endl;
     {
     vec< pair<double,double> > R, A;
     for ( int i = 0; i < events.isize( ); i++ )
          R.push( events[i].second, events[i].third );
     for ( int i = 0; i < model.isize( ); i++ )
          A.push( model[i].second, model[i].third);
     ostringstream out;
     vec<int> kmer_ids;
     vec<int> path;
     double clock = WallClockTime( );
     Align( R, A, out, kmer_ids, next, path );
     double elapsed = WallClockTime( ) - clock;
     double msec_per_event = 1000.0 * elapsed / R.size( );
     PRINT2_TO( sout, R.size( ), msec_per_event );

     // Now compute mean difference.

     /*
     double model_mean = 0, event_mean = 0;
     for ( int i = 0; i < kmer_ids.isize( ); i++ )
     {    double m = model[ kmer_ids[i] ].second;
          int v = path[i];
          int ei = (v-2) / A.size( );
          PRINT4( v, A.size( ), ei, events.size( ) ); // XXXXXXXXXXXXXXXXXXXXXXXXXXX
          double e = events[ei].second;
          model_mean += m;
          event_mean += e;    }
     model_mean /= kmer_ids.size( );
     event_mean /= kmer_ids.size( );
     PRINT2_TO( sout, model_mean, event_mean );
     */

     Bool print_align_output = False;
     if (print_align_output) sout << out.str( );
     // for ( int i = 0; i < kmer_ids.isize( ); i++ )
     //      sout << i << "\t" << kmers[ kmer_ids[i] ] << endl;

     // Print base calling comparison.

     sout << "\n";
     if (vc.compare_base_calling)
     {    int equals = 0, xcount = 0, movesum = 0;
          for ( int i = 1; i < path.isize( ) - 1; i++ )
          {    int e = ( path[i] - 2 ) / A.size( );
               movesum += move[e];
               String kmer = kmers[ kmer_ids[i-1] ];
               int p = BinPosition( kmers, kmer );
               int q = BinPosition( kmers, events[e].first );
               double obs = events[e].second;
               sout << movesum << "\t" << obs << "\t";
               {    int v = path[i];
                    String vn;
                    if ( v == 0 ) vn = "begin";
                    else if ( v == 1 ) vn = "end";
                    else
                    {    int i = (v-2) / A.size( );
                         int j = v - 2 - (i * A.size( ));
                         vn = "r" + ToString(i) + "." + kmer;    }
                    sout << vn;    }
               sout << "[" << ToString( 
                    Delta( events[e].second, events[e].third,
                    model[p].second /* model_level[e] */, model[p].third ) , 3 )
                    << "] " << '\t' << events[e].first 
                    << "[" << ToString( 
                         Delta( events[e].second, events[e].third,
                              model[q].second /* model_level[e] */, 
                              model[q].third ) , 3 )
                    << "]";
               sout << "\t";
               if ( p == q ) 
               {    sout << "equal";
                    equals++;    }
               vec<int> h;
               for ( int j = 0; j < ghits[q].isize( ); j++ )
               {    if ( Abs( ghits[q][j] - movesum ) > 10 + movesum/10 ) continue;
                    h.push_back( ghits[q][j] );    }
               sout << "\t" << printSeq(h);
               sout << endl;
               if ( e == 0 || events[e].first != events[e-1].first ) xcount++;    }
          sout << "\nequals = " << equals << endl;    }

     // Convert to fasta.  Not quite right.

     String read_file = "qqq." + ToString(id) + ".fasta";
     basevector q;
     for ( int i = 1; i < kmer_ids.isize( ); i++ )
     {    if ( kmer_ids[i] != kmer_ids[i-1] ) 
               q.push_back( as_char( kmers[ kmer_ids[i] ][4] ) );    }
     if (PRINT_BASES) q.Print( sout, id );
     calls[id] = q;

     // Assess.

     int start, stop;
     String GS;
     double score2 = Evaluate( q, read_file, id, nreads, "2", start, stop, 
          GS, GG, K, Glocs, sout, VISUAL, perfs );
     scores[id].second = score2;
     // Remove(read_file);

     sout << "\n" << Date( ) << ": done" << endl;    }

          }

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_Int_OrDefault_Doc(N, -1,
          "number of reads to process; all if unspecified");
     CommandArgument_Bool_OrDefault_Doc(VISUAL, False, 
          "print visual alignment of our basecalled read");
     CommandArgument_Bool_OrDefault_Doc(PRINT_BASES, False, 
          "print our basecalled read");
     CommandArgument_String_OrDefault_Doc(OUT_FASTB, "",
          "put our basecalled reads in this fastb file");
     CommandArgument_Bool_OrDefault_Doc(CONSENSUS, False, 
          "use ONT consensus instead of our own basecalls");
     CommandArgument_String_OrDefault_Doc(FC, "",
          "use only this flowcell");
     CommandArgument_String_OrDefault(DATA_DIR, "");
     EndCommandArguments;

     RunTime( );

     // Find reads.

     vec<String> rns;
     if ( DATA_DIR == "" )
         DATA_DIR = "/wga/scr4/vendor/on/2014-07-22/ecoli_r7/convert";
     vec<String> all = AllFiles( DATA_DIR );
     for ( int i = 0; i < all.isize( ); i++ )
     {    if ( all[i].Contains( ".fast5.fasta", -1 ) && !all[i].Contains( "comp" )
               && !all[i].Contains( "con" ) && all[i].Contains(FC) )
          {    rns.push_back( all[i].RevBefore( ".fast5.fasta" ) );    }
          if ( rns.isize( ) == N ) break;    }
     PRINT( rns.size( ) );

     // Define parameters.

     verb_control vc;
     vc.compare_base_calling = False;
     int nreads = -1;
     
     // Set up for evaluation.  Compress homopolymers.

     vecbasevector GG( "/wga/dev/references/Escherichia_coli/genome.fastb" );
     for ( int g = 0; g < (int) GG.size( ); g++ )
     {    basevector x;
          for ( int j = 0; j < GG[g].isize( ); j++ )
          {    if ( j >= 5 && GG[g][j] == GG[g][j-1] && GG[g][j] == GG[g][j-2]
                    && GG[g][j] == GG[g][j-3] && GG[g][j] == GG[g][j-4]
                    && GG[g][j] == GG[g][j-5] )
               {    continue;     }
               x.push_back( GG[g][j] );    }
          GG[g] = x;    }
     basevector Z;
     Z.SetToSubOf( GG[0], 10000, 500 );
     Z.Print( cout, "Z" );

     vecbasevector GGR(GG);
     for ( int g = 0; g < (int) GG.size( ); g++ )
          GGR[g].ReverseComplement( );
     GG.Append(GGR);
     const int K = 12;
     VecIntPairVec Glocs;
     CreateGlocs( GG, K, Glocs );

     // Do the work.

     int rcount = ( nreads < 0 ? rns.isize( ) : nreads );
     PRINT(rcount);
     vec< pair<double,double> > scores(rcount);
     vecbasevector calls(rcount);
     vec<ho_interval> perfs;
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int i = 0; i < rcount; i++ )
     {    ostringstream sout;
          Stuff( DATA_DIR, rns[i], i, sout, vc, GG, Glocs, scores,
               VISUAL, PRINT_BASES, calls, nreads, perfs, CONSENSUS );
          #pragma omp critical
          {    PRINT( rns[i] );
               cout << sout.str( );    }    }
     if ( OUT_FASTB != "" ) calls.WriteAll(OUT_FASTB);

     Sort(perfs);
     int start = perfs[0].Start( ), stop = perfs[0].Stop( );
     for ( int i = 0; i < perfs.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < perfs.isize( ); j++ )
          {    if ( stop - perfs[j].Start( ) < 5 ) break;
               stop = Max( stop, perfs[j].Stop( ) );    }
          cout << "INTERVAL: " << start << "-" << stop << endl;
          if ( j < perfs.isize( ) )
          {    start = perfs[j].Start( ), stop = perfs[j].Stop( );    }
          i = j - 1;    }

     // Report scores.

     cout << "\nSCORES:\n";
     vec<double> s1, s2;
     for ( int i = 0; i < rcount; i++ )
     {    s1.push_back( 100.0 * scores[i].first );
          s2.push_back( 100.0 * scores[i].second );    }
     for ( int i = 0; i < rcount; i++ )
     {    cout << "[" << i+1 << "] " << setprecision(3) << s1[i]
               << "% --> " << s2[i] << "%" << endl;    }
     Sort(s1), Sort(s2);
     cout << "MEDIAN: " << Median(s1) << "% --> " << Median(s2) << "%" << endl;
     cout << "MEAN: " << Mean(s1) << "% --> " << Mean(s2) << "%" << endl;    }
