///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Experimental code to use 10X data to extend from a CN1 edge in a DISCOVAR 
// de novo assembly.
//
// Problems.
// 1. Inefficient.
// 2. Have to force our way over assembly defects.
// 3. Some weaks accepted, making random choice, need to track better.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "FastIfstream.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "kmers/KmerRecord.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "paths/HyperBasevector.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/large/Lines.h"
#include "paths/long/large/tenx/TenxDirs.h"
#include "paths/long/large/tenx/TenxTools.h"
#include "random/Random.h"

template<int K> class walker10 {
     public:
     vec<int> path;
     int kmers;
     vec<uint32_t> bcs;
     vec<uint32_t> bx;
     vec<vec<kmer<K>>> reads;
     vec<vec<int>> tracks;
     vec<Bool> stale;
     vec<int64_t> kb_bci; 
     vec<vec<int>> ehits;

     int Seen( const int dist )
     {    int nbc = 0;
          for ( int b = 0; b < bx.isize( ); b++ )
          {    if ( tracks[b].nonempty( ) && kmers - tracks[b].back( ) <= dist ) 
                    nbc++;    }
          return nbc;    }

     vec<int> BarcodesSeen( const int dist )
     {    vec<int> s;
          for ( int b = 0; b < bx.isize( ); b++ )
          {    if ( tracks[b].nonempty( ) && kmers - tracks[b].back( ) <= dist ) 
                    s.push_back(b);    }
          return s;    }

};

template<int K> void Push( walker10<K>& W, const HyperBasevectorX& hb, 
     const int e, double& push_clock );

template<int K> void AddBarcodes( walker10<K>& W, const int x, 
     const HyperBasevectorX& hb, const vec<int>& inv, const String& tdir, 
     double& add_clock, double& tracks_clock );

template<int K> void FindHits( walker10<K>& W, const HyperBasevectorX& hb, 
     const int v, vec< vec< kmer<K> > >& hits, double& hits_clock, 
     const Bool use_stale = True );

void ClassifyDiffs( const align& a, const basevector& E1, const basevector& E2 );

int main( )
{    RunTime( );
     double clock = WallClockTime( );
     double setup_clock = -WallClockTime( );

     // Hardcoded directories.

     // int N = 1;
     int N = 3;
     String dir, odir, tdir;
     SetTenxDirs( N, dir, odir, tdir );

     // Verbosity control.

     const Bool print_tracks = False;
     const Bool display_final_tracks = False;

     // Set up data structure.

     const int K = 88;
     walker10<K> W;

     // Set up for examples.

     int seed;
     vec<vec<int>> force;

     // Define seed edge.

     // S=8648940,L59411,L41215,L36775,L80579,7988288,9920679,L48589,L24937,
     // 9873304,L6855
     // G=True D=3 LINENO=True COUNT=True SHOW_ALIGN=False
     /*
     seed = 8648940;
     force = {
          {731421,9722108,9111930},  // 9722108 should fold back into 9111930
          {830943,7988288,9920679},  // 7988288 should connect to dead end 9920679
          {2215772,9873304,9807636}, // e2 should fold back into e3(?)
          {10504599,2791068},        // uncaptured gap at frayed ends
          {8801991,8794969},         // uniquely closed by Moombat
          {9721554,10318553},        // uniquely closed by Moombat
          {142074,7353937},          // uncaptured gap
          {217266,9075254},          // just below enough support
          {9075254,7893892},         // uniquely closed by Moombat
          {10549307,7078142},        // closed by Moombat, line with two bubbles
          {8580518,9123459},         // uncaptured gap
          {9123459,8143338},         // uncaptured gap
          {1098130,9728128},         // Hop L=200 ==> 1 path
          {7371187,661500},          // Hop L=600 ==> 36 paths
          {8008795,1155581},         // uncaptured gap with 3 links
          {761235,5505961},          // uncaptured gap
          {410997,362759},           // Hop L=1000 ==> 1 path
          {5870568,8987388},         // uniquely closed by Moombat
          {8987388,8771582},         // uncaptured gap with 1 link
          {6380486,7806241},         // from one line to the next
          {8038907,8678454}          // nasty uncaptured gap
               };
     */

     // parallel haplotype, MUCH sloppier analysis
     /*
     seed = 7999966;
     force = {
          {336276,9775966},          // Hop L=100 ==> 1 path
          {10042442,9340828},        // looks like edge to fold back but didn't check
          {7259589,10504599,2791068},// uncaptured gap (?)
          {8801991,8794969},         // uniquely closed by Moombat
          {9733624,10318553},
          {9089802,9695101},         // Hop L=400 ==> 5 paths
          {9695101,9721642},         // easy hop
          {142074,7353937},          // uncaptured gap
          {9075254,7893892},         // uniquely closed by Moombat (not sure right)
          {7893892,7341260,9123459},
          {217266,7078139},          // from L38193 to L33755
          {8580518,2301231},         // from L33755 to L21685
          {1098130,8002201},         // from L21685 to L8150
          {7371187,7391121},         // from L8150 to L40502
          {8008795,1020842},         // from L40502 to L73659
          {761235,8452297},          // from L73659 to L26793
          {410997,962866},           // from L26793 to L4024
          {5870568,6383376},         // from L4024 to L97634
          {6380486,8521359},         // from L97634 to L9202
          {8038907,8704181}          // from L9202 to L2167
               };

     */

     seed = 8471468;
     force = { 
          {582806,6269133,8499969},
          {1904222,1958045,2011261}, // branch to -chr18, not inferrable
          {3845647,3865595},
          {9439918,5199755},
          {6579498,9017168},
          {7709052,7858322},
          {2916402,2874913}
          };

     // Goes with previous example:
     // Note that without second force, have
     // 54(70580). Dangerous branch, taking e1 = 7988288[184] over e2 = 9814143[132].
     // un1 = 45, un2 = 7
     // bun1.size( ) = 45, bun2.size( ) = 7
     // and wonder where the 7 barcodes come from

// xxx.hap1 xxx.hap2

/*
xxx.hap1, right:
52(77979). Simple branch, taking e1 = 9954408[110] over e2 = 9965965[2].
262(291322). Simple branch, taking e2 = 10395920[33] over e1 = 10394166[6].
358(342305). Simple branch, taking e1 = 10504591[64] over e2 = 10504595[0].
457(394860). Simple branch, taking e2 = 9778097[20] over e1 = 9764419[2].
685(588507). Simple branch, taking e1 = 6318169[13] over e2 = 6318170[0].
854(674961). Simple branch, taking e2 = 6476245[14] over e1 = 6476244[0].
1121(819100). Simple branch, taking e2 = 9711135[17] over e1 = 9103445[0].
1173(850897). Simple branch, taking e2 = 7078144[6] over e1 = 7078143[0].
1194(872192). Simple branch, taking e1 = 7078188[5] over e2 = 7078189[0].

xxx.hap1, now wrong:
1236(954420). Simple branch, taking e2 = 8680547[4] over e1 = 8678509[0].
1279(981353). Simple branch, taking e2 = 10025630[7] over e1 = 9343896[0].
1333(1003404). Simple branch, taking e2 = 10242407[9] over e1 = 9481506[0].
1394(1063177). Simple branch, taking e1 = 5270164[22] over e2 = 5270165[0].
1534(1212235). Simple branch, taking e2 = 10127925[23] over e1 = 10121249[0].
*/


     // Load counts and assembly.

     cout << "\n" << Date( ) << ": loading counts and assembly" << endl;
     HyperBasevectorX hb;
     vec<int> count, inv;
     vec<vec<vec<vec<int>>>> lines;
     vec< vec< pair<int,int> > > aligns;
     #pragma omp parallel sections
     {    
          #pragma omp section
          {    BinaryReader::readFile( odir + "/10X.count", &count );    }
          #pragma omp section
          {    BinaryReader::readFile( dir + "/a.hbx", &hb );    }
          #pragma omp section
          {    BinaryReader::readFile( dir + "/a.inv", &inv );    }
          #pragma omp section
          {    BinaryReader::readFile( dir + "/a.lines", &lines );    }
          #pragma omp section
          {    BinaryReader::readFile( dir + "/a.aligns", &aligns );    }    }
     vec< pair<int,int> > tol2( hb.E( ) );
     for ( int l = 0; l < lines.isize( ); l++ )
     for ( int i = 0; i < lines[l].isize( ); i++ )
     for ( int j = 0; j < lines[l][i].isize( ); j++ )
     for ( int k = 0; k < lines[l][i][j].isize( ); k++ )
          tol2[ lines[l][i][j][k] ] = make_pair( l, i );
     vec<int> lens;
     GetLineLengths( hb, lines, lens );
     vec<String> genome_names;
     fast_ifstream in( dir + "/../genome.names" );
     String line;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          genome_names.push_back(line);    }

     // Load alignments and barcodes and other 10X data structures.

     cout << Date( ) << ": loading barcodes" << endl;
     vec<double> bc_frac;
     vec< vec<int> > lhitsb;   // line to barcodes
     vec<vec<int>> lbc;        // barcode to lines
     #pragma omp parallel sections
     {    
          #pragma omp section
          {    BinaryReader::readFile( tdir + "/10X.bc", &W.bcs );    }
          #pragma omp section
          {    BinaryReader::readFile( tdir + "/10X.kb_bci", &W.kb_bci );    }
          #pragma omp section
          {    BinaryReader::readFile( odir + "/10X.ehits", &W.ehits );    }
          #pragma omp section
          {    BinaryReader::readFile( odir + "/10X.bc_frac", &bc_frac );    }
          #pragma omp section
          {    BinaryReader::readFile( odir + "/10X.lhitsb", &lhitsb );    }
          #pragma omp section
          {    BinaryReader::readFile( odir + "/10X.lbc", &lbc );    }    }
     setup_clock += WallClockTime( );

     // Identify initial barcodes.

     cout << Date( ) << ": identifying barcodes" << endl;
     for ( int i = 0; i < W.ehits[seed].isize( ); i++ )
          W.bx.push_back( W.bcs[ W.ehits[seed][i] ] );
     for ( int i = 0; i < W.ehits[ inv[seed] ].isize( ); i++ )
          W.bx.push_back( W.bcs[ W.ehits[ inv[seed] ][i] ] );
     UniqueSort(W.bx);
     cout << "using " << W.bx.size( ) << " barcodes initially" << endl;

     // Load reads.

     cout << Date( ) << ": loading reads now" << endl;
     W.reads.resize( W.bx.size( ) );
     for ( int i = 0; i < W.bx.isize( ); i++ )
     {    W.reads[i].ReadRange( 
               tdir + "/10X.kb", W.kb_bci[ W.bx[i] ], W.kb_bci[ W.bx[i] + 1 ] );    }

     // Initialize tracks.

     cout << Date( ) << ": initializing tracks" << endl;
     W.tracks.resize( W.bx.size( ) );
     for ( int l = 0; l <= hb.Bases(seed) - K; l++ )
     {    kmer<K> x;
          x.SetToSubOf( hb.EdgeObject(seed), l );
          kmer<K> y(x);
          y.ReverseComplement( );
          #pragma omp parallel for
          for ( int b = 0; b < W.bx.isize( ); b++ )
          {    if ( BinMember( W.reads[b], x ) ) W.tracks[b].push_back(l);
               else
               {    if ( BinMember( W.reads[b], x ) )
                         W.tracks[b].push_back(l);    }    }    }
     if (print_tracks)
     {    cout << Date( ) << ": initial tracks:" << endl;
          for ( int b = 0; b < W.bx.isize( ); b++ )
               cout << W.bx[b] << ": " << printSeq( W.tracks[b] ) << endl;    }
     W.stale.resize( W.bx.size( ), False );

     // Initialize path.

     W.kmers = 0;
     int weaks = 0;
     const int weaks_allowed = 300;

     // Timers.

     double push_clock = 0;
     double eval_solved_clock = 0;
     double add_clock = 0;
     double hits_clock = 0;
     double tracks_clock = 0;

     // Go forward.

     cout << Date( ) << ": begin walk\n" << endl;
     Push( W, hb, seed, push_clock );
     cout << "Starting at " << seed << "." << endl;
     int good_choices = 0, random_choices = 0, force_choices = 0;
     int captured_gaps = 0;
     vec< quad<int,int,int,int> > unsolved_bubbles, solved_bubbles;
     vec<int> unsolved_bubbles_kmers, solved_bubbles_kmers;
     while(1)
     {    int v = hb.ToRight( W.path.back( ) );

          // Force.

          Bool forced = False;
          for ( int i = 0; i < force.isize( ); i++ )
          {    if ( force[i][0] == W.path.back( ) )
               {    for ( int j = 1; j < force[i].isize( ); j++ )
                    {    Push( W, hb, force[i][j], push_clock );
                         cout << "Forcing " << force[i][j] << "." << endl;    }
                    forced = True;    }
               if (forced) break;    }
          if (forced) 
          {    force_choices++;
               continue;    }

          // Dead end.

          if ( hb.From(v).size( ) == 0 )
          {    cout << "Stopped." << endl;
               break;    }

          // Unique extension.

          if ( hb.From(v).size( ) == 1 )
          {    int e = hb.IFrom( v, 0 );
               Push( W, hb, e, push_clock );
               cout << "Taking unique " << e;
               if ( hb.Bases(e) == 0 ) 
               {    cout << " (gap)";
                    captured_gaps++;    }
               else cout << " (L" << tol2[e].first << ")";
               cout << "; " << W.Seen(10000)
                    << " barcodes seen in last 10 kb." << endl;

               // Analyzed barcodes.  Commented out because it doesn't help.

               /*
               cout << "Analyzing barcodes." << endl;
               vec<int> nbc;
               for ( int i = 0; i < W.ehits[e].isize( ); i++ )
                    nbc.push_back( W.bcs[ W.ehits[e][i] ] );
               for ( int i = 0; i < W.ehits[ inv[e] ].isize( ); i++ )
                    nbc.push_back( W.bcs[ W.ehits[ inv[e] ][i] ] );
               UniqueSort(nbc);
               for ( int i = 0; i < nbc.isize( ); i++ )
               {    int b = nbc[i];
                    if ( BinMember( W.bx, b ) ) continue;
                    vec< vec< kmer<K> > > reads(1);
                    {    reads[0].ReadRange( tdir + "/10X.kb", 
                              W.kb_bci[b], W.kb_bci[ b + 1 ] );    }
                    int supp1 = 0, supp2 = 0;
                    for ( int i = 0; i < solved_bubbles.isize( ); i++ )
                    {    int e1 = solved_bubbles[i].first; 
                         int e2 = solved_bubbles[i].second;
                         int u1 = solved_bubbles[i].third; 
                         int u2 = solved_bubbles[i].fourth;
                         int v = hb.ToLeft(e1);
                         vec< vec< kmer<K> > > hits(2);
                         for ( int j = 0; j < 2; j++ )
                         {    int e = hb.IFrom( v, j );
                              #pragma omp parallel for
                              for ( int l = 0; l <= hb.Bases(e) - K; l++ )
                              {    kmer<K> x;
                                   x.SetToSubOf( hb.EdgeObject(e), l );
                                   Bool found = False;
                                   int64_t p = BinPosition( reads[0], x );
                                   if ( p >= 0 )
                                   {    
                                        #pragma omp critical
                                        {    hits[j].push_back(x);    }
                                        found = True;    }
                                   if ( !found )
                                   {    kmer<K> y(x);
                                        y.ReverseComplement( );
                                        int64_t p = BinPosition( reads[0], y );
                                        if ( p >= 0 )
                                        {    
                                             #pragma omp critical
                                             {    hits[j].push_back(x);   }    }    }
                                                  }
                              Sort( hits[j] );    }
                         int u1_new = 0, u2_new = 0;
                         for ( int j = 0; j < hits[0].isize( ); j++ )
                              if ( !BinMember( hits[1], hits[0][j] ) ) u1_new++;
                         for ( int j = 0; j < hits[1].isize( ); j++ )
                              if ( !BinMember( hits[0], hits[1][j] ) ) u2_new++;
                         if ( u1 > u2 )
                         {    supp1 += u1_new;
                              supp2 += u2_new;    }
                         if ( u2 > u1 ) 
                         {    supp1 += u2_new;
                              supp2 += u1_new;    }    }
                    PRINT3( b, supp1, supp2 );    }
               */
               continue;    }

          // Simple bubble.

          if ( hb.From(v).size( ) == 2 && hb.From(v)[0] == hb.From(v)[1] )
          {    int e1 = hb.IFrom( v, 0 ), e2 = hb.IFrom( v, 1 );
               vec< vec< kmer<K> > > hits;
               FindHits( W, hb, v, hits, hits_clock );
               int u1 = 0, u2 = 0;
               for ( int j = 0; j < hits[0].isize( ); j++ )
                    if ( !BinMember( hits[1], hits[0][j] ) ) u1++;
               for ( int j = 0; j < hits[1].isize( ); j++ )
                    if ( !BinMember( hits[0], hits[1][j] ) ) u2++;

               if ( u1 > u2 && u1 >= 4 )
               {    Push( W, hb, e1, push_clock );
                    cout << "Simple branch, taking e1 = " << e1 << "[" << u1 
                         << "] over e2 = " << e2 << "[" << u2 << "]." << endl;
                    AddBarcodes( W, e1, hb, inv, tdir, add_clock, tracks_clock );
                    good_choices++;
                    solved_bubbles.push( e1, e2, u1, u2 );
                    solved_bubbles_kmers.push_back(W.kmers);
                    continue;    }
               if ( u2 > u1 && u2 >= 4 )
               {    Push( W, hb, e2, push_clock );
                    cout << "Simple branch, taking e2 = " << e2 << "[" << u2 
                         << "] over e1 = " << e1 << "[" << u1 << "]." << endl;
                    AddBarcodes( W, e2, hb, inv, tdir, add_clock, tracks_clock );
                    good_choices++;
                    solved_bubbles.push( e1, e2, u1, u2 );
                    solved_bubbles_kmers.push_back(W.kmers);
                    continue;    }
               // if ( u1 != u2 )
               {    cout << "Weak data at simple "
                         << "bubble, e1 = " << e1 << "[" << u1 << "], e2 = " << e2 
                         << "[" << u2 << "], giving up." << endl;
                    vec<int> nbc1, nbc2;
                    for ( int i = 0; i < W.ehits[e1].isize( ); i++ )
                         nbc1.push_back( W.bcs[ W.ehits[e1][i] ] );
                    for ( int i = 0; i < W.ehits[ inv[e1] ].isize( ); i++ )
                         nbc1.push_back( W.bcs[ W.ehits[ inv[e1] ][i] ] );
                    for ( int i = 0; i < W.ehits[e2].isize( ); i++ )
                         nbc2.push_back( W.bcs[ W.ehits[e2][i] ] );
                    for ( int i = 0; i < W.ehits[ inv[e2] ].isize( ); i++ )
                         nbc2.push_back( W.bcs[ W.ehits[ inv[e2] ][i] ] );

                    // For all barcodes that are assigned to e1 or e2, evaluate
                    // versus solved bubbles.  Try to recover.
                    //
                    // Note that this reads barcodes that are known.
                    // Need to think about it.

                    eval_solved_clock -= WallClockTime( );
                    int supp1 = 0, supp2 = 0;
                    for ( int pass = 1; pass <= 2; pass++ )
                    {    vec<uint32_t> bx;
                         if ( pass == 1 )
                         {    for ( int j = 0; j < nbc1.isize( ); j++ )
                                   bx.push_back( nbc1[j] );    }
                         else
                         {    for ( int j = 0; j < nbc2.isize( ); j++ )
                                   bx.push_back( nbc2[j] );    }
                         vec< vec< kmer<K> > > reads( bx.size( ) );
                         for ( int i = 0; i < bx.isize( ); i++ )
                         {    reads[i].ReadRange( tdir + "/10X.kb", 
                                   W.kb_bci[ bx[i] ], W.kb_bci[ bx[i] + 1 ] );    }
                         for ( int i = 0; i < solved_bubbles.isize( ); i++ )
                         {    int e1 = solved_bubbles[i].first; 
                              int e2 = solved_bubbles[i].second;
                              int u1 = solved_bubbles[i].third; 
                              int u2 = solved_bubbles[i].fourth;
                              int v = hb.ToLeft(e1);
                              vec< vec< kmer<K> > > hits(2);
                              for ( int j = 0; j < 2; j++ )
                              {    int e = hb.IFrom( v, j );
                                   #pragma omp parallel for
                                   for ( int l = 0; l <= hb.Bases(e) - K; l++ )
                                   {    kmer<K> x;
                                        x.SetToSubOf( hb.EdgeObject(e), l );
                                        Bool found = False;
                                        for ( int b = 0; b < bx.isize( ); b++ )
                                        {    int64_t p = BinPosition( reads[b], x );
                                             if ( p >= 0 )
                                             {    
                                                  #pragma omp critical
                                                  {    hits[j].push_back(x);    }
                                                  found = True;
                                                  break;    }    }    
                                        if ( !found )
                                        {    kmer<K> y(x);
                                             y.ReverseComplement( );
                                             for ( int b = 0; b < bx.isize( ); b++ )
                                             {    int64_t p 
                                                       = BinPosition( reads[b], y );
                                                  if ( p >= 0 )
                                                  {    
                                                       #pragma omp critical
                                                       {    hits[j].push_back(x);   }
                                                       found = True;
                                                       break;    }    }    }    }
                                   Sort( hits[j] );    }
                              int u1_new = 0, u2_new = 0;
                              for ( int j = 0; j < hits[0].isize( ); j++ )
                                   if ( !BinMember( hits[1], hits[0][j] ) ) 
                                        u1_new++;
                              for ( int j = 0; j < hits[1].isize( ); j++ )
                                   if ( !BinMember( hits[0], hits[1][j] ) ) 
                                        u2_new++;

                              if ( ( pass == 1 && u1 > u2 )
                                   || ( pass == 2 && u2 > u1 ) )
                              {    supp1 += u1_new;
                                   supp2 += u2_new;    }
                              else
                              {    supp1 += u2_new;
                                   supp2 += u1_new;    }

                              // if ( u1_new > 0 || u2_new > 0 )
                              //      PRINT6( e1, e2, u1, u2, u1_new, u2_new );     
                                   }    }
                    cout << "eval of solved bubbles: ";
                    PRINT2( supp1, supp2 );
                    Bool win1 = ( supp1 >= 5 && supp1 >= 5 * supp2 );
                    Bool win2 = ( supp2 >= 5 && supp2 >= 5 * supp1 );
                    eval_solved_clock += WallClockTime( );
                    if ( win1 || win2 )
                    {    cout << "Aha, saved!" << endl;
                         if (win1)
                         {    Push( W, hb, e1, push_clock );
                              cout << "Simple branch, taking e1 = " << e1 << "[" 
                                   << u1 << "] over e2 = " << e2 << "[" << u2 
                                   << "]." << endl;
                              AddBarcodes( W, e1, hb, inv, tdir,
                                   add_clock, tracks_clock );
                              good_choices++;
                              solved_bubbles.push( e1, e2, u1, u2 );
                              solved_bubbles_kmers.push_back(W.kmers);
                              continue;    }
                         if (win2)
                         {    Push( W, hb, e2, push_clock );
                              cout << "Simple branch, taking e2 = " << e2 << "[" 
                                   << u2 << "] over e1 = " << e1 << "[" << u1 
                                   << "]." << endl;
                              AddBarcodes( W, e2, hb, inv, tdir,
                                   add_clock, tracks_clock );
                              good_choices++;
                              solved_bubbles.push( e1, e2, u1, u2 );
                              solved_bubbles_kmers.push_back(W.kmers);
                              continue;    }    }

                    // OK, no luck.

                    unsolved_bubbles.push( e1, e2, u1, u2 );
                    unsolved_bubbles_kmers.push_back(W.kmers);
                    const basevector &E1 = hb.EdgeObject(e1); 
                    const basevector &E2 = hb.EdgeObject(e2);
                    alignment al;
                    SmithWatAffine( E1, E2, al );
                    align a = al;

                    ClassifyDiffs( a, E1, E2 );

                    UniqueSort(nbc1), UniqueSort(nbc2);
                    cout << "Total barcodes aligned uniquely to an edge in this "
                         << "bubble = " << nbc1.size( ) << ", "
                         << nbc2.size( ) << "." << endl;

                    cout << W.Seen(10000) 
                         << " barcodes seen within last 10000 bases" << endl;

                    if ( weaks < weaks_allowed )
                    {    ++weaks;
                         Push( W, hb, e1, push_clock );
                         cout << "Making random choice of " << e1 << " (" 
                              << hb.Kmers(e1) << "," << hb.Kmers(e2)
                              << " kmers)." << endl;
                         random_choices++;
                         continue;    }    }
               /*
               else
               {    cout << "Tie at simple branch, e1 = "
                         << e1 << "[" << u1 << "],  e2 = " << e2 << "["
                         << u2 << "], giving up." << endl;    }
               */
               break;    }

          // Dangerous branch.

          if ( hb.From(v).size( ) == 2 )
          {    int e1 = hb.IFrom( v, 0 ), e2 = hb.IFrom( v, 1 );

               vec< vec< kmer<K> > > hits;
               FindHits( W, hb, v, hits, hits_clock );
               int u1 = 0, u2 = 0;
               for ( int j = 0; j < hits[0].isize( ); j++ )
                    if ( !BinMember( hits[1], hits[0][j] ) ) u1++;
               for ( int j = 0; j < hits[1].isize( ); j++ )
                    if ( !BinMember( hits[0], hits[1][j] ) ) u2++;

               int un1 = 0, un2 = 0;
               vec<int> nbc1, nbc2;
               for ( int i = 0; i < W.ehits[e1].isize( ); i++ )
                    nbc1.push_back( W.bcs[ W.ehits[e1][i] ] );
               for ( int i = 0; i < W.ehits[ inv[e1] ].isize( ); i++ )
                    nbc1.push_back( W.bcs[ W.ehits[ inv[e1] ][i] ] );
               for ( int i = 0; i < W.ehits[e2].isize( ); i++ )
                    nbc2.push_back( W.bcs[ W.ehits[e2][i] ] );
               for ( int i = 0; i < W.ehits[ inv[e2] ].isize( ); i++ )
                    nbc2.push_back( W.bcs[ W.ehits[ inv[e2] ][i] ] );
               UniqueSort(nbc1), UniqueSort(nbc2);
               vec<int> bun1, bun2;
               for ( int i = 0; i < nbc1.isize( ); i++ )
               {    if ( BinMember( W.bx, nbc1[i] ) )
                    {    un1++;
                         bun1.push_back( nbc1[i] );    }    }
               for ( int i = 0; i < nbc2.isize( ); i++ )
               {    if ( BinMember( W.bx, nbc2[i] ) )
                    {    un2++;
                         bun2.push_back( nbc2[i] );    }    }
               UniqueSort(bun1), UniqueSort(bun2);

               if ( un1 >= 10 * Max( 1, un2 ) )
               {    Push( W, hb, e1, push_clock );
                    cout << "Strong enough branch, taking e1 = " << e1 << "[" << u1 
                         << "] over e2 = " << e2 << "[" << u2 << "]." << endl;
                    PRINT2( un1, un2 );
                    AddBarcodes( W, e1, hb, inv, tdir, add_clock, tracks_clock );
                    good_choices++;
                    continue;    }
               if ( un2 >= 10 * Max( 1, un1 ) )
               {    Push( W, hb, e2, push_clock );
                    cout << "Strong enough branch, taking e2 = " << e2 << "[" << u2 
                         << "] over e1 = " << e1 << "[" << u1 << "]." << endl;
                    PRINT2( un1, un2 );
                    AddBarcodes( W, e2, hb, inv, tdir, add_clock, tracks_clock );
                    good_choices++;
                    continue;    }
               if ( u1 > u2 )
               {    
                    // Commenting out as I don't think it makes sense:
                    // Push( W, hb, e1, push_clock );
                    cout << "Dangerous branch, could take e1 = " << e1 << "[" << u1 
                         << "] over e2 = " << e2 << "[" << u2 << "]." << endl;
                    PRINT2( un1, un2 );
                    cout << "losing barcodes = " << printSeq(bun2) << endl;    }
               if ( u2 > u1 )
               {    
                    // Commenting out as I don't think it makes sense:
                    // Push( W, hb, e2, push_clock );
                    cout << "Dangerous branch, could take e2 = " << e2 << "[" << u2 
                         << "] over e1 = " << e1 << "[" << u1 << "]." << endl;
                    PRINT2( un1, un2 );
                    cout << "losing barcodes = " << printSeq(bun1) << endl;    }
               if ( u1 == u2 )
               {    cout << "Tie at dangerous branch, e1 = " << e1 << "[" << u1 
                         << "],  e2 = " << e2 << "[" << u2 << "], giving up." 
                         << endl;
                    PRINT2( un1, un2 );    }
               int x1 = tol2[e1].first, y1 = tol2[e1].second;
               int x2 = tol2[e2].first, y2 = tol2[e2].second;
               if ( x1 == x2 && y1 == y2 )
               {    const vec<vec<vec<int>>>& L = lines[x1];
                    vec<int> s = {e1,e2};
                    vec<int> t;
                    for ( int j = 0; j < L[y1].isize( ); j++ )
                         t.push_back( L[y1][j][0] );
                    UniqueSort(s), UniqueSort(t);
                    if ( s == t )
                    {    if ( weaks < weaks_allowed )
                         {    ++weaks;
                              W.path.append( L[y1][0] );
                              for ( int j = 0; j < L[y1][0].isize( ); j++ )
                              {    Push( W, hb, L[y1][0][j], push_clock );    }
                              cout << "Making random choice of path " 
                                   << printSeq( L[y1][0] ) << "." << endl;
                              random_choices++;
                              continue;    }    }    }
               break;    }

          // Cell.

          int li = tol2[ W.path.back( ) ].first, lin = tol2[ W.path.back( ) ].second;
          if ( lin != lines[li].isize( ) - 1 )
          {    if ( weaks < weaks_allowed )
               {    ++weaks;
                    const vec<vec<vec<int>>>& L = lines[li];
                    W.path.append( L[lin+1][0] );
                    for ( int j = 0; j < L[lin+1][0].isize( ); j++ )
                    {    Push( W, hb, L[lin+1][0][j], push_clock );    }
                    cout << "Making random choice of path " 
                         << printSeq( L[lin+1][0] ) << "." << endl;
                    random_choices++;
                    continue;    }    }

          // Something else.

          vec<int> count( hb.From(v).size( ), 0 );
          cout << W.path.size( ) + 1 << ". More complex branch, giving up." << endl;
          break;    }

     // Study edges near last vertex.

     /*
     cout << "\n" << Date( ) << ": Studying edges near last vertex." << endl;
     int vf = hb.ToRight( W.path.back( ) );
     const int max_nears = 1000;
     vec< pair<int,int> > nears;
     for ( int i = 0; i < (int) hb.From(vf).size( ); i++ )
          nears.push( 1, hb.IFrom(vf,i) );
     while(1)
     {    int ns = nears.size( );
          if ( ns >= max_nears ) break;
          for ( int j = 0; j < ns; j++ )
          {    int e = nears[j].second;
               int v =  hb.ToRight(e);
               for ( int i = 0; i < (int) hb.From(v).size( ); i++ )
               {    Bool known = False;
                    int f = hb.IFrom(v,i);
                    for ( int l = 0; l < nears.isize( ); l++ )
                    {    if ( nears[l].second == f )
                         {    known = True;
                              break;   }   }
                    if (known) continue;
                    nears.push( nears[j].first + 1, f );    }    }
          if ( nears.isize( ) == ns ) break;    }
     PRINT( nears.size( ) );
     cout << Date( ) << ": finding kmers" << endl;
     vecbasevector stuff;
     for ( int i = 0; i < nears.isize( ); i++ )
          stuff.push_back( hb.EdgeObject( nears[i].second ) );
     vec< triple<kmer<K>,int,int> > kmers_plus;
     MakeKmerLookup0( stuff, kmers_plus );
     cout << Date( ) << ": searching kmers" << endl;
     vec<int> active= W.BarcodesSeen(10000);
     cout << active.size( ) << " active barcodes" << endl;
     Ofstream( out, "local.fasta" ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     int acount = 0; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     for ( int i = 0; i < active.isize( ); i++ )
     {    int b = active[i];
          basevector p; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          for ( int j = 0; j < W.reads[b].isize( ); j++ )
          {    kmer<K> x = W.reads[b][j];
               x.GetBasevector(p);
               p.Print( out, acount++ ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               for ( int pass = 1; pass <= 2; pass++ )
               {    if ( pass == 2 ) x.ReverseComplement( );
                    int64_t low = LowerBound1( kmers_plus, x );
                    int64_t high = UpperBound1( kmers_plus, x );
                    for ( int64_t m = low; m < high; m++ )
                    {    int nid = kmers_plus[m].second;
                         int pos = kmers_plus[m].third;
                         int e = nears[nid].second;
                         int d = nears[nid].first;
                         PRINT4( e, pos, d, b );    }    }    }    }
     cout << "\n" << Date( ) << ": done\n" << endl;
     */

     // Show view forward.  Start by finding the last 5kb+ line seen.

     int L = -1;
     vec<int> ls, xls;
     for ( int i = 0; i < W.path.isize( ); i++ )
          xls.push_back( tol2[ W.path[i] ].first );
     UniqueSort(xls);
     for ( int i = W.path.isize( ) - 1; i >= 0; i-- )
     {    int l = tol2[ W.path[i] ].first;
          if ( lens[l] >= 5000 )
          {    ls.push_back(l);
               break;    }    }
     if ( ls.nonempty( ) )
     {    cout << "\nLines seen forward from L" << ls[0] << ":\n";

          // Go through the lines.

          vec<String> reports( ls.size( ) );
          #pragma omp parallel for
          for ( int li = 0; li < ls.isize( ); li++ )
          {    int L = ls[li];

               // Find the lines.

               const vec<int>& xb = lhitsb[L];
               vec<int> M;
               for ( int i = 0; i < xb.isize( ); i++ )
                    M.append( lbc[ xb[i] ] );
               Sort(M);
               vec< pair<int,int> > lhn;
               for ( int i = 0; i < M.isize( ); i++ )
               {    int j = M.NextDiff(i);
                    lhn.push( j - i, M[i] );
                    i = j - 1;    }
     
               // Now do the math.

               ostringstream out;
               out << "\n";
               int prints = 0;
               for ( int i = 0; i < lhn.isize( ); i++ )
               {    int l = lhn[i].second, l2 = -1;
                    if ( lens[l] < 1000 ) continue;

                    // Check for line already visited.

                    if ( BinMember( xls, l ) ) continue;
                    if ( l >= 0 && BinMember( xls, l-1 ) && lines[l].front( )[0][0] 
                         == inv[ lines[l-1].back( )[0][0] ] )
                    {    continue;    }
                    if ( l < lines.isize( ) - 1 
                         && BinMember( xls, l+1 ) && lines[l].front( )[0][0] 
                         == inv[ lines[l+1].back( )[0][0] ] )
                    {    continue;    }

                    // Proceed.

                    Bool pair = False;
                    if ( i < lhn.isize( ) - 1 )
                    {    l2 = lhn[i+1].second;
                         if ( lines[l].front( )[0][0] 
                              == inv[ lines[l2].back( )[0][0] ] )
                         {    pair = True;    }    }
                    vec<String> chrs;
                    for ( int j = 0; j < lines[l].isize( ); j++ )
                    {    for ( int r = 0; r < lines[l][j].isize( ); r++ )
                         for ( int s = 0; s < lines[l][j][r].isize( ); s++ )
                         {    int e = lines[l][j][r][s];
                              for ( int m = 0; m < aligns[e].isize( ); m++ )
                              {    chrs.push_back( 
                                        genome_names[aligns[e][m].first] );    }
                              e = inv[e];
                              for ( int m = 0; m < aligns[e].isize( ); m++ )
                              {    chrs.push_back( 
                                        genome_names[aligns[e][m].first] );
                                       }    }    }
                    UniqueSort(chrs);
                    int obs = lhn[i].first, nbc = xb.size( );
                    double expect = nbc * bc_frac[l];
                    double devs = (obs-expect) / sqrt(expect);
                    if ( devs < 25 ) continue;
                    out << "[#" << ++prints << ", " << obs << " barcodes, " 
                         << expect << " expect, " << devs << " devs";
                    out << "] L" << l;
                    if ( !pair ) out << " (l=" << lens[l] << ")";
                    else
                    {    out << "/" << l2 << " (l=" << lens[l] << ")";
                         i++;    }
                    if ( chrs.nonempty( ) ) out << " chrs = " << printSeq(chrs);
                    out << "\n";    }
               reports[li] = out.str( );    }

          // Print reports.

          for ( int i = 0; i < reports.isize( ); i++ )
               cout << reports[i];    }

     // Show edges that are linked to by multiple active barcodes.

     /*
     const int min_see = 5;
     cout << "\nEdges linked to by at least " << min_see << " barcodes:\n";
     vec<vec<int>> see( hb.E( ) );
     vec< pair<int,int> > places;
     BinaryReader::readFile( odir + "/10X.aligns", &places );
     for ( int i = 0; i < active.isize( ); i++ )
     {    vec<int> es;
          int b = W.bx[ active[i] ];
          for ( int64_t rid = W.kb_bci[b]; rid < W.kb_bci[b+1]; rid++ )
          {    int e = places[rid].first;
               if ( e >= 0 ) 
               {    if ( tol2[ inv[e] ].first < tol2[e].first ) e = inv[e];
                    es.push_back(e);    }    }
          UniqueSort(es);
          for ( int j = 0; j < es.isize( ); j++ )
               see[ es[j] ].push_back(b);    }
     vec<triple<int,int,int>> sees;
     for ( int e = 0; e < hb.E( ); e++ )
     {    if ( see[e].isize( ) >= min_see )
               sees.push( tol2[e].first, e, see[e].size( ) );    }
     Sort(sees);
     for ( int i = 0; i < sees.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < sees.isize( ); j++ )
               if ( sees[j].first != sees[i].first ) break;
          for ( int k = i; k < j; k++ )
          {    int l = sees[k].first, e = sees[k].second, n = sees[k].third;
               int l2 = tol2[ inv[e] ].first;
               if ( BinMember( xls, l ) | BinMember( xls, l2 ) ) continue;
               cout << e << "[L=" << l << "/" << l2 << ", bc=" << n << "]\n";    }
          i = j - 1;    }
     */

     // Reassess simple bubbles.

     cout << "\nReassessment of simple bubbles:\n" << endl;
     double reassess_clock = -WallClockTime( );
     vec< quad<int,int,int,int> > resolved_bubbles;
     vec<int> resolved_bubbles_kmers;
     vec<Bool> solved( unsolved_bubbles.size( ), False );
     for ( int i = 0; i < unsolved_bubbles.isize( ); i++ )
     {    int e1 = unsolved_bubbles[i].first, e2 = unsolved_bubbles[i].second;
          int u1 = unsolved_bubbles[i].third, u2 = unsolved_bubbles[i].fourth;
          int v = hb.ToLeft(e1);
          vec< vec< kmer<K> > > hits;
          FindHits( W, hb, v, hits, hits_clock, False );
          int u1_new = 0, u2_new = 0;
          for ( int j = 0; j < hits[0].isize( ); j++ )
               if ( !BinMember( hits[1], hits[0][j] ) ) u1_new++;
          for ( int j = 0; j < hits[1].isize( ); j++ )
               if ( !BinMember( hits[0], hits[1][j] ) ) u2_new++;

          vec<int> nbc1, nbc2;
          for ( int i = 0; i < W.ehits[e1].isize( ); i++ )
               nbc1.push_back( W.bcs[ W.ehits[e1][i] ] );
          for ( int i = 0; i < W.ehits[ inv[e1] ].isize( ); i++ )
               nbc1.push_back( W.bcs[ W.ehits[ inv[e1] ][i] ] );
          for ( int i = 0; i < W.ehits[e2].isize( ); i++ )
               nbc2.push_back( W.bcs[ W.ehits[e2][i] ] );
          for ( int i = 0; i < W.ehits[ inv[e2] ].isize( ); i++ )
               nbc2.push_back( W.bcs[ W.ehits[ inv[e2] ][i] ] );
          int supp1 = 0, supp2 = 0;
          for ( int pass = 1; pass <= 2; pass++ )
          {    vec<uint32_t> bx;
               if ( pass == 1 )
               {    for ( int j = 0; j < nbc1.isize( ); j++ )
                         bx.push_back( nbc1[j] );    }
               else
               {    for ( int j = 0; j < nbc2.isize( ); j++ )
                         bx.push_back( nbc2[j] );    }
               // cout << "using " << ( pass == 1 ? "nbc1" : "nbc2" )
               //      << endl;
               vec< vec< kmer<K> > > reads( bx.size( ) );
               for ( int i = 0; i < bx.isize( ); i++ )
               {    reads[i].ReadRange( tdir + "/10X.kb", 
                         W.kb_bci[ bx[i] ], W.kb_bci[ bx[i] + 1 ] );    }
               for ( int z = 0; z < solved_bubbles.isize( ); z++ )
               {    if ( Abs( solved_bubbles_kmers[z] - unsolved_bubbles_kmers[i] )
                         > 100000 )
                    {    continue;    }
                    int e1 = solved_bubbles[z].first, e2 = solved_bubbles[z].second;
                    int u1 = solved_bubbles[z].third, u2 = solved_bubbles[z].fourth;
                    int v = hb.ToLeft(e1);
                    vec< vec< kmer<K> > > hits(2);
                    for ( int j = 0; j < 2; j++ )
                    {    int e = hb.IFrom( v, j );
                         #pragma omp parallel for
                         for ( int l = 0; l <= hb.Bases(e) - K; l++ )
                         {    kmer<K> x;
                              x.SetToSubOf( hb.EdgeObject(e), l );
                              Bool found = False;
                              for ( int b = 0; b < bx.isize( ); b++ )
                              {    int64_t p = BinPosition( reads[b], x );
                                   if ( p >= 0 )
                                   {    
                                        #pragma omp critical
                                        {    hits[j].push_back(x);    }
                                        found = True;
                                        break;    }    }    
                              if ( !found )
                              {    kmer<K> y(x);
                                   y.ReverseComplement( );
                                   for ( int b = 0; b < bx.isize( ); b++ )
                                   {    int64_t p 
                                             = BinPosition( reads[b], y );
                                        if ( p >= 0 )
                                        {    
                                             #pragma omp critical
                                             {    hits[j].push_back(x);   }
                                             found = True;
                                             break;    }    }    }    }
                         Sort( hits[j] );    }
                    int u1_new = 0, u2_new = 0;
                    for ( int j = 0; j < hits[0].isize( ); j++ )
                         if ( !BinMember( hits[1], hits[0][j] ) ) 
                              u1_new++;
                    for ( int j = 0; j < hits[1].isize( ); j++ )
                         if ( !BinMember( hits[0], hits[1][j] ) ) 
                              u2_new++;

                    if ( ( pass == 1 && u1 > u2 )
                         || ( pass == 2 && u2 > u1 ) )
                    {    supp1 += u1_new;
                         supp2 += u2_new;    }
                    else
                    {    supp1 += u2_new;
                         supp2 += u1_new;    }    }    }

          if ( u1_new > 0 || u2_new > 0 || supp1 > 0 || supp2 > 0 )
          {    // cout << "e1/2 = " << e1 << "/" << e2 << ", u1/2 = " << u1 << "/"
               //      << u2 << ", u1/2_new = " << u1_new << "/" << u2_new << ", "
               //      << "supp1/2 = " << supp1 << "/" << supp2;
               Bool s1 = ( supp1 >= 5 * Max( supp2, 1 ) );
               Bool s2 = ( supp2 >= 5 * Max( supp1, 1 ) );
               if ( s1 || s2 ) 
               {    resolved_bubbles.push( e1, e2, supp1, supp2 );
                    resolved_bubbles_kmers.push_back( unsolved_bubbles_kmers[i] );
                    // cout << " [resolved!]";
                    for ( int i = 0; i < W.path.isize( ); i++ )
                    {    if ( W.path[i] == e1 && s2 ) W.path[i] = e2;
                         if ( W.path[i] == e2 && s1 ) W.path[i] = e1;    }
                    random_choices--;
                    good_choices++;
                    solved[i] = True;    }
               else
               {    cout << "e1/2 = " << e1 << "/" << e2 << ", u1/2 = " << u1 << "/"
                         << u2 << ", u1/2_new = " << u1_new << "/" << u2_new << ", "
                         << "supp1/2 = " << supp1 << "/" << supp2 << endl;    }
               // cout << endl;    
                    }    }
     solved_bubbles.append(resolved_bubbles);
     solved_bubbles_kmers.append(resolved_bubbles_kmers);
     EraseIf( unsolved_bubbles, solved );
     EraseIf( unsolved_bubbles_kmers, solved );
     reassess_clock += WallClockTime( );

     // Print stats.

     cout << "\n";
     PRINT2( good_choices, random_choices );
     cout << PERCENT_RATIO( 3, random_choices, random_choices + good_choices )
          << " random" << endl;
     cout << captured_gaps << " captured gaps and " << force_choices
           << " forced choices" << endl;
     int nbc2 = 0, nbc10 = 0;
     for ( int b = 0; b < W.bx.isize( ); b++ )
     {    if ( W.tracks[b].nonempty( ) && W.kmers - W.tracks[b].back( ) <= 2000 ) 
               nbc2++;
          if ( W.tracks[b].nonempty( ) && W.kmers - W.tracks[b].back( ) <= 10000 ) 
               nbc10++;    }
     cout << ToStringAddCommas(W.kmers) << " bases in path" << endl;
     cout << W.Seen(2000) << " barcodes seen within last 2000 bases" << endl;
     cout << W.Seen(10000) << " barcodes seen within last 10000 bases" << endl;

     // Compute final tracks.

     vec<vec<int>> tracksf( W.bx.isize( ) );
     if ( print_tracks || display_final_tracks )
     {    cout << Date( ) << ": computing final tracks" << endl;
          vec< vec< int64_t > > reads_ids( W.bx.size( ) );
          for ( int i = 0; i < W.bx.isize( ); i++ )
          {    reads_ids[i].ReadRange( tdir + "/10X.kb_ids", 
                    W.kb_bci[ W.bx[i] ], W.kb_bci[ W.bx[i] + 1 ] );    }
          int kmers = 0;
          for ( int i = 0; i < W.path.isize( ); i++ )
          {    int e = W.path[i];
               int start = ( i == 0 ? 0 : hb.K( ) - K );
               for ( int l = start; l <= hb.Bases(e) - K; l++ )
               {    kmer<K> x;
                    x.SetToSubOf( hb.EdgeObject(e), l );
                    kmer<K> y(x);
                    y.ReverseComplement( );
                    #pragma omp parallel for
                    for ( int b = 0; b < W.bx.isize( ); b++ )
                    {    int64_t p = BinPosition( W.reads[b], x );
                         if ( p >= 0 && count[reads_ids[b][p]] < 100 )
                              tracksf[b].push_back( kmers + l );
                         else
                         {    int64_t p = BinPosition( W.reads[b], y );
                              if ( p >= 0 && count[reads_ids[b][p]] < 100 )
                                   tracksf[b].push_back( kmers + l );    }    }    }
               kmers += hb.Kmers(e);    }    }
     if (print_tracks)
     {    cout << Date( ) << ": final tracks:" << endl;
          for ( int b = 0; b < W.bx.isize( ); b++ )
               cout << W.bx[b] << ": " << printSeq( tracksf[b] ) << endl;
          cout << endl;    }

     // Plot the tracks.  Sort and delete tracks having few points.

     if (display_final_tracks)
     {    
          // Display heuristics.

          const int min_points   = 10;
          const double min_dist  = 0.7;
          const int color_search = 3;
          const double max_color = 2.3;
          const double vert_mult = 3.1;

          // Create the display.

          Ofstream( out, "xxx.data" )
          {    vec<int> M( W.bx.isize( ), 0 );
               for ( int b = 0; b < tracksf.isize( ); b++ )
                    if ( tracksf[b].nonempty( ) ) M[b] = Median( tracksf[b] );
               SortSync( M, tracksf );
               vec<Bool> to_delete( tracksf.isize( ), False );
               for ( int b = 0; b < tracksf.isize( ); b++ )
                    if ( tracksf[b].isize( ) < min_points ) to_delete[b] = True;
               EraseIf( tracksf, to_delete );
               vec<vec<double>> c( tracksf.size( ), vec<double>(3) );
               for ( int b = 0; b < tracksf.isize( ); b++ )
               {    get_color:
                    c[b][0] = ( randomx( ) % 1000 ) / 1000.0;
                    c[b][1] = ( randomx( ) % 1000 ) / 1000.0;
                    c[b][2] = ( randomx( ) % 1000 ) / 1000.0;
                    if ( Sum(c[b]) > max_color ) goto get_color;
                    double s = 0;
                    for ( int u = 1; u <= color_search; u++ )
                    {    if ( b - u < 0 ) break;
                         for ( int j = 0; j < 3; j++ )
                              s += ( c[b][j] - c[b-u][j] ) * ( c[b][j] - c[b-u][j] );
                         if ( sqrt(s) < min_dist ) goto get_color;    }
                    for ( int j = 0; j < tracksf[b].isize( ); j++ )
                    {    out << tracksf[b][j]/1000.0 << " " << b << " " << c[b][0] 
                              << " " << c[b][1] << " " << c[b][2] 
                              << "\n";    }    }    }
          int yh = int( round( vert_mult * tracksf.size( ) ) );
          SystemSucceed( "PlotPoints IN=xxx.data "
               "OUT=/wga/dev/jaffe/BroadCRD/xxx.png COLOR_BY_POINT=True "
               "POINTSIZE=1 X_AXIS_OFFSET=0 YHEIGHT=" + ToString(yh) 
               + " TITLE=\"10X molecules across assembled region, as "
               + "function of distance in kb\"" + " TITLE_FONTSIZE=11.5 XMIN=0"
               + " X_AXIS_EXTEND=0 Y_AXIS_EXTEND=0 NH=True > /dev/null 2>&1" );    }
     
     // Done.

     cout << ConvertTime(setup_clock) << " used in setup" << endl;
     cout << ConvertTime(eval_solved_clock) << " used evaluating solved "
          << "bubbles" << endl;
     cout << ConvertTime(push_clock) << " used pushing" << endl;
     cout << ConvertTime(add_clock) << " used adding barcodes" << endl;
     cout << ConvertTime(add_clock) << " used updating tracks (inside adding) " 
          << endl;
     cout << ConvertTime(hits_clock) << " used finding hits" << endl;
     cout << ConvertTime(reassess_clock) 
          << " used reassessing (overlaps hits clock)" << endl;
     cout << TimeSince(clock) << " used in total" << endl;
     cout << "\n" << Date( ) << ": done" << endl << endl;
     Scram(0);    }

// =================================================================================

template<int K> void Push( walker10<K>& W, const HyperBasevectorX& hb, 
     const int e, double& push_clock )
{    
     push_clock -= WallClockTime( );

     // Update path.

     W.path.push_back(e);

     // Update tracks.

     for ( int l = hb.K( ) - K; l <= hb.Bases(e) - K; l++ )
     {    kmer<K> x;
          x.SetToSubOf( hb.EdgeObject(e), l );
          kmer<K> y(x);
          y.ReverseComplement( );
          #pragma omp parallel for
          for ( int b = 0; b < W.bx.isize( ); b++ )
          {    if ( W.stale[b] ) continue;
               int64_t p = BinPosition( W.reads[b], x );
               if ( p >= 0
                    // && count[reads_ids[b][p]] < 100 
                         )
               {    W.tracks[b].push_back( W.kmers + l );    }
               else
               {    int64_t p = BinPosition( W.reads[b], y );
                    if ( p >= 0 
                         // && count[reads_ids[b][p]] < 100 
                              )
                    {    W.tracks[b].push_back( W.kmers + l );    }    }    }    }

     // Update kmers.

     int n = 0;
     if ( hb.Bases(e) > 0 ) n = hb.Kmers(e);
     cout << W.path.size( ) << "(" << W.kmers << "). ";    
     W.kmers += n;    

     // Mark stale barcodes.

     const int max_gap = 50000;
     for ( int b = 0; b < W.bx.isize( ); b++ )
     {    if ( W.tracks[b].empty( ) || W.kmers - W.tracks[b].back( ) > max_gap )
               W.stale[b] = True;    }
     push_clock += WallClockTime( );    }

template<int K> void AddBarcodes( walker10<K>& W, const int x, 
     const HyperBasevectorX& hb, const vec<int>& inv, const String& tdir, 
     double& add_clock, double& tracks_clock )
{    double clock = WallClockTime( );
     add_clock -= WallClockTime( );
     vec<int> nbc;
     for ( int i = 0; i < W.ehits[x].isize( ); i++ )
     {    int b = W.bcs[ W.ehits[x][i] ];
          if ( BinMember( W.bx, b ) ) W.stale[ BinPosition( W.bx, b ) ] = False;
          else nbc.push_back(b);    }
     for ( int i = 0; i < W.ehits[ inv[x] ].isize( ); i++ )
     {    int b = W.bcs[ W.ehits[ inv[x] ][i] ];
          if ( BinMember( W.bx, b ) ) W.stale[ BinPosition( W.bx, b ) ] = False;
          else nbc.push_back(b);    }
     if ( nbc.empty( ) ) 
     {    add_clock += WallClockTime( );
          return;    }
     UniqueSort(nbc);
     vec<uint32_t> bx_orig(W.bx);
     W.bx.append(nbc);
     W.reads.resize( W.bx.size( ) );
     W.tracks.resize( W.bx.size( ) );
     W.stale.resize( W.bx.size( ), False );
     SortSync( W.bx, W.reads, W.tracks, W.stale );

     // Load reads.
     
     for ( int i = 0; i < W.bx.isize( ); i++ )
     {    if ( W.stale[i] ) continue;
          if ( BinMember( bx_orig, W.bx[i] ) ) continue;
          W.reads[i].ReadRange( 
               tdir + "/10X.kb", W.kb_bci[ W.bx[i] ], W.kb_bci[ W.bx[i] + 1 ] );    }
     cout << "Added " << nbc.size( ) << " barcodes (total non-stale barcodes = " 
          << W.bx.size( ) - Sum(W.stale) << "), time used = " 
          << TimeSince(clock) << "." << endl;    
     // cout << "added barcodes = " << printSeq(nbc) << endl;

     // Update tracks.

     tracks_clock -= WallClockTime( );
     vec<uint32_t> bxu;
     for ( int b = 0; b < W.bx.isize( ); b++ )
     {    if ( W.stale[b] ) continue;
          if ( BinMember( bx_orig, W.bx[b] ) ) continue;
          bxu.push_back(b);    }
     int np = W.path.size( );
     vec< vec< kmer<K> > > xs(np), ys(np);

     const int kmers_back = 50000;
     int total_kmers = 0;
     for ( int i = 0; i < np; i++ )
          total_kmers += hb.Kmers( W.path[i] );
     int pi, kmers = 0;
     for ( pi = 0; pi < np; pi++ )
     {    int nk = hb.Kmers( W.path[pi] );
          if ( total_kmers - kmers - nk <= kmers_back ) break;
          kmers += nk;    }

     #pragma omp parallel for
     for ( int i = pi; i < np; i++ )
     {    int e = W.path[i];
          if ( hb.Bases(e) == 0 ) continue;
          int start = ( i == 0 ? 0 : hb.K( ) - K ), stop = hb.Bases(e) - K;
          xs[i].resize( stop - start + 1 ), ys[i].resize( stop - start + 1 );
          for ( int l = start; l <= stop; l++ )
          {    xs[i][l-start].SetToSubOf( hb.EdgeObject(e), l );
               ys[i][l-start] = xs[i][l-start];
               ys[i][l-start].ReverseComplement( );    }    }
     #pragma omp parallel for
     for ( int ib = 0; ib < bxu.isize( ); ib++ )
     {    int b = bxu[ib], kmers = 0;
          for ( int i = 0; i < pi; i++ )
               kmers += hb.Kmers( W.path[i] );
          for ( int i = pi; i < np; i++ )
          {    int e = W.path[i];
               if ( hb.Bases(e) == 0 ) continue;
               int start = ( i == 0 ? 0 : hb.K( ) - K );
               for ( int l = start; l <= hb.Bases(e) - K; l++ )
               {    const kmer<K> &x = xs[i][l-start], &y = ys[i][l-start];
                    int64_t p = BinPosition( W.reads[b], x );
                    if ( p >= 0
                         /* && count[reads_ids[b][p]] < 100 */ )
                    {    W.tracks[b].push_back( kmers + l );    }
                    else
                    {    int64_t p = BinPosition( W.reads[b], y );
                         if ( p >= 0
                              /* && count[reads_ids[b][p]] < 100 */ )
                         {    W.tracks[b].push_back( kmers + l );    }    }    }
               kmers += hb.Kmers(e);    }    }
     tracks_clock += WallClockTime( );
     add_clock += WallClockTime( );    }

template<int K> void FindHits( walker10<K>& W, const HyperBasevectorX& hb, 
     const int v, vec< vec< kmer<K> > >& hits, double& hits_clock, 
     const Bool use_stale )
{
     hits_clock -= WallClockTime( );
     hits.clear( );
     hits.resize(2);
     for ( int j = 0; j < 2; j++ )
     {    int e = hb.IFrom( v, j );
          #pragma omp parallel for
          for ( int l = 0; l <= hb.Bases(e) - K; l++ )
          {    kmer<K> x;
               x.SetToSubOf( hb.EdgeObject(e), l );
               Bool found = False;
               for ( int b = 0; b < W.bx.isize( ); b++ )
               {    if ( use_stale && W.stale[b] ) continue;
                    int64_t p = BinPosition( W.reads[b], x );
                    if ( p >= 0 )
                    {    
                         #pragma omp critical
                         {    hits[j].push_back(x);    }
                         found = True;
                         break;    }    }    
               if ( !found )
               {    kmer<K> y(x);
                    y.ReverseComplement( );
                    for ( int b = 0; b < W.bx.isize( ); b++ )
                    {    if ( use_stale && W.stale[b] ) continue;
                         int64_t p = BinPosition( W.reads[b], y );
                         if ( p >= 0 )
                         {    
                              #pragma omp critical
                              {    hits[j].push_back(x);    }
                              found = True;
                              break;    }    }    }    }
          Sort( hits[j] );    }
     hits_clock += WallClockTime( );    }

void ClassifyDiffs( const align& a, const basevector& E1, const basevector& E2 )
{    cout << "differences: ";
     int subs = 0;
     Bool first = True;
     int p1 = a.pos1( ), p2 = a.pos2( );
     for ( int j = 0; j < a.Nblocks( ); j++ ) 
     {    if ( a.Gaps(j) > 0 )  
          {    if ( !first ) cout << ", ";
               first = False;
               Bool same = True;
               int n = a.Gaps(j);
               for ( int l = 1; l < n; l++ )
                    if ( E2[p2+l] != E2[p2+l-1] ) same = False;
               int add = 0;
               if (same)
               {    for ( int l = -1; ; l-- )
                    {    if ( p2+l < 0 ) break;
                         if ( E2[p2+l] != E2[p2] ) break;
                         add++;    }
                    for ( int l = n; ; l++ )
                    {    if ( p2+l >= E2.isize( ) ) break;
                         if ( E2[p2+l] != E2[p2] ) break;
                         add++;    }    }
               if ( same && add > 0 )
               {    char base = as_base( E2[p2] );
                    cout << "homopolymer indel " << base << "^" 
                         << add << " vs " << base << "^" << add + n;    }
               else cout << "indel of size " << n << endl;
               p2 += n;    }
          if ( a.Gaps(j) < 0 ) 
          {    if ( !first ) cout << ", ";
               first = False;
               Bool same = True;
               int n = -a.Gaps(j);
               for ( int l = 1; l < n; l++ )
                    if ( E1[p1+l] != E1[p1+l-1] ) same = False;
               int add = 0;
               if (same)
               {    for ( int l = -1; ; l-- )
                    {    if ( p1+l < 0 ) break;
                         if ( E1[p1+l] != E2[p1] ) break;
                         add++;    }
                    for ( int l = n; ; l++ )
                    {    if ( p1+l >= E1.isize( ) ) break;
                         if ( E1[p1+l] != E2[p1] ) break;
                         add++;    }    }
               if ( same && add > 0 )
               {    char base = as_base( E1[p1] );
                    cout << "homopolymer indel " << base << "^" 
                         << add << " vs " << base << "^" << add + n;    }
               else cout << "indel of size " << n << endl;
               p1 += n;    }
          for ( int x = 0; x < a.Lengths(j); x++ ) 
          {    if ( E1[p1] != E2[p2] ) subs++;
               ++p1;
               ++p2;    }    }
     if ( subs > 0 && !first ) cout << ", ";
     if ( subs == 1 ) cout << "substitution";
     if ( subs > 1 ) cout << subs << " substitutions";
     cout << endl;    }
