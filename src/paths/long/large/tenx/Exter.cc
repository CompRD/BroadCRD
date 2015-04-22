///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Experimental code to use 10X data to extend lines in a DISCOVAR de novo assembly.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "FastIfstream.h"
#include "MainTools.h"
#include "feudal/PQVec.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "paths/HyperBasevector.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/large/Lines.h"
#include "paths/long/large/tenx/TenxDirs.h"
#include "paths/long/large/tenx/TenxTools.h"

int main( )
{    RunTime( );

     // Control.

     const Bool verbose1 = False;
     const Bool verbose2 = False;
     const Bool show_gaps = False;
     Bool show_nearby_lines = False;
     const Bool avoid_ends = False;

     // Define seed and ancillary data.

     // 3         0      10,000 **REF GAP**
     // 3    11,998     261,734 line[5976]   2999427..6660799   len=243511  cov=1.00
     // 3   263,028     360,600 line[18587]  6258032..303622    len=98184   cov=1.02
     // ----------------------------------------------------------------------------
     // 3   361,625     573,863 line[7535]   75344..1359938     len=212242  cov=1.02
     // 3   574,547     577,285 line[65599]  1018410..887520    len=2539    cov=1.18
     // 3   577,183     580,287 line[62902]  9670913..9670913   len=2905    cov=1.15
     // 3   580,934     581,895 line[97954]  364326..399247     len=1063    cov=1.05
     // 3   582,095     583,430 line[94129]  9671027..9671027   len=1136    cov=1.26
     // 3   584,910     588,383 line[60529]  8646111..8646111   len=3274    cov=1.07
     // 3   588,357     592,591 line[56925]  9839651..9839651   len=4035    cov=1.21
     // 3   592,392     595,884 line[60409]  9871553..9871553   len=3293    cov=1.12
     // 3   597,609     600,552 line[64021]  9847773..9847773   len=2744    cov=1.17
     // 3   600,402     607,147 line[50936]  9800311..9800311   len=6546    cov=1.10
     // 3   608,222     611,745 line[60248]  9726491..7936805   len=3325    cov=1.09
     // 3   611,573     697,983 line[20555]  9701140..8794013   len=86269   cov=1.05
     //
     // Basically missing from PNAS, just found this:
     //
     // scaffold 1111 (l = 558244):
     // ...
     // -- (48 +/- 12) --> 50887 (l = 2608)
     // -- (862 +/- 74) --> 28120 (l = 1077)
     // scaffold 2882 (l = 6413):
     // 10942 (l = 1291)
     // -- (1050 +/- 35) --> 14784 (l = 1371)
     // -- (998 +/- 29) --> 286011 (l = 1703)

     // #0
     // UNRESOLVED FROM L18587 to L7535
     // int seed = 393288;
     // const int depth = 25;
     // const int len_goal = 5000;
     // answer should be (?)
     // 393288,{7951594,346614},303622,263024,223593,185530,149944 
     // 
     // FIRST ALTERNATIVE
     // 95473,73114,73113,95472 .. 128934,128935,127252,127251,126201,7344438
     // single read says intermediate is:
     // 8611123,107781,137870,137869,171299,7919370,117906,7341063,8616543..
     // 10593530,97677,74936,74935,97676,119755,7917739,163165,129712,129713,
     // 7894077,185119,8613701,8795921,236549,7346337,197905,7927561,162175
     // 
     //
     // (7344438,10531545)
     // ... 75344

     // #1
     // int seed = 1359938;
     // const int depth = 25; // keeps going
     // const int len_goal = 6000;
     // [1] 1359938,1290221,1221523,1153257,1085577,1018410,952800,
     //     887520,822781,758540,694659,9670913
     // gaps >= 100:
     // 244 378
     // 459 566
     // 7088 7188
     // 7269 7642
     // missing edges:

     // #2
     int seed = 9670913;
     const int depth = 15; // reduced; keeps going
     const int len_goal = 6248; // note rigged
     show_nearby_lines = True;
     // [1,2] 9670913,271178,243029,204278,9658733,280493,9008805,364326,
     //       {412782,412783},399247,352079,352078,9671027
     // missing edges:

     // #3
     // int seed = 9671027;
     // const int depth = 25; // keeps going
     // const int len_goal = 6000;
     // [1] 9671027,565977,7355365,9713507,1148615,1148616,1148617,1216896,
     //     9155886,1425082,1494543,8646111
     // missing edges: 1148616 1148617 1425082 1494543
     // [#1, 32 barcodes, 116 reads] L94130/94129 (l=1136)
     // [#2, 22 barcodes, 2372 reads] L7536/7535 (l=212242)
     // [#3, 21 barcodes, 182 reads] L62903/62902 (l=2905)
     // [#4, 21 barcodes, 270 reads] L56926/56925 (l=4035)
     // [#5, 19 barcodes, 164 reads] L60530/60529 (l=3274)
     // [#6, 16 barcodes, 66 reads] L77403/77402 (l=1658)
     // [#7, 16 barcodes, 162 reads] L65599/65598 (l=2539)
     // [#8, 16 barcodes, 140 reads] L64022/64021 (l=2744)
     // [#9, 15 barcodes, 250 reads] L50937/50936 (l=6546)
     // [#10, 14 barcodes, 90 reads] L60410/60409 (l=3293)

     // #4
     // int seed = 8646111;
     // const int depth = 25;
     // const int len_goal = 6000;
     // [1] 8646111,9175857,9839651
     // missing edges:
     // [#1, 123 barcodes, 900 reads] L60530/60529 (l=3274)
     // [#2, 78 barcodes, 766 reads] L56926/56925 (l=4035)
     // [#3, 65 barcodes, 904 reads] L50937/50936 (l=6546)
     // [#4, 63 barcodes, 5326 reads] L7536/7535 (l=212242)
     // [#5, 54 barcodes, 452 reads] L64022/64021 (l=2744)
     // [#6, 52 barcodes, 486 reads] L65599/65598 (l=2539)
     // [#7, 48 barcodes, 336 reads] L60410/60409 (l=3293)
     // [#8, 45 barcodes, 356 reads] L62903/62902 (l=2905)
     // [#9, 45 barcodes, 2656 reads] L20555/20554 (l=86269)
     // [#10, 43 barcodes, 202 reads] L77403/77402 (l=1658)

     // #5
     // int seed = 9839651;
     // const int depth = 25;
     // const int len_goal = 6000;
     // [1] 9839651,9871553
     // missing edges:

     // #6
     // int seed = 9871553;
     // const int depth = 25;
     // const int len_goal = 6000;
     // [1] 9871553,2185248,2185247,2129311,2072233,2014225,1954990,
     //     1954989,1954991,9847773
     // missing edges:

     // #7
     // int seed = 9847773;
     // const int depth = 25; // keeps going
     // const int len_goal = 6000;
     // [1] 9847773,1581145,1514245,9800311
     // missing edges:

     // #8
     // int seed = 9847773;
     // const int depth = 25; // keeps going
     // const int len_goal = 6000;
     // [1] 9847773,1581145,1514245,9800311
     // missing edges:

     // #9
     // int seed = 9800311;
     // const int depth = 25;
     // const int len_goal = 10000; // note larger
     // [1] 9800311,1100950,1100949,9137367,7417579,902416,837799,773727,9726491
     // missing edges:

     // #10
     // int seed = 7936805;
     // const int depth = 15; // note smaller, keeps going
     // const int len_goal = 6000;
     // [1] 7936805,7931796,8623688,10596564,156826,8619778,8618206,143078,9701140
     // missing edges:

     // #11 = 117 kb edge = one chromosomal locus, not in GRCh38 = voom.13
     // hard to sort out
     // int seed = 10022646;
     // const int depth = 25;
     // const int len_goal = 6000;
     // show_nearby_lines = True;

     // Voom: 42, 44, 46, 48, 50, 52
     // Voom.42 = L15640
     // Voom.44 = L17076

     // #12 = Voom.44 = L17076, right edge
     // int seed = 9368699;
     // int seed = 9925309; // L55506
     // int seed = 9283658; // L56613
     // const int depth = 25;
     // const int len_goal = 6000;
     // show_nearby_lines = True;
     // line links plus graph seem to imply order
     // L17076,L64892 x 3,L49184,L54506,L55506,L73206,L56613,gap,L64026

     // #13 = Voom.46 = L17670
     // int seed = 9187203;
     // int seed = 8142111; // L51817
     // const int depth = 50;
     // const int len_goal = 12000;
     // show_nearby_lines = True;
     // order L17670,L51817,L58966
     // from L51817 with these parameters get 14 paths, takes 8.5 minutes
     // (so working but very slow)

     // #14 = Voom.50 = L19020
     // int seed = 9012225;
     // const int depth = 25;
     // const int len_goal = 6000;
     // show_nearby_lines = True;
     // not a bad illustration, 8 solutions, continues on through L41723
     // presumably three het sites make 8 solutions

     // Hardcoded directories.

     int N = 1;
     String dir, odir, tdir;
     SetTenxDirs( N, dir, odir, tdir );

     // Load counts and assembly.

     cout << "\n" << Date( ) << ": loading counts and assembly" << endl;
     HyperBasevectorX hb;
     vec<int> count, inv, tol, lens;
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
     GetTol( hb, lines, tol );
     GetLineLengths( hb, lines, lens );
     vec<String> genome_names;
     fast_ifstream in( dir + "/../genome.names" );
     String line;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          genome_names.push_back(line);    }

     // Find extensions.

     cout << Date( ) << ": finding extensions" << endl;
     vec<vec<int>> paths;
     vec<int> p = {seed};
     paths.push_back(p);
     Bool stable = False;
     for ( int d = 1; d <= depth; d++ )
     {    vec<vec<int>> paths2;
          for ( int i = 0; i < paths.isize( ); i++ )
          {    if ( hb.Cat( paths[i] ).isize( ) >= len_goal )
                    paths2.push_back( paths[i] );
               else
               {    int v = hb.ToRight( paths[i].back( ) );
                    for ( int j = 0; j < (int) hb.From(v).size( ); j++ )
                    {    vec<int> p = paths[i];
                         p.push_back( hb.IFrom( v, j ) );
                         paths2.push_back(p);    }    }    }
          if ( paths == paths2 )
          {    stable = True;
               break;    }
          paths = paths2;    }
     vec<vec<int>> pathsx;
     for ( int i = 0; i < paths.isize( ); i++ )
     {    if ( hb.Cat( paths[i] ).isize( ) >= len_goal )
               pathsx.push_back( paths[i] );    }
     paths = pathsx;
     cout << "found " << paths.size( ) << " paths, ";
     if (stable) cout << "stable" << endl;
     else cout << "unstable" << endl;
     vec<basevector> bpaths;
     for ( int i = 0; i < paths.isize( ); i++ )
          bpaths.push_back( hb.Cat( paths[i] ) );

     // Load alignments and barcodes.

     cout << Date( ) << ": loading alignments" << endl;
     vec< pair<int,int> > places;
     BinaryReader::readFile( odir + "/10X.aligns", &places );
     cout << Date( ) << ": loading barcodes" << endl;
     vec<uint32_t> bcs;
     BinaryReader::readFile( tdir + "/10X.bc", &bcs );
     vec<int64_t> bci;
     BinaryReader::readFile( tdir + "/10X.bci", &bci );
     int nbc_total = bci.size( ) - 1;
     cout << "total barcodes = " << ToStringAddCommas(nbc_total) << endl;
     int64_t total_bases = places.size( ) * 88;
     cout << "total bases = " << ToStringAddCommas(total_bases) << endl;
     cout << Date( ) << ": loading line hits" << endl;
     vec<int> lhits;
     BinaryReader::readFile( odir + "/10X.line_hits", &lhits );
     vec<double> bc_frac;
     BinaryReader::readFile( odir + "/10X.bc_frac", &bc_frac );

     // Identify the barcodes.

     vec<uint32_t> bx;
     cout << Date( ) << ": identifying barcodes" << endl;
     for ( int64_t id = 0; id < places.jsize( ); id++ )
     {    int e = places[id].first, epos = places[id].second;
          if ( e == seed || e == inv[seed] ) 
          {    if ( hb.Bases(e) - epos > 20000 ) continue;
               // hardcoded!!
               if ( !avoid_ends || ( epos >= hb.K( ) 
                    && epos <= hb.EdgeObject(e).isize( ) - hb.K( ) - 88 ) )
               {    bx.push_back( bcs[id] );    }    }    }
     UniqueSort(bx);
     cout << "using " << bx.size( ) << " barcodes" << endl;

     // Identify the reads.

     cout << Date( ) << ": identifying reads" << endl;
     vec<vec<int64_t>> hits( bx.size( ) );
     #pragma omp parallel for
     for ( int64_t id = 0; id < places.jsize( ); id++ )
     {    int p = BinPosition( bx, bcs[id] );
          if ( p < 0 ) continue;
          #pragma omp critical
          {    hits[p].push_back(id);    }    }
     for ( int i = 0; i < bx.isize( ); i++ )
          Sort( hits[i] );

     // Load the reads.

     cout << Date( ) << ": loading reads" << endl;
     vec<int64_t> all;
     vec<int> alli;
     for ( int i = 0; i < bx.isize( ); i++ )
     for ( int j = 0; j < hits[i].isize( ); j++ )
     {    int64_t pid = hits[i][j]/2;
          // BUG: paired reads are getting pushed twice ****************************
          all.push_back( 2*pid, 2*pid+1 );
          alli.push_back(i);    }
     cout << "total reads in barcode set = " << ToStringAddCommas( all.size( ) )
          << endl;
     vecbasevector bases;
     VecPQVec quals;
     bases.Read( tdir + "/10X.fastb", all );
     quals.Read( tdir + "/10X.qualp", all );
     if (verbose1)
     {    for ( int64_t id = 0; id < places.jsize( ); id++ )
          {    int e = places[id].first, epos = places[id].second;
               if ( e == seed || e == inv[seed] ) 
               {    // hardcoded!!
                    if ( !avoid_ends || ( epos >= hb.K( ) 
                         && epos <= hb.EdgeObject(e).isize( ) - hb.K( ) - 88 ) )
                    {    int m = BinPosition( all, id );
                         int p = BinPosition( bx, bcs[id] );
                         cout << BaseAlpha(p) << " from " << bases[m].ToString( ) 
                              << endl;    }    }    }    }

     // Show lines having most hits.

     if (show_nearby_lines)
     {    FindTenxNhood( hb, inv, lines, lens, places, lhits, bc_frac,
               all, bcs, bx.size( ), nbc_total, genome_names, aligns );    }

     // Align the reads to the paths.

     cout << Date( ) << ": aligning the reads to the paths" << endl;
     vecbasevector stuff(bases);
     stuff.Append(bases);
     for ( int j = bases.size( ); j < (int) stuff.size( ); j++ )
          stuff[j].ReverseComplement( );
     for ( int i = 0; i < bpaths.isize( ); i++ )
          stuff.push_back( bpaths[i] );
     const int K = 88;
     vec< triple<kmer<K>,int,int> > kmers_plus;
     MakeKmerLookup0( stuff, kmers_plus );
     vec<vec< pair<int,int> >> sets( paths.size( ) );
     for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
     {    int64_t j, m;
          for ( j = i + 1; j < kmers_plus.jsize( ); j++ )
               if ( kmers_plus[j].first != kmers_plus[i].first ) break;
          for ( m = i; m < j; m++ )
               if ( kmers_plus[m].second >= 2 * (int) bases.size( ) ) break;
          for ( int64_t k1 = i; k1 < m; k1++ )
          for ( int64_t k2 = m; k2 < j; k2++ )
          {    int rid = kmers_plus[k1].second;
               Bool rc = False;
               if ( rid >= (int) bases.size( ) )
               {    rid -= bases.size( );
                    rc = True;    }
               int pid = kmers_plus[k2].second - 2 * (int) bases.size( ); 
               int pos = kmers_plus[k2].third;
               // cout << "read " << ( !rc ? "+" : "-" ) << rid << " placed at "
               //      << pid << "." << pos << endl;    
               sets[pid].push( pos, rid );    }
          i = j - 1;    }

     // Study the read sets.

     cout << Date( ) << ": studying the read sets" << endl;
     for ( int i = 0; i < paths.isize( ); i++ )
          Sort( sets[i] );
     vec<Bool> beaten( paths.size( ), False );
     #pragma omp parallel for
     for ( int i1 = 0; i1 < paths.isize( ); i1++ )
     for ( int i2 = i1+1; i2 < paths.isize( ); i2++ )
     {    ostringstream* out = NULL;
          if (verbose2) out = new ostringstream;
          if (verbose2)
          {    *out << "\ncomparing path " << i1 << "[l=" << bpaths[i1].size( ) 
                    << "]" << " to path " << i2 << "[l=" << bpaths[i2].size( ) 
                    << "]" << endl;
               *out << "p1 = " << printSeq( paths[i1] ) << endl;
               *out << "p2 = " << printSeq( paths[i2] ) << endl;    }
          const vec<pair<int,int>> &s1 = sets[i1], &s2 = sets[i2];
          vec<int> bc1, bc2;
          int n1 = 0, n2 = 0;
          for ( int j1 = 0; j1 < s1.isize( ); j1++ )
          {    int rpos1 = s1[j1].first, rid1 = s1[j1].second;
               if ( rpos1 > bpaths[i2].isize( ) - K ) break;
               if ( count[ all[rid1] ] > 100 ) continue;
               Bool found = False;
               for ( int j2 = 0; j2 < s2.isize( ); j2++ )
               {    if ( s2[j2].second == rid1 )
                    {    found = True;
                         break;    }    }
               if (found) continue;
               if (verbose2)
               {    *out << BaseAlpha( BinPosition( bx, bcs[ all[rid1] ] ) ) << ".";
                    *out << rid1 << "." << rpos1 << " unique to first path" << endl;
                    *out << bases[rid1].ToString( ) << endl;    }
               if ( rpos1 <= bpaths[i2].isize( ) - K - 200 ) 
               {    n1++;
                    bc1.push_back( bcs[ all[rid1] ] );    }    }
          for ( int j2 = 0; j2 < s2.isize( ); j2++ )
          {    int rpos2 = s2[j2].first, rid2 = s2[j2].second;
               if ( rpos2 > bpaths[i1].isize( ) - K ) break;
               if ( count[ all[rid2] ] > 100 ) continue;
               Bool found = False;
               for ( int j1 = 0; j1 < s1.isize( ); j1++ )
               {    if ( s1[j1].second == rid2 )
                    {    found = True;
                         break;    }    }
               if (found) continue;
               if (verbose2)
               {    *out << BaseAlpha( BinPosition( bx, bcs[ all[rid2] ] ) ) << ".";
                    *out << rid2 << "." << rpos2 << " unique to second path" << endl;
                    *out << bases[rid2].ToString( ) << endl;    }
               if ( rpos2 <= bpaths[i1].isize( ) - K - 200 ) 
               {    n2++;
                    bc2.push_back( bcs[ all[rid2] ] );    }    }
          UniqueSort(bc1), UniqueSort(bc2);
          n1 = bc1.size( ), n2 = bc2.size( );
          if ( n1 >= 4 && ( n2 == 0 || n1 > 10 * n2 ) )
          {    if (verbose2) *out << i1 << " beats " << i2 << endl;
               beaten[i2] = True;    }
          if ( n2 >= 4 && ( n1 == 0 || n2 > 10 * n1 ) )
          {    if (verbose2) *out << i2 << " beats " << i1 << endl;
               beaten[i1] = True;    }
          if (verbose2)
          {
               #pragma omp critical
               {    cout << out->str( );    }    
               delete out;    }    }
     cout << "\nunbeaten paths:\n";
     int pcount = 0;
     int survivors = paths.isize( ) - Sum(beaten);
     for ( int i = 0; i < paths.isize( ); i++ )
     {    if ( beaten[i] ) continue;
          const vec<int>& p = paths[i];
          cout << "[" << ++pcount << "=" << i << "] " << printSeq(p) << endl;    
          if ( survivors <= 2 )
          {    int len = hb.Cat(p).size( );
               cout << "len = " << len << endl;
               if (show_gaps)
               {    cout << "gaps >= 100:\n";
                    if ( sets[i][0].first >= 100 ) 
                         cout << "* " << sets[i][0].first << endl;
                    for ( int j = 1; j < sets[i].isize( ); j++ )
                    {    if ( sets[i][j].first - sets[i][j-1].first >= 100 )
                         {    cout << sets[i][j-1].first << " " << sets[i][j].first
                                   << endl;    }    }
                    if ( len - sets[i].back( ).first >= 100 )
                         cout << sets[i].back( ).first << " *" << endl;    }
               cout << "missing edges:";
               int start = 0;
               vec<ho_interval> epos;
               for ( int j = 0; j < p.isize( ); j++ )
               {    epos.push( start, start + hb.Bases( p[j] ) );
                    start += hb.Kmers( p[j] );    }
               for ( int j = 0; j < p.isize( ); j++ )
               {    Bool found = False;
                    for ( int l = 0; l < sets[i].isize( ); l++ )
                    {    ho_interval h( sets[i][l].first, sets[i][l].first + 88 );
                         if ( Subset( h, epos[j] ) ) 
                         {    found = True;
                              break;    }    }
                    if ( !found ) cout << " " << p[j];    }
               cout << endl;    }    }

     // Done.

     cout << Date( ) << ": done\n" << endl;
     Scram(0);    }
