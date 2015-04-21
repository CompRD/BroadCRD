///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Voom.  Take human assemblies and look for 10 kb stretches that are devoid
// of 400 base perfect matches to the reference.  Also check secondary paths.
//
// Related to NewHuman.  This is different in that it can find stretches within
// a line.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "FastIfstream.h"
#include "FeudalMimic.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/KmerBaseBroker.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/LongReadsToPaths.h"
#include "paths/long/fosmid/Ebv.cc"
#include "paths/long/large/Lines.h"
#include "paths/long/MakeKmerStuff.h"

template<int KG> void MarkReferenceKmers( const vecbasevector& x, 
     const vec< triple<kmer<KG>,int,int> >& genome_kmers_plus, vecbitvector& ref )
{    Mimic( x, ref );
     for ( int i = 0; i < (int) x.size( ); i++ )
     for ( int j = 0; j <= x[i].isize( ) - KG; j++ )
     {    kmer<KG> y;
          y.SetToSubOf( x[i], j );
          int64_t low = LowerBound1(genome_kmers_plus, y);
          int64_t high = UpperBound1(genome_kmers_plus, y);
          if ( low < high ) ref[i].Set( j, True );
          else
          {    y.ReverseComplement( );
               int64_t low = LowerBound1(genome_kmers_plus, y);
               int64_t high = UpperBound1(genome_kmers_plus, y);
               if ( low < high ) ref[i].Set( j, True );    }    }    }

void ExtractReferenceFree( const vecbasevector& x, const vecbitvector& ref, 
     const int KG, const int min_len, vec<vec<basevector>>& novel )
{    novel.clear( );
     for ( int i = 0; i < (int) x.size( ); i++ )
     for ( int j = 0; j <= x[i].isize( ) - KG; j++ )
     {    if ( ref[i][j] ) continue;
          int start_i = i, start_j = j;
          int stop_i = start_i, stop_j = start_j + KG;
          while(1)
          {    if ( ref[i][j] ) break;
               stop_i = i, stop_j = j + KG;
               if ( j < x[i].isize( ) - KG ) j++;
               else 
               {    j = 0;
                    i++;
                    if ( i == (int) x.size( ) ) break;    }    }
          vec<basevector> n;
          if ( start_i == stop_i )
          {    basevector b( x[start_i], start_j, stop_j - start_j );
               n.push_back(b);    }
          else
          {    basevector b1( x[start_i], start_j, x[start_i].isize( ) - start_j );
               n.push_back(b1);
               for ( int l = start_i + 1; l < stop_i; l++ )
                    n.push_back( x[l] );
               basevector b2( x[stop_i], 0, stop_j );
               n.push_back(b2);    }
          int nsize = 0;
          for ( int l = 0; l < n.isize( ); l++ )
               nsize += n[l].size( );
          if ( nsize >= min_len ) novel.push_back(n);
          if ( i == (int) x.size( ) ) break;    }    }

template<int KE> Bool IsEbv( const vec<vec<basevector>>& novel, 
     const vec< triple<kmer<KE>,int,int> >& ebv_kmers_plus )
{    for ( int u = 0; u < novel.isize( ); u++ )
     for ( int i = 0; i < novel[u].isize( ); i++ )
     for ( int j = 0; j <= novel[u][i].isize( ) - KE; j++ )
     {    kmer<KE> x;
          x.SetToSubOf( novel[u][i], j );
          int64_t low = LowerBound1( ebv_kmers_plus, x );
          int64_t high = UpperBound1( ebv_kmers_plus, x );
          if ( low < high ) return True;    }
     return False;    }

int main( )
{    RunTime( );
     double clock = WallClockTime( );

     // Controls.

     int nsamples = 25;
     // int nsamples = 0; // default, to use all samples
     String OUT_INSTANCE = "new_human_exp";
     // String OUT_INSTANCE = "new_human"; // default




     // String dir = "/wga/scr4/jaffe/GapToy/51400.newchem/a.final";
     // String dir = "/wga/scr4/jaffe/GapToy/51476.F3/a.final";
     // String dir = "/wga/scr4/jaffe/GapToy/51472.F2/a.final";
     // String dir = "/wga/scr4/jaffe/GapToy/50858.HCC1143+BL/a.final";

     String root = "/wga/scr4/jaffe/GapToy";

     // Define samples.

     // vec<String> I = { "51400.newchem", "51476.F3" };

     vec<String> I = { "51400.newchem", "humans/HG00096", "humans/HG00268",
          "humans/HG00419", "humans/HG00759", "humans/HG01051", "humans/HG01112",
          "humans/HG01500", "humans/HG01565", "humans/HG01583", "humans/HG01595",
          "humans/HG01879", "humans/HG02568", "humans/HG02922", "humans/HG03006",
          "humans/HG03052", "humans/HG03642", "humans/HG03742", "humans/NA18525",
          "humans/NA18939", "humans/NA19017", "humans/NA19625", "humans/NA19648",
          "humans/NA20502", "humans/NA20845" };
     if ( nsamples > 0 ) I.resize(nsamples);
     int ns = I.size( );

     for ( int i = 0; i < ns; i++ )
          ForceAssert( IsDirectory( root + "/" + I[i] ) );

     // Assembly data structures.

     vec<HyperBasevector> HB(ns);
     vec<vec<vec<vec<vec<int>>>>> LINES(ns);
     vec<vec<int>> INV(ns);
     vec<vec<vec<covcount>>> COVS(ns);
     vec<vec<int>> COUNTS(ns);
     vec< vec< vec< pair<int,int> > > > HITS(ns);
     vec<vec<int>> LENS(ns);
     vec<vec<String>> SUBSAM_NAMES(ns);

     // Load assemblies.

     cout << Date( ) << ": loading assemblies" << endl;
     for ( int s = 0; s < I.isize( ); s++ )
     {    String dir = root + "/" + I[s] + "/a.final";
          BinaryReader::readFile( dir + "/a.lines", &LINES[s] );
          BinaryReader::readFile( dir + "/a.hbv", &HB[s] );
          BinaryReader::readFile( dir + "/a.inv", &INV[s] );
          BinaryReader::readFile( dir + "/a.aligns", &HITS[s] );
          BinaryReader::readFile( dir + "/a.covs", &COVS[s] );
          if ( IsRegularFile( dir + "/../subsam.names" ) )
               BinaryReader::readFile( dir + "/../subsam.names", &SUBSAM_NAMES[s] );
          else SUBSAM_NAMES[s].resize( 1, "C" );
          COUNTS[s].resize( HB[s].E( ) );
          fast_ifstream in( dir + "/a.counts" );
          String line;
          getline( in, line );
          for ( int e = 0; e < HB[s].E( ); e++ )
          {    getline( in, line );
               COUNTS[s][e] = line.After( " " ).Int( );    }
          GetLineLengths( HB[s], LINES[s], LENS[s] );    }

     // Get sample names.

     vec<String> subsam_names;
     for ( int s = 0; s < I.isize( ); s++ )
     {    String head;
          if ( I[s].Contains( ".newchem" ) )
               head = "NA12878";
          else if ( I[s].Contains( "humans/" ) ) head = I[s].After( "humans/" );
          else
          {    cout << "Confused by subsample names." << endl;
               Scram(1);    }
          subsam_names.push_back(head);    }

     // Hash EBV.

     cout << Date( ) << ": hashing EBV" << endl;
     const int KE = 60;
     vec< triple<kmer<KE>,int,int> > ebv_kmers_plus;
     {    vecbasevector ebv;
          ebv.push_back( EBV( ) );
          ebv.push_back( EBV( ) );
          ebv.back( ).ReverseComplement( );
          MakeKmerLookup0( ebv, ebv_kmers_plus );    }

     // Find long novel sequences.

     vec< triple< int, int, vec<basevector> > > novels;
     vec< triple< int, int, vec<String> > > novelms;
     {    
          // Build 400-mer lookup table for human.

          cout << Date( ) << ": building genome lookup table" << endl;
          vecbasevector genome( "/wga/scr4/bigrefs/grch38/genome.fastb" );
          const int KG = 400;
          vec< triple<kmer<KG>,int,int> > genome_kmers_plus;
          MakeKmerLookup0( genome, genome_kmers_plus );

          // Find long novel sequences.

          cout << Date( ) << ": finding long novel sequences" << endl;
          int lcount = 0;
          for ( int s = 0; s < I.isize( ); s++ )
          {    
               // Create lookup table for novel stuff we've already seen.

               vecbasevector gn;
               for ( int i = 0; i < novels.isize( ); i++ )
               for ( int j = 0; j < novels[i].third.isize( ); j++ )
                    gn.push_back( novels[i].third[j] );
               vec< triple<kmer<KG>,int,int> > gn_kmers_plus;
               MakeKmerLookup0( gn, gn_kmers_plus );

               const int min_len = 10000;
               #pragma omp parallel for schedule(dynamic, 1000)
               for ( int li = 0; li < LINES[s].isize( ); li++ )
               {    if ( LENS[s][li] < min_len ) continue;
                    const vec<vec<vec<int>>>& L = LINES[s][li];
                    int e1 = L.front( )[0][0], e2 = L.back( )[0][0];
                    if ( make_pair( INV[s][e2], INV[s][e1] ) < make_pair( e1, e2 ) )
                         continue;

                    #pragma omp critical
                    {    if ( lcount % 100 == 0 )
                              DPRINT2( lcount, LINES[s].size( ) ); // XXXXXXXXXXXXXX
                         lcount++;    }

                    // Split line into contigs.

                    vec<vec<vec<vec<int>>>> tigs;
                    MakeTigs( L, tigs );

                    // Get sequence for each contig, always choosing the first path.

                    vecbasevector x;
                    for ( int c = 0; c < tigs.isize( ); c++ )
                    {    const vec<vec<vec<int>>>& L = tigs[c];
                         vec<int> e;
                         for ( int j = 0; j < L.isize( ); j++ )
                         for ( int k = 0; k < L[j][0].isize( ); k++ )
                              e.push_back( L[j][0][k] );
                         x.push_back( HB[s].Cat(e) );    }

                    // Find the 400-mers that match the reference.  This could
                    // be much faster, since we're only interested in 10 kb
                    // stretches that are devoid of matches.

                    vecbitvector ref1, ref2;
                    MarkReferenceKmers( x, genome_kmers_plus, ref1 );
                    MarkReferenceKmers( x, gn_kmers_plus, ref2 );
                    vecbitvector ref(ref1);
                    for ( int j = 0; j < (int) ref.size( ); j++ )
                    for ( int l = 0; l < (int) ref[j].size( ); l++ )
                         if ( ref2[j][l] ) ref[j].Set( l, True );

                    // Find the ranges that are reference-free.

                    vec<vec<basevector>> novel;
                    ExtractReferenceFree( x, ref, KG, min_len, novel );

                    // Check for EBV.
               
                    if ( IsEbv( novel, ebv_kmers_plus ) ) continue;

                    // Align novel to the common segments from L, and mark the
                    // parts of novel that are not aligned as "lower case".

                    vec<vec<String>> novelm( novel.size( ) );
                    for ( int i = 0; i < novel.isize( ); i++ )
                    for ( int j = 0; j < novel[i].isize( ); j++ )
                    {    String s = novel[i][j].ToString( );
                         for ( int l = 0; l < s.isize( ); l++ )
                              s[l] = tolower( s[l] );
                         novelm[i].push_back(s);    }
                    vec<basevector> commons;
                    for ( int i = 0; i < L.isize( ); i += 2 )
                         commons.push_back( HB[s].EdgeObject( L[i][0][0] ) );
                    vecbasevector all;
                    for ( int i = 0; i < commons.isize( ); i++ )
                         all.push_back( commons[i] );
                    vec< pair<int,int> > origin;
                    for ( int i = 0; i < novel.isize( ); i++ )
                    for ( int j = 0; j < novel[i].isize( ); j++ )
                    {    all.push_back( novel[i][j] );
                         origin.push( i, j );    }
                    const int KX = 200;
                    vec< triple<kmer<KX>,int,int> > kmers_plus;
                    MakeKmerLookup3( all, kmers_plus );
                    for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
                    {    int64_t j;
                         for ( j = i + 1; j < kmers_plus.jsize( ); j++ )
                         {    if ( kmers_plus[j].first != kmers_plus[i].first ) 
                                   break;    }
                         for ( int64_t l1 = i; l1 < j; l1++ )
                         for ( int64_t l2 = i; l2 < j; l2++ )
                         {    int id1 = kmers_plus[l1].second;
                              int id2 = kmers_plus[l2].second - commons.isize( );
                              int pos2 = kmers_plus[l2].third;
                              if ( id1 >= commons.isize( ) || id2 < 0 ) continue;
                              int m1 = origin[id2].first, m2 = origin[id2].second;
                              for ( int r = pos2; r < pos2 + KX; r++ )
                              {    novelm[m1][m2][r] 
                                        = toupper( novelm[m1][m2][r] );    }    }
                         i = j - 1;    }

                    // Save.

                    #pragma omp critical
                    {    for ( int l = 0; l < (int) novel.size( ); l++ )
                         {    novels.push( s, li, novel[l] );
                              novelms.push( s, li, novelm[l] );    }    }    }

               // Check second paths.  Note that we're checking versus all novels,
               // across samples, may not make sense.

               int depth = 10;
               for ( int di = 1; di < depth; di++ )
               {    cout << Date( ) << ": checking second paths, depth = "
                         << di << endl;
                    vecbasevector gn;
                    for ( int i = 0; i < novels.isize( ); i++ )
                    for ( int j = 0; j < novels[i].third.isize( ); j++ )
                         gn.push_back( novels[i].third[j] );
                    vec< triple<kmer<KG>,int,int> > gn_kmers_plus;
                    MakeKmerLookup0( gn, gn_kmers_plus );
                    #pragma omp parallel for schedule(dynamic, 100)
                    for ( int li = 0; li < LINES[s].isize( ); li++ )
                    {    if ( LENS[s][li] < min_len ) continue;
                         const vec<vec<vec<int>>>& L = LINES[s][li];
                         int e1 = L.front( )[0][0], e2 = L.back( )[0][0];
                         if ( make_pair( INV[s][e2], INV[s][e1] ) 
                              < make_pair( e1, e2 ) )
                         {    continue;    }
                         for ( int i = 1; i < L.isize( ); i += 2 )
                         {    if ( L[i].isize( ) <= di ) continue;
                              Bool has_gap = False;
                              for ( int j = 0; j < L[i][di].isize( ); j++ )
                              {    if ( HB[s].Bases( L[i][di][j] ) == 0 ) 
                                        has_gap = True;    }
                              if (has_gap) continue;
                              basevector b = HB[s].Cat( L[i][di] );
                              if ( b.isize( ) < min_len ) continue;
                              vecbasevector x;
                              x.push_back(b);
                              vecbitvector ref1, ref2;
                              MarkReferenceKmers( x, genome_kmers_plus, ref1 );
                              MarkReferenceKmers( x, gn_kmers_plus, ref2 );
                              vecbitvector ref(ref1);
                              for ( int j = 0; j < (int) ref.size( ); j++ )
                              for ( int l = 0; l < (int) ref[j].size( ); l++ )
                                   if ( ref2[j][l] ) ref[j].Set( l, True );
                              vec<vec<basevector>> novel;
                              ExtractReferenceFree( x, ref, KG, min_len, novel );
                              if ( novel.empty( ) ) continue;
                              if ( IsEbv( novel, ebv_kmers_plus ) ) continue;
                              #pragma omp critical
                              {    cout << "found";
                                   for ( int j = 0; j < novel.isize( ); j++ )
                                   for ( int l = 0; l < novel[j].isize( ); l++ )
                                        cout << " " << novel[j][l].size( );
                                   cout << endl;
                                   for ( int l = 0; l < (int) novel.size( ); l++ )
                                   {    novels.push( s, li, novel[l] );    
                                        vec<String> x;
                                        for ( int j = 0; j < novel[l].isize( ); 
                                             j++ )
                                        {    x.push_back( 
                                                  novel[l][j].ToString( ) );    }
                                        novelms.push( s, li, x );
                                        }    }    }    }    }    }    }

     // Sorting novels.

     cout << Date( ) << ": sorting novels" << endl;
     ParallelSort(novels), ParallelSort(novelms);

     // Print the novels.

     int part = 1;
     int last_li = -1;
     vec<int> totals( 1000000, 0 );
     int64_t all = 0;
     for ( int i = 0; i < novelms.isize( ); i++ )
     {    int s = novelms[i].first, li = novelms[i].second;
          if ( li != last_li ) part = 1;
          else part++;
          double cov = 2 * COVS[s][0][ LINES[s][li][0][0][0] ].Cov( );
          cout << ">part_" << part << "_of_line_" << li << " (sample="
               << subsam_names[s] << ",cov=" << cov << ")\n";
          int nc = 0;
          for ( int j = 0; j < (int) novelms[i].third.size( ); j++ )
          {    if ( j > 0 )
               {    for ( int l = 0; l < 100; l++ )
                    {    if ( nc > 0 && nc % 80 == 0 ) cout << "\n";
                         nc++;
                         cout << "N";    }    }
               for ( int l = 0; l < novelms[i].third[j].isize( ); l++ )
               {    if ( nc > 0 && nc % 80 == 0 ) cout << "\n";
                    nc++;
                    cout << novelms[i].third[j][l];    }    }
          cout << "\n";
          for ( int j = 0; j < (int) novelms[i].third.size( ); j++ )
          {    int n = novelms[i].third[j].size( );
               all += n;
               int c = int(round(cov));
               totals[c] += n;    }
          last_li = li;    }
     
     // Print summary.

     cout << "\n";
     vec<vec<String>> rows;
     vec<String> row = { "CN", "BASES" };
     rows.push_back(row);
     for ( int j = 0; j < totals.isize( ); j++ )
     {    if ( totals[j] > 0 )
          {    row.clear( );
               row.push_back( ToString(j), ToStringAddCommas(totals[j]) );
               rows.push_back(row);    }    }
     row.clear( );
     row.push_back( "total", ToStringAddCommas(all) );
     rows.push_back(row);
     PrintTabular( cout, rows, 2, "lr" );
     cout << "\n" << Date( ) << ": done, time used = " << TimeSince(clock)
          << endl;    }
