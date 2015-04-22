///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// NewHuman.  Starting from some human assemblies, build a new assembly by
// harvesting their lines that have size at least 10 kb and do not have any
// 400 base perfect match with GRCh38.  This merged human assembly is minimal:
// it has just enough info to be viewable with NhoodInfo.  It is a K=200 unipath
// graph.
//
// See also Voom.
//
// The 400 base match condition is clearly too stringent.
//
// Coverage values to be taken with a grain of salt because:
// If multiple edges ei in one of the original assemblies map to a given edge e in 
// the combined assembly, then the coverage assigned to e is the MAXIMUM of the
// coverages of the ei's.
//
// This needs a terabyte box to run.  Run time is roughly 15 minutes per sample,
// if you have a lot of samples.
//
// In principle we could try to translate the relevant paths from the originating
// assemblies.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "FastIfstream.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/KmerBaseBroker.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/LongReadsToPaths.h"
#include "paths/long/fosmid/Ebv.cc"
#include "paths/long/large/Lines.h"

int main( )
{    RunTime( );

     // Controls.

     int nsamples = 1; // zero to use all samples; the case of one sample is special
     String OUT_INSTANCE = "new_human";

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

     // Find long unaligned lines.

     vec<vec<Bool>> UN(ns);
     {    // Build 400-mer lookup table for human.

          cout << Date( ) << ": building genome lookup table" << endl;
          vecbasevector genome( "/wga/scr4/bigrefs/grch38/genome.fastb" );
          const int KG = 400;
          vec< triple<kmer<KG>,int,int> > genome_kmers_plus;
          MakeKmerLookup0( genome, genome_kmers_plus );

          cout << Date( ) << ": finding long unaligned lines" << endl;
          for ( int s = 0; s < I.isize( ); s++ )
          {    UN[s].resize( HB[s].E( ), False );
               const int min_len = 10000;
               for ( int i = 0; i < LINES[s].isize( ); i++ )
               {    if ( LENS[s][i] < min_len ) continue;
                    const vec<vec<vec<int>>>& L = LINES[s][i];
                    Bool aligned = False;
                    Bool ebv = False;
                    for ( int j = 0; j < L.isize( ); j ++ )
                    for ( int k = 0; k < L[j].isize( ); k++ )
                    for ( int l = 0; l < L[j][k].isize( ); l++ )
                    {    int e = L[j][k][l];
                         if ( HITS[s][e].nonempty( ) 
                              || HITS[s][ INV[s][e] ].nonempty( ) )
                         {    aligned = True;    }
                         if ( !aligned && HB[s].Bases(e) > 0 )
                         {    kmer<KE> x;
                              x.SetToSubOf( HB[s].EdgeObject(e), 0 );
                              int64_t low = LowerBound1( ebv_kmers_plus, x );
                              int64_t high = UpperBound1( ebv_kmers_plus, x );
                              if ( low < high ) ebv = True;    }    }
                    if (aligned) continue;
                    if (ebv) continue;

                    Bool known = False;
                    for ( int j = 0; j < L.isize( ); j++ )
                    for ( int k = 0; k < L[j].isize( ); k++ )
                    for ( int l = 0; l < L[j][k].isize( ); l++ )
                    {    int e = L[j][k][l];
                         basevector E = HB[s].EdgeObject(e);
                         for ( int pass = 1; pass <= 2; pass++ )
                         {    if ( pass == 2 ) E.ReverseComplement( );
                              for ( int m = 0; m <= E.isize( ) - KG; m++ )
                              {    kmer<KG> x;
                                   x.SetToSubOf( E, m );
                                   int64_t low = LowerBound1(genome_kmers_plus, x);
                                   int64_t high = UpperBound1(genome_kmers_plus, x);
                                   if ( low < high ) known = True;    }    }    }
                    /*
                    if (known) // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    {    cout << "line " << I[s] << "." << i // XXXXXXXXXXXXXXXXXXXX
                              << " is known" << endl;    } // XXXXXXXXXXXXXXXXXXXXXX
                    */
                    if (known) continue;

                    for ( int j = 0; j < L.isize( ); j++ )
                    for ( int k = 0; k < L[j].isize( ); k++ )
                    for ( int l = 0; l < L[j][k].isize( ); l++ )
                    {    int e = L[j][k][l];
                         if ( HB[s].Bases(e) == 0 ) continue;
                         UN[s][e] = True;    }    }    }    }

     // Combine the unaligned stuff.

     cout << Date( ) << ": combining unaligned stuff" << endl;
     vecbasevector x;
     for ( int s = 0; s < ns; s++ )
     {    vec<int> to_left, to_right;
          HB[s].ToLeft(to_left), HB[s].ToRight(to_right);
          for ( int e = 0; e < HB[s].E( ); e++ )
          {    if ( !UN[s][e] ) continue;
               x.push_back( HB[s].EdgeObject(e) );
               int v = to_left[e], w = to_right[e];
               for ( int i = 0; i < HB[s].To(v).isize( ); i++ )
               {    int e1 = HB[s].ITo( v, i );
                    if ( !UN[s][e1] ) continue;
                    x.push_back( HB[s].Cat( e1, e ) );    }
               for ( int i = 0; i < HB[s].From(w).isize( ); i++ )
               {    int e2 = HB[s].IFrom( w, i );
                    if ( !UN[s][e2] ) continue;
                    x.push_back( HB[s].Cat( e, e2 ) );    }    }    }

     // Form merged assembly.

     cout << Date( ) << ": forming into new assembly" << endl;
     HyperBasevector hb;
     const int K1 = 200;
     {    HyperKmerPath h;
          vecKmerPath paths, paths_rc;
          unsigned const COVERAGE = 5u;
          LongReadsToPaths( x, K1, COVERAGE, 0, False,
               &hb, &h, &paths, &paths_rc );    }

     // Get ancillary data structures.

     vec<int> inv;
     hb.Involution(inv);
     vec<vec<vec<vec<int>>>> lines;
     int MAX_CELL_PATHS = 50;
     int MAX_DEPTH = 10;
     FindLines( hb, inv, lines, MAX_CELL_PATHS, MAX_DEPTH );

     // Set up to compute coverage.  Traverse samples.

     cout << Date( ) << ": start coverage computation" << endl;
     vec<vec<covcount>> covs(ns);
     for ( int s = 0; s < ns; s++ )
          covs[s].resize( hb.E( ), -2 );
     for ( int s = 0; s < ns; s++ )
     {    PRINT(s);

          // Build lookup table.

          cout << Date( ) << ": building all" << endl;
          vecbasevector all;
          for ( int e = 0; e < hb.E( ); e++ )
               all.push_back( hb.EdgeObject(e) );
          for ( int e = 0; e < HB[s].E( ); e++ )
               all.push_back( HB[s].EdgeObject(e) );
          cout << Date( ) << ": building lookup table" << endl;
          vec< triple<kmer<K1>,int,int> > kmers_plus;
          MakeKmerLookup0( all, kmers_plus );
     
          // Compute coverage.

          cout << Date( ) << ": computing coverage" << endl;
          const int64_t batches = 100;
          vec<int64_t> bstart(batches+1);
          int64_t top = kmers_plus.size( );
          for ( int64_t i = 0; i <= batches; i++ )
               bstart[i] = ( top * i ) / batches;

          /*
          for ( int64_t i = 0; i < batches - 1; i++ )
          {    int64_t& t = bstart[i+1];
               while ( t > 0 && t < top 
                    && kmers_plus[t-1].first == kmers_plus[t].first )
               {    t++;    }
               for ( int64_t j = i+1; j < batches - 1; j++ )
               {    if ( bstart[j] == bstart[j+1] && bstart[j+1] < top )
                         bstart[j+1]++;    }    }
          */

          for ( int64_t i = 1; i < batches; i++ )
          {    int64_t& t = bstart[i];
               while( t > bstart[i-1] 
                    && kmers_plus[t].first == kmers_plus[t-1].first )
               {    t--;    }    }

          for ( int64_t i = 1; i < batches; i++ )
          {    if ( bstart[i] > 0 && kmers_plus[ bstart[i] ].first
                    == kmers_plus[ bstart[i] - 1 ].first )
               {    cout << "Oops, problem 1." << endl;
                    Scram(0);    }
               if ( bstart[i] < bstart[i-1] )
               {    cout << "Oops, problem 2." << endl;
                    Scram(0);    }    }
          if ( bstart[0] != 0 )
          {    cout << "Oops, problem 3." << endl;
               Scram(0);    }
          if ( bstart.back( ) != top )
          {    cout << "Oops, problem 4." << endl;
               Scram(0);    }
          for ( int64_t i = 1; i < batches; i++ )
          {    int64_t n = bstart[i] - bstart[i-1];
               int64_t expect = top / batches;
               double ratio = double(n)/double(expect);
               double eps = 0.00001;
               if ( ratio < 1 - eps || ratio > 1 + eps )
               {    PRINT(ratio);
                    cout << "Oops, problem 5." << endl;
                    Scram(0);    }    }

          #pragma omp parallel for schedule(dynamic, 1)
          for ( int64_t bi = 0; bi < batches; bi++ )
          {    
               #pragma omp critical
               {    cout << "starting batch " << bi+1 << endl;    }
               for ( int64_t i = bstart[bi]; i < bstart[bi+1]; i++ )
               {    int64_t j, m;
                    for ( j = i + 1; j < bstart[bi+1]; j++ )
                         if ( kmers_plus[j].first != kmers_plus[i].first ) break;
                    for ( m = i; m < j; m++ )
                         if ( kmers_plus[m].second >= hb.E( ) ) break;
                    for ( int64_t k1 = i; k1 < m; k1++ )
                    for ( int64_t k2 = m; k2 < j; k2++ )
                    {    int e = kmers_plus[k1].second; 
                         int f = kmers_plus[k2].second - hb.E( );
     
                         /*
                         if ( covs[s][e].Def( )
                              && covs[s][e].Cov( ) != COVS[s][0][f].Cov( ) )
                         {    cout << "Alt cov values for edge " << f << ": "
                                   << covs[s][e].Cov( ) << ", "
                                   << COVS[s][0][f].Cov( ) << endl;    }
                         */
     
                         #pragma omp critical
                         {    if ( !covs[s][e].Def( ) 
                                   || covs[s][e].Cov( ) < COVS[s][0][f].Cov( ) )
                              {    covs[s][e].Set( 
                                        COVS[s][0][f].Cov( ) );    }    }    }
                    i = j - 1;    }    }

          // If no source assembly edge shared a kmer with a combined assembly
          // edge, and it's at least 1kb, then its coverage is declared to be zero.

          for ( int e = 0; e < hb.E( ); e++ )
          {    if ( covs[s][e].Cov( ) == -2 && hb.Bases(e) >= 1000 ) 
                    covs[s][e].Set(0);    }    }
     
     // Write files.

     cout << Date( ) << ": writing files" << endl;
     String out_dir = "/wga/scr4/jaffe/GapToy/" + OUT_INSTANCE;
     Mkdir777(out_dir);
     BinaryWriter::writeFile( out_dir + "/a.hbv", hb );
     HyperBasevectorX hbx(hb);
     BinaryWriter::writeFile( out_dir + "/a.hbx", hbx );
     vecbasevector edges;
     for ( int e = 0; e < hb.E( ); e++ )
          edges.push_back( hb.EdgeObject(e) );
     edges.WriteAll( out_dir + "/a.fastb" );
     BinaryWriter::writeFile( out_dir + "/a.inv", inv );
     BinaryWriter::writeFile( out_dir + "/a.covs", covs );
     BinaryWriter::writeFile( out_dir + "/a.lines", lines );
     BinaryWriter::writeFile( out_dir + "/subsam.names", subsam_names );
     cout << Date( ) << ": done" << endl;
     Scram(0);     }
