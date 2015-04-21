///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// The maize dataset for DISCOVAR de novo also contains sorghum as a contaminant.
// To split the assembly we use the maize and sorghum reference sequences.  For a 
// given line in the combined assembly, we compute its number of 80-mers that are 
// present in each of the references, denoted here m and s.
//
// Note that the assembly size may be too low because we use a 1 kb cutoff and
// a substantial fraction of the genome may be in smaller lines.
//
// The same approach ought to work for other mixed datasets.

#include "Basevector.h"
#include "MainTools.h"
#include "paths/HyperBasevector.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/Lines.h"

int main( )
{    RunTime( );

     // Load the assembly.

     String final_dir = "/wga/scr4/jaffe/GapToy/50158.maize/a.final";
     HyperBasevector hb;
     BinaryReader::readFile( final_dir + "/a.hbv", &hb );
     vec<int> inv;
     BinaryReader::readFile( final_dir + "/a.inv", &inv );
     vec<vec<vec<vec<int>>>> lines;
     BinaryReader::readFile( final_dir + "/a.lines", &lines );
     vec<int> tol, lens;
     GetTol( hb, lines, tol );
     GetLineLengths( hb, lines, lens );

     // Flag edges that are in the 'first path' of a line.

     vec<Bool> first( hb.E( ), False );
     for ( int l = 0; l < lines.isize( ); l++ )
     for ( int j = 0; j < lines[l].isize( ); j++ )
     for ( int k = 0; k < lines[l][j][0].isize( ); k++ )
          first[ lines[l][j][0][k] ] = True;

     // Load reference sequences.

     cout << Date( ) << ": loading genome" << endl;
     vecbasevector genome;
     genome.ReadAll( "/wga/scr4/bigrefs/maize/genome.fastb" );
     int nmaize = genome.size( );
     genome.ReadAll( "/wga/scr4/bigrefs/sorghum/genome.fastb", True );
     int nsorghum = genome.size( ) - nmaize;

     // Add in assembly.

     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
          genome.push_back( hb.EdgeObject(e) );

     // Build kmers.

     const int K = 80;
     vec< triple<kmer<K>,int,int> > kmers_plus;
     cout << Date( ) << ": making lookup table" << endl;
     MakeKmerLookup( genome, kmers_plus );

     // Find shared kmers.

     cout << Date( ) << ": finding shared kmers" << endl;
     vec<int> maize_count( lines.size( ), 0 ), sorghum_count( lines.size( ), 0 );
     for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
     {    int64_t j, m, n;
          for ( j = i + 1; j < kmers_plus.jsize( ); j++ )
               if ( kmers_plus[j].first != kmers_plus[i].first ) break;
          for ( m = i; m < j; m++ )
               if ( kmers_plus[m].second >= nmaize ) break;
          for ( n = m; n < j; n++ )
               if ( kmers_plus[n].second >= nmaize + nsorghum ) break;
          int maize = m - i, sorghum = n - m;
          for ( int64_t l = n; l < j; l++ )
          {    int e = kmers_plus[l].second - nmaize - nsorghum;
               if ( !first[e] ) continue;
               if ( maize > 0 ) maize_count[tol[e]]++;
               if ( sorghum > 0 ) sorghum_count[tol[e]]++;    }
          i = j - 1;    }
     for ( int l = 0; l < lens.isize( ); l++ )
     {    if ( lens[l] < 1000 ) continue;
          PRINT4( l, lens[l], maize_count[l], sorghum_count[l] );    }

     // Classify lines.

     vec<char> nature( lines.size( ), '?' );
     for ( int l = 0; l < lines.isize( ); l++ )
     {    if ( maize_count[l] == 0 && sorghum_count[l] == 0 ) continue;
          if ( maize_count[l] >= 10 * sorghum_count[l] ) nature[l] = 'M';
          if ( sorghum_count[l] >= 10 * maize_count[l] ) nature[l] = 'S';    }
     int m_tot = 0, s_tot = 0; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     for ( int l = 0; l < lines.isize( ); l++ ) // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     {    if ( nature[l] == 'M' ) m_tot++; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          if ( nature[l] == 'S' ) s_tot++;    } // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     PRINT2( m_tot, s_tot ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     // Compute stats.

     int64_t n_maize = 0, n_sorghum = 0, n_unknown = 0;
     int64_t gc_maize = 0, gc_sorghum = 0, gc_unknown = 0;
     vec<int> mtigs, stigs, utigs;
     vec<int> msuper, ssuper, usuper;
     for ( int l = 0; l < lines.isize( ); l++ )
     {    if ( lens[l] < 1000 ) continue;
          if ( nature[l] == 'M' ) n_maize += lens[l];
          else if ( nature[l] == 'S' ) n_sorghum += lens[l];
          else n_unknown += lens[l];

          // Split line into contigs.

          const vec<vec<vec<int>>>& L = lines[l];
          vec<vec<vec<vec<int>>>> contigs;
          MakeTigs( L, contigs );

          // Get stats from contigs.

          int tlen = 0;
          for ( int c = 0; c < contigs.isize( ); c++ )
          {    const vec<vec<vec<int>>>& L = contigs[c];
               vec<int> e;
               for ( int j = 0; j < L.isize( ); j++ )
               for ( int k = 0; k < L[j][0].isize( ); k++ )
                    e.push_back( L[j][0][k] );
               basevector b = hb.Cat(e);
               tlen += b.size( );
               int gc = 0;
               for ( int j = 0; j < b.isize( ); j++ )
                    if ( as_base( b[j] ) == 'G' || as_base( b[j] ) == 'C' ) gc++;
               if ( nature[l] == 'M' ) 
               {    mtigs.push_back( b.size( ) );
                    gc_maize += gc;    }
               else if ( nature[l] == 'S' ) 
               {    stigs.push_back( b.size( ) );
                    gc_sorghum += gc;    }
               else 
               {    utigs.push_back( b.size( ) );
                    gc_unknown += gc;    }    }

          // Get scaffold stats.

          if ( nature[l] == 'M' ) msuper.push_back(tlen);
          else if ( nature[l] == 'S' ) ssuper.push_back(tlen);
          else usuper.push_back(tlen);    }

     // Print summary stats.

     Sort(mtigs), Sort(stigs), Sort(utigs), Sort(msuper), Sort(ssuper), Sort(usuper);
     cout << "\nmaize\n";
     cout << "assembly size = " << ToStringAddCommas(n_maize/2) << endl;
     cout << "gc content = " << PERCENT_RATIO( 3, gc_maize, n_maize ) << endl;
     if ( mtigs.nonempty( ) ) cout << "N50 contig = " << N50(mtigs) << endl;
     if ( msuper.nonempty( ) ) cout << "N50 scaffold = " << N50(msuper) << endl;
     cout << "\nsorghum\n";
     cout << "assembly size = " << ToStringAddCommas(n_sorghum/2) << endl;
     cout << "gc content = " << PERCENT_RATIO( 3, gc_sorghum, n_sorghum ) << endl;
     if ( stigs.nonempty( ) ) cout << "N50 contig = " << N50(stigs) << endl;
     if ( ssuper.nonempty( ) ) cout << "N50 scaffold = " << N50(ssuper) << endl;
     cout << "\nunknown\n";
     cout << "assembly size = " << ToStringAddCommas(n_unknown/2) << endl;
     cout << "gc content = " << PERCENT_RATIO( 3, gc_unknown, n_unknown ) << endl;
     if ( utigs.nonempty( ) ) cout << "N50 contig = " << N50(utigs) << endl;
     if ( usuper.nonempty( ) ) cout << "N50 scaffold = " << N50(usuper) << endl;
     Scram(0);    }
