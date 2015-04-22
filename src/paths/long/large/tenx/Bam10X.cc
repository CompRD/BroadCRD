///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Parse a 10X BAM file.  Hardcoded for the moment.  For N=1, took 
// 4.6 hours on crd26.

// Generate:
// 1. parallel:
//    10X.fastb   bases (with reads groups by barcode)
//    10X.qualp   packed quals
//    10X.bc      canonicalized barcodes
// 2. 10X.bci     barcode start points
// 3. more stuff not yet described.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "Qualvector.h"
#include "feudal/PQVec.h"
#include "kmers/KmerRecord.h"
#include "paths/long/large/tenx/TenxDirs.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_Int_Doc(N, "1 or 2 or 3 or 4 or 6");
     EndCommandArguments;

     // Set up directories etc.

     String dir, odir, tdir, bam;
     SetTenxDirs( N, dir, odir, tdir );
     if ( N == 1 ) bam = tdir + "/../NA12878_WGS_reads_pos_sorted.bam";
     if ( N == 2 ) bam = tdir + "/../HCC1143_BL_possorted_bam.bam";
     if ( N == 3 ) bam = tdir + "/../HCC1143_possorted_bam.bam";
     if ( N == 4 ) bam = tdir + "/../S0117_XDP_possorted_bam.bam";
     if ( N == 5 ) bam = tdir + "/../S0118_XDP_possorted_bam.bam";
     if ( N == 6 ) bam = tdir + "/../S0119_XDP_possorted_bam.bam";
     Mkdir777(dir);

     // Define data structures.

     vec<basevector> bases;
     vec<qualvector> quals;
     vec<String> names;
     vec<uint32_t> bcs;

     // Read bam.

     fast_pipe_ifstream in( "samtools view " + bam );
     String line;
     vec<String> x;
     double clock = WallClockTime( );
     cout << Date( ) << ": reading" << endl;
     int64_t secondaries = 0, supplementaries = 0, missing = 0, nlines = 0;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          if ( line.Contains( "@", 0 ) ) continue;
          nlines++;
          Tokenize( line, x );
          ForceAssert( x.size( ) >= 11 );
          Bool rc = x[1].Int( ) & 0x10, first = x[1].Int( ) & 0x40;
          Bool secondary = ( x[1].Int( ) & 0x100 ) != 0;
          Bool supplementary = ( x[1].Int( ) & 0x800 ) != 0;
          if (secondary) 
          {    secondaries++;
               continue;    }
          if (supplementary) 
          {    supplementaries++;
               continue;    }
          String seq = x[9], qual = x[10], bc;
          String name = x[0];
          if (first) name += "/1";
          else name += "/2";
          for ( int i = 0; i < seq.isize( ); i++ )
          {    if ( seq[i] == 'N' ) 
               {    if ( !rc ) seq[i] = 'A';
                    else seq[i] = 'T';    }    }
          for ( int i = 0; i < x.isize( ); i++ )
          {    if ( x[i].Contains( "BX:Z:1-", 0 ) ) bc = x[i].After( "BX:Z:1-" );
               if ( x[i].Contains( "BX:Z:", 0 ) && x[i].Contains( "-1", -1 ) )
                    bc = x[i].Between( "BX:Z:", "-1" );    }
          if ( bc == "" )
          {    missing++;
               static Bool mspoken = False;
               if ( !mspoken && nlines >= 10000 )
               {    mspoken = True;
                    cout << "\nInitial look shows "
                         << PERCENT_RATIO( 3, missing, nlines )
                         << " of bam entries as missing barcode.\n" << endl;    }
               continue;    }
          basevector b(seq);
          qualvector q( seq.size( ) );
          ForceAssertEq( seq.size( ), qual.size( ) );
          for ( int j = 0; j < (int) seq.size( ); j++ )
               q[j] = (int) qual[j] - 33;
          if (rc)
          {    b.ReverseComplement( );
               q.ReverseMe( );    }
          ForceAssertEq( bc.isize( ), 14 );
          uint32_t m = 0;
          for ( int j = 0; j < bc.isize( ); j++ )
          {    m *= 4;
               m += as_char( bc[j] );    }
          bases.push_back(b), quals.push_back(q); 
          bcs.push_back(m), names.push_back(name);    }
     cout << nlines << " bam entries" << endl;
     cout << missing << " = " << PERCENT_RATIO( 3, missing, nlines )
          << " of bam entries were missing a barcode" << endl;
     cout << secondaries << " = " << PERCENT_RATIO( 3, secondaries, nlines )
          << " of bam entries were marked as secondary" << endl;
     cout << supplementaries << " = " << PERCENT_RATIO( 3, supplementaries, nlines )
          << " of bam entries were marked as supplementary" << endl;
     cout << PERCENT_RATIO( 3, bcs.jsize( ), nlines )
          << " of bam entries were used" << endl;
     cout << bcs.size( ) << " passing bam entries" << endl << endl;
     if ( bcs.size( ) == 0 )
     {    cout << "Uh oh, nothing to do." << endl;
          Scram(0);    }

     // Sort.

     cout << Date( ) << ": sorting" << endl;
     for ( int64_t i = 0; i < names.jsize( ); i++ )
          names[i] = ToString( bcs[i] ) + "." + names[i];
     ParallelSortSync( names, bcs, bases, quals );

     // Fix names and count orphans.

     int multis = 0;
     for ( int64_t i = 0; i < names.jsize( ) - 1; i++ )
          if ( names[i] == names[i+1] ) multis++;
     cout << Date( ) << ": fixing names" << endl;
     for ( int64_t i = 0; i < names.jsize( ); i++ )
          names[i] = names[i].After( "." ).RevBefore( "/" );
     int orphans = 0;
     int bads = 0;
     for ( int64_t i = 0; i < names.jsize( ); i++ )
     {    Bool paired = False;
          if ( i % 2 == 1 && names[i] == names[i-1] ) 
          {    paired = True;
               if ( bcs[i] != bcs[i-1] ) bads = True;    }
          else if ( i % 2 == 0 && i < names.jsize( ) - 1 && names[i] == names[i+1] ) 
          {    paired = True;
               if ( bcs[i] != bcs[i+1] ) bads = True;    }
          if ( !paired ) orphans++;    }
     PRINT3( multis, orphans, bads );
     cout << PERCENT_RATIO( 3, orphans, bases.jsize( ) )
          << " of reads are orphaned" << endl;
     if ( multis > 0 || orphans > 0 || bads > 0 )
          cout << "WARNING: there should be no multis or bads or orphans" << endl;

     // Reduce barcodes.

     {    int bc = 0;
          vec<uint32_t> bcs2( bcs.size( ) );
          for ( int64_t i = 0; i < bcs.jsize( ); i++ )
          {    if ( i > 0 && bcs[i] != bcs[i-1] ) bc++;
               bcs2[i] = bc;    }
          bcs = bcs2;    }

     // Find barcode start points.

     vec<int64_t> bci( bcs.back( ) + 2 );
     bci.back( ) = bcs.size( );
     for ( int64_t j = bcs.jsize( ) - 1; j >= 0; j-- )
          bci[ bcs[j] ] = j;

     // Compress quality scores.

     cout << Date( ) << ": compressing quals" << endl;
     VecPQVec pquals;
     convertAssignParallel( quals.begin( ), quals.end( ), pquals );

     // Pack bases.

     cout << Date( ) << ": packing bases" << endl;
     vecbasevector basesf;
     for ( int64_t i = 0; i < bases.jsize( ); i++ )
          basesf.push_back( bases[i] );

     // Write.

     cout << Date( ) << ": writing" << endl;
     basesf.WriteAll( tdir + "/10X.fastb" );
     pquals.WriteAll( tdir + "/10X.qualp" );
     BinaryWriter::writeFile( tdir + "/10X.bc", bcs );
     BinaryWriter::writeFile( tdir + "/10X.bci", bci );

     // Build alt indices.

     cout << Date( ) << ": building alt indices" << endl;
     const int K = 88;
     vec< pair<kmer<K>,int64_t> > kb;
     vec<int64_t> kb_bci( bci.size( ), 0 );
     kmer<K> xx;
     for ( int64_t i = 0; i < (int64_t) basesf.size( ); i++ )
     {    kb_bci[ bcs[i] + 1 ] = kb.size( ) + 1;
          if ( basesf[i].isize( ) == K )
          {    xx.SetToSubOf( basesf[i], 0 );
               kb.push( xx, i );    }    }
     #pragma omp parallel for
     for ( int i = 0; i < bci.isize( ) - 1; i++ )
          sort( kb.begin( ) + kb_bci[i], kb.begin( ) + kb_bci[i+1] );
     vec< kmer<K> > kb1( kb.size( ) );
     vec< int64_t > kb2( kb.size( ) );
     #pragma omp parallel for
     for ( int64_t i = 0; i < kb.jsize( ); i++ )
     {    kb1[i] = kb[i].first;
          kb2[i] = kb[i].second;    }
     BinaryWriter::writeFile( tdir + "/10X.kb", kb1 );
     BinaryWriter::writeFile( tdir + "/10X.kb_ids", kb2 );
     BinaryWriter::writeFile( tdir + "/10X.kb_bci", kb_bci );

     // Done.

     cout << TimeSince(clock) << " used" << endl;
     cout << Date( ) << ": done" << endl;
     Scram(0);    }
