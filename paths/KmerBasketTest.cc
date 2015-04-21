///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// KmerBasket.  David's private testcode.

// Note should really use trimmed jumps.

// Run getem first.

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "Qualvector.h"
#include "ReadPairing.h"
#include "paths/HyperBasevector.h"
#include "paths/KmerBasket.h"

int main( )
{
     // Setup.

     RunTime( );

     // Define assembly directory.

     String dir = "/wga/scr4/ALLPATHS/bacteria/Escherichia_coli_29193/"
          "3031H.2[15.3%]+302GJ.2[17.6%]";

     // Unibase and its reverse complement.

     int s1 = 223, s2 = 2160;

     // Get the read ids incident upon these by grepping some files.

     /* (read getem) */

     // Load the lists of read ids.

     vec<int> xf1, xf2, xj1, xj2;
     fast_ifstream inf1( "listf1" );
     fast_ifstream inf2( "listf2" );
     fast_ifstream inj1( "listj1" );
     fast_ifstream inj2( "listj2" );
     String line;
     while(1)
     {    getline( inf1, line );
          if ( inf1.fail( ) ) break;
          xf1.push_back( line.Int( ) );    }
     while(1)
     {    getline( inf2, line );
          if ( inf2.fail( ) ) break;
          xf2.push_back( line.Int( ) );    }
     while(1)
     {    getline( inj1, line );
          if ( inj1.fail( ) ) break;
          xj1.push_back( line.Int( ) );    }
     while(1)
     {    getline( inj2, line );
          if ( inj2.fail( ) ) break;
          xj2.push_back( line.Int( ) );    }

     // Find their partners.

     vec<read_pairing> pf1, pf2, pj1, pj2;
     ReadsToPairs( dir, xf1, pf1, "frag_reads_orig" );
     ReadsToPairs( dir, xf2, pf2, "frag_reads_orig" );
     ReadsToPairs( dir, xj1, pj1, "jump_reads_orig" );
     ReadsToPairs( dir, xj2, pj2, "jump_reads_orig" );
     for ( int j = 0; j < pf1.isize( ); j++ )
          if ( !pf1[j].Dead( ) ) xf2.push_back( pf1[j].Partner( xf1[j] ) );
     for ( int j = 0; j < pf2.isize( ); j++ )
          if ( !pf2[j].Dead( ) ) xf1.push_back( pf2[j].Partner( xf2[j] ) );
     for ( int j = 0; j < pj1.isize( ); j++ )
          if ( !pj1[j].Dead( ) ) xj2.push_back( pj1[j].Partner( xj1[j] ) );
     for ( int j = 0; j < pj2.isize( ); j++ )
          if ( !pj2[j].Dead( ) ) xj1.push_back( pj2[j].Partner( xj2[j] ) );
     UniqueSort(xf1), UniqueSort(xf2), UniqueSort(xj1), UniqueSort(xj2);

     // Dump read ids.

     cout << "\nusing forward fragment reads:" << endl;
     for ( int i = 0; i < xf1.isize( ); i++ )
          cout << xf1[i] << endl;
     cout << "\nusing reverse fragment reads:" << endl;
     for ( int i = 0; i < xf2.isize( ); i++ )
          cout << xf2[i] << endl;
     cout << "\nusing forward jump reads:" << endl;
     for ( int i = 0; i < xj1.isize( ); i++ )
          cout << xj1[i] << endl;
     cout << "\nusing reverse jump reads:" << endl;
     for ( int i = 0; i < xj2.isize( ); i++ )
          cout << xj2[i] << endl;
     cout << endl;

     // Load the reads.

     vecbasevector Rf1, Rf2, Rj1, Rj2;
     vecqualvector Qf1, Qf2, Qj1, Qj2;
     Rf1.Read( dir + "/frag_reads_orig.fastb", xf1 );
     Rf2.Read( dir + "/frag_reads_orig.fastb", xf2 );
     Qf1.Read( dir + "/frag_reads_orig.qualb", xf1 );
     Qf2.Read( dir + "/frag_reads_orig.qualb", xf2 );
     for ( size_t i = 0; i < Rf2.size( ); i++ )
     {    Rf2[i].ReverseComplement( );
          Qf2[i].ReverseMe( );    }
     Rj1.Read( dir + "/jump_reads_orig.fastb", xj1 );
     Rj2.Read( dir + "/jump_reads_orig.fastb", xj2 );
     Qj1.Read( dir + "/jump_reads_orig.qualb", xj1 );
     Qj2.Read( dir + "/jump_reads_orig.qualb", xj2 );
     for ( size_t i = 0; i < Rj2.size( ); i++ )
     {    Rj2[i].ReverseComplement( );
          Qj2[i].ReverseMe( );    }

     // Build common set of reads.

     vecbasevector R;
     vecqualvector Q;
     R.Append(Rf1), R.Append(Rf2), R.Append(Rj1), R.Append(Rj2);
     Q.Append(Qf1), Q.Append(Qf2), Q.Append(Qj1), Q.Append(Qj2);
     
     // Assemble.

     HyperBasevector hb;
     KmerBasket( R, Q, hb, 40, 9, 2, cout );

     // Print HyperKmerPath.

     cout << hb;
     {    Ofstream( out, "xxx.fasta" );
          out << hb;    }

     // Align it.

     cout << "\n" << String( 80, '=' ) << "\n";
     System( "cd \"" + dir + "\"; QueryLookupTable K=12 MM=12 MC=0.15 "
          + "SEQS=/wga/dev/jaffe/BroadCRD/xxx.fasta L=../genome.lookup VISUAL=True "
          + "QUERY_NAMING=from_record QUIET=True NH=True" );    }
