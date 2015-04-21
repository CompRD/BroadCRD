///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// DiffHbv: Find some of the differences between two HyperBasevectors.  They
// must have the same K.  We assume that they are in the same orientation.
//
// Can also read SupportedHyperBasevectors.

#include "MainTools.h"
#include "paths/HyperBasevector.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/SupportedHyperBasevector.h"

template<int K> void DiffHbvCore( const HyperBasevector& hb1, 
     const HyperBasevector& hb2, const Bool PRINT_BASES )
{    for ( int pass = 1; pass <= 2; pass++ )
     {    const HyperBasevector& h1 = ( pass == 1 ? hb1 : hb2 );
          const HyperBasevector& h2 = ( pass == 1 ? hb2 : hb1 );
          int id1 = ( pass == 1 ? 1 : 2 );
          int id2 = ( pass == 1 ? 2 : 1 );
          vec< triple<kmer<K>,int,int> > kmers_plus;
          vecbasevector edges1;
          for ( int e = 0; e < h1.EdgeObjectCount( ); e++ )
               edges1.push_back( h1.EdgeObject(e) );
          MakeKmerLookup0( edges1, kmers_plus );
          vec< kmer<K> > kmers( kmers_plus.size( ) );
          for ( size_t i = 0; i < kmers.size( ); i++ )
               kmers[i] = kmers_plus[i].first;
          for ( int e = 0; e < h2.EdgeObjectCount( ); e++ )
          {    const basevector& E = h2.EdgeObject(e);
               vec<int> missing;
               for ( int j = 0; j <= E.isize( ) - K; j++ )
               {    kmer<K> x;
                    x.SetToSubOf( E, j );
                    if ( !BinMember( kmers, x ) ) missing.push_back(j);    }
               for ( int i = 0; i < missing.isize( ); i++ )
               {    int j;
                    for ( j = i + 1; j < missing.isize( ); j++ )
                         if ( missing[j] != missing[j-1] + 1 ) break;
                    int start = missing[i], stop = missing[j-1] + 1;
                    if (PRINT_BASES) cout << "\n";
                    cout << "hbv " << id2 << ", kmers " << e << "." << start
                         << "-" << stop << " are missing from hbv " << id1 << endl;
                    if (PRINT_BASES)
                    {    for ( int k = start; k < stop + K - 1; k++ )
                         {    cout << as_base( E[k] );
                              if ( ( k - start ) % 80 == 79 && k != stop+K-2 ) 
                                   cout << "\n";    }
                         cout << "\n";    }
                    i = j - 1;   }    }    }    }

int main(int argc, char *argv[])
{    
     RunTime( ); 
     
     BeginCommandArguments;
     CommandArgument_String_OrDefault(HB1, "");
     CommandArgument_String_OrDefault(HB2, "");
     CommandArgument_String_OrDefault(SHB1, "");
     CommandArgument_String_OrDefault(SHB2, "");
     CommandArgument_Bool_OrDefault(PRINT_BASES, False);
     EndCommandArguments;

     // Load HyperBasevectors.

     ForceAssert( HB1 == "" || SHB1 == "" );
     ForceAssert( HB2 == "" || SHB2 == "" );
     HyperBasevector hb1, hb2;
     if ( HB1 != "" ) BinaryReader::readFile( HB1, &hb1 );
     else
     {    SupportedHyperBasevector shb1;
          BinaryReader::readFile( SHB1, &shb1 );
          hb1 = shb1;    }
     if ( HB2 != "" ) BinaryReader::readFile( HB2, &hb2 );
     else
     {    SupportedHyperBasevector shb2;
          BinaryReader::readFile( SHB2, &shb2 );
          hb2 = shb2;    }
     ForceAssertEq( hb1.K( ), hb2.K( ) );

     // Compare them.

     int K = hb1.K( );
     if ( K == 60 ) DiffHbvCore<60>( hb1, hb2, PRINT_BASES );
     else if ( K == 80 ) DiffHbvCore<80>( hb1, hb2, PRINT_BASES );
     else if ( K == 200 ) DiffHbvCore<200>( hb1, hb2, PRINT_BASES );
     else
     {    cout << "That K is unsupported." << endl;
          Scram(1);    }    }
