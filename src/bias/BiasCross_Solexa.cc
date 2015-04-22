/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/// BiasCross_Solexa.  Given two lanes of Solexa data from the same compute bias by
/// (a) comparing half the reads from each lane to the other half from that lane;
/// (b) comparing half the reads from each lane to half from the other lane.
/// The bias numbers thusly computed are comparable to each other but not 
/// comparable to the bias numbers computed using StartBias.

#include "Basevector.h"
#include "MainTools.h"
#include "PackAlign.h"
#include "bias/StartBias.h"
#include "lookup/PerfectCount.h"
#include "math/Functions.h"
#include "solexa/SolexaTools.h"
#include "solexa/SolexaMetrics.h"
#include "solexa/SolexaTools.h"

int main( int argc, char *argv[] ) 
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(HEAD1, "flowcell.lane for first lane");
     CommandArgument_String_Doc(HEAD2, "flowcell.lane for second lane");
     CommandArgument_Int_OrDefault_Doc(N, 27, "length to truncate reads to");
     CommandArgument_String_OrDefault_Doc(REF, "", "alternative reference");
     EndCommandArguments;

     // Load metrics.  Quit if there is no reference or genome is big.

     solexa_metric_db db( HEAD1 + ".metrics" );
     Bool have_ref = IsRegularFile( HEAD1 + ".reference.lookup" );
     if ( !have_ref ) exit(0);
     ForceAssert( db.Defined( "genome_size" ) );
     if ( db.Value( "genome_size" ).Int( ) > 100 * 1000 * 1000 ) exit(0);

     // Load bases.

     vecbasevector bases1( HEAD1 + ".fastb" ), bases2( HEAD2 + ".fastb" );

     // Shrink to desired size. 

     for ( size_t i = 0; i < bases1.size( ); i++ )
          bases1[i].resize(N);
     for ( size_t i = 0; i < bases2.size( ); i++ )
          bases2[i].resize(N);

     // Load more data.

     String refname = ( REF == "" ? HEAD1 + ".reference" : REF );
     vecbasevector ref( refname + ".fastb" );
     vecbasevector rref = ref;
     ReverseComplement(rref);

     // Find piles.

     vec<placement_mark> places1, places2;
     PerfectPick( bases1, refname + ".lookup", FW_OR_RC, places1 );
     PerfectPick( bases2, refname + ".lookup", FW_OR_RC, places2 );

     // Report bias.

     PRINT2( places1.size( ), places2.size( ) );
     int P = Min( places1.size( ), places2.size( ) ) / 2;
     RandomShuffle(places1), RandomShuffle(places2);

     vec<placement_mark> places11, places12, places21, places22;
     for ( int i = 0; i < P; i++ )
     {    places11.push_back( places1[i] ), places12.push_back( places1[P+i] );
          places21.push_back( places2[i] ), places22.push_back( places2[P+i] );    }
     
     cout << "\ndirect:\n";
     cout << "11 vs 12: ";
     cout << StartBiasRelRaw( ref, rref, places11, places12 ) << "\n";
     cout << "21 vs 22: ";
     cout << StartBiasRelRaw( ref, rref, places21, places22 ) << "\n";
     cout << "\ncross:\n";
     cout << "11 vs 21: ";
     cout << StartBiasRelRaw( ref, rref, places11, places21 ) << "\n";
     cout << "11 vs 22: ";
     cout << StartBiasRelRaw( ref, rref, places11, places22 ) << "\n";
     cout << "12 vs 21: ";
     cout << StartBiasRelRaw( ref, rref, places12, places21 ) << "\n";
     cout << "12 vs 22: ";
     cout << StartBiasRelRaw( ref, rref, places12, places22 ) << "\n";
     cout << "\nfull cross:\n";
     static vec<placement_mark> places1x, places2x;
     int reps = 20;
     double cross = 0.0;
     for ( int i = 0; i < reps; i++ )
     {    places1x = places1, places2x = places2;
          RandomShuffle(places1x), RandomShuffle(places2x);
          places1x.resize(2*P), places2x.resize(2*P);
          cross += StartBiasRelRaw( ref, rref, places1x, places2x ).Double( );    }
     cout << ToString( cross / double(reps), 1 ) << "\n";    }
