/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/// SimulateBias.  Generate simulated perfect reads at random from a genome.
/// Bias in three ways:
/// (a) duplicate reads at rate DUP;
/// (b) reduce genome coverage to COV;
/// (c) apply given distribution to first base [via BASE1="{*%,*%,*%,*%}"].
/// What bias score do we get?

#include "Basevector.h"
#include "MainTools.h"
#include "ParseSet.h"
#include "bias/StartBias.h"
#include "lookup/PerfectCount.h"
#include "random/FindRandomReads.h"
#include "random/Random.h"

int main( int argc, char *argv[] ) 
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(GENOME, 
                                "Files GENOME.fastb and GENOME.lookup must exist.");
     CommandArgument_Int(NREADS);
     CommandArgument_Double_OrDefault_Doc(DUP, 0.0, "Fraction to be duplicated.");
     CommandArgument_Double_OrDefault_Doc(COV, 1.0, "Fraction of genome covered.");
     CommandArgument_Int_OrDefault(READLENGTH, 27);
     CommandArgument_Bool_OrDefault_Doc(CROSS, False, "Compute cross-bias instead.");
     CommandArgument_String_OrDefault_Doc(FIRST_BASE, "", 
          "First base distribution." );
     EndCommandArguments;

     // Check arguments.

     ForceAssert( COV > 0.0 && COV <= 1.0 );
     vec<double> first_base;
     if ( FIRST_BASE != "" )
     {    ParseDoubleSet( FIRST_BASE, first_base, False );
          ForceAssertEq( first_base.isize( ), 4 );    }

     // Generate the reads.

     vecbasevector genome( GENOME + ".fastb" ), reads;
     vecbasevector rgenome = genome;
     ReverseComplement(rgenome);
     FindRandomReads F;
     vecbasevector sgenome(genome);
     if ( COV < 1.0 )
     {    for ( size_t i = 0; i < genome.size( ); i++ )
          {    sgenome[i].resize( 
                    int( round( double( genome[i].size( ) ) * COV ) ) );    }    }
     if ( FIRST_BASE == "" ) F.Reads( NREADS, sgenome, reads, READLENGTH );
     else F.Reads( NREADS, first_base, sgenome, reads, READLENGTH );
     int extra = int( round( DUP * double(NREADS) ) );
     for ( int i = 0; i < extra; i++ )
          reads.push_back_reserve( reads[ randomx( ) % reads.size( ) ] );

     // Compute bias.

     if ( !CROSS )
     {    vec<placement_mark> places;
          PerfectPick( reads, GENOME + ".lookup", FW_OR_RC, places );
          cout << StartBias( reads, genome, rgenome, places, True, "", True, True )
               << "\n";    }
     else
     {    vecbasevector reads1, reads2;
          reads1.reserve(reads.size()*5/8);
          reads2.reserve(reads.size()*5/8);
          for ( size_t i = 0; i < reads.size( ); i++ )
          {    if ( randomx( ) % 2 == 0 ) reads1.push_back( reads[i] );
               else reads2.push_back( reads[i] );    }
          vec<placement_mark> places1, places2;
          PerfectPick( reads1, GENOME + ".lookup", FW_OR_RC, places1 );
          PerfectPick( reads2, GENOME + ".lookup", FW_OR_RC, places2 );
          cout << StartBiasRelRaw( genome, rgenome, places1, places2 ) 
               << "\n";    }    }
                                                        
