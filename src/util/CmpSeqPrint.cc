// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


// Print binary aligns produced by CmpSeq.  Use the same arguments as those
// given to CmpSeq.

#include "Alignment.h"
#include "Basevector.h"
#include "FetchReads.h"
#include "system/ParsedArgs.h"
#include "Quality.h"
#include "Qualvector.h"
#include "system/RunTime.h"
#include "String.h"
#include "system/System.h"
#include "Vec.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String_OrDefault(BASE, "");

     // Data for first sequences:

     CommandArgument_String(FILE1);
     CommandArgument_UnsignedInt_OrDefault(COUNT1, 0);
     CommandArgument_String_OrDefault(NAME1, "");
     CommandArgument_String_OrDefault(QUAL1, "");

     // Data for second sequences:

     CommandArgument_String(FILE2);
     CommandArgument_UnsignedInt_OrDefault(COUNT2, 0);
     CommandArgument_String_OrDefault(NAME2, "");
     CommandArgument_String_OrDefault(QUAL2, "");

     // Other options:

     CommandArgument_String_OrDefault(OUTPUT, "");
     CommandArgument_String_OrDefault(OMIT_IF_PROPER_OVERLAP_EXISTS, "False");
     CommandArgument_Bool_OrDefault(NO_BREAK, False);
     CommandArgument_UnsignedInt_OrDefault(K, 24);
     CommandArgument_UnsignedInt_OrDefault(MIN_OVERLAP, 50);
     CommandArgument_Double_OrDefault(MIN_COVERAGE1, 0.0);
     CommandArgument_Bool_OrDefault(USE_ONLY_SEMIPROPER_ALIGNMENTS, False);
     CommandArgument_Bool_OrDefault(OMIT_PROPER_ALIGNERS, False);
     CommandArgument_UnsignedInt_OrDefault(MIN_SIZE1, 0);
     CommandArgument_Double_OrDefault(MAX_ERROR_RATE, 1.0);
     CommandArgument_UnsignedInt_OrDefault(MAXCLIQ, 1000);
     CommandArgument_UnsignedInt_OrDefault(MAX_ALIGNS, 10000);
     CommandArgument_UnsignedInt_OrDefault(MIN_MUTMER, 0);
     CommandArgument_UnsignedInt_OrDefault(BLOCK_SIZE, 500);
     CommandArgument_Bool_OrDefault(AVOID_PROMISCUOUS_KMERS, False);
     CommandArgument_Bool_OrDefault(TRIM_ALIGNMENT_ENDS, False);
     CommandArgument_UnsignedInt_OrDefault(MIN_PERFECT_MATCH, 0);
     CommandArgument_Double_OrDefault(MAX_QUAL_SCORE, 0.0);
     CommandArgument_Bool_OrDefault(PERFECT_ALIGNMENTS, False);
     CommandArgument_Bool_OrDefault(SMITH_WATERMAN, False);

     // Output options:

     CommandArgument_Bool_OrDefault(OUT_PROPER_MATCHES, True);
     CommandArgument_Bool_OrDefault(OUT_IMPROPER_MATCHES, False);
     CommandArgument_Bool_OrDefault(OUT_MISSING1, False);
     CommandArgument_Bool_OrDefault(OUT_COVERAGE2, False);
     CommandArgument_Bool_OrDefault(OUT_READS1, False);
     CommandArgument_Bool_OrDefault(CALL_IMPROPER_PROPER, False);
     CommandArgument_Bool_OrDefault(LT, False);
     CommandArgument_Bool_OrDefault(NE, False);
     CommandArgument_Bool_OrDefault(VISUAL_COVERAGE, False);
     CommandArgument_String_OrDefault(BINARY_ALIGNMENTS_FILE, "");

     EndCommandArguments;

     ForceAssert( BINARY_ALIGNMENTS_FILE != "" );
     ForceAssert( QUAL1 != "" && QUAL2 != "" );

     String PREBASE = PRE + "/" + BASE + "/";
     String contigs1_file = PREBASE + FILE1, contigs2_file = PREBASE + FILE2;

     ostream *outp;
     String output = PRE + "/" + OUTPUT;
     outp = ( OUTPUT.empty() 
	      ? (ostream*) &cout 
	      : (ostream*) new ofstream( output.c_str() ) );
     ostream &out = *outp;

     // Read in the contigs, including their id's.

     vecbasevector EE, EErc;
     int N0, N1;
     {    vecbasevector cons2;

          int breaker = !NO_BREAK ? 20 : 0;
          if ( !contigs1_file.Contains( ".fastb", -1 ) )
               FetchReads( EE, 0, contigs1_file, breaker );
          else EE.ReadAll(contigs1_file);

          if ( COUNT1 > 0 )
          {    if ( EE.size( ) < COUNT1 )
                    FatalErr( "COUNT1 value is too large" );
               EE.resize(COUNT1);    }
          N0 = EE.size( );

          if ( !contigs2_file.Contains( ".fastb", -1 ) )
               FetchReads( cons2, 0, contigs2_file, breaker );
          else cons2.ReadAll(contigs2_file);

          if ( COUNT2 > 0 )
          {    if ( cons2.size( ) < COUNT2 )
                    FatalErr( "COUNT2 value is too large" );
               cons2.resize(COUNT2);    }
          N1 = cons2.size( );
          for ( size_t i = 0; i < cons2.size( ); i++ )
               EE.push_back( cons2[i] );    }
     EErc.resize( EE.size( ) );
     for ( size_t i = 0; i < EE.size( ) ; i++ )
     {    EErc[i].Setsize( EE[i].size( ) );
          EErc[i].ReverseComplement( EE[i] );    }
     vec<int> EE1_length(N0), EE2_length(N1);
     for ( int i = 0; i < N0; i++ )
          EE1_length[i] = EE[i].size( );
     for ( int i = 0; i < N1; i++ )
          EE2_length[i] = EE[ N0 + i ].size( );

     vecqualvector Q, Qrc;
     if ( QUAL1 != "" ) ReadQualityScores( PREBASE + QUAL1, EE1_length, Q, Qrc );
     if ( QUAL2 != "" ) 
          ReadQualityScores( PREBASE + QUAL2, EE2_length, Q, Qrc, False, True );

     vec<String> ids1, ids2;
     String line;
     if ( contigs1_file.Contains( ".fasta", -1 ) )
     {    Ifstream( c1, contigs1_file );
          while(1)
          {    getline( c1, line );
               if ( !c1 ) break;
               if ( !line.Contains( ">", 0 ) ) continue;
               while(1)
               {    if ( line.Contains( " ", -1 ) ) line.erase(line.size( ) - 1, 1);
                    else break;    }
               ids1.push_back( line.After( ">" ) );    }    }
     else if ( NAME1 != "" )
     {    READ( PREBASE + NAME1, vec<String>, ids );
          ids1.resize(N0);
          for ( int i = 0; i < N0; i++ )
               ids1[i] = ids[i];    }
     else
     {    ids1.resize(N0);
          for ( int i = 0; i < N0; i++ )
               ids1[i] = ToString(i);    }
     if ( contigs2_file.Contains( ".fasta", -1 ) )
     {    Ifstream( c2, contigs2_file );
          while(1)
          {    getline( c2, line );
               if ( !c2 ) break;
               if ( !line.Contains( ">", 0 ) ) continue;
               while(1)
               {    if ( line.Contains( " ", -1 ) ) line.erase(line.size( ) - 1, 1);
                    else break;    }
               ids2.push_back( line.After( ">" ) );    }    }
     else if ( NAME2 != "" )
     {    READ( PREBASE + NAME2, vec<String>, ids );
          ids2.resize(N1);
          for ( int i = 0; i < N1; i++ )
               ids2[i] = ids[i];    }
     else
     {    ids2.resize(N1);
          for ( int i = 0; i < N1; i++ )
               ids2[i] = ToString(i);    }
     AssertEq( (int) ids1.size( ), N0 );
     AssertEq( (int) ids2.size( ), N1 );

     // Here's where the work is done:

     READ( PRE + "/" + BINARY_ALIGNMENTS_FILE, vec<alignment_plus>, bin_aligns );
     for ( unsigned int i = 0; i < bin_aligns.size( ); i++ )
     {    alignment_plus& ap = bin_aligns[i];
          int id1 = ap.Id1( ), id2 = ap.Id2( );
          out << "\n[" << ids1[id1] << " vs " 
               << ( ap.Rc2( ) ? "rc to " : "" ) << ids2[id2] << "]\n";
          ap.Print( True, out, EE[id1], EE[id2 + N0], Q[id1], Q[id2 + N0] );  }  }
