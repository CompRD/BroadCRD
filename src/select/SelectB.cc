///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// SelectB: select entries from both a fastb and a qualb file.

#include <algorithm>

#include "Basevector.h"
#include "MainTools.h"
#include "ParseSet.h"
#include "Qualvector.h"
#include "feudal/PQVec.h"

int main( int argc, char *argv[] )
{
     RunTime();

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(IN_HEAD);
     CommandArgument_String(OUT_HEAD);
     CommandArgument_String_OrDefault(IDS, "");
     // overrides IDS:
     CommandArgument_String_OrDefault(PIDS, "");
     CommandArgument_Bool_OrDefault( RC, False );
     CommandArgument_Bool_OrDefault( PRE_OUT, False );
     CommandArgument_Bool_OrDefault( USE_QUAL, True );
     CommandArgument_Bool_OrDefault( PRESERVE_ORDER, False );
     CommandArgument_String_OrDefault(NAMES, "");
     CommandArgument_Bool_OrDefault(NAMES_PARALLEL, False);
     EndCommandArguments;

     String input = PRE + "/" + IN_HEAD, output = OUT_HEAD;
     if (PRE_OUT) output = PRE + "/" + OUT_HEAD;

     ForceAssert( IDS != "" || PIDS != "" );
     if ( PIDS != "" )
     {    vec<int> pids;
          ParseIntSet( PIDS, pids );
          vec<int> ids;
          for ( int i = 0; i < pids.isize( ); i++ )
               ids.push_back( 2*pids[i], 2*pids[i]+1 );
          ostringstream out;
          out << printSeq(ids);
          IDS = "{" + out.str( ) + "}";    }

     vec<int> ids, ids_unsorted;

     vecString names;
     if ( NAMES != "" ) names.ReadAll( PRE + "/" + NAMES );

     int64_t total_reads = MastervecFileObjectCount( input + ".fastb" );
     ParseIntSet( IDS, ids_unsorted, false, 0, total_reads );
     ids = ids_unsorted;
     Sort(ids);

     vecbasevector EE;
     EE.Read( input + ".fastb", ids, 0 );

     if (RC)
     {    for ( size_t i = 0; i < EE.size(); ++i )
               EE[i].ReverseComplement();    }

     Ofstream( fastout, output + ".fasta" );
     for ( int i = 0; i < (int) EE.size( ); i++ )
     {    int j = i;
          if (PRESERVE_ORDER) j = BinPosition( ids, ids_unsorted[i] );
          if ( names.size( ) == 0 )
               EE[j].Print( fastout, "sequence_" + ToString( ids[j] ) );
          else if (NAMES_PARALLEL) EE[j].Print( fastout, names[j] );
          else EE[j].Print( fastout, names[ ids[j] ] );    }

     if (USE_QUAL)
     {    if ( !IsRegularFile( input + ".qualb" ) 
               && !IsRegularFile( input + ".qualp" ) )
          {    cout << "Warning: " << input << ".qualb doesn't exist.\n";
               return 0;    }
          vecqualvector Q;
          if ( IsRegularFile( input + ".qualp" ) )
          {    VecPQVec qualsp;
               qualsp.Read( input + ".qualp", ids, 0 );
               Q.resize( ids.size( ) );
               for ( int i = 0; i < ids.isize( ); i++ )
                    qualsp[i].unpack( &Q[i] );    }
          else Q.Read( input + ".qualb", ids, 0 );
          if (RC)
          for ( vecqvec::size_type i = 0; i < Q.size(); ++i )
	       Q[i].ReverseMe( );
          if ( !PRESERVE_ORDER ) Q.WriteAll( output + ".qualb" );
          else
          {    vecqualvector Qor;
               Qor.reserve( Q.size( ) );
               for ( vecqvec::size_type i = 0; i < Q.size( ); i++ )
               {    int j = BinPosition( ids, ids_unsorted[i] );
                    Qor.push_back( Q[j] );    }
               Qor.WriteAll( output + ".qualb" );    }    }    }
