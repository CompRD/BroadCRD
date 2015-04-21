///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Fastavector.h"
#include "ParseSet.h"

/**
 * SelectFastaByIds
 *
 * FASTA_IN: input fasta
 * FASTA_OUT: output fasta
 * SELECT: ids of selected entries (parsed with ParseIntSet)
 */
int main( int argc, char** argv )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( FASTA_IN );
  CommandArgument_String( FASTA_OUT );
  CommandArgument_String( SELECT );
  EndCommandArguments;

  vec<fastavector> fastas;
  vec<String> names;
  LoadFromFastaFile( FASTA_IN, fastas, names );

  vec<int> select;
  ParseIntSet( SELECT, select );

  ofstream out( FASTA_OUT.c_str( ) );
  for (size_t ii=0; ii<select.size( ); ii++) {
    int fasta_id = select[ii];
    const fastavector &fasta = fastas[fasta_id];
    const String &name = names[fasta_id];
    fasta.Print( out, name );
  }
  out.close( );
  
}
