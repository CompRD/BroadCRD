///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"

#include "Fastavector.h"
#include "Superb.h"
#include "paths/SaveScaffoldGraph.h"

/**
 * ReferenceAsAssembly
 *
 * Save a reference genome as an allpath assembly (one contig per
 * super).
 *
 * Input:
 *  <REF_HEAD>.fasta
 *
 * Output:
 *  <OUT_HEAD>.superb
 *  <OUT_HEAD>.assembly.fasta
 *  <OUT_HEAD>.contigs.fasta
 *  <OUT_HEAD>.contigs.fast{b,amb}
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( REF_HEAD );
  CommandArgument_String( OUT_HEAD );
  EndCommandArguments;
  
  String ref_bases = REF_HEAD + ".fasta";
  
  vec<fastavector> ref;
  LoadFromFastaFile( ref_bases, ref );
  
  vec<superb> supers( ref.size( ) );
  for (size_t ii=0; ii<ref.size( ); ii++)
    supers[ii].PlaceFirstTig( ii, ref[ii].size( ) );
  
  SaveScaffoldAssembly( OUT_HEAD, supers, ref, 0, true );
  
}
