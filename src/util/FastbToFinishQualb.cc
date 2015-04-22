///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "MainTools.h"
#include "Qualvector.h"

/**
 * FastbToFinishQualb
 *
 * Load a fastb and generate its finished-grade qualb analogue.
 *
 * BASE_DIR: where all files are
 * HEAD: head of fastb file (output is saved as <BASE_DIR>/<HEAD>.qualb)
 * QUAL: quality score for bases
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( BASE_DIR );
  CommandArgument_String( HEAD );
  CommandArgument_UnsignedInt_OrDefault( QUAL, 50 );
  EndCommandArguments;

  String bases_file = BASE_DIR + "/" + HEAD + ".fastb";
  String out_file = BASE_DIR + "/" + HEAD + ".qualb";
  
  vecbasevector bases;
  bases.ReadAll( bases_file );

  vecqualvector quals;

  quals.reserve( bases.size( ) );
  
  for (size_t ii=0; ii<bases.size( ); ii++) {
    quals.push_back(qualvector(bases[ii].size(), int(QUAL)));
  }
  
  quals.WriteAll( out_file );

}
