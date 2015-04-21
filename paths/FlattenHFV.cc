///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"
#include "Fastavector.h"
#include "paths/HyperFastavector.h"
#include "paths/FlattenHyperFastavector.h"

/**
 * FlattenHFV
 *
 * A wrapper around FlattenHyperFastavector. 
 *
 * HYPER: full path name (an HyperFastavector)
 * HEAD_OUT: full path name of head of output (allpaths assembly) 
 */

int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( HYPER );
  CommandArgument_String( HEAD_OUT );
  EndCommandArguments;

  HyperFastavector theHFV( HYPER );
  
  FlattenHyperFastavector( cout, theHFV, HEAD_OUT );  
}
