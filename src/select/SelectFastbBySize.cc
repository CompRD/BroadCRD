///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"

/**
 * SelectFastbBySize
 *
 * FASTB_IN: input fastb
 * FASTB_OUT: output fastb
 * MIN_LEN: select entries >= MIN_LEN
 */
int main( int argc, char** argv )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( FASTB_IN );
  CommandArgument_String( FASTB_OUT );
  CommandArgument_Int( MIN_LEN );
  EndCommandArguments;

  vecbvec bases_in( FASTB_IN );

  vecbvec bases_out;
  bases_out.reserve( bases_in.size( ) );
  for (size_t ii=0; ii<bases_in.size( ); ii++) {
    if( (int)bases_in[ii].size( ) < MIN_LEN ) continue;
    bases_out.push_back( bases_in[ii] );
  }

  bases_out.WriteAll( FASTB_OUT );

  cout << "Entries in:  " << ToStringAddCommas( bases_in.size( ) ) << "\n"
       << "Entries out: " << ToStringAddCommas( bases_out.size( ) ) << "\n"
       << endl;

}
