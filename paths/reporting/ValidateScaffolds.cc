///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Fastavector.h"
#include "Superb.h"
#include "Vec.h"
#include "VecUtilities.h"

/**
 * ValidateScaffolds
 *
 * Check that the lengths reported by a superb file match the lengths of the
 * associated vec<fastavector>. It loads files
 *    <SCAFFOLDS><FASTA_EXTENSTION>.fasta
 *    <SCAFFOLDS>.superb
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( SCAFFOLDS );
  CommandArgument_String_OrDefault( FASTA_EXTENSION, ".contigs" )
  EndCommandArguments;

  String contigs_file = SCAFFOLDS + FASTA_EXTENSION + ".fasta";
  String scaffolds_file = SCAFFOLDS + ".superb";
     
  vec<fastavector> contigs;
  LoadFromFastaFile( contigs_file, contigs );
     
  vec<superb> scaffolds;
  ReadSuperbs( scaffolds_file, scaffolds );
  
  int n_tigs = 0;
  for (int ii=0; ii<scaffolds.isize( ); ii++) {
    n_tigs += scaffolds[ii].Ntigs( );
    for (int jj=0; jj<scaffolds[ii].Ntigs( ); jj++) {
      int tig = scaffolds[ii].Tig( jj );
      if ( tig >= (int)contigs.size( ) ) {
	cout << "s" << ii << "_" << jj << "." << scaffolds[ii].Ntigs( )
	     << " = " << scaffolds[ii].Len( jj ) << " bp,"
	     << "   invalid idx!\n";
	continue;
      }
      if ( (int)contigs[tig].size( ) != scaffolds[ii].Len( jj ) ) {
	cout << "s" << ii << "_" << jj << "." << scaffolds[ii].Ntigs( )
	     << " = " << scaffolds[ii].Len( jj ) << " bp,"
	     << "   c" << tig << " = " << contigs[tig].size( ) << " bp\n";
      }
    }
  }

  cout << "\n"
       << "n scaffolds:            " << scaffolds.size( ) << "\n"
       << "n contigs in fasta:     " << contigs.size( ) << "\n"
       << "n contigs in scaffolds: " << n_tigs << "\n"
       << endl;
  
}
