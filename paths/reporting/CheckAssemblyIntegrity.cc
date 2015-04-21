/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"

#include "Basevector.h"
#include "Fastavector.h"
#include "Superb.h"
#include "efasta/EfastaTools.h"
#include "paths/AssemblyIntegrity.h"

/**
 * CheckAssemblyIntegrity
 *
 * Formally validate an allpaths assembly integrity (check files are
 * in sync, etc.)
 */ 
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( HEAD );
  EndCommandArguments;
  
  String fasta_file = HEAD + ".contigs.fasta";
  String efasta_file = HEAD + ".contigs.efasta";
  String fastb_file = HEAD + ".contigs.fastb";
  String supers_file = HEAD + ".superb";
  
  vec<fastavector> *fastas = 0;
  vec<fastavector> loaded_fastas;
  if ( IsRegularFile( fasta_file ) ) {
    cout << Date( ) << ": loading fastas" << endl;
    LoadFromFastaFile( fasta_file, loaded_fastas );
    fastas = &loaded_fastas;
  }

  vecbvec *fastbs = 0;
  vecbvec loaded_fastbs;
  if ( IsRegularFile( fastb_file ) ) {
    cout << Date( ) << ": loading fastbs" << endl;
    loaded_fastbs.ReadAll( fastb_file );
    fastbs = &loaded_fastbs;
  }
  
  vecbvec *efastas = 0;
  vecbvec loaded_efastas;
  if ( IsRegularFile( efasta_file ) ) {
    cout << Date( ) << ": loading flattened efastas" << endl;
    LoadEfastaFlat( efasta_file, loaded_efastas );
    efastas = &loaded_efastas;
  }

  cout << Date( ) << ": loading superbs" << endl;
  vec<superb> supers;
  ReadSuperbs( supers_file, supers );

  cout << Date( ) << ": done loading\n" << endl;

  AssemblyIntegrity( &cout, &supers, fastas, fastbs, efastas );

  cout << Date( ) << ": done" << endl;  

}

