/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2012) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

const char *DOC =
  "Takes an IN_HEAD.contigs.efasta and converts it to a fastg file.";


#include "CoreTools.h"
#include "MainTools.h"
#include "system/System.h"

#include "efasta/EfastaTools.h"
#include "fastg/FastgTools.h"

int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandDoc(DOC);
  CommandArgument_String(IN);

  EndCommandArguments;

  // Load data from efasta 
  
  VecEFasta econtigs;
  vec<String> headers;
  if ( IsRegularFile( IN + ".efasta" ) ) {        
    cout << Date() << " : Loading efasta..." << endl;
    LoadEfastaIntoStrings( IN + ".efasta", econtigs, headers );

  } else{
    FatalErr(IN + ".efasta file not found");
  }

  cout << Date() << " : Writing out fastg ..." << endl;

//  fastg_meta FGM;
  Ofstream(out_g, IN + ".fastg");
//  out_g << FGM.GetFileHeader( IN ) << "\n";
  out_g << fastg_meta::GetFileHeader( IN ) << "\n";
  for (size_t it = 0; it < econtigs.size(); it++) {
    basefastg bfg( econtigs[it] );
    headfastg hfg( headers[it] );
    recfastg  rfg( hfg, bfg );
    rfg.Print( out_g );
  }
//  out_g << FGM.GetFileFooter() << "\n";
  out_g << fastg_meta::GetFileFooter() << "\n";
  
  cout << Date( ) << " : Done." << endl;

}
 

