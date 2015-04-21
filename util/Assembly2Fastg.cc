/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2012) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

const char *DOC =
  "Takes an IN.contigs.efasta, IN.superb and converts it to a fastg file.";


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
  if ( IsRegularFile( IN + ".contigs.efasta" ) ) {        
    cout << Date() << " : Loading efasta..." << endl;
    LoadEfastaIntoStrings( IN + ".contigs.efasta", econtigs, headers );

  } else FatalErr(IN + ".contigs.efasta file not found");

  
  // Load superb
  

  
  vec<superb> scaffolds;
  if ( IsRegularFile( IN + ".superb" ) ) { 
    cout << Date() << " : Loading superb..." << endl; 
    ReadSuperbs(IN + ".superb", scaffolds);
  } else  FatalErr(IN + ".superb file not found");
  PRINT( scaffolds.size() );


  size_t n_scaffolds = scaffolds.size();
  size_t n_contigs = 0;
  for (size_t i = 0; i < scaffolds.size( ); i++)
    n_contigs += scaffolds[i].Ntigs( );
  cout << Date() << " : Assembly contains " << n_contigs << " contigs in " 
       << n_scaffolds << " scaffolds." << endl;
  

  cout << Date() << " : Writing out fastg ..." << endl;
  
//  fastg_meta FGM;
  Ofstream(out_g, IN + ".fastg");
//  out_g << FGM.GetFileHeader( IN ) << "\n";
  out_g << fastg_meta::GetFileHeader( IN ) << "\n";
  
  
  for (size_t si = 0; si < scaffolds.size( ); si++){
    basefastg bfg;
    for ( int ti = 0; ti < scaffolds[si].Ntigs(); ti++ ) {
      int tigId = scaffolds[si].Tig(ti);
      basefastg tfg( econtigs.at(tigId) );
      bfg += tfg;
      if ( ti < scaffolds[si].Ntigs() -1 ){
	basefastg gfg( scaffolds[si].Gap(ti), scaffolds[si].Dev(ti) );
	bfg += gfg;
      }
    }
    recfastg rfg( ToString(si), bfg );
    rfg.Print( out_g );    
  }
//  out_g << FGM.GetFileFooter() << "\n";
  out_g << fastg_meta::GetFileFooter() << "\n";

  cout << Date( ) << " : Done." << endl;
  
}

 

