/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2012) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

const char *DOC =
  "Takes an IN_HEAD.contigs.<fastg|efasta|fastb|fasta> file and it's associated IN_HEAD.superb or"
  "IN_HEAD.assembly.fastg"
  "file and generates the following new files:\n"
  "OUT_HEAD.contigs.{fastg,efasta,fasta,fastb}\n"
  "OUT_HEAD.assembly.{fastg,efasta,fasta}\n"
  "OUT_HEAD.superb";


#include "CoreTools.h"
#include "MainTools.h"
#include "system/System.h"
#include "Fastavector.h"
#include "efasta/EfastaTools.h"
#include "fastg/FastgTools.h"
#include "paths/AssemblyCleanupTools.h"

int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandDoc(DOC);
  CommandArgument_String(IN_HEAD);
  CommandArgument_String(OUT_HEAD);
  CommandArgument_String_OrDefault_Doc(CONTIG_PREFIX, "contig_",
    "Add this string to the contig number to create the contig name");
  CommandArgument_Bool_OrDefault_Doc(REORDER, False,
    "If True, then reorder contigs to match scaffolds - removing any orphan contigs in the process.");

  EndCommandArguments;

  // Load data from efasta or fastb or fasta
  
  VecEFasta econtigs;
  vec<recfastg> fcontigs;
  vec<recfastg> fscaffolds;
  vec<superb> scaffolds;
  
  if ( ! IsRegularFile( IN_HEAD + ".assembly.fastg" ) ) {


    if ( ! IsRegularFile( IN_HEAD + ".superb" ) )
      FatalErr("Missing file: " + IN_HEAD + ".superb" );
    cout << Date() << " : Loading superb..." << endl; 
    ReadSuperbs(IN_HEAD + ".superb", scaffolds);
    

    if ( IsRegularFile( IN_HEAD + ".contigs.fastg" ) ) {               // Use fastg
      
      cout << Date() << " : Loading .contigs.fastg ..." << endl;
      LoadFastg( IN_HEAD + ".contigs.fastg", fcontigs );
      econtigs.resize( fcontigs.size() );
      efasta tmp;
      for ( size_t i = 0; i < fcontigs.size(); i++ ){
	ForceAssert( fcontigs[i].IsGapless() );
	tmp.clear();
	fcontigs[i].AsEfasta( tmp ,fastg_meta::MODE_2);
	econtigs[i] = tmp;
      }
      

    } else if ( IsRegularFile( IN_HEAD + ".contigs.efasta" ) ) {       // Use efasta
      
      cout << Date() << " : Loading efasta..." << endl;
      LoadEfastaIntoStrings( IN_HEAD + ".contigs.efasta", econtigs );
      fcontigs.resize( econtigs.size() );
      for ( size_t i = 0; i < econtigs.size(); i++ )
	fcontigs[i] = recfastg( ToString(i), econtigs[i] );

    } else if (IsRegularFile( IN_HEAD + ".contigs.fasta" ) ) {         // Use fasta
      
      cout << Date() << " : Loading fasta..." << endl;
      vec<fastavector> fastas;
      LoadFromFastaFile(  IN_HEAD + ".contigs.fasta", fastas );
      econtigs.resize( fastas.size() );
      fcontigs.resize( fastas.size() );
      for ( size_t i = 0; i < fastas.size(); i++ ){
	econtigs[i] = efasta( fastas[i] );
	fcontigs[i] = recfastg( ToString(i), econtigs[i] );
      }
    } else if (IsRegularFile( IN_HEAD + ".contigs.fastb" ) ) {        // Use fastb
      
      cout << Date() << " : Loading fastb..." << endl;
      vecbasevector bases( IN_HEAD + ".contigs.fastb");
      econtigs.resize( bases.size() );
      fcontigs.resize( bases.size() );
      for ( size_t i = 0; i < bases.size( ); i++ ){
	econtigs[i] = fastavector(bases[i]);
	fcontigs[i] = recfastg( ToString(i), econtigs[i] );
      }
    } else 
      FatalErr("Could not find a fastg, efasta, fastb or fasta file: " + IN_HEAD + ".contigs.<efasta|fastb|fasta>" );

    fscaffolds.resize( scaffolds.size() );
    for ( size_t is = 0; is < scaffolds.size(); is++ ){
      basefastg fbases( scaffolds[is], fcontigs );
      fscaffolds[is] = recfastg( ToString( is ), fbases );
    }
  
  } else {
    
    cout << Date() << " : Loading " << IN_HEAD << ".assembly.fastg ..." << endl;
    LoadFastg( IN_HEAD + ".assembly.fastg", fscaffolds );    
    PRINT( fscaffolds.size() );

    int lastTid = -1;
    for ( size_t i = 0; i < fscaffolds.size(); i++ ){
      superb s_c;
      vec<recfastg> fcontigs_c;
      fscaffolds[i].AsScaffold( s_c, fcontigs_c, lastTid );
      scaffolds.push_back( s_c );
      for ( size_t ifc = 0; ifc < fcontigs_c.size(); ifc++ ){
	fcontigs.push_back( fcontigs_c[ifc] );
	ForceAssert( fcontigs_c[ifc].IsGapless() );
	efasta ef;
	fcontigs_c[ifc].AsEfasta( ef,fastg_meta::MODE_2 );
	econtigs.push_back( ef );
	ForceAssertEq( fcontigs.back().Length1(), econtigs.back().Length1() );
	ForceAssertEq( fcontigs.back().MinLength(fastg_meta::MODE_2), econtigs.back().MinLength() );
	ForceAssertEq( fcontigs.back().MaxLength(fastg_meta::MODE_2), econtigs.back().MaxLength() );
      }
	
    }
    PRINT( lastTid );
  } 

  ForceAssertEq( scaffolds.size(), fscaffolds.size() );

  size_t n_scaffolds = scaffolds.size();
  size_t n_contigs = 0;
  for (size_t i = 0; i < scaffolds.size( ); i++)
    n_contigs += scaffolds[i].Ntigs( );
  cout << Date() << " : Assembly contains " << n_contigs << " contigs in " 
       << n_scaffolds << " scaffolds." << endl; 


  // Sanity Check
  if (n_contigs != fcontigs.size())
    cout << Date() << " : WARNING: superb and efasta contig count inconsistency found" << endl;

  ForceAssertEq( fcontigs.size(), econtigs.size() );
  for ( size_t i = 0; i < fcontigs.size(); i++ ){
    if ( econtigs[i].size() == 0u ){
      cout << "Empty Efasta string\n";
      PRINT( fcontigs[i].bases() );
      PRINT( econtigs[i] );
      FatalErr("Empty efasta string");
    }
    ForceAssert( fcontigs[i].IsGapless() );
    ForceAssertEq( fcontigs[i].Length1(), econtigs[i].Length1() );
    ForceAssertEq( fcontigs[i].MinLength(fastg_meta::MODE_2), econtigs[i].MinLength() );
    ForceAssertEq( fcontigs[i].MaxLength(fastg_meta::MODE_2), econtigs[i].MaxLength() );
    
  }

  ForceAssertEq( fcontigs.size(), econtigs.size() );


  if (REORDER) {

    // Renumber contigs according to scaffold position and reorder fastg accordingly

    vec<uint32_t> used(fcontigs.size(),0);
    vec<recfastg> fcontigs_reorder;

    cout << Date() << " : Re-ordering contigs to match scaffold" << endl;

    fcontigs_reorder.resize(n_contigs);
    size_t new_index = 0;
    for (size_t scaffold_index = 0; scaffold_index < n_scaffolds; scaffold_index++) {
      for (size_t tig_index = 0; tig_index < static_cast<size_t>(scaffolds[scaffold_index].Ntigs( ) ); tig_index++) {
	size_t old_index = scaffolds[scaffold_index].Tig(tig_index);
	scaffolds[scaffold_index].SetTig(tig_index, new_index);
	fcontigs_reorder[new_index] = fcontigs[old_index];
	used[old_index]++;
	new_index++;
      }
    }
    
    fcontigs = fcontigs_reorder;

    for ( size_t is = 0; is < fscaffolds.size(); is++ ){
      basefastg fbases( scaffolds[is], fcontigs );
      fscaffolds[is] = recfastg( ToString( is ), fbases );
    }
  }


  // Check integrity of scaffolds
  Assembly assembly_in( scaffolds, fcontigs );
  assembly_in.check_integrity();


  
  // Write output
  
  cout << Date() << " : Writing new files to : " 
       << OUT_HEAD << ".*" << endl;
  
  vec<FastaVec> flattened_fasta(fcontigs.size());
  vecbasevector flattened_fastb(fcontigs.size());
  VecEFasta     flattened_efasta(fcontigs.size());
  
  {
    
    Ofstream(out_e, OUT_HEAD + ".contigs.efasta");
    Ofstream(out_a, OUT_HEAD + ".contigs.fasta");
    Ofstream(out_g, OUT_HEAD + ".contigs.fastg");


//    fastg_meta FGM;
//    out_g << FGM.GetFileHeader( OUT_HEAD + ".contigs.fastg" ) << "\n";
    out_g << fastg_meta::GetFileHeader( OUT_HEAD + ".contigs.fastg" ) << "\n";

    efasta tmp;
    for (size_t it = 0; it < fcontigs.size(); it++) {
      fcontigs[it].AsFasta(flattened_fasta[it]);
      fcontigs[it].AsFastb(flattened_fastb[it]);
      tmp.clear();
      fcontigs[it].AsEfasta(tmp,fastg_meta::MODE_2);
      flattened_efasta[it] = tmp;
      flattened_efasta[it].Print(out_e, CONTIG_PREFIX + ToString(it));
      flattened_fasta[it].Print(out_a, CONTIG_PREFIX + ToString(it));
      fcontigs[it].Print(out_g, CONTIG_PREFIX + ToString(it));
    }
//    out_g << FGM.GetFileFooter() << "\n";
    out_g << fastg_meta::GetFileFooter() << "\n";

    out_a.close();
    out_e.close();
    out_g.close();


    flattened_fastb.WriteAll( OUT_HEAD + ".contigs.fastb" );  
    
    
    WriteFastg( OUT_HEAD + ".assembly.fastg", fscaffolds );
    WriteScaffoldedEFasta( OUT_HEAD + ".assembly.efasta", flattened_efasta, scaffolds );
    WriteScaffoldedFasta( OUT_HEAD + ".assembly.fasta", flattened_fasta, scaffolds );
//    WriteFastaEfasta(OUT_HEAD + ".assembly.fasta",OUT_HEAD + ".assembly.efasta",fscaffolds);    
    WriteSuperbs( OUT_HEAD + ".superb", scaffolds );
    WriteSummary( OUT_HEAD + ".summary", scaffolds );
  }
  
  //Assembly assembly( scaffolds, flattened_efasta );
  Assembly assembly( scaffolds, fcontigs );
  assembly.check_integrity();

  cout << Date( ) << " : Done." << endl;
}

