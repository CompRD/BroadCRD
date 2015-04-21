/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"
#include "paths/Alignlet.h"
#include "util/RunCommand.h"
#include "util/SearchFastb2Core.h"
#include <omp.h>
// MakeDepend: library OMP

/**
 * MapECJumps
 *
 * Map error corrected jump reads onto given contigs (only perfect
 * alignments are found).
 *
 * CONTIGS_HEAD: head of contigs
 * RUN_DIR: where jump/long_jump reads are (fastb needed)
 * CONTIGS_TAIL: extension of contigs file (eg, ".fastb")
 * JUMPS: head of jump (ec) reads
 * LONG_JUMPS: head of long (ec) jumps reads 
 * K: kmer size (needed by SearchFastb2Core)
 * NUM_THREADS: needed by GetAlignsFast (if HEAD_REF is given)
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( CONTIGS_HEAD );
  CommandArgument_String( RUN_DIR );
  CommandArgument_String_OrDefault( CONTIGS_TAIL, ".fastb" );
  CommandArgument_String_OrDefault( JUMPS, "" );
  CommandArgument_String_OrDefault( LONG_JUMPS, "" );
  CommandArgument_Int_OrDefault( K, 96 );
  CommandArgument_UnsignedInt_OrDefault( NUM_THREADS, 0 );
  EndCommandArguments;

  // Dir and file names.
  String contigsF = CONTIGS_HEAD + CONTIGS_TAIL;
  String jumpF =  JUMPS == "" ? "" : RUN_DIR + "/" + JUMPS + ".fastb";
  String JumpF = LONG_JUMPS == "" ? "" : RUN_DIR + "/" + LONG_JUMPS + ".fastb";
  
  String jalignsF = CONTIGS_HEAD + ".jump.aligns";
  String JalignsF = CONTIGS_HEAD + ".Jump.aligns";
  String jidxF = CONTIGS_HEAD + ".jump.idx";
  String JidxF = CONTIGS_HEAD + ".Jump.idx";
 
  vec<String> needed;
  needed.push_back( contigsF );
  if ( jumpF != "" ) needed.push_back( jumpF );
  if ( JumpF != "" ) needed.push_back( JumpF );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;
  
  if ( jumpF == "" && JumpF == "" ) {
    cout << "Fatal error: neither jump nor long_jump reads given.\n" << endl;
    return 1;
  }
  
  // Thread control.
  NUM_THREADS = configNumThreads( NUM_THREADS );
  omp_set_num_threads( NUM_THREADS );

  // Load.
  cout << Date( ) << ": loading contigs" << endl;
  vecbvec contigs( contigsF );

  // Align reads.
  if ( jumpF != "" ) {
    cout << Date( ) << ": aligning jump reads" << endl;
    vec< triple<int64_t,int64_t,int> > jaligns;
    SearchFastb2( jumpF, contigsF, K, &jaligns );

    cout << Date( ) << ": saving " << jaligns.size( ) << " aligns" << endl;
    vecbvec jumps( jumpF );
    size_t n_jumps = jumps.size( );
    vec<alignlet> jalignlets;
    SF2AlignsToAlignlets( jumps, contigs, jaligns, jalignlets );
    BinaryWriter::writeFile( jalignsF, jalignlets );

    cout << Date( ) << ": saving index" << endl;
    vec<int> idx( n_jumps, -2 );
    vec<int> placs( n_jumps, 0 );
    for (size_t ii=0; ii<jaligns.size( ); ii++) {
      int id = jaligns[ii].first;
      idx[id] = ( idx[id] == -2 ? ii : -1 );
      placs[id] += 1;
    }
    BinaryWriter::writeFile( jidxF, idx );    
    
    size_t n_uniq = 0;
    size_t n_mult = 0;
    for (size_t ii=0; ii<n_jumps; ii++) {
      if ( placs[ii] == 1 ) n_uniq++;
      else if ( placs[ii] > 1 ) n_mult++;
    }
    size_t n_unplaced = n_jumps - n_uniq - n_mult;
    double ratio_unplaced = SafeQuotient( n_unplaced, n_jumps );
    double ratio_uniq = SafeQuotient( n_uniq, n_jumps );
    double ratio_mult = SafeQuotient( n_mult, n_jumps );

    vec< vec<String> > table;
    
    table.push_back( MkVec( String( "  unplaced" ),
			    ToStringAddCommas( n_unplaced ),
			    ToString( 100. * ratio_unplaced, 1 ) ) );
		     
    table.push_back( MkVec( String( "unique" ),
			    ToStringAddCommas( n_uniq ),
			    ToString( 100. * ratio_uniq, 1 ) ) );

    table.push_back( MkVec( String( "multiple" ),
			    ToStringAddCommas( n_mult ),
			    ToString( 100. * ratio_mult, 1 ) ) );

    cout << "\nPLACEMENT OF JUMP READS\n\n";
    PrintTabular( cout, table, 3, "rrr" );
    cout << endl;
  }

  if ( JumpF != "" ) {
    cout << Date( ) << ": aligning long jump reads" << endl;
    vec< triple<int64_t,int64_t,int> > Jaligns;
    SearchFastb2( JumpF, contigsF, K, &Jaligns );

    cout << Date( ) << ": saving " << Jaligns.size( ) << " aligns" << endl;
    vecbvec Jumps( JumpF );
    size_t n_Jumps = Jumps.size( );
    vec<alignlet> Jalignlets;
    SF2AlignsToAlignlets( Jumps, contigs, Jaligns, Jalignlets );
    BinaryWriter::writeFile( JalignsF, Jalignlets );
    
    cout << Date( ) << ": saving index" << endl;
    vec<int> idx( n_Jumps, -2 );
    vec<int> placs( n_Jumps, 0 );
    for (size_t ii=0; ii<Jaligns.size( ); ii++) {
      int id = Jaligns[ii].first;
      idx[id] = ( idx[id] == -2 ? ii : -1 );
      placs[id] += 1;
    }
    BinaryWriter::writeFile( jidxF, idx );
    
    size_t n_uniq = 0;
    size_t n_mult = 0;
    for (size_t ii=0; ii<n_Jumps; ii++) {
      if ( placs[ii] == 1 ) n_uniq++;
      else if ( placs[ii] > 1 ) n_mult++;
    }
    size_t n_unplaced = n_Jumps - n_uniq - n_mult;
    double ratio_unplaced = SafeQuotient( n_unplaced, n_Jumps );
    double ratio_uniq = SafeQuotient( n_uniq, n_Jumps );
    double ratio_mult = SafeQuotient( n_mult, n_Jumps );

    vec< vec<String> > table;
    
    table.push_back( MkVec( String( "  unplaced" ),
			    ToStringAddCommas( n_unplaced ),
			    ToString( 100. * ratio_unplaced, 1 ) ) );
		     
    table.push_back( MkVec( String( "unique" ),
			    ToStringAddCommas( n_uniq ),
			    ToString( 100. * ratio_uniq, 1 ) ) );

    table.push_back( MkVec( String( "multiple" ),
			    ToStringAddCommas( n_mult ),
			    ToString( 100. * ratio_mult, 1 ) ) );

    cout << "\nPLACEMENT OF LONG JUMP READS\n\n";
    PrintTabular( cout, table, 3, "rrr" );
    cout << endl;
  }
  
  // Done.
  cout << Date( ) << ": MapECJumps done" << endl;
  
}

