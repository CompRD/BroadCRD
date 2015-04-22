///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"
#include "util/RunCommand.h"

/**
 * TestContigMapping
 *
 * Consistency test used by the path-based assisted assembly package.
 * Make sure that mapping works (where new contigs are mapped back to
 * some original contigs via a .orig.ids map file).
 *
 * FASTB_ORIG: full path name to original fastb file
 * HEAD_NEW: it loads <HEAD_NEW>.{fastb,orig_ids}
 * CONTIGS: if  not empty, load <HEAD_NEW>.<CONTIGS>.fastb 
 */
int main( int argc, char *argv[] )
{
  RunTime( );
 
  BeginCommandArguments;
  CommandArgument_String( FASTB_ORIG );
  CommandArgument_String( HEAD_NEW );
  CommandArgument_String_OrDefault( CONTIGS, "" );
  CommandArgument_Bool_OrDefault( VERBOSE, False );
  EndCommandArguments;
  
  // Dir and file names.
  String str_contigs = ( CONTIGS == "" ? "" : "." + CONTIGS );
  String new_fastb_file = HEAD_NEW + str_contigs + ".fastb";
  String map_file = HEAD_NEW + ".orig.ids";
  
  vec<String> needed;
  needed.push_back( FASTB_ORIG );
  needed.push_back( new_fastb_file );
  needed.push_back( map_file );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;

  // Load.
  cout << Date( ) << ": loading contigs" << endl;
  vecbvec orig_bases( FASTB_ORIG );
  vecbvec new_bases( new_fastb_file );

  cout << Date( ) << ": loading map " << endl;
  READ( map_file, vec<int>, orig_ids );
  
  // Files out of sync.
  if ( orig_ids.size( ) != new_bases.size( ) ) {
    cout << "The files orig.ids and contigs.fastb have different sizes\n"
	 << " size of orig.ids:      " << orig_ids.size( ) << "\n"
	 << " size of contigs.fastb: " << new_bases.size( ) << "\n"
	 << "\n"
	 << "TestContigMapping failed, exiting now.\n" << endl;
    return 1;
  }

  // Test mapping.
  int dotter = 100;
  cout << Date( ) << ": testing " << new_bases.size( ) << " new contigs";
  if ( ! VERBOSE ) cout << " (. = " << dotter << " contigs)";
  else cout << "\n";
  cout << endl;
  int n_passed = 0;
  for (int ii=0; ii<(int)new_bases.size( ); ii++) {
    if ( ( ! VERBOSE ) && ii % dotter == 0 ) Dot( cout, ii / dotter );

    int orig_id = orig_ids[ii] < 0 ? - orig_ids[ii] - 1 : orig_ids[ii];
    bool orig_rc = orig_ids[ii] < 0;

    if ( VERBOSE )
      cout << "c" << ii << " vs u" << orig_id << ( orig_rc ? "-\t" : "+\t" );
    
    if ( size_t( orig_id ) >= orig_bases.size( ) ) {
      if ( VERBOSE ) cout << "FAIL (orig_id out of bound)\n";
      continue;
    }

    bvec &newb = new_bases[ii];
    bvec origb_rc;
    if ( orig_rc ) {
      origb_rc = orig_bases[orig_id];
      origb_rc.ReverseComplement( );
    }
    bvec &origb = orig_rc ? origb_rc : orig_bases[orig_id];
    if ( newb != origb ) {
      if ( VERBOSE ) cout << "FAIL (bvecs do not match)\n";
      continue;
    }

    if ( VERBOSE ) cout << "ok\n";
    n_passed++;
  }
  if ( ! VERBOSE ) cout << endl;

  // Print result.
  double ratio = 1.0 - SafeQuotient( n_passed, (int)new_bases.size( ) );
  
  cout << "\n"
       << "RESULTS\n"
       << "\n"
       << "n_tested:  " << ToStringAddCommas( new_bases.size( ) ) << "\n"
       << "n_passed:  " << ToStringAddCommas( n_passed ) << "\n"
       << "% failed:  " << ToString( 100.0 * ratio, 4 ) << "\n"
       << "\n";
  
  // Done.
  cout << Date( ) << ": done" << endl;

}
