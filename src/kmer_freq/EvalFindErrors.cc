/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/*****************************************************************************
 *
 * EvalFindErrors: A tool to test the effectiveness of error correction as in
 * FindErrors.
 *
 * In order to run EvalFindErrors, you must have both the original and the
 * error-corrected versions of a set of reads, and they must be the same size
 * (i.e., uncorrectable reads have not yet been removed.)  You must have
 * aligned both of these read sets to a reference genome , which should be
 * high-quality (or the evaluation won't be very informative.)  The only reads
 * considered for error-correction effectivess are those which align uniquely
 * to the same location both before and after the error correction.
 *
 *
 * INPUT FILES:
 *
 * <GENOME>.fastb
 * <DIR>/<READS_ORIG>.fastb
 * <DIR>/<READS_EDIT>.fastb
 * <DIR>/<READS_ORIG>.qltout
 * <DIR>/<READS_EDIT>.qltout
 *
 *
 * Output goes to stdout.
 *
 *
 *
 * Josh Burton
 * September 2009
 *
 ****************************************************************************/


#include "Basevector.h"
#include "MainTools.h"
#include "feudal/BinaryStream.h"
#include "lookup/LookAlign.h"

int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( GENOME );
  CommandArgument_String( DIR );
  CommandArgument_String_OrDefault( READS_ORIG, "reads_orig" );
  CommandArgument_String_OrDefault( READS_EDIT, "reads_edit" );
  EndCommandArguments;
  
  // Load genome.
  cout << Date( ) << ": Loading genome..." << endl;
  vecbasevector genome( GENOME + ".fastb" );
  
  // Load reads.
  cout << Date( ) << ": Loading reads..." << endl;
  vecbasevector reads_orig( DIR + "/" + READS_ORIG + ".fastb" );
  vecbasevector reads_edit( DIR + "/" + READS_EDIT + ".fastb" );
  ForceAssertEq( reads_orig.size(), reads_edit.size() );
  ForceAssertEq( reads_orig.sumSizes( ), reads_edit.sumSizes( ) );
  size_t n_reads = reads_orig.size();
  
  vec<Bool> keep;
  BinaryReader::readFile( DIR + "/" + READS_EDIT + ".keep", &keep );
  
  
  
  // Load alignments - of both original and edited reads.
  cout << Date( ) << ": Loading alignments..." << endl;
  
  vec<look_align> aligns_orig, aligns_edit;
  vec< vec<int> > aligns_orig_index, aligns_edit_index;
  LoadLookAligns( DIR + "/" + READS_ORIG + ".qltout", aligns_orig, aligns_orig_index, n_reads );
  LoadLookAligns( DIR + "/" + READS_EDIT + ".qltout", aligns_edit, aligns_edit_index, n_reads );
  
  
  // Loop over all reads to examine their alignments.
  cout << Date( ) << ": Main algorithm!  Examining reads' alignments before and after error correction" << endl;
  int n_reads_examined = 0;
  int n_errors_orig_total = 0, n_errors_edit_total = 0;
  int n_reads_improved = 0, n_reads_worsened = 0;
  
  for ( size_t i = 0; i < n_reads; i++ ) {
    
    // Filter the reads.  We only want to examine reads whose alignments
    // we can be confident of.  In order to be acceptable:
    // 1. The read must not be marked for deletion by FindErrors.
    // 2. The read must align EXACTLY once, in both original and edited form.
    // 3. The read's alignment must be in the same place in both the original
    //    and edited form.
    // 4. The read's alignment must be indel-free - hence all errors are
    //    mismatches.
    if ( !keep[i] ) continue;
    if ( !aligns_orig_index[i].solo( ) ||
	 !aligns_edit_index[i].solo( ) ) continue;
    const look_align & la_orig = aligns_orig[ aligns_orig_index[i][0] ];
    const look_align & la_edit = aligns_edit[ aligns_edit_index[i][0] ];
    if ( la_orig.TargetId( )      != la_edit.TargetId( ) ) continue;
    if ( la_orig.StartOnTarget( ) != la_edit.StartOnTarget( ) ) continue;
    if ( la_orig.IsQueryRC( )     != la_edit.IsQueryRC( ) ) continue;
    if ( la_orig.indels > 0 || la_edit.indels > 0 ) continue;
    
    // Find the sequence on the reference where this read aligns.
    int read_length = reads_edit[i].size( );
    basevector target( genome[ la_orig.TargetId( ) ], la_orig.StartOnTarget( ), read_length );
    if ( la_orig.IsQueryRC( ) ) target.ReverseComplement( );
    
    
    // Manually count the number of errors in this alignment, both before
    // and after the read's error correction.  Because FindErrors only corrects
    // mutations and not indels, we can assume that any change in the number of 
    // errors is directly due to a correction made by FindErrors.
    int n_errors_orig = 0, n_errors_edit = 0;
    for ( int j = 0; j < read_length; j++ ) {
      if ( target[j] != reads_orig[i][j] ) n_errors_orig++;
      if ( target[j] != reads_edit[i][j] ) n_errors_edit++;
    }
    
    n_errors_orig_total += n_errors_orig;
    n_errors_edit_total += n_errors_edit;
    if ( n_errors_orig > n_errors_edit ) n_reads_improved++;
    if ( n_errors_orig < n_errors_edit ) n_reads_worsened++;
    
    
    n_reads_examined++;
  }
  
  double pct_errors_fixed = 100.0 * ( 1 - ( (double) n_errors_edit_total ) / n_errors_orig_total );
  double pct_reads_improved = 100.0 * ( (double) n_reads_improved ) / n_reads_examined;
  double pct_reads_worsened = 100.0 * ( (double) n_reads_worsened ) / n_reads_examined;
  
  
  // Report!
  cout << Date( ) << ": Examined " << n_reads_examined << " reads out of " << n_reads << " for errors." << endl;
  cout << "\tThese reads contained a total of " << n_errors_orig_total << " errors before error correction, and " << n_errors_edit_total << " after, for a reduction of " << pct_errors_fixed << "%." << endl;
  cout << "\t" << n_reads_improved << " (" << pct_reads_improved << "%) reads were improved, and " << n_reads_worsened << " (" << pct_reads_worsened << "%) reads were worsened." << endl;
  
  
  
  cout << Date( ) << ": Done!" << endl;
}
