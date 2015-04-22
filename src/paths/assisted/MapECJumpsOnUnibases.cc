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

void ReportMissingAligns( const size_t n_reads,
			  const String &descrip,
			  const vec< triple<int64_t,int64_t,int> > &aligns,
			  ostream &out );

/**
 * MapECJumpsOnUnibases
 *
 * Map error corrected jump reads onto given unibases (only perfect
 * alignments are found). This code has a large overlap with
 * MapECJumps, the main difference being that in here we align each
 * read to all unibases (hence, each read will own exactly two aligns,
 * one fw and one rc).
 *
 * INPUT:
 *   <RUN_DIR>/<READS>.unibases.k<K>
 *   <RUN_DIR>/<JUMPS>.fastb
 *   <RUN_DIR>/<LONG_JUMPS>.fastb
 *
 * OUTPUT:
 *   <RUN_DIR>/<JUMPS>.aligns         (if JUMPS != "" )
 *   <RUN_DIR>/<LONG_JUMPS>.aligns    (if LONG_JUMPS != "" )
 *
 * NB: at least one between JUMPS and LONG_JUMPS must be given (both
 *   is fine).
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( RUN_DIR );
  CommandArgument_String( READS );
  CommandArgument_String_OrDefault( JUMPS, "" );
  CommandArgument_String_OrDefault( LONG_JUMPS, "" );
  CommandArgument_Int_OrDefault( K, 96 );
  CommandArgument_UnsignedInt_OrDefault( NUM_THREADS, 0 );
  EndCommandArguments;

  // Check args.
  if ( JUMPS == "" && LONG_JUMPS == "" ) {
    cout << "Fatal error: neither jump nor long_jump reads given.\n" << endl;
    return 1;
  }

  // Thread control.
  NUM_THREADS = configNumThreads( NUM_THREADS );
  omp_set_num_threads( NUM_THREADS );

  // Dir and file names.
  String strK = ToString( K );

  String unibases_file = RUN_DIR + "/" + READS + ".unibases.k" + strK;

  String jumps_head = RUN_DIR + "/" + JUMPS;
  String jumps_file = jumps_head + ".fastb";
  String jaligns_file = jumps_head + ".aligns";

  String Jumps_head = RUN_DIR + "/" + LONG_JUMPS;
  String Jumps_file = Jumps_head + ".fastb";
  String Jaligns_file = Jumps_head + ".aligns";

  vec<String> needed;
  needed.push_back( unibases_file );
  if ( JUMPS != "" ) needed.push_back( jumps_file );
  if ( LONG_JUMPS != "" ) needed.push_back( Jumps_file );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;
  
  // Align jumps.
  if ( JUMPS != "" ) {
    const size_t n_reads = MastervecFileObjectCount( jumps_file );
    const String descrip = "jump";

    cout << Date( ) << ": aligning jump reads" << endl;
    vec< triple<int64_t,int64_t,int> > jaligns;
    SearchFastb2( jumps_file, unibases_file, K, &jaligns );
    
    cout << Date( ) << ": saving aligns" << endl;
    BinaryWriter::writeFile( jaligns_file, jaligns );
    
    ReportMissingAligns( n_reads, descrip, jaligns, cout );
  }

  // Align long jumps.
  if ( LONG_JUMPS != "" ) {
    const size_t n_reads = MastervecFileObjectCount( Jumps_file );
    const String descrip = "long jump";
    
    cout << Date( ) << ": aligning long jump reads" << endl;
    vec< triple<int64_t,int64_t,int> > Jaligns;
    SearchFastb2( Jumps_file, unibases_file, K, &Jaligns );

    cout << Date( ) << ": saving " << Jaligns.size( ) << " aligns" << endl;
    BinaryWriter::writeFile( Jaligns_file, Jaligns );

    ReportMissingAligns( n_reads, descrip, Jaligns, cout );
  }
  
  // Done.
  cout << Date( ) << ": MapECJumpsOnUnibases done" << endl;
  
}

/**
 * ReportMissingAligns
 *
 * Report number of reads with no aligns found.
 */
void ReportMissingAligns( const size_t n_reads,
			  const String &descrip,
			  const vec< triple<int64_t,int64_t,int> > &aligns,
			  ostream &out )
{
  const size_t n_aligns = aligns.size( );

  vec<int> used( n_reads, 0 );
  for (size_t ii=0; ii<n_aligns; ii++)
    ( used[aligns[ii].first] )++;
  
  size_t n_missing = 0;
  size_t n_bad_other = 0;
  for (size_t ii=0; ii<used.size( ); ii++) {
    if ( used[ii] != 2 ) {
      if ( used[ii] == 0 ) n_missing++;
      else n_bad_other++;
    }
  }

  if ( n_missing < 1 && n_bad_other < 1 ) return;

  out << "\n"
      << "WARNING! Problems found aligning " + descrip << " reads:\n"
      << "\n"
      << "reads in input:  " << ToStringAddCommas( n_reads ) << "\n"
      << "no align found:  " << ToStringAddCommas( n_missing ) << "\n"
      << "other problem:   " << ToStringAddCommas( n_bad_other ) << "\n"
      << endl;
}

