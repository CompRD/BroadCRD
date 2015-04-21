/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Basevector.h"
#include "LocsHandler.h"
#include "PairsHandler.h"
#include "SeqInterval.h"

/**
 * FragmentsInfo
 *
 * Load the given locs file and look at how fragments are placed on
 * the contigs.
 *
 * ID: id of contig
 * BEGIN: start on contig
 * END: end on contig
 * MAX_STRETCH: only show uniquely placed and not-too-streched fragments
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String( SUBDIR );
  CommandArgument_Int( ID );
  CommandArgument_Int( BEGIN );
  CommandArgument_Int( END );
  CommandArgument_Double_OrDefault( MAX_STRETCH, 5.0 );
  EndCommandArguments;
  
  // Dir and file names.
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/" + SUBDIR;

  String reads_file = run_dir + "/reads.fastb";
  String pairs_file = run_dir + "/reads.pairto";

  String locs_file = sub_dir + "/mergedcontigs_orig.locs";
  String contigs_file = sub_dir + "/mergedcontigs.fastb";

  // Load.
  int n_reads = MastervecFileObjectCount( reads_file );
  int n_contigs = MastervecFileObjectCount( contigs_file );

  cout << Date( ) << ": loading pairs" << endl;
  phandler pairs( n_reads, pairs_file );

  cout << Date( ) << ": loading locs" << endl;
  lhandler locs( n_reads, n_contigs, locs_file );

  cout << Date( ) << ": done loading\n" << endl;

  // Initial selection of ids of pairs.
  vec<int> select;
  int fpos = locs.FirstLoc( ID );
  if ( fpos > -1 ) {
    for (int ii=fpos; ii<locs.Size( ); ii++) {
      if ( locs[ii].Contig( ) != ID ) break;
      if ( locs[ii].StartOnContig( ) > END ) break;
      int pair_id = pairs.GetPairId( locs[ii].ReadId( ) );
      if ( pair_id < 0 ) continue;
      int begin = locs[ii].StartOnContig( );
      int end = 1 + locs[ii].StopOnContig( );
      bool begin_in = ( BEGIN <= begin && begin < END );
      bool end_in = ( BEGIN <= end && end < END );
      if ( ! ( begin_in || end_in ) ) continue;
      select.push_back( pair_id );
    }
  }

  sort( select.begin( ), select.end( ) );
  select.erase( unique( select.begin( ), select.end( ) ), select.end( ) );

  // Refine selection (store pair_id in interval_id).
  vec<seq_interval> fragments;
  fragments.reserve( select.size( ) );
  for (int ii=0; ii<select.isize( ); ii++) {
    int pair_id = select[ii];
    const read_pairing &pair = pairs[pair_id];

    int id1 = pair.id1;
    int id2 = pair.id2;
    int sep = pair.sep;
    int sd = pair.sd;

    const read_location *loc1 = locs.GetPlacement( id1 );
    const read_location *loc2 = locs.GetPlacement( id2 );
    if ( ! ( loc1 || loc2 ) ) continue;

    int cg1 = loc1->Contig( );
    int cg2 = loc2->Contig( );
    if ( cg1 != cg2 ) continue;

    bool fw1 = loc1->OrientationOnContig( ) == ForwardOr;
    bool fw2 = loc2->OrientationOnContig( ) == ForwardOr;
    if ( fw1 == fw2 ) continue;

    int osep = 0;
    if ( fw1 ) osep = loc2->StartOnContig( ) - loc1->StopOnContig( ) - 1;
    else osep = loc1->StartOnContig( ) - loc2->StartOnContig( ) - 1;

    double stretch = SafeQuotient( osep - sep, sd );
    if ( stretch > MAX_STRETCH ) continue;

    int beg = fw1 ? loc1->StartOnContig( ) : loc2->StartOnContig( ) ;
    int end = fw1 ? 1 + loc2->StopOnContig( ) : 1 + loc1->StopOnContig( );
    fragments.push_back( seq_interval( pair_id, cg1, beg, end ) );
  }

  sort( fragments.begin( ), fragments.end( ) );

  // Print output.
  for (int ii=0; ii<fragments.isize( ); ii++) {
    int pair_id = fragments[ii].IntervalId( );
    int id1 = pairs[pair_id].id1;
    int id2 = pairs[pair_id].id2;
    const read_location *loc1 = locs.GetPlacement( id1 );
    const read_location *loc2 = locs.GetPlacement( id2 );

    bool fw1 = loc1->OrientationOnContig( ) == ForwardOr;
    const read_location *fwloc = fw1 ? loc1 : loc2;
    const read_location *rcloc = fw1 ? loc2 : loc1;

    int sep = pairs[pair_id].sep;
    int sd = pairs[pair_id].sd;
    int osep = rcloc->StartOnContig( ) - fwloc->StopOnContig( ) - 1;
    double stretch = SafeQuotient( osep - sep, sd );

    cout << "p" << pair_id
	 << "   c" << ID
	 << " [" << fragments[ii].Begin( )
	 << ", " << fragments[ii].End( )
	 << ")_" << loc1->LengthOfContig( )
	 << "   r" << fwloc->ReadId( ) << "_fw"
	 << " r" << rcloc->ReadId( ) << "_rc"
	 << "   gsep (" << sep << " +/- " << sd << ")"
	 << "   osep " << osep
	 << "   stretch " << ToString( stretch, 2 )
	 << "\n";
  }
  cout << endl;

  // Done.
  cout << Date( ) << ": done" <<  endl;
  return 0;
}
