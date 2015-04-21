///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "VecUtilities.h"
#include "graph/Digraph.h"
#include "lookup/LookAlign.h"
#include "paths/PdfEntry.h"
#include "paths/reporting/ReftigUtils.h"
#include "util/RunCommand.h"

/**
 * ShowCandidateGaps
 *
 * Align unibases to a reference, and show gaps between unibases
 * longer than MIN_LEN bases, with CN <= MAX_CN (skip gaps >= MAX_GAP
 * bases). Output is sent to cout (some files are cached to disk).
 *
 * WARNING: it will only show the aligns of unibases aligned fw on the
 * reference.
 *
 * REF_HEAD: it loads <REF_HEAD>.lookup
 * READS: it loads ../<READS>.unibases.k<K>
 * MAX_CN: show gaps between unibases with CN <= MAX_CN
 * MIN_LEN: min length (in bases) of unibases to show
 * MAX_GAP: skip gaps >= MAX_GAP
 * FORCE: do not used cached alignments
 */
int main(int argc, char *argv[])
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_Int( K );
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String( REF_HEAD );
  CommandArgument_String_OrDefault( READS, "extended" );
  CommandArgument_Int_OrDefault( MAX_CN, 1 );
  CommandArgument_Int_OrDefault( MIN_LEN, 1000 );
  CommandArgument_Int_OrDefault( MAX_GAP, 5000 );
  CommandArgument_Bool_OrDefault( FORCE, False );
  EndCommandArguments;

  // Dir and file names.
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String out_dir = run_dir + "/" + READS + ".ShowCandidateGaps";

  String strK = ToString( K );
  String lookup_file = REF_HEAD + ".lookup";
  String head_in = run_dir + "/" + READS;
  String unibases_file = head_in + ".unibases.k" + strK;
  String copynum_file = head_in + ".unipaths.predicted_count.k" + strK;

  String aligns_file = out_dir + "/unibases.qltout";
  String sel_unibases_file = out_dir + "/sel_unibases.fastb";
  
  Mkpath( out_dir );
  
  // These files must exist.
  vec<String> needed;
  needed.push_back( lookup_file );
  needed.push_back( unibases_file );
  needed.push_back( copynum_file );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;
  
  // Select unibases.
  if ( FORCE || ! IsRegularFile( sel_unibases_file ) ) {
    cout << Date( ) << ": loading copy numbers" << endl;
    vec<int> copynum;
    LoadCopyNumbers( copynum_file, copynum );
    
    cout << Date( ) << ": loading unibases" << endl;
    vecbvec unibases( unibases_file );
    
    cout << Date( ) << ": selecting unibases" << endl;
    size_t n_sel = 0;
    vecbvec sel_unibases( unibases.size( ) );
    for (size_t ii=0; ii<unibases.size( ); ii++) {
      if ( unibases[ii].size( ) < size_t( MIN_LEN ) ) continue;
      if ( copynum[ii] > MAX_CN ) continue;
      sel_unibases[ii] = unibases[ii];
      n_sel++;
    }

    cout << Date( ) << ": saving " << n_sel << " selected unibases" << endl;
    sel_unibases.WriteAll( sel_unibases_file );
  }

  // Get aligns (aligns are cached to disk), and sort them.
  vec<look_align> all_aligns;
  GetAlignsFast( K, sel_unibases_file, lookup_file, aligns_file, all_aligns,
		 ! FORCE, out_dir );
  
  // Remove rc aligns, sort fw aligns.
  cout << Date( ) << ": sorting fw aligns" << endl;
  vec<look_align> aligns;
  aligns.reserve( all_aligns.size( ) );
  for (size_t ii=0; ii<all_aligns.size( ); ii++)
    if ( all_aligns[ii].IsQueryFW( ) )
      aligns.push_back( all_aligns[ii] );

  order_lookalign_TargetBegin sorter;
  sort( aligns.begin( ), aligns.end( ), sorter );
  
  // Show gaps.
  vec< vec<String> > table;

  vec<String> line = MkVec( ToString( "uni_left" ),
			    ToString( "gap" ),
			    ToString( "t_id" ),
			    ToString( "t_range" ),
			    ToString( "uni_right" ) );
  table.push_back( line );

  for (size_t ii=1; ii<aligns.size( ); ii++) {
    const look_align &left_al = aligns[ii-1];
    const look_align &right_al = aligns[ii];
    
    if ( left_al.target_id != right_al.target_id ) continue;

    int left_end = left_al.a.Pos2( );
    int right_beg = right_al.a.pos2( );
    bool overlap = ( right_beg < left_end );
    if ( right_beg - left_end > MAX_GAP ) continue;

    int len_al1 = left_al.a.Pos2( ) - left_al.a.pos2( );
    int len_al2 = left_al.query_length;
    if ( len_al1 < len_al2 ) continue;

    len_al1 = right_al.a.Pos2( ) - right_al.a.pos2( );
    len_al2 = right_al.query_length;
    if ( len_al1 < len_al2 ) continue;

    line.clear( );
    line.push_back( ToString( left_al.query_id )
		    + " [" + ToString( left_al.query_length ) + "]" );
    line.push_back( ToString( right_beg - left_end ) );
    line.push_back( ToString( left_al.target_id ) );
    line.push_back( "[" + ToString( overlap ? right_beg : left_end )
		    + ", " + ToString( overlap ? left_end : right_beg )
		    + ")" );
    line.push_back( ToString( right_al.query_id )
		    + " [" + ToString( right_al.query_length ) + "]" );
    table.push_back( line );
  }
  
  cout << "\n"
       << "LEGEND\n"
       << "\n"
       << " uni_left: id of left unibase [length of unibase]\n"
       << " gap:      size of gap (overlap if <0)\n"
       << " t_id:     id of target\n"
       << " t_range:  range of gap on target\n"
       << " un_right: id of right unibase [length of unibase]\n"
       << "\n";

  cout << "GAPS ON TARGET\n\n";
  PrintTabular( cout, table, 3, "rrrrr" );
  cout << endl;
  
  // Done.
  cout << Date( ) << ": done" << endl;
  
}
