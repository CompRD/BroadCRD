///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "paths/ReadLoc.h"
#include "feudal/BinaryStream.h"

/**
 * BuildLocalDataset
 *
 * Select a range of contigs from a given scaffold, and collect the
 * ids of all the reads involved and, optionally, their partners. Save
 * these to files.
 *
 * Special mode SINGLE_GAP: select all reads in the two contigs
 * adjacent a gap, plus all the reads that could fall into the gap. In
 * other words, do not select unassembled reads that would fall either
 * before the left or after the right contig.
 * 
 * HEAD_OUT: full path name of output files (head)
 * SCAFFOLD_ID: id of selected scaffold
 * CGPOS_BEGIN: position of first selected contig in scaffold
 * CGPOS_COUNT: how many contigs to select
 * ADD_MATES: save also ids of partners (defaults to true!)
 * SINGLE_GAP: special case for a single gap between two contigs
 */
int main(int argc, char *argv[])
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( ALIGNS_IN )
  CommandArgument_String( SCAFFOLDS_IN )
  CommandArgument_String( HEAD_OUT );
  CommandArgument_Int( SCAFFOLD_ID );
  CommandArgument_Int( CGPOS_BEGIN );
  CommandArgument_Int( CGPOS_COUNT );
  CommandArgument_Bool_OrDefault( ADD_MATES, True );
  CommandArgument_Bool_OrDefault( SINGLE_GAP, False );
  EndCommandArguments;

  // Check some arguments.
  if ( SINGLE_GAP && CGPOS_COUNT != 2 ) {
    cout << "Fatal error: CGPOS_COUNT must be 2, if SINGLE_GAP=True.\n" << endl;
    return 0;
  }

  // Dir and file names.
  String run_dir = Dirname( ALIGNS_IN );
  String readlocs_head = ALIGNS_IN.SafeBefore( ".readlocs" );

  String out_frag_file = HEAD_OUT + ".frag.ids";
  String out_jump_file = HEAD_OUT + ".jump.ids";
  String out_long_jump_file = HEAD_OUT + ".long_jump.ids";

  // Make sure path to HEAD_OUT exists.
  if ( HEAD_OUT.Contains( "/" ) ) {
    String base_dir = HEAD_OUT.SafeBeforeLast( "/" );
    Mkpath( base_dir );
  }
  
  // Load.
  cout << Date( ) << ": loading supers" << endl;
  vec<superb> supers;
  ReadSuperbs( SCAFFOLDS_IN, supers );
    
  // Select contigs (sorted).
  longlong tot_len = 0;
  vec<int> selected_tigs;
  {
    const superb &sup = supers[SCAFFOLD_ID];
    int begin = CGPOS_BEGIN;
    int end = CGPOS_BEGIN + CGPOS_COUNT;
    if ( begin < 0 || end > sup.Ntigs( ) ) {
      cout << "\tFatal error: scaffold " << SCAFFOLD_ID
	   << " has only " << sup.Ntigs( )
	   << " contigs, while you selected\nCGPOS_BEGIN=" << CGPOS_BEGIN
	   << ", CGPOS_COUNT=" << CGPOS_COUNT
	   << "\nExit.\n" << endl;
      return 1;
    }
    selected_tigs.reserve( CGPOS_COUNT );
    for (int ii=begin; ii<end; ii++) {
      selected_tigs.push_back( sup.Tig( ii ) );
      tot_len += sup.Len( ii );
    }
    sort( selected_tigs.begin( ), selected_tigs.end( ) );
  }
  
  // Parse contigs.
  read_locs_on_disk locs_file( readlocs_head, run_dir );
  
  cout << Date( ) << ": parsing "
       << CGPOS_COUNT << " contigs from s"
       << SCAFFOLD_ID << " ("
       << ToStringAddCommas( tot_len ) << " bases)\n"
       << endl;
  
  vec<int> frag_ids;
  vec<int> jump_ids;
  vec<int> long_jump_ids;

  // Standard mode.
  if ( ! SINGLE_GAP ) {
    for (int it=0; it<selected_tigs.isize( ); it++ ) {
      Dot( cout, it );
      int tig_id = selected_tigs[it];
      vec<read_loc> locs;
      locs_file.LoadContig( tig_id, locs );
      for (int jj=0; jj<locs.isize( ); jj++) {
	const read_loc& rl = locs[jj];
	if ( rl.Frag( ) ) {
	  frag_ids.push_back( rl.ReadId( ) );
	  if ( ADD_MATES ) frag_ids.push_back( rl.PartnerReadId( ) );
	}
	if ( rl.Jump( ) ) {
	  jump_ids.push_back( rl.ReadId( ) );
	  if ( ADD_MATES ) jump_ids.push_back( rl.PartnerReadId( ) );
	}
	if ( rl.LongJump( ) ) {
	  long_jump_ids.push_back( rl.ReadId( ) );
	  if ( ADD_MATES ) long_jump_ids.push_back( rl.PartnerReadId( ) );
	}
      }
    }
  }

  // SINGLE_GAP mode.
  if ( SINGLE_GAP ) {
    for (int pass=0; pass<2; pass++) {
      Dot( cout, pass );
      int tig_id = selected_tigs[pass];
      vec<read_loc> locs;
      locs_file.LoadContig( tig_id, locs );
      for (int jj=0; jj<locs.isize( ); jj++) {
	const read_loc& rl = locs[jj];
	bool pick_mate = ( pass == 0 && rl.Fw( ) ) || ( pass == 1 && rl.Rc( ) );
	if ( rl.Frag( ) ) {
	  frag_ids.push_back( rl.ReadId( ) );
	  if ( pick_mate ) frag_ids.push_back( rl.PartnerReadId( ) );
	}
	if ( rl.Jump( ) ) {
	  jump_ids.push_back( rl.ReadId( ) );
	  if ( pick_mate ) jump_ids.push_back( rl.PartnerReadId( ) );
	}
	if ( rl.LongJump( ) ) {
	  long_jump_ids.push_back( rl.ReadId( ) );
	  if ( pick_mate ) long_jump_ids.push_back( rl.PartnerReadId( ) );
	}
      }
    }
  }
  
  cout << endl;

  // Sort and unique.
  cout << Date( ) << ": sorting" << endl;
  UniqueSort( frag_ids );
  UniqueSort( jump_ids );
  UniqueSort( long_jump_ids );

  // Overall stats.
  cout << "\nOVERALL STATS\n\n"
       << "      frag reads : " << frag_ids.isize( ) << "\n"
       << "      jump reads : " << jump_ids.isize( ) << "\n"
       << " long jump reads : " << long_jump_ids.isize( ) << "\n"
       << endl;
  
  // Sort and save.
  cout << Date( ) << ": saving" << endl;
  BinaryWriter::writeFile( out_frag_file.c_str(), frag_ids );
  BinaryWriter::writeFile( out_jump_file.c_str(), jump_ids );
  BinaryWriter::writeFile( out_long_jump_file.c_str(), long_jump_ids );

  // Done.
  cout << Date( ) << ": done" << endl;

}
