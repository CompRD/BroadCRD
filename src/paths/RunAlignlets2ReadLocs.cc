/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "PairsManager.h"
#include "VecUtilities.h"
#include "lookup/LookAlign.h"
#include "paths/Alignlets2ReadLocs.h"
#include "util/RunCommand.h"

/**
 * RunAlignlets2ReadLocs
 *
 * Generate read_locs from given alignlets and indexes.
 *
 * INPUT
 *   <CONTIGS_HEAD>.fastb
 *   <RUN_DIR>/<FRAGS>.pairs            (if FRAGS is given)
 *   <RUN_DIR>/<JUMPS>.pairs            (if JUMPS is given)
 *   <RUN_DIR>/<LONG_JUMPS>.pairs       (if LONG_JUMPS is given)
 *   <CONTIGS_HEAD>.frag.{aligns,idx}   (if FRAGS is given)
 *   <CONTIGS_HEAD>.jumps.{aligns,idx}  (if JUMPS is given)
 *   <CONTIGS_HEAD>.Jumps.{aligns,idx}  (if LONG_JUMPS is given)
 *
 * OUTPUT:
 *   <CONTIGS_HEAD>.readlocs
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( CONTIGS_HEAD );
  CommandArgument_String( RUN_DIR );
  CommandArgument_String_OrDefault( FRAGS, "" );
  CommandArgument_String_OrDefault( JUMPS, "" );
  CommandArgument_String_OrDefault( LONG_JUMPS, "" );
  EndCommandArguments;
  
  // Dir and file names.
  String contigs_F = CONTIGS_HEAD + ".fastb";
  
  String frag_pairs_F = RUN_DIR + "/" + FRAGS + ".pairs";
  String frag_aligns_F = CONTIGS_HEAD + ".frag.aligns";
  String frag_idx_F = CONTIGS_HEAD + ".frag.idx";

  String jump_pairs_F = RUN_DIR + "/" + JUMPS + ".pairs";
  String jump_aligns_F = CONTIGS_HEAD + ".jump.aligns";
  String jump_idx_F = CONTIGS_HEAD + ".jump.idx";

  String Jump_pairs_F = RUN_DIR + "/" + LONG_JUMPS + ".pairs";
  String Jump_aligns_F = CONTIGS_HEAD + ".Jump.aligns";
  String Jump_idx_F = CONTIGS_HEAD + ".Jump.idx";
  
  bool use_frags = ( FRAGS != "" );
  bool use_jumps = ( JUMPS != "" );
  bool use_Jumps = ( LONG_JUMPS != "" );
  if ( ! ( use_frags || use_jumps || use_Jumps ) ) {
    cout << "Fatal error. No aligns found!\n" << endl;
    return 1;
  }
  
  vec<String> needed;
  needed.push_back( contigs_F );
  if ( use_frags ) needed.push_back( frag_pairs_F );
  if ( use_frags ) needed.push_back( frag_aligns_F );
  if ( use_frags ) needed.push_back( frag_idx_F );
  if ( use_jumps ) needed.push_back( jump_pairs_F );
  if ( use_jumps ) needed.push_back( jump_aligns_F );
  if ( use_jumps ) needed.push_back( jump_idx_F );
  if ( use_Jumps ) needed.push_back( Jump_pairs_F );
  if ( use_Jumps ) needed.push_back( Jump_aligns_F );
  if ( use_Jumps ) needed.push_back( Jump_idx_F );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;
  
  // Create empty alignlets datasets.
  int n_tigs = MastervecFileObjectCount( contigs_F );
  SAlignletsData frags;
  SAlignletsData jumps;
  SAlignletsData Jumps;

  PairsManager fpairs;
  vec<alignlet> faligns;
  vec<int> fidx;

  PairsManager jpairs;
  vec<alignlet> jaligns;
  vec<int> jidx;

  PairsManager Jpairs;
  vec<alignlet> Jaligns;
  vec<int> Jidx;

  // Fill alignlets datasets.
  if ( use_frags ) {
    cout << Date( ) << ": loading frag data" << endl;
    fpairs.Read( frag_pairs_F );
    BinaryReader::readFile( frag_aligns_F, &faligns );
    BinaryReader::readFile( frag_idx_F, &fidx );
    frags.Set( &fpairs, &faligns, &fidx );
  }
  
  if ( use_jumps ) {
    cout << Date( ) << ": loading jump data" << endl;
    jpairs.Read( jump_pairs_F );
    BinaryReader::readFile( jump_aligns_F, &jaligns );
    BinaryReader::readFile( jump_idx_F, &jidx );
    jumps.Set( &jpairs, &jaligns, &jidx );
  }

  if ( use_Jumps ) {
    cout << Date( ) << ": loading long jump data" << endl;
    Jpairs.Read( Jump_pairs_F );
    BinaryReader::readFile( Jump_aligns_F, &Jaligns );
    BinaryReader::readFile( Jump_idx_F, &Jidx );
    Jumps.Set( &Jpairs, &Jaligns, &Jidx );
  }
  
  // Generate readlocs (logging within).
  vec<read_loc> locs;
  Alignlets2ReadLocs( frags, jumps, Jumps, locs, &cout );

  // Save.
  cout << Date( ) << ": saving readlocs" << endl;
  WriteReadLocs( CONTIGS_HEAD, n_tigs, locs );
  
  // Done.
  cout << Date( ) << ": RunAlignlets2ReadLocs done" << endl;
  return 0;

}

