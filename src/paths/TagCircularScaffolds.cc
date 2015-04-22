///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"

#include "PairsManager.h"
#include "PrintAlignment.h"
#include "Superb.h"
#include "Vec.h"
#include "VecUtilities.h"
#include "efasta/EfastaTools.h"
#include "math/Functions.h"
#include "paths/AssemblyCleanupTools.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "paths/IdentifyCircularScaffolds.h"

/**
 * TagCircularScaffolds
 *
 * Run IdentifyCircularScaffolds to identify and tag circular
 * scaffolds. Optionally, clip one end of circular scaffolds, if their
 * ends overlap perfectly by an amount comparable with the overlap
 * estimated by the links. It saves a full allpaths assembly in
 * output.
 *
 * ASSEMBLY_IN: head of input assembly name
 * ASSEMBLY_OUT: head of output assembly name
 * MIN_LINKS: required to decide if scaffold is circular
 * *.READS: which reads to use (at least one must be given)
 * NUM_THREADS: if negative, use all available
 * CLIP: if true, clip one end of perfectly overlapping circular scaffolds
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String( SUBDIR );
  CommandArgument_String( ASSEMBLY_IN );
  CommandArgument_String( ASSEMBLY_OUT );
  CommandArgument_Int_OrDefault( MIN_LINKS, 14 );
  CommandArgument_String_OrDefault( FRAG_READS, "" );
  CommandArgument_String_OrDefault( JUMP_READS, "" );
  CommandArgument_String_OrDefault( LONG_JUMP_READS, "" );
  CommandArgument_UnsignedInt_OrDefault( NUM_THREADS, 0 );
  CommandArgument_Bool_OrDefault( CLIP, False );
  EndCommandArguments;

  // Constants.
  int MIN_TO_CLIP = 95;   // clip if perfect overlap is >= MIN_TO_CLIP

  // Validate args.
  ForceAssert( FRAG_READS != "" || JUMP_READS != "" || LONG_JUMP_READS != "" );

  // Thread control.
  NUM_THREADS = configNumThreads( NUM_THREADS );

  // Dir and file names.
  String data_dir = PRE + "/" + DATA;
  String run_dir = data_dir + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;

  String head = sub_dir + "/" + ASSEMBLY_IN;
  String contigs_file = head + ".contigs.fasta";
  String efasta_file = head + ".contigs.efasta";
  String supers_file = head + ".superb";
  String readlocs_file = head + ".readlocs";

  String head_out = sub_dir + "/" + ASSEMBLY_OUT;
  String out_rings_file = head_out + ".rings";
  String out_readlocs_file = head_out + ".readlocs";
  
  // Load.
  read_locs_on_disk locs_file( head, run_dir );

  cout << Date( ) << ": loading contigs fasta" << endl;
  vec<fastavector> contigs;
  LoadFromFastaFile( contigs_file, contigs );
  
  VecEFasta efastas;
  if ( IsRegularFile( efasta_file ) ) {
    cout << Date( ) << ": loading contigs efasta" << endl;
    LoadEfastaIntoStrings( efasta_file, efastas );
  }
  else {
    cout << Date( ) << ": generating contigs efasta" << endl;
    efastas.reserve( contigs.size( ) );
    for (size_t ii=0; ii<contigs.size( ); ii++)
      efastas.push_back( efasta( contigs[ii] ) );
  }

  cout << Date( ) << ": loading supers" << endl;
  vec<superb> supers;
  ReadSuperbs( supers_file, supers );

  // Import pairs
  vec<String> pairs_files( 3, "" );
  vec<PairsManager> pairs;
  if ( FRAG_READS != "" ) pairs_files[0] = FRAG_READS;
  if ( JUMP_READS != "" ) pairs_files[1] = JUMP_READS;
  if ( LONG_JUMP_READS != "" ) pairs_files[2] = LONG_JUMP_READS;
  for (int ii=0; ii<3; ii++) {
    if ( pairs_files[ii] == "" ) {
      pairs.push_back( PairsManager( ) );
      continue;
    }
    cout << Date( ) << ": loading " << pairs_files[ii] << " pairs file" << endl;
    String full_pairs_file = run_dir + "/" + pairs_files[ii] + ".pairs";
    if ( ! IsRegularFile( full_pairs_file ) )
      FatalErr( "Unable to find pairs file for: " + pairs_files[ii] );
    pairs.push_back( PairsManager( full_pairs_file ) );
  }
  
  cout << Date( ) << ": done loading\n" << endl;
  
  // Run code.
  vec<CRing> rings;
  IdentifyCircularScaffolds( pairs, contigs, supers, locs_file, rings,
			     &cout, MIN_LINKS );
  
  // Test for perfect overlaps.
  int nclip = 0;
  vec<int> clippers( contigs.size( ), 0 );
  if ( CLIP ) {
    cout << "\n" << Date( ) << ": testing for perfect overlap" << endl;
    for (size_t ring_id=0; ring_id<rings.size( ); ring_id++) {
      int super_id = rings[ring_id].SuperId( );
      int gap = rings[ring_id].GapSize( );
      int dev = rings[ring_id].GapDev( );
      int nlinks = rings[ring_id].NLinks( );

      int right_cid = supers[super_id].Tig( 0 );
      int left_cid = supers[super_id].Tig( supers[super_id].Ntigs( ) - 1 );
      int min_overlap = Max( gap - 3 * dev, MIN_TO_CLIP );
      int max_overlap = gap + 3 * dev;
      if ( (int)contigs[left_cid].size( ) < min_overlap ) continue;
      if ( (int)contigs[right_cid].size( ) < min_overlap ) continue;
      
      const bvec lb = contigs[left_cid].ToBasevector( );
      const bvec rb = contigs[right_cid].ToBasevector( );
      
      int offset = Min( (int)lb.size( ) + gap, (int)lb.size( ) );
      int band = 3 * dev;
      int errors = 0;
      align al;
      SmithWatBandedA2<unsigned int>( lb, rb, offset, band, al, errors );
      
      // Align is not perfect.
      if ( errors > 0 ) continue;
      
      // Align is too short.
      int al_len = al.Pos2( ) - al.pos2( );
      if ( al_len < MIN_TO_CLIP ) continue;
      
      // Align would clip away the whole contig (this should not happen).
      if ( al.pos1( ) < 1 ) continue;
      
      // Make sure there are no ambiguous bases in the clipped portion.
      bool amb_found = false;
      size_t clip_size = contigs[left_cid].size( ) - al.pos1( );
      const efasta &theE = efastas[left_cid];
      for (size_t ii=theE.size( ) - clip_size; ii<=theE.size( ); ii++) {
	if ( theE[ii] == '{' || theE[ii] == '}' ) {
	  amb_found = true;
	  break;
	}
      }
      if ( amb_found ) continue;

      // Ok: store amount to keep, and adjust super, and CRing (gap=0, dev=1).
      nclip++;
      clippers[left_cid] = al.pos1( );
      rings[ring_id] = CRing( super_id, 0, 1, nlinks );
      supers[super_id].SetLen( supers[super_id].Ntigs( )-1, al.pos1( ) );
      cout << " s" << super_id << ": perfect match of " << al_len << " bases\n";
    }
    
    // Report clipping events, and resize contigs.
    if ( nclip < 1 )
      cout << Date( ) << ": no valid perfect overlap found\n" << endl;
    else {
      cout << Date( ) << ": clipping " << nclip << " contig(s)\n" << endl;
      for (size_t cid=0; cid<clippers.size( ); cid++) {
	if ( clippers[cid] > 0 ) {
	  size_t clip_amt = contigs[cid].size( ) - clippers[cid];
	  efastas[cid].resize( efastas[cid].size( ) - clip_amt );
	  contigs[cid].resize( clippers[cid] );
	}
      }
    }
  }
  
  // Save rings.
  cout << Date( ) << ": saving rings file" << endl;
  ofstream out( out_rings_file.c_str( ) );
  if ( rings.size( ) < 1 )
    out << "No circular scaffolds found.\n" << endl;
  else {
    vec< vec<String> > table;
    vec<String> line;
    rings[0].GenerateLegend( line );
    table.push_back( line );
    for (int ii=0; ii<rings.isize( ); ii++) {
      int sid = rings[ii].SuperId( );
      int cid = supers[sid].Tig( supers[sid].Ntigs( ) - 1 );
      rings[ii].PrintableInfo( line );
      table.push_back( line );
    }
    PrintTabular( out, table, 3, "rrrr" );
    out.close( );
  }
  
  // Save readlocs.
  if ( CLIP && nclip > 0 ) {
    cout << Date( ) << ": saving readlocs... " << flush;
    vec<read_loc> locs;
    vec<read_loc> all_locs;
    int n_removed = 0;
    all_locs.reserve( locs_file.getLocsSize( ) );
    for (size_t cid=0; cid<contigs.size( ); cid++) {
      locs_file.LoadContig( cid, locs );
      if ( clippers[cid] < 1 )
	copy( locs.begin( ), locs.end( ), back_inserter( all_locs ) );
      else {
	for (int loc_id=0; loc_id<locs.isize( ); loc_id++) {
	  const read_loc &rloc = locs[loc_id];
	  if ( rloc.Stop( ) > clippers[cid] ) {
	    n_removed++;
	    continue;
	  }
	  all_locs.push_back( rloc );
	}
      }
    }
    double ratio = 10000. * SafeQuotient( n_removed, locs_file.getLocsSize( ) );
    cout << "(" << n_removed
	 << " readlocs removed, ie " << ToString( ratio, 1 )
	 << "%%)" << endl;
    WriteReadLocs( head_out, contigs.size( ), all_locs );
  }
  else {
    cout << Date( ) << ": copying readlocs" << endl;
    Cp( readlocs_file, out_readlocs_file );
  }

  // Save assembly (logging within).
  Assembly theAssembly( supers, efastas );
  theAssembly.check_integrity( );
  theAssembly.WriteAll( head_out );
  
  // Done.
  cout << Date( ) << ": done" << endl;
  
}

