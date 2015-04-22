///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Alignment.h"
#include "Basevector.h"
#include "ParseSet.h"
#include "VecUtilities.h"
#include "lookup/LookAlign.h"
#include "util/RunCommand.h"
// MakeDepend: dependency DisplayLocs
// MakeDepend: dependency EvalScaffolds
// MakeDepend: dependency MakeLookupTable

/**
 * class CClosure
 *
 * A closure is an interval on a POST_PATCH contig. Remark: begin_ is
 * always <= end_ (the member variable overlap_ is used to tag
 * closures coming from overlapping contigs). A closure can optionally
 * store the alignment of the interval onto a reference.
 */
class CClosure {

public: 
  
  CClosure( ) : 
    contig_id_ ( -1 ),
    begin_ ( 0 ),
    end_ ( 0 ),
    overlap_ ( false ),
    al_ ( ),
    ref_id_ ( 0 ) { }
  
  CClosure( int contig_id, int begin, int end, bool overlap,
	    align al, int ref_id ) :
    contig_id_ ( contig_id ),
    begin_ ( begin ),
    end_ ( end ),
    overlap_ ( overlap ),
    al_ ( al ),
    ref_id_ ( ref_id ) { ForceAssertLe( begin_, end_); }

  int ContigId( ) const { return contig_id_; }
  int Begin( ) const { return begin_; }
  int End( ) const { return end_; }
  bool Overlap( ) const { return overlap_; }
  align Align( ) const { return al_; }
  int RefId( ) const { return ( ref_id_ < 0 ? - ref_id_ - 1 : ref_id_ ); }
  int Gap( ) const { return ( overlap_ ? -1 : +1 ) * ( end_ - begin_ ); }
  bool RcRef( ) const { return ref_id_ < 0; }
  bool IsAligned( ) const { return al_.Nblocks( ) > 0; }
  
  void PrintInfo( ostream &out,
		  const vec<String> &pile_ups,
		  const vecbvec *contigs = 0,
		  const vecbvec *genome = 0 ) { 
    
    // Bases of contig and (possibly rc-ed) reference, or null.
    const bvec *contig = contigs ? & ( (*contigs)[contig_id_] ) : 0;
    const bvec *ref = 0;
    bvec rc_ref;
    if ( genome ) {
      if ( this->RcRef( ) ) {
	rc_ref = (*genome)[this->RefId( )];
	rc_ref.ReverseComplement( );
      }
      ref = & ( this->RcRef( ) ? rc_ref : (*genome)[this->RefId( )] );
    }
    
    // Summary line.
    out << "tig" << contig_id_
	<< " [" << begin_ << ", " << end_ << ")"
	<< "  gap = " << this->Gap( );
    if ( this->IsAligned( ) ) {
      out << "  aligns " << ( this->RcRef( ) ? "rc of " : "" )
	  << "ref" << this->RefId( )
	  << " [" << al_.pos2( ) << ", " << al_.Pos2( ) << ")";
      if ( contig && ref ) {
	vector<int> errors = al_.MutationsGap1Gap2( *contig, *ref );
	out << "  mut = " << errors[0]
	    << "  ins = " << errors[2]
	    << "  del = " << errors[1];
      }
    }
    out << "\n\n";

    // Print pile ups (without reference).
    if ( ! ( this->IsAligned( ) && contigs && ref ) ) {
      vec< vec<String> > table;
      for (int ii=begin_; ii<end_; ii++) {
	String sTig = String( as_base( (*contig)[ii] ) );	
	table.push_back( MkVec( ToString( ii ), sTig, pile_ups[ii] ) );
      }
      PrintTabular( out, table, 2, "rl" );
      return;
    }

    // Print pile ups (with reference).
    vec< vec<String> > table;
    String sIndel = "-";
    String sMismatch = "*";
    String s0 = "";

    int pos1 = al_.pos1( );
    int pos2 = al_.pos2( );
    for (int block_id=0; block_id<al_.Nblocks( ); block_id++) {

      // Move on gap.
      int gap = al_.Gaps( block_id );
      if ( gap < 0 ) {
	for (int ii=0; ii<-gap; ii++) {
	  String sPos = ToString( pos1 );
	  String sTig = String( as_base( (*contig)[pos1] ) );
	  table.push_back( MkVec( sPos, sIndel, s0, sTig, pile_ups[pos1] ) );
	  pos1++;
	}
      }
      if ( gap > 0 ) {
	for (int ii=0; ii<gap; ii++) {
	  String sRef = String( as_base( (*ref)[pos2] ) );
	  table.push_back( MkVec( s0, sIndel, sRef, s0, s0 ) );
	  pos2++;
	}
      }
      
      // Move on length.
      int length = al_.Lengths( block_id );      
      for (int ii=0; ii<length;ii++) {
	String sPos = ToString( pos1 );
	String sTig = String( as_base( (*contig)[pos1] ) );
	String sRef = String( as_base( (*ref)[pos2] ) );
	String sDecor = ( sTig == sRef ? s0 : sMismatch );
	table.push_back( MkVec( sPos, sDecor, sRef, sTig, pile_ups[pos1] ) );
	pos1++;
	pos2++;
      }
    }

    PrintTabular( out, table, 2, "rllll" );
    
  }
    
    
private:
  
  int contig_id_; // id of POST_PATCH contig
  int begin_;     // begin of closure interval
  int end_;       // end of closure interval
  bool overlap_;  // if gap between PRE_PATCH contigs is <=0 (contigs overlap)

  align al_;      // align of region of POST_PATCH contig on reference
  int ref_id_;    // reference id (<0 if reference is rc)
  
};

/**
 * EvalPatching
 *
 * Evaluation code for patched intervals.
 *
 * PRE_PATCH: assembly head before patching
 * POST_PATCH: assembly head after patching
 * USE_REF: align POST_PATCH contigs against reference
 * FORCE: do not use cached data
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_String_OrDefault( PRE_PATCH, "linear_scaffolds0.clean.patched" );
  CommandArgument_String_OrDefault( POST_PATCH, "linear_scaffolds0.clean.patched.pacbio" );
  CommandArgument_Bool_OrDefault( USE_REF, False );
  CommandArgument_Bool_OrDefault( FORCE, False );
  EndCommandArguments;
  
  // Dir and file names.
  String pdr = "PRE=" + PRE + " DATA=" + DATA + " RUN=" + RUN;
  String pdrs = pdr + " SUBDIR=" + SUBDIR;

  String data_dir = PRE + "/" + DATA;
  String run_dir = data_dir + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
  String pre_patch_head = sub_dir + "/" + PRE_PATCH;
  String post_patch_head = sub_dir + "/" + POST_PATCH;
  String pre_patch_eval_dir = pre_patch_head + ".EvalOnPostPatched";
  String post_patch_eval_dir = post_patch_head + ".Eval";
  String post_patch_display_dir = post_patch_head + ".DisplayLocs";

  String genome_fastb = data_dir + "/genome.fastb";
  String genome_lookup = data_dir + "/genome.lookup";
  String post_patch_contigs_lookup = post_patch_head + ".contigs.lookup";
  String post_patch_contigs_fastb = post_patch_head + ".contigs.fastb";
  String post_patch_aligns_file = post_patch_eval_dir + "/aligns.qlt";
  String pre_patch_aligns_file = pre_patch_eval_dir + "/aligns.qlt";
  
  // No genome.lookup found: fatal error.
  if ( USE_REF && ! IsRegularFile( genome_lookup ) ) {
    cout << "genome.lookup file missing. Leaving now.\n" << endl;
    return 1;
  }

  // Eval POST_PATCH assembly.
  if ( USE_REF && ( FORCE || ! IsRegularFile( post_patch_aligns_file ) ) ) {
    String theCommand
      = "EvalScaffolds LOOKUP=" + genome_lookup
      + " SCAFFOLDS=" + post_patch_head
      + " FORCE=" + ( FORCE ? "True" : "False" );
    if ( ! RunCommandBool( theCommand ) ) {
      cout << "EvalScaffolds failed. Leaving now.\n" << endl;
      return 1;
    }
  }

  // Generate lookup table of POST_PATCH contigs.
  if ( FORCE || ! IsRegularFile( post_patch_contigs_lookup ) ) {
    String theCommand
      = "MakeLookupTable LOOKUP_ONLY=True SOURCE=" + post_patch_head
      + ".contigs.fastb OUT_HEAD=" + post_patch_head
      + ".contigs";
    if ( ! RunCommandBool( theCommand ) ) {
      cout << "MakeLookupTable failed. Leaving now.\n" << endl;
      return 1;
    }
  }

  // Eval PRE_PATCH against POST_PATCH.
  if ( FORCE || ! IsRegularFile( pre_patch_aligns_file ) ) {
    String theCommand
      = "EvalScaffolds LOOKUP=" + post_patch_contigs_lookup
      + " SCAFFOLDS=" + pre_patch_head
      + " OUT_DIR=" + pre_patch_eval_dir
      + " FORCE=" + ( FORCE ? "True" : "False" );
    if ( ! RunCommandBool( theCommand ) ) {
      cout << "EvalScaffolds failed. Leaving now.\n" << endl;
      return 1;
    }
  }
  
  // Load aligns of PRE_PATCH contigs against POST_PATCH contigs.
  cout << Date( ) << ": loading aligns of pre patch contigs" << flush;
  vec<look_align_plus> pre_hits;
  LoadLookAlignPlus( pre_patch_aligns_file, pre_hits );
  
  // PRE_PATCH contigs should be subsumed in POST_PATCH contigs.
  int n_removed_pre = 0;
  for (int ii=0; ii<pre_hits.isize( ); ii++) {
    const look_align_plus &lap = pre_hits[ii];
    if ( lap.a.pos1( ) != 0 || lap.a.Pos1( ) != (int)lap.query_length ) {
      pre_hits.erase( pre_hits.begin( ) + ii );
      n_removed_pre++;
      ii--;
    }
  }
  if ( n_removed_pre > 0 )
    cout << " (" << n_removed_pre << " aligns discarded)";
  cout << endl;
  
  // Sort pre_hits by start on POST_PATCH contigs.
  order_lookalign_TargetBegin sorter;
  sort( pre_hits.begin( ), pre_hits.end( ) );

  // Load aligns of POST_PATCH contigs against reference.
  vec<look_align_plus> post_hits;
  if ( USE_REF ) {
    cout << Date( ) << ": loading aligns of post patch contigs" << endl;
    LoadLookAlignPlus( post_patch_aligns_file, post_hits );
  }  
  
  // Load bases.
  cout << Date( ) << ": loading bases of contigs" << endl;
  vecbvec contigs( post_patch_contigs_fastb );

  vecbvec full_genome;
  if ( USE_REF ) {
    cout << Date( ) << ": loading bases of contigs" << endl;
    full_genome.ReadAll( genome_fastb );
  }
  const vecbvec *genome = USE_REF ? &full_genome : 0;

  cout << Date( ) << ": done loading\n" << endl;
  
  // Map POST_PATCH contigs to post_hits.
  vec< vec<size_t> > post_contig2hit;
  if ( USE_REF ) {
    int ntigs = MastervecFileObjectCount( post_patch_contigs_fastb );
    post_contig2hit.resize( ntigs );
    for (size_t ii=0; ii<post_hits.size( ); ii++) {
      const look_align_plus &lap = post_hits[ii];
      post_contig2hit[lap.query_id].push_back( ii );
    }
  }
  
  // Identify all closures, as intervals on POST_PATCH contigs.
  vec<CClosure> closures;
  closures.reserve( pre_hits.size( ) );
  for (int ii=1; ii<pre_hits.isize( ); ii++) {
    const look_align_plus &left = pre_hits[ii-1];
    const look_align_plus &right = pre_hits[ii];

    // Some filtering.
    if ( left.target_id != right.target_id ) continue;
    if ( left.rc1 || right.rc1 ) continue;

    // The window of closure on POST_PATCH.
    int tid = left.target_id;
    int begin = left.a.Pos2( );
    int end = right.a.pos2( );
    bool overlap = ( begin >= end );
    if ( overlap ) swap( begin, end );

    // Loop over all aligns for this POST_PATCH contig.
    int local_count = 0;
    if ( USE_REF ) {
      for (int jj=0; jj<post_contig2hit[tid].isize( ); jj++) {
	const look_align_plus &post_hit = post_hits[ post_contig2hit[tid][jj] ];
	align al = post_hit.a;
	
	// Reverse al (if needed), so contig is always fw.
	if ( post_hit.rc1 ) {
	  int contig_len = post_hit.query_length;
	  int target_len = post_hit.target_length;
	  al.ReverseThis( contig_len, target_len );
	}
	
	// Portion containing window of closure is not aligned.
	if ( al.pos1( ) > begin ) continue;
	if ( end > al.Pos1( ) ) continue;
	
	// Trim align and instantiate closure.
	TrimAlignmentFront( al, begin - al.pos1( ) );
	TrimAlignmentBack( al, al.Pos1( ) - end );
	int ref_id = post_hit.target_id;
	if ( post_hit.rc1 ) ref_id = - ref_id - 1;
	
	closures.push_back( CClosure( tid, begin, end, overlap, al, ref_id ) );
	local_count++;
      }
    }
    
    // Default closure.
    if ( local_count < 1 )
      closures.push_back( CClosure( tid, begin, end, overlap, align( ), 0 ) );
  }
  
  // Run DisplayLocs.
  Mkpath( post_patch_display_dir );
  String display_head = post_patch_display_dir + "/tig";
  for (int ii=0; ii<closures.isize( ); ii++) {
    int contig_id = closures[ii].ContigId( );
    String display_file = display_head + ToString( contig_id ) + ".display";
    if ( FORCE || ! IsRegularFile( display_file ) ) {
      cout << Date( ) << ": displaying locs of tig" << contig_id << endl;
      String theCommand
	= "DisplayLocs " + pdrs
	+ " ASSEMBLY=" + POST_PATCH
	+ " TIGS=" + ToString( contig_id )
	+ " SHOW_LOCS=False SHOW_PILEUPS=True >&" + display_file;
      if ( ! RunCommandBool( theCommand ) ) {
	cout << "DisplayLocs failed. Leaving now.\n" << endl;
	return 1;
      }
    }
  }
  
  // Print closures.
  String separator;
  for (int ii=0; ii<80; ii++)
    separator += "=";

  cout << separator << "\n";
  int cached_id = -1;
  vec<String> pile_ups;
  for (int ii=0; ii<closures.isize( ); ii++) {
    int contig_id = closures[ii].ContigId( );

    // Cache pile_ups.
    if ( cached_id != contig_id ) {
      cout << Date( ) << ": caching pile-ups of tig" << contig_id << endl;
      cached_id = contig_id;
      pile_ups.clear( );
      pile_ups.resize( contigs[contig_id].size( ) );

      String display_file = display_head + ToString( contig_id ) + ".display";
      ifstream in( display_file.c_str( ) );
      String line;
      while ( in ) {
	getline( in, line );
	if ( ! in ) break;
	if ( ! line.Contains( " " ) ) continue;
	if ( line.Before( " " ) != "0" ) continue;
	pile_ups[0] = line.After( " " );
	for (int ii=1; ii<contigs[contig_id].isize( ); ii++) {
	  getline( in, line );
	  ForceAssert( in );
	  ForceAssertEq( line.Before( " " ), ToString( ii ) );
	  pile_ups[ii] = line.After( " " );
	}
      }
      in.close( );
    }
    
    cout << "\n";
    closures[ii].PrintInfo( cout, pile_ups, &contigs, genome );
    cout << "\n" << separator << "\n";
  }

}

