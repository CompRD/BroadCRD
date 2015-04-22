///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>
#include "CoreTools.h"
#include "MainTools.h"
#include "Superb.h"
#include "math/HoInterval.h"
#include "paths/ReadLoc.h"

/**
 */


/**
 * ProbeInversions
 *
 
*/


int main( int argc, char *argv[] ){
  
  RunTime( );
  
  BeginCommandArguments;

  CommandArgument_String_OrDefault( PRE, "" );
  CommandArgument_String( SUPERB );
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");

  CommandArgument_UnsignedInt_OrDefault_Doc(COVERAGE_THRESH, 10, 
    "Threshold for base coverage in the suspect region to report");
  CommandArgument_Int_OrDefault_Doc(LENGTH_THRESH, 100, 
    "Length threshold for reporting suspect region");
  CommandArgument_Bool_OrDefault( VERBOSE, False );
  EndCommandArguments;

  // Set number of threads.
  NUM_THREADS = configNumThreads( NUM_THREADS 	);
  omp_set_num_threads( NUM_THREADS );


  if ( PRE != "" ) PRE += "/";

  String superb_head = SUPERB.SafeAfterLast("/").SafeBefore(".super");
  PRINT(superb_head);
  String superb_dir = SUPERB.SafeBeforeLast("/");
  if ( superb_dir == "" ) superb_dir = "./";
  PRINT( superb_dir);

  String superb_file = PRE + SUPERB;
  String contigs_file = PRE + superb_head + ".contigs.fastb";
  PRINT( contigs_file );
  ForceAssert( IsRegularFile( contigs_file ) );
  ForceAssert( IsRegularFile( superb_file ) );

  cout << Date() << ": preparing readloc reading" << endl;
  read_locs_on_disk locs_file( superb_head, superb_dir );

  cout << Date() << ": reading contigs file" << endl;
  vecbasevector contigs( contigs_file );

  cout << Date() << ": reading superbs" << endl;
  vec<superb> supers;
  ReadSuperbs( superb_file, supers);

  const int n_contigs = contigs.size( );
  
  // Maps.
  cout << Date() << ": computing maps" << endl;
  vec<int> super_len( supers.size( ), 0 );
  vec<int> start_on_super( n_contigs );
  vec<int> to_super( n_contigs );
  vec<int> tig_len( n_contigs );
  for (int ii=0; ii<supers.isize( ); ii++) {
    int pos = 0;
    super_len[ii] = supers[ii].TrueLength( );
    for (int jj=0; jj<supers[ii].Ntigs( ); jj++) {
      start_on_super[ supers[ii].Tig( jj ) ] = pos;
      to_super[ supers[ii].Tig( jj ) ] = ii;
      tig_len[ supers[ii].Tig( jj ) ] = supers[ii].Len( jj );
      pos += supers[ii].Len( jj );
      if ( jj < supers[ii].Ntigs( ) - 1 ) pos += supers[ii].Gap( jj );
    }
  }

  vec< vec<ho_interval> > tig_alregs( n_contigs );

  Bool expectEarlyFw = True;
  size_t n_earlyFw = 0, n_earlyRc = 0;
  cout << Date() << ": calculating raw coverages" << endl;
  #pragma omp parallel for
  for (int super_id=0; super_id<supers.isize( ); super_id++) {
    const superb &sup = supers[super_id];
    const int slen = super_len[super_id];
    vec<NormalDistribution> nds;

    // Loop over contigs in super.
    for (int cgpos=0; cgpos<sup.Ntigs( ); cgpos++) {
      int tig1 = sup.Tig( cgpos );
      
      // Loop over all locs in contig.
      vec<read_loc> locs;
      #pragma omp critical
      locs_file.LoadContig( tig1, locs );
      for (int loc_id=0; loc_id<locs.isize( ); loc_id++) {
	const read_loc &rloc = locs[loc_id];
	int tig2 = rloc.PartnerContigId( );

	if ( super_id != to_super.at(tig2) ) continue;	
	if ( rloc.ReadClass() != 1 ) continue;
	if ( ! rloc.PartnerPlaced() )  continue;
	
	int start1 = rloc.Start() < 0 ? 0 : rloc.Start();
	int end1 = rloc.Start() + rloc.ReadLength() -1 < tig_len[tig1] ? 
	  rloc.Start() + rloc.ReadLength() -1 : tig_len[tig1]; 

	start1 += start_on_super[tig1];
	end1 += start_on_super[tig1];
	
	int start2 = rloc.PartnerStart() < 0 ? 0 : rloc.PartnerStart();
	int end2 = rloc.PartnerStart() + rloc.PartnerReadLength() -1 < tig_len[tig2] ? 
	  rloc.PartnerStart() + rloc.PartnerReadLength() -1 : tig_len[tig2]; 
	
	start2 += start_on_super[tig2];
	end2 += start_on_super[tig2];
	

	Bool earlyFw = False, earlyRc = False;
	if ( start1 < start2 && rloc.Fw() || start2 < start1 && rloc.PartnerFw() )
	  earlyFw = True;
	if ( start1 < start2 && rloc.Rc() || start2 < start1 && rloc.PartnerRc() )
	  earlyRc = True;

	// Only keep fw-fw and rc-rc pairs in the same super.
	if ( earlyFw ) n_earlyFw++;
	if ( earlyRc ) n_earlyRc++;

	if ( rloc.Fw() && rloc.PartnerRc() ) continue;
	if ( rloc.Rc() && rloc.PartnerFw() ) continue;


	if ( expectEarlyFw ){
	  if ( earlyFw ){
	    if ( start1 < start2 ) tig_alregs[tig2].push_back( ho_interval( start2 - start_on_super[tig2], end2 - start_on_super[tig2] ) );
	    else if ( start2 < start1 ) tig_alregs[tig1].push_back( ho_interval( start1 - start_on_super[tig1], end1 - start_on_super[tig1] ) );
	  }
	  if ( earlyRc ){
	    if ( start1 < start2 ) tig_alregs[tig1].push_back( ho_interval( start1 - start_on_super[tig1], end1 - start_on_super[tig1] ) );
	    else if ( start2 < start1 ) tig_alregs[tig2].push_back( ho_interval( start2 - start_on_super[tig2], end2 - start_on_super[tig2] ) );
	  }
	}
	else{
	  if ( earlyRc ){
	    if ( start1 < start2 ) tig_alregs[tig2].push_back( ho_interval( start2 - start_on_super[tig2], end2 - start_on_super[tig2] ) );
	    else if ( start2 < start1 ) tig_alregs[tig1].push_back( ho_interval( start1 - start_on_super[tig1], end1 - start_on_super[tig1] ) );
	  }
	  if ( earlyFw ){
	    if ( start1 < start2 ) tig_alregs[tig1].push_back( ho_interval( start1 - start_on_super[tig1], end1 - start_on_super[tig1] ) );
	    else if ( start2 < start1 ) tig_alregs[tig2].push_back( ho_interval( start2 - start_on_super[tig2], end2 - start_on_super[tig2] ) );
	  }
	}
	  
      }
    }
  }

  cout << Date() << ": computing suspect segments" << endl;
  vec< vec<ho_interval> > tig_covs( n_contigs );
  for ( int tid = 0; tid < n_contigs; tid++ ){
    vec<ho_interval> cov;
    ExtractGivenCoverage( tig_len[tid], COVERAGE_THRESH, tig_alregs[tid], cov );
    for ( size_t ic = 0; ic < cov.size(); ic++ ){
      if ( cov[ic].Stop() - cov[ic].Start() > LENGTH_THRESH )
 	tig_covs[tid].push_back( cov[ic] );
    }
  }
  
  cout << Date() << ": printing suspect segments" << endl;
  for ( size_t sid = 0; sid < supers.size(); sid++){
    for ( int i = 0; i < supers[sid].Ntigs(); i++ ){
      int tid = supers[sid].Tig( i );
      if ( tig_covs[tid].size() == 0 ) continue;
      cout << endl; PRINT3( sid, tid, tig_len[tid] );
      for ( size_t ic = 0; ic < tig_covs[tid].size(); ic++ ){
	int start = tig_covs[tid][ic].Start(); int end = tig_covs[tid][ic].Stop() -1;
	int len = end - start + 1;
	cout << start << "-" << end << "(l=" << len << ")\t";
      }
      cout << "\n";
    }
  }
  
  PRINT2( n_earlyFw, n_earlyRc );
  cout << Date( ) << ": IdentifyInversions done" << endl;
  
}
