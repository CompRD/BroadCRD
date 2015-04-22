///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/* ReadBadCoverage.  Find positions in the reference that are not covered by any one read extending +/- RADIUS to left and right */

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS


#include "MainTools.h"
#include "Superb.h"
#include "paths/ReadLoc.h"
#include "system/LockedData.h"
#include "system/ParsedArgs.h"
#include "system/SpinLockedData.h"
#include "VecUtilities.h"
#include "math/HoInterval.h"
#include "lookup/LookAlign.h"

#include <omp.h>


void align2errorless( const basevector& b, const basevector& ref,  const align& a, vec<ho_interval>& v_hos, const int type = 0 ){
  // type: 0 - errorless, 1 - indelless, 2 - all
  v_hos.clear();
  if ( type == 2 ){
    v_hos.push_back( ho_interval( a.pos2(), a.Pos2() ) );
    //ForceAssertEq( b.isize(), a.Pos2() - a.pos2() );
    return;
  }
  int p1 = a.pos1( ), p2 = a.pos2( );
  ForceAssertGe( p1, 0 );
  ForceAssertLt( p2, ref.isize() );
  Bool is_open = False;
  int start = -1, stop = -1;
  for ( int j = 0; j < a.Nblocks( ); j++ ){    
    if ( a.Gaps(j) > 0 )  p2 += a.Gaps(j);    
    if ( a.Gaps(j) < 0 )  p1 -= a.Gaps(j);
    if ( type == 0 ){
      is_open = False;
      for ( int x = 0; x < a.Lengths(j); x++ ){    
	if ( ref[p2] == b[p1] ){
	  if ( ! is_open ) start = p2, is_open = True;
	}
	else{
	  if ( is_open ){
	    stop = p2;
	    is_open = False;
	    v_hos.push_back( ho_interval( start, stop ) );
	  }
	}
	++p1; ++p2;    
      }
      if ( is_open ){
	stop = p2; 
	is_open = False;
	v_hos.push_back( ho_interval( start, stop ) ); 
      }
    }else{ // indelless
      v_hos.push_back( ho_interval( p2, p2 + a.Lengths(j) ) );
      p1 += a.Lengths(j);
      p2 += a.Lengths(j); 
    }
  }   
  ForceAssertLe( p1, b.isize() );
  ForceAssertLe( p2, ref.isize() );
  ForceAssertEq( p1, a.Pos1() );
  ForceAssertEq( p2, a.Pos2() );
  return;
}


int main(int argc, char *argv[]){
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String(READ_LOCS);
  CommandArgument_String(REF_FILE);
  CommandArgument_String_OrDefault_Doc(LIB_TYPES, "frag", 
	 "can be {frag,jump,long_jump} or a subset");
  CommandArgument_Int_OrDefault_Doc(COV_TYPE, "0", 
	 "can be 0 - errorless, 1 - gapless, 2 - any intervals in read alignments");
  CommandArgument_Int_OrDefault_Doc(EDGE_SKIP, "300", 
	 "skip edges of reference contigs");

  CommandArgument_String_OrDefault(FRAG_HEAD, "frag_reads_filt_cpd" );
  CommandArgument_String_OrDefault(JUMP_HEAD, "jump_reads_filt_cpd" );
  CommandArgument_String_OrDefault(LONG_JUMP_HEAD, "long_jump_reads" );

  CommandArgument_Int_OrDefault_Doc(RADIUS, 10,
	 "radius around reference position to compute coverage");
  CommandArgument_Int_OrDefault_Doc(COV_THRESH, -1,
	 "segment coverage threshold");
  CommandArgument_Double_OrDefault_Doc(COV_THRESH_FRAC, -0.1,
	 "segment coverage threshold as fraction of average coverage");
  CommandArgument_String_OrDefault_Doc(MASKING_ALIGNS_FILE, "",
	 "if specified load lookaligns that define masked regions "
		     "(does computations outside of those regions)");
  CommandArgument_Bool_OrDefault_Doc(MASK_COMPLEMENT, False,
	 "if set True then lookaligns in MASKING_ALIGNS_FILE define complement set masking intervals");
  
  CommandArgument_Int_OrDefault_Doc( NUM_THREADS, 16,
	 "number of threads to be used" );
  EndCommandArguments;
  
  NUM_THREADS = configNumThreads(NUM_THREADS);
  omp_set_num_threads(NUM_THREADS);

  
  if ( COV_THRESH < 0 && COV_THRESH_FRAC < 0 ){
    FatalErr("SPECIFY COV_THRESH or COV_THESH_FRAC");
  }

  // find out which classes of reads to use and read them in


 // read in reads

  vecbasevector frag_reads;
  vecbasevector jump_reads;
  vecbasevector long_jump_reads;

  String reads_dir = READ_LOCS.SafeBeforeLast("/");
  DPRINT( reads_dir );
  
  // read reference
  
  ForceAssert( IsRegularFile(REF_FILE) );
  vecbasevector reftigs( REF_FILE ); 
  DPRINT( reftigs.size() );
  
  vec< vec<ho_interval> > tig2read_hos( reftigs.size() );

  int64_t sumAliLen = 0, sumReadLen = 0;

  if ( READ_LOCS.Contains(".readloc") ){
    vec<Bool> lib_types( 3, False );
    vec<String> lib_types_s;
    ParseStringSet( LIB_TYPES, lib_types_s );
    for ( int i = 0; i < lib_types_s.isize( ); i++ ){    
      if ( lib_types_s[i] == "frag" ){
	lib_types[0] = True;
	frag_reads.ReadAll( reads_dir + "/" + FRAG_HEAD + ".fastb" );
      }
      else if ( lib_types_s[i] == "jump" ){
	lib_types[1] = True;
	jump_reads.ReadAll( reads_dir + "/" + JUMP_HEAD + ".fastb" );
      }
      else if ( lib_types_s[i] == "long_jump" ){
	lib_types[2] = True;
	long_jump_reads.ReadAll( reads_dir + "/" + LONG_JUMP_HEAD + ".fastb" );
      }
      else{    
	cout << "Illegal LIB_TYPES." << endl;
	return 1;    
      }    
    }

  
    ForceAssert( IsRegularFile( READ_LOCS ) );
    String read_locs_dir = READ_LOCS.SafeBeforeLast("/");
    String read_locs_head = (READ_LOCS.SafeAfterLast("/")).SafeBefore(".");
    DPRINT2( read_locs_dir, read_locs_head );
    read_locs_on_disk locs_file( READ_LOCS.SafeBeforeLast("."), read_locs_dir );
    
    cout << Date() << ": finding read alignment intervals that don't have errors" << endl;
    cout << Date() << ": going through reference contigs" << endl;
    int ntigs_done = 0;
    #pragma omp parallel for
    for ( size_t rtig = 0; rtig < reftigs.size( ); rtig++ ){
      
      #pragma omp critical
      {
	ntigs_done++;
	Dot(cout, 100.0 * (double)ntigs_done/(double)reftigs.size() );	
      }
      
      vec<read_loc> locs;
      #pragma omp critical
      { locs_file.LoadContig( rtig, locs ); }
      
      if ( locs.size() < 1 ) continue;
      
      vec<ho_interval> v_hos;
      for ( size_t rli = 0; rli < locs.size(); rli++ ){
	read_loc rl = locs[rli];
	if ( rl.Rc() ) continue;
	uint64_t rid = rl.ReadId();
	basevector b;
	if ( rl.ReadClass( ) == 0 ) b = frag_reads[rid];
	if ( rl.ReadClass( ) == 1 ) b = jump_reads[rid];
	if ( rl.ReadClass( ) == 2 ) b = long_jump_reads[rid];
	align a;
	rl.GetAlign( a, b, reftigs[rtig] );
	// look_align la;
	//       la.ResetFromAlign( a, b, reftigs[rtig] );
	//       ForceAssertEq( b.isize(), la.Pos2() - la.pos2() );
	vec<ho_interval> v_hos_loc;
	align2errorless( b, reftigs[rtig], a, v_hos_loc, COV_TYPE );
	v_hos.append( v_hos_loc ); 
	sumAliLen += a.Pos1() - a.pos1();
	sumReadLen += b.isize();
      } 
      #pragma omp critical
      { tig2read_hos[rtig] = v_hos; }
      // mark reference position not well covered by reads
    }
  }else{
    String reads_file = READ_LOCS.SafeBeforeLast(".") + ".fastb";
    PRINT( reads_file );
    ForceAssert( IsRegularFile( reads_file ) );
    vecbasevector reads( reads_file );
    cout << Date() << ": loading read alignments from " << READ_LOCS << endl;
    vec<look_align> aligns;
    LoadLookAligns( READ_LOCS, aligns );
    for ( size_t il = 0; il < aligns.size(); il++ ){
      look_align la = aligns[il];
      uint64_t rid = la.QueryId();
      int tid = la.TargetId();
      align a = la.a;
      vec<ho_interval> v_hos_loc;
      align2errorless( reads[rid], reftigs[tid], a, v_hos_loc, COV_TYPE );
      tig2read_hos[tid].append( v_hos_loc ); 
      sumAliLen += a.Pos1() - a.pos1();
      sumReadLen += reads[rid].isize();
    }
  }
  
  if ( sumReadLen > 0 )
    DPRINT3( sumReadLen, sumAliLen, 100 * sumAliLen / sumReadLen );



  vec< vec<int> > tigCov( reftigs.size() );
  for ( size_t rtig = 0; rtig < reftigs.size(); rtig++ )
    tigCov[rtig].resize( reftigs[rtig].size(), 0 );
  #pragma omp parallel for
  for ( size_t rtig = 0; rtig < reftigs.size( ); rtig++ ){
    for ( size_t hoi = 0; hoi < tig2read_hos[rtig].size(); hoi++ ){
      ho_interval& ho = tig2read_hos[rtig][hoi];
      if ( ho.Length() < 2 * RADIUS + 1 ) continue;
      for ( int p2 = ho.Start() + RADIUS; p2 < ho.Stop() - RADIUS; p2++ )
	tigCov[rtig][p2]++;
    }
  }

  // compute average coverage
  double aveCov = 0; 
  int64_t sumPos = 0;
  for ( size_t rtig = 0; rtig < reftigs.size( ); rtig++ )
    for ( size_t p = EDGE_SKIP; p < reftigs[rtig].size() - EDGE_SKIP; p++ ){
      sumPos++;
      aveCov += tigCov[rtig][p];
    }      
  aveCov = aveCov / sumPos;
  DPRINT2( aveCov, sumPos );
  
  int THRESH = COV_THRESH;
  if ( COV_THRESH_FRAC >= 0 )
    THRESH = floor( aveCov * COV_THRESH_FRAC );

  DPRINT( THRESH );

  
  int nRefUncov = 0;
  for ( size_t rtig = 0; rtig < reftigs.size( ); rtig++ )
    for ( size_t p = EDGE_SKIP; p < reftigs[rtig].size() - EDGE_SKIP; p++ )
      if ( tigCov[rtig][p] < THRESH ) nRefUncov++;

  DPRINT( nRefUncov );

 
  //cout << " done!" << endl;
  //cout << Date() << ": finished going through contigs" << endl;


  // mask regions if requested

  vec< vec<ho_interval> > reftig2masked_hos( reftigs.size() );

  if ( MASKING_ALIGNS_FILE.nonempty() ){
    cout << Date() << ": loading masking alignments from " << MASKING_ALIGNS_FILE << endl;
    ForceAssert( IsRegularFile( MASKING_ALIGNS_FILE ) );
    vec<look_align> aligns;
    LoadLookAligns( MASKING_ALIGNS_FILE, aligns );
    DPRINT( aligns.size() );
    for ( size_t li = 0; li < aligns.size(); li++ ){
      look_align& la = aligns[li];
      align a = la.a;
      int tid = la.TargetId();
      reftig2masked_hos.at(tid).push_back( ho_interval( la.StartOnTarget(), la.EndOnTarget() ) );
    }
  }

  
  vec< vec<int> > reftig2masked( reftigs.size() );
  for ( size_t i = 0; i < reftigs.size(); i++ )
    reftig2masked[i].resize( reftigs[i].size(), 0 );
  for ( size_t tid = 0; tid < reftig2masked_hos.size(); tid++ ){
    vec<ho_interval> masked;
    if ( ! MASK_COMPLEMENT )
      ExtractGivenCoverage( reftigs.size(), 1, reftig2masked_hos[tid], masked );
    else
      Uncovered( reftigs.size(), reftig2masked_hos[tid], masked );
    
    for ( size_t ih = 0; ih < masked.size(); ih++ ){
      for ( int p = masked[ih].Start(); p < masked[ih].Stop(); p++ )
	reftig2masked[tid][p] = 1;
    }

  }
  


  // compute and report uncovered

  size_t nuncov = 0, nall = 0, nunmasked = 0;
  vec<int> cov_data;
  for ( size_t rtig = 0; rtig < tigCov.size(); rtig++ ){
    nall += tigCov[rtig].size();
    for ( size_t p2 = EDGE_SKIP; p2 < tigCov[rtig].size() - EDGE_SKIP; p2++ ){
      if ( ! reftig2masked[rtig][p2] ){
	nunmasked++;
	cov_data.push_back( tigCov[rtig][p2] );
	if ( tigCov[rtig][p2] < THRESH ) nuncov++;
      }
    }
  }

  double mean_cov = Mean( cov_data );
  double sigma_cov = StdDev( cov_data, mean_cov );

  DPRINT5( nuncov, nunmasked, nall, mean_cov, sigma_cov );
  
  

  vec<int> covLens; 
  for ( size_t rtig = 0; rtig < tigCov.size(); rtig++ ){
    for ( size_t p2 = EDGE_SKIP; p2 < tigCov[rtig].size() - EDGE_SKIP; p2++ ){
      if ( tigCov[rtig][p2] >= THRESH ){
	int newlen = 0;
	while( p2 < tigCov[rtig].size() - EDGE_SKIP && tigCov[rtig][p2] >= THRESH ){
	  newlen++; 
	  p2++;
	}
	covLens.push_back(newlen);
      }
    }
  }

  ReverseSort( covLens );
  int n25 = -1, n50 = -1, n75 = -1, n95 = -1;
  if ( sumPos > 0 ){
    int64_t sumLen = 0;
    for ( size_t ic = 0; ic < covLens.size(); ic++ ){
      sumLen += covLens[ic];
      if ( n25 < 0 && (double)sumLen / sumPos >= 0.25 ) n25 = covLens[ic];
      if ( n50 < 0 && (double)sumLen / sumPos >= 0.50 ) n50 = covLens[ic];
      if ( n75 < 0 && (double)sumLen / sumPos >= 0.75 ) n75 = covLens[ic];
      if ( n95 < 0 && (double)sumLen / sumPos >= 0.95 ) n95 = covLens[ic];
    }
    DPRINT4( n25, n50, n75, n95 );
    PRINT3( covLens.size(), sumLen, sumPos );
  }



  
  // compute number of regions with missing coverage
  int REG_RADIUS=48;
  int nreg_miss = 0;
  vec< vec<ho_interval> > tig2uncov_hos( reftigs.size() );
  for ( size_t rtig = 0; rtig < tigCov.size(); rtig++ ){
    for ( int p2 = 0; p2 < tigCov[rtig].isize(); p2++ ){
      if ( tigCov[rtig][p2] < THRESH &&  ! reftig2masked[rtig][p2] ){
	int start = p2 > REG_RADIUS ? p2 - REG_RADIUS : 0;
	int stop  = p2 < tigCov[rtig].isize() - REG_RADIUS ? p2 : tigCov[rtig].isize();
	tig2uncov_hos[rtig].push_back( ho_interval(start, stop) );
      }
    }
  }

  
  
  for ( size_t rtig = 0; rtig < tigCov.size(); rtig++ ){
    vec<ho_interval> cov;
    Uncovered( tigCov[rtig].isize(), tig2uncov_hos[rtig], cov );
    vec<ho_interval> uncov;
    Uncovered( tigCov[rtig].isize(), cov, uncov );
    nreg_miss += uncov.isize();
    
  }
  DPRINT( nreg_miss );
  cout << Date() << ": Done with ReadBadCoverage!" << endl;
}    
    
