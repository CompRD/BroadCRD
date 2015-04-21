///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "FetchReads.h"
#include "Vec.h"
#include "VecUtilities.h"
#include "lookup/LookAlign.h"
#include "FastaFileset.h"
#include "paths/ReadLoc.h"
#include "PairsManager.h"

#include <map>

#include <omp.h>
// MakeDepend: library OMP
// MakeDepend: dependency QueryLookupTable
// MakeDepend: dependency MakeLookupTable
// MakeDepend: dependency SAM2CRDDump


class pdata{

public:
  size_t id1_, id2_;
  bool fw1_, fw2_;
  short read_class_;

  pdata( uint64_t id1, bool fw1, uint64_t id2, bool fw2 ){
    id1_ = id1; id2_ = id2;
    fw1_ = fw1; fw2_ = fw2;
  }
  
};

Bool cmp_target( const look_align& la1, const look_align& la2 ){    
  if ( la1.target_id < la2.target_id ) return True;
  if ( la1.target_id > la2.target_id ) return False;
  if ( la1.pos2( ) < la2.pos2( ) ) return True;
  if ( la1.pos2( ) > la2.pos2( ) ) return False;
  if ( la1.query_id < la2.query_id ) return True;
  return False;    
}

bool cmp_errors( const look_align& la1, const look_align& la2 ){    
  if ( la1.Errors() < la2.Errors() )
    return true;
  else return false; 
}

/**
 * ProbeMissingRefDuplications
 * Suppose you observed a duplication (inversion) in your assembly that is not clearly
 * identified by direct alignment to a reference. This code checks if new sequencing data
 * (uses jumping read pairs for human and mouse) suggest valid polimorphic difference 
 * between new assembly and the old reference. It looks at the number of broken jumping-pair
 * alignments.
 
*/


int main( int argc, char *argv[] ){
  
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String( SUBDIR );
  
  CommandArgument_String( INVERSIONS );
  CommandArgument_String( GENOME_NAME );
  CommandArgument_String_OrDefault( GENOME, "genome" );
  CommandArgument_String_OrDefault( ADJUST_SUPER, "" );
  CommandArgument_String_OrDefault( ADJUST_RANGE, "" );
  CommandArgument_Bool_OrDefault( MAPPED_PAIRS_ONLY, False );
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  CommandArgument_Bool_OrDefault( REBUILD_LOOKUP, False );
  CommandArgument_Bool_OrDefault( VERBOSE, False );
  EndCommandArguments;

  // Thread control (OMP used in FixScaffoldsCore)
  
  NUM_THREADS = configNumThreads(NUM_THREADS);
  omp_set_num_threads( NUM_THREADS );

  String data_dir = PRE + "/" + DATA + "/";
  String run_dir  = data_dir + RUN + "/";
  String sub_dir   = run_dir + SUBDIR + "/";
  String work_dir  = sub_dir + "/inversions/";

  if ( ! IsRegularFile(work_dir) )
    Mkdir777( work_dir );
  

  // Read in inversion coordinates
  
  String inversions_file = work_dir + INVERSIONS;
  String fileOut = work_dir + INVERSIONS.SafeBeforeLast(".") + ".confirmations";
  PRINT( inversions_file );
  vec< vec<int> > inversions_in;
  {
    ForceAssert( IsRegularFile( inversions_file ) );
    ifstream ifin( inversions_file.c_str() );
    String line;
    while( getline(ifin,line) ){
      vec<String> fields;
      Tokenize( line, '\t', fields ); //tig1 \t start1 \t end1 \t tig2 \t start2 \t end2\n";
      vec<int> values;
      for ( size_t i = 0; i < fields.size(); i++ ){
	values.push_back( atoi(fields[i].c_str()) );
      }
      inversions_in.push_back( values );
    }
  }

  //String genome_file = data_dir + GENOME + ".fastb";
  //vecbasevector reference( genome_file );

  String sheared_jump_lanes = "";
  if ( GENOME_NAME == "mouse" )
    sheared_jump_lanes = "\"{61NDRAAXX.{1,2,3,4,5,6,7},61NGGAAXX.{2,3,4,5,6},61NCCAAXX.{1,2,3,4}}\"";
  else if ( GENOME_NAME == "human" )
    sheared_jump_lanes = "\"2025JABXX.{1,2,3,4,5,6,7,8}\"";
  else
    FatalErr("GENOME = " + GENOME_NAME + " is not supported");

  String picard_head = "/seq";

  String extract_args = " PICARD_HEAD=" + picard_head + " TMP=" + work_dir + " WRITE_NAMES=False WRITE_ALIGNS=True UNMAPPED=False ";
  
  String libinfo = "/wga/scr2/dexter/libinfo/dexter_libs";
  

  PRINT( inversions_in.size() );
  ofstream fout( fileOut.c_str() );
  int extensionSize = 10000;
  for ( size_t i = 0; i < inversions_in.size(); i++ ){
    if ( inversions_in[i].size() < 1 ){
      fout << endl;
      continue;
    }
    int super = ADJUST_SUPER == "" ? inversions_in[i][0] : inversions_in[i][0] + atoi(ADJUST_SUPER.c_str());
    int r0 = ADJUST_RANGE == "" ? inversions_in[i][1] : inversions_in[i][1] + atoi(ADJUST_RANGE.c_str());
    int r1 = ADJUST_RANGE == "" ? inversions_in[i][2] : inversions_in[i][2] + atoi(ADJUST_RANGE.c_str());
    
    int ext0 = r0 > extensionSize ? extensionSize : r0;
    int ext1 = extensionSize;
  

    String super_range = "chr" + ToString(super) + ":" + ToString(r0 - ext0) + "-" + ToString(r1 + ext1);

    String target_head = work_dir + "jumps_r" + ToString(i);

    // Get read data associated with inversions
    PRINT2( super_range, target_head );
    SystemSucceed("ExtractFromBAM.pl " + extract_args
		  + ARG(TARGET_HEAD, target_head)
		  + ARG(RANGE, super_range)		   
		  + ARG(LANES, sheared_jump_lanes)
		  + ARG(VERSION, "newest" )	
		  + ARG(REQUIRE_LIBINFO, "True")
		  + ARG(WRITE_ALIGNS, "True")
		  + ARG(WRITE_NAMES, "True")
		  + ARG(WRITE_PAIRS, "True")
		  + ARG(MAPPED_PAIRS_ONLY, MAPPED_PAIRS_ONLY)
		  + ARG(LIBINFO,libinfo)			
		  + ">&" +  work_dir + "ExtractFromBAM.r" + ToString(i) + ".out" );
    
    
    String aligns_file = target_head + ".qltout";
    vec<look_align> aligns;
    LoadLookAligns( aligns_file, aligns );
    

    PairsManager pairs( target_head + ".pairs" );
    size_t n_reads = pairs.nReads();
    
    vec<size_t> aligns_index( n_reads, n_reads );
    for ( size_t i = 0; i < aligns.size(); i++ ){
      size_t qid = aligns[i].QueryId();
      ForceAssertLt( qid, n_reads );
      ForceAssertEq( aligns_index[qid], n_reads );
      aligns_index.at(qid) = i;
    }

    int nPaired = 0, nUnpaired = 0;
    int nunstretched = 0, nproper= 0, nstretched = 0, nhalf = 0; 
    vec<bool> idUsed( pairs.nReads(), false );
    for ( size_t id1 = 0; id1 < pairs.nReads(); id1++ ){
      if ( idUsed.at(id1) ) continue;
      else idUsed.at(id1) = true;

      if ( pairs.isPaired(id1) ){
	longlong id2 = pairs.getPartnerID( id1 );
	longlong pid = pairs.getPairID( id1 );
	
	if ( idUsed.at(id2) ) continue;
	else idUsed.at(id2) = true;

	if ( aligns_index[id1] < n_reads && aligns_index[id2] < n_reads ){
	  const look_align& la1 = aligns[ aligns_index[id1] ];
	  const look_align& la2 = aligns[ aligns_index[id2] ];
	  if ( ( la1.StartOnTarget() >= r0 && la1.EndOnTarget() <= r1 && (la2.StartOnTarget() > r1 || la2.EndOnTarget() < r0 ) ) ||
	       ( la2.StartOnTarget() >= r0 && la2.EndOnTarget() <= r1 && (la1.StartOnTarget() > r1 || la1.EndOnTarget() < r0 ) )
	       ){
	    nPaired++;
	    //both aligned
	    int sep = pairs.sep(pid);
	    int sd = pairs.sd(pid);
	    if ( abs(la1.StartOnTarget() - la2.StartOnTarget()) < sep + 3 * sd ){
	      if ( la1.Fw1() != la2.Fw1() ) nproper++;
	      nunstretched++;
	    }else nstretched++;
	  }
	}else if ( aligns_index[id1] < n_reads ||  aligns_index[id2] < n_reads ){
	  const look_align& la = aligns_index[id1] < n_reads ? aligns[ aligns_index[id1] ] : aligns[ aligns_index[id2] ];
	  if ( la.StartOnTarget() >= r0 && la.EndOnTarget() <= r1 ){
	    nPaired++;
	    nhalf++;
	  }
	}
      }else{
	const look_align& la = aligns[ aligns_index[id1] ];
	if ( la.StartOnTarget() >= r0 && la.EndOnTarget() <= r1 ){
	  nUnpaired++;
	  nhalf++;
	}
      }
    }
    PRINT2( pairs.nReads(), aligns.size() );
    PRINT6( nunstretched, nproper, nstretched, nhalf, nPaired, nUnpaired );
    PRINT6_TO( fout, nunstretched, nproper, nstretched, nhalf, nPaired, nUnpaired );
  }

  cout << Date( ) << ": done" << endl;

}


