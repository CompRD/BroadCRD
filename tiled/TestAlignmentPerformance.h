/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


#include "TaskTimer.h"
#include "random/Random.h"
#include "tiled/AlignsOnConsensus.h"

// A function which measures the time taken to align a small subset of reads in a 
// contig to the contig using AlignsOnConsensus, and if the runtime is greater than 
// the threshold, return 'true'.  

// This should identify those contigs for which a full alignment via AlignsOnConsensus 
// would take an unreasonable amount of time.
//
// Current parameters were established from testing on the Horse assembly so
// that only the worst-performing contig returned false.
//
//

bool TestAlignmentPerformance(const int contig_id,
			      const vecbasevector &consensus_bases,
			      const vecqualvector &consensus_quals,
			      const vecbasevector &read_bases,
			      const vecqualvector &read_quals,
			      const vec<int> &first_locs,
			      const vec<read_location> &locs)
{

  // parameters
  int numReadsInSubset(200);
  int minNumReadsInContig(1000);
  float sysRunTimeThreshold(0.1);



  int ii(first_locs[contig_id]);
  for (; ii<(int)locs.size( ); ii++)
    if ( locs[ii].Contig( ) != contig_id )
      break;
      
  int numreads(ii-first_locs[contig_id]);
  if ( numreads > minNumReadsInContig )
  {
    TaskTimer t;
    t.Start();
	
    vec<t_align> dummy_aligns;
    AlignsOnConsensus( dummy_aligns,
		       contig_id,
		       consensus_bases,
		       consensus_quals,
		       read_bases,
		       read_quals,
		       first_locs,
		       locs,
		       false,
		       128,
		       1,
		       false,
		       0,
		       0,
		       true,
		       numReadsInSubset);
    
    t.Stop();
    
    if ( t.SysSecs() > sysRunTimeThreshold )
      return true;
    
  }
  return false;
}
