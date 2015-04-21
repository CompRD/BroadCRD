///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "ParallelVecUtilities.h"
#include "efasta/EfastaTools.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/ultra/ConsensusScoreModel.h"
#include "paths/long/Friends.h"
#include "paths/long/ultra/LongErrorModel.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/ultra/MakeCorrectedReads.h"

template<int K> void MakeCorrectedReadsCore( const vecbasevector& reads,
     const unsigned int RANDOM_SEED, const String& RID, const Bool CHEAT_MODEL,
     const unsigned int NUM_THREADS, 
     const double ERR_SUB, const double ERR_INS, const double ERR_DEL,
     const long_heuristics& heur, const long_logging_control& log_control, 
     const long_logging& logc, const ref_data& ref, VecEFasta& corrected, 
     vec<int>& cid, String const& TMP )
{
     double clock = WallClockTime( );
     int N = reads.size( );

     // Find friends.

     Mkdir777(TMP);
     String friendsFile(TMP+"/frag_reads_orig.friends");
     FindFriends(reads,friendsFile,heur.USE_BROKEN_FF);
     IAndOsVec F(friendsFile);

     // Align the reads to the genome.  Define kmer alignments to read 0. 
     // (Using genome!!)  This is not used currently and is here just as a
     // reminder of how to do this.

     /*
     vec< vec< pair<int,int> > > a;
     BuildKmerAlignsUsingGenome( K, reads, g, a );
     */
     REPORT_TIME( clock, "used in setup" );

     // Estimate error model.

     double mclock = WallClockTime( );
     double psubs = -1.0, pinserts = -1.0, pdels = -1.0;
     if ( !CHEAT_MODEL )
     {    DefineErrorModel<K>( RANDOM_SEED, heur.READ_SAMPLE_SIZE, reads, F, 
               psubs, pinserts, pdels, heur, log_control, logc );    }
     else
     {    psubs = ERR_SUB;
          pinserts = ERR_INS;
          pdels = ERR_DEL;    }
     ConsensusScoreModel error_model; // The model to score consensus sequences. 
     error_model.Init(pdels, pinserts, psubs);
     REPORT_TIME( mclock, "used creating error model" );

     // Go through the reads.
     
     vec<int> rid;
     if ( RID != "" ) ParseIntSet( RID, rid, true, 0, N );
     else rid = vec<int>( N, vec<int>::IDENTITY );
     BuildCorrectedReads<K>( reads, F, rid, corrected, cid, 
          error_model, heur, log_control, logc, ref, NUM_THREADS );

     // Sort reads, to fix random order arising from parallel loop.

     double sclock = WallClockTime( );
     ParallelSortSync( cid, corrected );
     REPORT_TIME( sclock, "used sorting corrected reads" );    }

void MakeCorrectedReads( const long_sim& sim, const vecbasevector& reads,
     const String& RID, const unsigned int NUM_THREADS, 
     const long_heuristics& heur, const long_logging_control& log_control, 
     const long_logging& logc, const ref_data& ref, VecEFasta& corrected, 
     vec<int>& cid, String const& TMP )
{
     if ( !heur.CORRECT )
     {    double pclock = WallClockTime( );
          for ( int i = 0; i < (int) reads.size( ); i++ )
          {    corrected.push_back( efasta( reads[i].ToString( ) ) );
               cid.push_back(i);    }
          REPORT_TIME( pclock, "used pushing reads" );
          return;    }

     int K = heur.K;

     if ( K == 16 )
     {    MakeCorrectedReadsCore<16>( reads, sim.RANDOM_SEED, RID, sim.CHEAT_MODEL, 
               NUM_THREADS, sim.ERR_SUB, sim.ERR_INS, sim.ERR_DEL, heur, 
               log_control, logc, ref, corrected, cid, TMP );    }
     else if ( K == 20 )
     {    MakeCorrectedReadsCore<20>( reads, sim.RANDOM_SEED, RID, sim.CHEAT_MODEL, 
               NUM_THREADS, sim.ERR_SUB, sim.ERR_INS, sim.ERR_DEL, heur, 
               log_control, logc, ref, corrected, cid, TMP );    }
     else if ( K == 24 )
     {    MakeCorrectedReadsCore<24>( reads, sim.RANDOM_SEED, RID, sim.CHEAT_MODEL, 
               NUM_THREADS, sim.ERR_SUB, sim.ERR_INS, sim.ERR_DEL, heur, 
               log_control, logc, ref, corrected, cid, TMP );    }
     else if ( K == 28 )
     {    MakeCorrectedReadsCore<28>( reads, sim.RANDOM_SEED, RID, sim.CHEAT_MODEL, 
               NUM_THREADS, sim.ERR_SUB, sim.ERR_INS, sim.ERR_DEL, heur, 
               log_control, logc, ref, corrected, cid, TMP );    }
     else ForceAssert( 0 == 1 );    }
