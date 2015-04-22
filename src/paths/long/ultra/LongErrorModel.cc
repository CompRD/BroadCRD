///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Alignment.h"
#include "CoreTools.h"
#include "VecUtilities.h"
#include "lookup/LookAlign.h"
#include "pairwise_aligners/LocalAlign.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "paths/LongReadTools.h"
#include "paths/long/ultra/ConsensusScoreModel.h"
#include "paths/long/Friends.h"
#include "paths/long/Logging.h"
#include "paths/long/ultra/GetFriendsAndAlignsInitial.h"
#include "paths/long/ultra/LongErrorModel.h"
#include "paths/long/LongProtoTools.h"

void error_model::GetAlignData( int& subs, int& inserts, int& dels, longlong& alen ){

  int64_t ng = gangs_->size( );
  vec<int> subsx(ng, 0), insertsx(ng, 0), delsx(ng, 0);
  vec<longlong> alenx(ng, 0);
  #pragma omp parallel for
  for ( size_t gi = 0; gi < (*gangs_).size(); gi++ ){
    const vecbasevector& reads = (*gangs_)[gi];
    //#pragma omp critical 
    //{ PRINT( reads.size() ); }
    if ( reads[0].size() == 0 ) continue;
    vec< vec<edit0> > edits; 
    edits.clear_and_resize( reads[0].size() +1 );
    vec<int> coverage( reads[0].size() +1, 0 );
    int locdels = 0, locins = 0, locsubs = 0;
    align a;
    look_align la;
    for ( size_t i = 1; i < reads.size( ); i++ ){
      if ( reads[i].size() == 0 ) continue; 
      //LocalAlign( reads[i], reads[0], a, 1, -1, -1 );
      int errors1;
      int offset = ( reads[i].isize() - reads[0].isize() ) / 2;
      int bandwidth = Max( reads[0].isize(), reads[i].isize() ) / 2;
      SmithWatBandedA( reads[i], reads[0], offset, bandwidth, a, errors1, 0, 1, 1 );

      int errors = a.Errors(reads[i],reads[0]);
      la.ResetFromAlign( a, reads[i], reads[0]  );
      double errorRate = la.ErrorRate();
      if ( errorRate > 0.2 ) continue; // skip bad alignments

      int galen = a.Pos1() - a.pos1() + a.Pos2() - a.pos2() + la.Indels();
      
      int gamin = round( (double)galen * 0.1 ) -1;
      gamin = Max( gamin, 0 );
      int gamax = galen - gamin;
      gamax = Min( gamax, galen -1 );
      //#pragma omp critical
      //{ PRINT3( galen, gamin, gamax ); }
      int gapos = 0;
      int p1 = a.pos1( ), p2 = a.pos2( );
      for ( int j = 0; j < a.Nblocks( ); j++ ) {    
	if ( a.Gaps(j) > 0 ){    
	  if ( gapos >= gamin && gapos <= gamax ){
	    locdels++;
	    edits[p2].push( DELETION, a.Gaps(j) );
	    for ( int i = 0; i <= a.Gaps(j); i++ )
	      coverage[p2-i]++;
	  }
	  p2 += a.Gaps(j);
	  gapos += a.Gaps(j);
	}
	if ( a.Gaps(j) < 0 ){ 
	  if ( gapos >= gamin && gapos <= gamax ){
	    locins++;
	    basevector b( reads[i], p1, -a.Gaps(j) );
	    edits[p2].push( INSERTION, b.ToString( ) );
	  }
	  p1 -= a.Gaps(j);
	  gapos -= a.Gaps(j);
	}
	for ( int x = 0; x < a.Lengths(j); x++ ) {
	  if ( gapos >= gamin && gapos <= gamax ){
	    if ( reads[0][p2] != reads[i][p1] ){
	      edits[p2].push( SUBSTITUTION, (char) as_base( reads[i][p1] ) );
	      locsubs++;
	    }
	    coverage[p2]++;
	  }
	  p1++, p2++;    
	  gapos++;
	}    
      }    

    }

    /*
    #pragma omp critical
    { //PRINT3( locsubs, locdels, locins );
       int checkdels = 0, checkins = 0, checksubs = 0;
       for ( size_t ei = 0; ei < edits.size(); ei++ ){
	 int s = 0, i = 0, d = 0;
	 for ( size_t pi = 0; pi < edits[ei].size(); pi++ ){
	   if ( edits[ei][pi].etype == INSERTION ){ 
	     checkins++; i++;
	   }else if (  edits[ei][pi].etype == DELETION ){
	     checkdels++; d++;
	   }else if (  edits[ei][pi].etype == SUBSTITUTION ){
	     checksubs++; s++;
	   }
	 }
	 //if ( s + d + i > 0.2 * coverage[ei] )
	 //PRINT5( ei, s, i, d, coverage[ei] );
       }
       //PRINT3( checksubs, checkdels, checkins );
    }
    */
   
    for ( size_t ei = 0; ei < edits.size(); ei++ ){
      if ( coverage[ei] < 10 ) continue;
      alenx[gi]++;
      if ( edits[ei].size() > 0 ){
	Sort( edits[ei] );
	for ( size_t pi = 0; pi < edits[ei].size(); pi++ ){
	  size_t pn = pi + 1;
	  for ( pn = pi +1; pn < edits[ei].size(); pn++){
	    if ( edits[ei][pi].etype == INSERTION && edits[ei][pn].etype != INSERTION )
	      break;
	    else if ( edits[ei][pi].etype == DELETION && edits[ei][pn].etype != DELETION )
	      break;
	    else if ( edits[ei][pi].etype == SUBSTITUTION && 
		      ( edits[ei][pi].etype != SUBSTITUTION || edits[ei][pi].seq != edits[ei][pn].seq ) )
	      break;
	  }
	  {   
	    if ( (double)(pn - pi) > 0.5 * (double)coverage[ei] ){
	      if ( edits[ei][pi].etype == INSERTION ){ 
		delsx[gi]++;
	      }else if (  edits[ei][pi].etype == DELETION ){
		insertsx[gi]++;
	      }else if (  edits[ei][pi].etype == SUBSTITUTION ){
		subsx[gi]++;
	      }
	      else FatalErr("What the heck" );
	    }
	  }
	  pi = pn -1;
	}
      }
    }
  }

  for ( int64_t gi = 0; gi < ng; gi++ )
  {    alen += alenx[gi];
       subs += subsx[gi];
       inserts += insertsx[gi];
       dels += delsx[gi];    }

}

template<int K> void DefineErrorModel( const unsigned int RANDOM_SEED, 
     const int READ_SAMPLE_SIZE, const vecbasevector& reads, const IAndOsVec& F,
     double& psubs, double& pinserts, double& pdels,
     const long_heuristics& heur, const long_logging_control& log_control,
     const long_logging& logc )
{
     int N = reads.size( );
     vec<vecbasevector>  gangs;
     srand(RANDOM_SEED);
     vec<int> rids = vec<int>( N, vec<int>::IDENTITY );
     random_shuffle( rids.begin(), rids.end() );
     int sample_size = N;
     if ( READ_SAMPLE_SIZE > 0 ) sample_size = Min( READ_SAMPLE_SIZE, N );
     else
     {    // lowest error rate we care 10-3
          // estimation error 10% therefore trying to get 100 events
          int lenThresh = 100 / 0.001, locLen = 0;
          for ( int rxi = 0; rxi < rids.isize(); rxi++ )
          {    locLen += reads[ rids[rxi] ].size();
	       sample_size = rxi + 1;
	       if ( locLen > lenThresh ) break;    }    }
     DPRINT( sample_size );
     #pragma omp parallel for
     for ( int rxi = 0; rxi < sample_size; rxi++ )
     {    int xi = rids[rxi];
          vecbasevector gang;
          vec< vec< pair<int,int> > > a;
          long_logging_control log_control_local(log_control);
          long_logging logc_local( "" );
          long_heuristics heur_local(heur);
          heur_local.COLUMN_FILTER = False;
          GetFriendsAndAlignsInitial<K>( F, reads, xi, gang, 0, a, cout, 
               ConsensusScoreModel(0.1, 0.1, 0.1), heur_local, log_control_local,
               logc_local );
          #pragma omp critical
          {    gangs.push_back(gang);    }    }
     Sort(gangs);
     DPRINT( gangs.size( ) );
     cout << Date() << ": building error model" << endl;
     error_model LEM( gangs );
     int subs = 0, inserts = 0, dels = 0;
     longlong alen = 0;
     LEM.GetAlignData( subs, inserts, dels, alen );
     psubs    = (double)subs/(double)alen;
     pinserts = (double)inserts/(double)alen;
     pdels    = (double)dels/(double)alen;
     cout << Date( ) << ": estimated err rates, " << setiosflags(ios::fixed) 
          << setprecision(2) << "del = " << 100.0 * pdels << "%, ins = " 
          << 100.0 * pinserts << "%, sub = " << 100.0 * psubs << "%" 
          << resetiosflags(ios::fixed) << endl;    }

template void DefineErrorModel<16>( const unsigned int RANDOM_SEED, 
     const int READ_SAMPLE_SIZE, const vecbasevector& reads, const IAndOsVec& F,
     double& psubs, double& pinserts, double& pdels, const long_heuristics& heur, 
     const long_logging_control& log_control, const long_logging& logc );

template void DefineErrorModel<20>( const unsigned int RANDOM_SEED, 
     const int READ_SAMPLE_SIZE, const vecbasevector& reads, const IAndOsVec& F,
     double& psubs, double& pinserts, double& pdels, const long_heuristics& heur, 
     const long_logging_control& log_control, const long_logging& logc );

template void DefineErrorModel<24>( const unsigned int RANDOM_SEED, 
     const int READ_SAMPLE_SIZE, const vecbasevector& reads, const IAndOsVec& F,
     double& psubs, double& pinserts, double& pdels, const long_heuristics& heur, 
     const long_logging_control& log_control, const long_logging& logc );

template void DefineErrorModel<28>( const unsigned int RANDOM_SEED, 
     const int READ_SAMPLE_SIZE, const vecbasevector& reads, const IAndOsVec& F,
     double& psubs, double& pinserts, double& pdels, const long_heuristics& heur, 
     const long_logging_control& log_control, const long_logging& logc );
