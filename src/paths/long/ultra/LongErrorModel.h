///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef LONG_ERROR_MODEL_H
#define LONG_ERROR_MODEL_H

#include "Basevector.h"
#include "CoreTools.h"
#include "polymorphism/Edit.h"
#include "paths/long/Friends.h"
#include "paths/long/Logging.h"
#include "paths/long/LongProtoTools.h"

class error_model {

     public:

     error_model( ) { }
     error_model( const vec<vecbasevector>& gangs ){
       gangs_ = &gangs;
       alive_.resize( gangs.size() );
       for ( size_t i = 0; i < gangs.size(); i++ )
	 alive_[i] = vec<Bool>( gangs[i].size(), True );
     }

     Bool Alive( const int gid, const int rid ) const { return alive_[gid][rid]; }
     void GetAlignData( int& subs, int& inerts, int& dels, longlong& alen );
     
     
     private:

     const vec<vecbasevector>* gangs_;
     vec< vec<Bool> > alive_;

};

template<int K> void DefineErrorModel( const unsigned int RANDOM_SEED, 
     const int READ_SAMPLE_SIZE, const vecbasevector& reads, const IAndOsVec& F,
     double& psubs, double& pinserts, double& pdels, const long_heuristics& heur, 
     const long_logging_control& log_control, const long_logging& logc );

#endif
