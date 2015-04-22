///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef C_NHOOD_EVAL_H
#define C_NHOOD_EVAL_H

#include "SeqInterval.h"
#include "TokenizeString.h"
#include "paths/HyperBasevector.h"
#include "system/System.h"

/**
 * class CNhoodEval
 *
 * Class to parse the output files from a local assembly, and to
 * evaluate a local assembly. For every edge in the nhood, there may
 * be two different types of placement: the cloud placement, and the
 * true placement:
  *   - The cloud placement is the one predicted by the linking
  *     information, is a system that puts the begin of the seed at
  *     base zero.
  *   - The true placement is only found if RunAllPathsLG was run with
  *     evalulation on, and it is given by the alignment of the edges
  *     onto the reference.
  * Alternately, you can run the method Eval( ), to realign the edges
  * onto the reference.
  */
 class CNhoodEval {
   
 public:
   
   CNhoodEval( ) : seed_id_ ( -1 ) { }
   
   CNhoodEval( const String &log_file ) { this->ParseLog( log_file ); }

   // Parse log file to capture core info.
   void ParseLog( const String &log_file );

   // Id of seed (as an entry of edge_id_), and id of unipath seed.
   int SeedId( ) const { return seed_id_; }
   int SeedUnipathId( ) const { return seed_id_ < 0 ? -1 : edge_id_[seed_id_]; }

   // Number of edges in this nhood.
   int Size( ) const { return edge_id_.isize( ); }

   // Constant accessors for the edge with local id ii.
   int EdgeId( int ii ) const { return edge_id_[ii]; }

   int CloudBegin( int ii ) const { return begin_[ii]; }
   int CloudEnd( int ii ) const { return end_[ii]; }
   int CloudLength( int ii ) const { return end_[ii] - begin_[ii]; }
   int CloudCopyNumber( int ii ) const { return cn_[ii]; }
   
   bool TrueLocsExist( ) const { return true_target_id_.size( ) > 0; }
   int TrueTargetId( int ii ) const { return true_target_id_[ii]; }
   int TrueBegin( int ii ) const { return true_begin_[ii]; }
   int TrueEnd( int ii ) const { return true_end_[ii]; }
   bool TrueRc( int ii ) const { return true_rc_[ii]; }
   int TrueCopyNumber( int ii ) const { return true_cn_[ii]; }
   
   // Begin and end of nhood.
   pair<int,int> CloudWin( ) const;
   pair<int,int> TrueWin( ) const;

   // Begin and end on reference (based on signed alignment of seed).
   pair<int,int> ImpliedOnRef( int signed_pos ) const;
   
   // Info for this nhood (mimics the output from the nhood log file).
   void PrintInfo( ostream &out ) const;
   
   // Evaluate local assembly. It will generate files parallel to head_.
   String Eval( const vecbvec &bases,
		const String &unibases_fastb,
		const triple<int64_t,int64_t,int> &align ) const;
   
   // How much of the target is covered by the given aligns.
   double Coverage( int offset, int tlen , const String &qltout ) const;
   
   
 private:
   
   String head_;       // head of in/out file names
   int seed_id_;       // edge_id_[seed_id_] = unipath id of seed
   vec<int> edge_id_;  // ids of edges
   
   vec<int> begin_;    // cloud locs (always filled)
   vec<int> end_;
   vec<int> cn_;
   
   vec<int> true_target_id_;   // true locs (only filled if eval code was run)
   vec<int> true_begin_; 
   vec<int> true_end_;
   vec<bool> true_rc_;
   vec<int> true_cn_;
   
 };

#endif
