///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef BTL__C_MONO_REP__H
#define BTL__C_MONO_REP__H

#include "Basevector.h"
#include "ReadLocation.h"
#include "String.h"
#include "Vec.h"
#include "lookup/LookAlign.h"
#include "tiled/TAlign.h"

#include "Corrector.h"
#include "BaseProb.h"
#include "ReadScore.h"
#include "math/Functions.h"
#include "tiled/CharBaseUtil.h"
#include "tiled/CSnp.h"
#include "tiled/ScHelp.h"
#include "tiled/Tiling.h"
#include "tiled/TilingPrinter.h"
#include "tiled/VecTiling.h"

/**
 * class CConsensus
 *
 * Given an initial consensus draft, and alignments of reads onto the
 * draft, it will generate fixes (as in class Corrector) for both bases
 * and quality scores.
 *
 * n_haplotypes_: only allowed are 1 or 2
 * win_size_: size of (non-overlapping) windows
 * max_pads_ratio_: discard reads with too many pads wrt local consensus
 * good_score_: read chunks are scored with as in ReadScore.h
 * verbose_: toggle verbose mode when logging
 */
class CConsensus {

public:
  
  CConsensus( int n_haplotypes, ostream *log = 0 );
  
  ~CConsensus( );

  void SetTilingAnalyzer( tiling_analyzer *analyzer,
			  const vec<Bool> *chimeras = 0 );
  
  void SetVerbose( bool verbose );
  
  void SetDefaults( );
  
  // Generate a set of corrections for the given contig.
  void GenerateFixes( vec<corrector> &fixes );
  
  
private:
  
  void ClearAll( );
  
  void ClearTempData( );

  void BuildWindows( );
  
  // Decide the length of the next window.
  int NextWinSize( int start );

  // Cluster reads and select candidates to lead consensus.
  void SelectLeadReads( int win_id );
  
  // Use selected lead reads to drive consensus.
  void LeadReadsToConsensus( );

  // Consensus at the given position in the rectangle.
  void ToConsensusCore( int pos, int look_at, bool is_poor );

  // Assign a quality score to pads by looking at the flanking bases.
  bool ReadBaseQual( int rPos, int bPos, char &base, int &qual );
  
  // Do not take action if old and new bases match.
  void AddFix( corrector &corr );
  
  // Print current rectangle.
  void PrintCurrentWin( ) const;

  // Weak test for good stretch (of read).
  bool IsGoodWeak ( const vec<char> &read, const vec<char> &consensus ) const;

  // Strong test for good stretch (strictly stronger than IsGoodWeak).
  bool IsGoodStrong( const vec<char> &read, const vec<char> &consensus ) const;
  
  
private:
  
  int n_haplotypes_;       // only allowed: 1 or 2
  int win_size_;           // minimum size of windows
  double max_pads_ratio_;  // to decide if a chunk of read is too padded
  double good_score_;      // to decide if a chunk of read is good qual
  bool verbose_;           // if verbose log
  ostream *log_;           // log stream (may be null)
  
  tiling_analyzer *analyzer_;  // input draft consensus (as tiling)
  const vec<Bool> *chimeras_;  // optional pointer for chimeric reads
  
  pair<int,int> master_win_;   // overall master window
  vec< pair<int,int> > win_;   // partition of master in non overlapping chunks
  vec<corrector> fixes_;       // set of corrections

  pair<int,int> winsel_;       // current selected window (from win_)      
  vec< vec<char> > rbases_;    // current rectangle (bases)
  vec< vec<int> > rquals_;     // current rectangle (quals)
  vec< int > read_pos_;        // current rectangle (read position)
  vec< scHelp > helpers_;      // score helpers (helpers_[0] is consensus)
  vec< vec<int> > helper_pos_; // reads partitioned and sorted by score
  vec< float > cscores_;       // cumulative score for helper_pos_[ii]

};

/**
 * PrintErrorCsvLine
 *
 * It generates a csv one-liner parallel to the one dumped by
 * EvalOnRef.
 */
inline void PrintErrorCsvLine( const String &csvline_file,
			       const String &name,
			       const String &ref )
{
  ofstream out( csvline_file.c_str( ) );

  if ( ref != "" ) out << ",,,," << name << ",ref,,,,,,FAILED\n";
  out << "0,,,," << name << ",,,,,,,FAILED\n";

  out.close( );
}

/**
 * SAligner
 *
 * Structure needed mainly to sort alignments of monomers on sequences.
 */
struct SAligner {
  
  SAligner( int offset, int len, int matches, float score ) :
    offset_ ( offset ),
    len_ ( len ),
    matches_ ( matches ),
    score_ ( score ) { }

  friend bool operator< ( const SAligner &left, const SAligner &right ) {
    if ( left.matches_ > right.matches_ ) return true;
    if ( left.matches_ < right.matches_ ) return false;
    if ( left.score_ < right.score_ ) return true;
    if ( left.score_ > right.score_ ) return false;
    return left.len_ > right.len_;
  }
  
  friend ostream &operator<< ( ostream &out, const SAligner &al ) {
    out << al.offset_ << "\t"
	<< al.len_ << "\t"
	<< al.matches_ << "\t"
	<< ToString( al.score_, 2 );
    return out;
  }

  int offset_;       // offset in monomer space
  int len_;          // length of alignment
  int matches_;      // number of matches
  float score_;      // alignment score (lower is better)
  
};



/**
 * class CMonoRep
 *
 * Class to represent a sequence as a list of monomers ids, for the
 * validation of TALEN products, as:
 *
 *    n0, id0,  n1, id1, ...   nN, id1N, n(N+1)
 *   gap,  id, gap,  id, .... gap,   id,    gap
 *
 * For example:
 *
 *    32, 1, 0, 3, 0, 2, 0, 2, ... 0, 3, 27
 *
 * n_i's are called "gaps", and they are the number of bases of
 * separation between monomers. They must start and finish a
 * sequence. In the example, there are 32 bases before the first
 * monomer in the sequence (which has id 1); there are 27 bases after
 * the last monomer in the sequence (id 3); and there are no "extra"
 * bases between any two adjacent monomers.
 *
 * id_i's are the numeric ids of the monomer in the sequence.
 *
 * WARNING: rc1 in chain_ is used to store the orientation of the
 *          target sequence, not the orientation of the aligning
 *          monomer (aligns are "Reversed", if rc1 is true).
 */
class CMonoRep {

public:
  
  CMonoRep( ) : name_( "" ), log_ ( 0 ) { }
  
  CMonoRep( const String &name, ostream *log = 0 );

  void Clear( );

  // Set this from a string representing the reference in monomer space.
  void SetFromReference( const String &ref,
			 const vecbvec &primers,
			 const vecbvec &parts,
			 const vecString &parts_ids );
  
  // Set this from a list of aligns of monomers on sequence(s).
  void SetFromHits( const vec<look_align_plus> &hits,
		    const vecbvec &parts,
		    const bvec &bases,
		    const qvec &quals,
		    const int &tid,
		    const bool rc = false );
  
  // Add an initial gap...
  void AddInitialGap( int gap0 );

  // ... and a whole new segment.
  void AddSegment( const bvec &bases,
		   const qvec &quals,
		   const look_align_plus &hit,
		   const float &score,
		   const int &gap,
		   const int &id );

  // Set this to the CMonoRep generated by merging left and right.
  bool Merge( const CMonoRep &left, const CMonoRep &right );

  // Refine consensus (bases and quals only! See .cc for details).
  bool Refine( const CMonoRep &fw1,
	       const CMonoRep &fw2,
	       const CMonoRep &rc,
	       const String &out_head );

  // Const accessors.
  const String &Name( ) const { return name_; }
  const bvec &Bases( ) const { return bases_; }
  const qvec &Quals( ) const { return quals_; }
  bool IsEmpty( ) const { return ids_.size( ) < 1; }
  bool IsAwful( ) const { return scores_.size()>0 && Min(scores_)>BAD_QUAL; }
  size_t NMonomers( ) const { return ids_.size( ); }
  look_align_plus Chain( int ii ) const { return chain_[ii]; }
  float Score( int ii ) const { return scores_[ii]; }
  int Gap( int ii ) const{ return gaps_[ii]; }
  int Id( int ii ) const { return ids_[ii]; }
  
  // Save bases and quals only.
  void SaveBasesAndQuals( const String &head ) const;

  // Block length, primer included (inferred from chain_).
  int BlockLength( int ii ) const;

  // Evaluate against reference (print on visual align of consensus on ref).
  void EvalOnRef( const CMonoRep &ref,
		  const vecString &parts_ids,
		  const String &name,
		  const String &csv_file,
		  ostream &visual ) const;
  
  // Simple minded Print.
  void Print( const vecString &parts_ids, ostream &out ) const;

  // Verbose print.
  void VerbosePrint( const vecString &parts_ids,
		     const vecbvec &parts,
		     ostream &out ) const;
  
  
private:

  const String name_;           // tag name (eg, fw1, rc, consensus, ref)
  mutable ostream *log_;        // optional log gile

  int PRIMER_LEN = 4;           // length of primers
  int FIN_QUAL = 50;            // finished grade quality (qual cap)
  float PERFECT_QUAL = 1.0;     // these define good vs poor aligns. Aligns'
  float GOOD_QUAL = 10.0;       //   scores computed with ScoreAlignment, and
  float BAD_QUAL = 100.0;       //   scores > BAD_QUAL are called "awful"
  
  bvec bases_;                  // bases of sequence
  qvec quals_;                  // quals of sequence
  vec<look_align_plus> chain_;  // chain of aligns (for one target only)
  vec<float> scores_;           // score of aligns (in sync with chain_)
  vec<int> gaps_;               // gaps (has size ids_.size( )+1)
  vec<int> ids_;                // ids of monomers

};

/**
 * Offset
 *
 * Find offset between two CMonoReps (in monomer space).
 *
 * Return false if the two objects do not align. NOTE: if left == fw1
 * and right == fw1, then Offset will narrow the search band around 6.
 * Optionally, return number of matches/mismatches (if pointers are
 * not null).
 */
bool Offset( const CMonoRep &left,
	     const CMonoRep &right,
	     int &offset,
	     ostream *log = 0,
	     int *matches = 0, 
	     int *mismatches = 0 );

/**
 * OffsetBp
 *
 * It converts offset_ms (the offset in monomer space) into the offset
 * in base space.
 */
int OffsetBp( const CMonoRep &left, 
	      const CMonoRep &right, 
	      const int offset_ms );

#endif
