// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

#include "Alignment.h"
#include "pairwise_aligners/AlignTwoBasevectors.h"
#include "math/Functions.h"
#include "ScoreAlignment.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "STLExtensions.h"
#include "system/System.h"
#include "tiled/AugmentAligns.h"

/*
 * augment_aligns
 * Constructor
 */
augment_aligns::augment_aligns( ) :
  hits_ ( 0 ),
  all_aligns_ ( 0 ),
  bases_ ( 0 ),
  quals_ ( 0 ),
  indiv_ ( 0 ),
  log_ ( 0 ),
  new_aligns_log_ ( 0 ),
  copy_only_ ( false ),
  max_new_ ( 3 ),
  max_discrepancy_ ( 30 ),
  max_error_rate_ ( 0.4 ),
  devnull_( "/dev/null" )
{
  this->Setup( );
}



/*
 * augment_aligns
 * Constructor
 */
augment_aligns::augment_aligns( const vec<look_align_plus> *hits,
				const vec_alignment_plus *all_aligns,
				const vecbasevector *bases,
				const vecqualvector *quals,
				const vec<int> *indiv,
				ostream *log,
				ostream *new_aligns_log ) :
  hits_ ( hits ),
  all_aligns_ ( all_aligns ),
  bases_ ( bases ),
  quals_ ( quals ),
  indiv_ ( indiv ),
  log_ ( log ),
  new_aligns_log_ ( new_aligns_log ),
  copy_only_ ( false ),
  max_new_ ( 3 ),
  max_discrepancy_ ( 30 ),
  max_error_rate_ ( 0.4 ),
  devnull_( "/dev/null" )
{
  this->Setup( );
}



/*
 * augment_aligns
 * SetPointers
 *
 * See the comment in the Constructor re hits being sorted.
 */
void augment_aligns::SetPointers( const vec<look_align_plus> *hits,
				  const vec_alignment_plus *all_aligns,
				  const vecbasevector *bases,
				  const vecqualvector *quals,
				  const vec<int> *indiv,
				  ostream *log,
				  ostream *new_aligns_log )
{
  hits_ = hits;
  all_aligns_ = all_aligns;
  bases_ = bases;
  quals_ = quals;
  indiv_ = indiv;
  log_ = log;
  new_aligns_log_ = new_aligns_log;

  this->Setup( );
}



/*
 * augment_aligns
 * Generate
 *
 * Generate/copy alignments on the full chromosome.
 */
void augment_aligns::Generate( )
{
  int begin = 0;
  int end = (int)hits_->size( );
  
  this->Generate( begin, end );
}



/*
 * augment_aligns
 * Generate
 *
 * First copy in sub_aligns all and only the alignments from all_aligns
 * for which at least one read belongs to the set of reads in [begin_hit,
 * end_hit). Then, search for and add the new alignments, provided
 * copy_only_ is not true.
 */
void augment_aligns::Generate( int begin_hit, int end_hit )
{
  this->CopyAligns( begin_hit, end_hit );
  this->GenerateNewAligns( begin_hit, end_hit );
}



/*
 * augment_aligns
 * SortAndSave
 */
void augment_aligns::SortAndSave( const String &aligns_file )
{
  ostream &log = ( log_ ) ? *log_ : devnull_;

  log << Date( ) << ": sorting " << sub_aligns_.size( ) << " aligns" << endl;
  sort( sub_aligns_.begin( ), sub_aligns_.end( ) );

  log << Date( ) << ": saving onto " << aligns_file << endl;
  WriteAlignsWithIndex( aligns_file, sub_aligns_ );

  log << Date( ) << ": done" << endl;
}



/*
 * augment_aligns
 * Setup
 */
void augment_aligns::Setup( )
{
  if ( hits_ ) {
    order_lookalign_TargetLookAlignOffset sorter;
    ForceAssert( is_sorted( hits_->begin( ), hits_->end( ), sorter ) );
  }
}



/*
 * augment_aligns
 * CopyAligns
 *
 * Copy pertinent alignments from all_aligns to sub_aligns_. Memory
 * for sub_aligns_ is reserved in here.
 */
void augment_aligns::CopyAligns( int begin_hit, int end_hit )
{
  ostream &log = ( log_ ) ? *log_ : devnull_;
  
  // Generate the vector of read ids in the given interval.
  vec<int> read_ids;
  read_ids.reserve( end_hit - begin_hit );
  for (int ii=begin_hit; ii<end_hit; ii++)
    read_ids.push_back( (*hits_)[ii].query_id );
  sort( read_ids.begin( ), read_ids.end( ) );
			
  // Count pertinent aligns.
  int dotter = 100000;
  log << Date( ) << ": count aligns to be copied - one dot corresponds to "
      << dotter << " alignments ("
      << all_aligns_->GetNumberAlignments( ) << " aligns total)"
      << endl;
  
  alignment_plus al;
  int n_killed = 0;
  int n_aligns = (int)all_aligns_->GetNumberAlignments( );
  vec<Bool> kill_align(n_aligns, False );
  for (int ii=0; ii<n_aligns; ii++) {
    if ( ii % dotter == 0 && ii > 0 )
      Dot( log, ii / dotter );
    
    int id1 = all_aligns_->GetAlignmentId1( ii );
    if ( binary_search( read_ids.begin( ), read_ids.end( ), id1 ) )
      continue;
    
    int id2 = all_aligns_->GetAlignmentId2( ii );
    if ( binary_search( read_ids.begin( ), read_ids.end( ), id2 ) )
      continue;
    
    n_killed++;
    kill_align[ii] = True;
  }
  log << "\n";
  
  // Reserve memory (with extra space for the new alignments which will
  //  be generated later on).
  int n_alive = (int)kill_align.size( ) - n_killed;
  int max_added = 2 * ( end_hit - begin_hit ) * max_new_;
  sub_aligns_.reserve( n_alive + max_added );
  
  // Copy alive alignments to sub_aligns_.
  log << Date( ) << ": copying " << n_alive << " aligns" << endl;
  for (int ii=0; ii<n_aligns; ii++) {
    if ( kill_align[ii] == False ) {
      all_aligns_->GetAlignment( al, ii );
      sub_aligns_.push_back( al );
    }
  }

  log << Date( ) << ": copying done" << endl;
}



/*
 * augment_aligns
 * GenerateNewAligns
 *
 * If two reads land close to each other onto the master sequence, and
 * if an alignment between them does not exist already in all_aligns_,
 * then the two reads are tested for alignment. If an align is found, then
 * it will be added to sub_aligns_. Notice that memory for sub_aligns_ is
 * reserved inside CopyAligns( ).
 */
void augment_aligns::GenerateNewAligns( int begin_hit, int end_hit )
{
  if ( copy_only_ )
    return;
  
  ostream &log = ( log_ ) ? *log_ : devnull_;
  ostream &new_log = ( new_aligns_log_ ) ? *new_aligns_log_ : devnull_;

  int dotter = 10000;
  log << Date( ) << ": generate new aligns - one dot corresponds to "
      << dotter << " hits ("
      << end_hit - begin_hit << " hits total)"
      << endl;
  
  int skip_max = 3;
  
  // Loop over the hits. Only check for overlap with hits on the right. 
  for (int hit_id=begin_hit; hit_id<end_hit; hit_id++) {
    if ( hit_id % dotter == 0 && hit_id > 0 )
      Dot( log, hit_id / dotter );
    
    // If read individual is not 0, skip it.
    if ( indiv_ ) {
      int the_id = (*hits_)[hit_id].query_id;
      int the_ind =  (*indiv_)[the_id];
      if ( the_ind != 0 ) {
	new_log << "Skipping hit " << hit_id
		<< " i.e. " << the_id
		<< ": it belongs to individual " << the_ind
		<< "\n";
	continue;
      }
    }

    // Build vector of hit_ids that will be tested for overlap with hit_id.
    vec<int> to_test;
    vec<int> not_to_test;
    int skipped = 0;
    for (int next_hit=hit_id+1; next_hit<end_hit; next_hit++) {
      if ( skipped > skip_max )
	break;
      if ( this->IsInAllAligns( hit_id, next_hit ) ) {
	not_to_test.push_back( next_hit );
	continue;
      }
      if ( this->Overlap( hit_id, next_hit ) > 0 ) {
	to_test.push_back( next_hit );
	continue;
      }
      skipped++;
      if ( not_to_test.size( ) < 1 )
	to_test.push_back( next_hit );
    }

    // Test for align and eventually add to the heap.
    int n_tested = 0;
    for (int ii=0; ii<(int)to_test.size( ); ii++) {
      if ( n_tested >= max_new_ )
	break;
      this->AddAlign( hit_id, to_test[ii] );
      n_tested++;
    }
  }
  log << "\n";

  log << Date( ) << ": "
      << sub_aligns_.size( ) << " aligns at the end"
      << endl;
  
}



/*
 * augment_aligns
 * AddAlign
 *
 * Test for align the two reads at hit1 and hit2, and eventually
 * add the two alignments (read1-read2 and read2-read1) to sub_aligns_.
 * Notice that the aligns are push_back'ed (it is assumed memory is
 * reserved before calling AddAlign). It returns true iff the aligns
 * have been added. Adding aligns will leave sub_aligns_ unsorted.
 */
bool augment_aligns::AddAlign( int hit1, int hit2 )
{
  ostream &log = ( new_aligns_log_ ) ? *new_aligns_log_ : devnull_;

  int read1 = (*hits_)[hit1].query_id;
  int read2 = (*hits_)[hit2].query_id;

  log << "Testing hits [" << hit1
      << ", " << hit2
      << "] (reads [" << read1
      << ", " << read2
      << "]); implied overlap: " << this->Overlap( hit1, hit2 )
      << "; ";
  
  // Try with SmithWatBandedA
  static align al;
  float al_score = -1.0;
  float al_error = 100.0;

  al_score = this->AlignBandedSW( hit1, hit2, al );
  if ( al_score >= 0 )
    al_error = this->ErrorRate( hit1, hit2, al );
  
  // If SmithWatBandedA fails, try with AlignTwoBasevectors.
  if ( al_score < 0 || al_error > max_error_rate_ ) {
    log << "sw failed, trying atb... ";
    al_score = this->AlignATB( hit1, hit2, al );
    if ( al_score >= 0)
      al_error = this->ErrorRate( hit1, hit2, al );
  }
  
  // Align not found.
  if ( al_score < 0 || al_error > max_error_rate_ ) {
    log << "align not found\n";
    return false;
  }
  
  // Success. Add the aligns to sub_aligns_.
  static alignment tempalign;
  static alignment_plus swap_ap;
  int len1 = (*hits_)[hit1].query_length;
  int len2 = (*hits_)[hit2].query_length;
  orientation orient1 = (*hits_)[hit1].rc1 ? ReverseOr : ForwardOr;
  orientation orient2 = (*hits_)[hit2].rc1 ? ReverseOr : ForwardOr;
  Bool RC2 = ( orient1 == orient2 ) ? False : True;
  tempalign.Set( al );

  alignment_plus the_ap( read1, read2, len1, len2, RC2, tempalign, al_score );
  swap_ap.SetToSwapOf( the_ap, len1, len2 );
  
  sub_aligns_.push_back( the_ap );
  sub_aligns_.push_back( swap_ap );
  
  log << "observed overlap: " << the_ap.a.Pos1( ) - the_ap.a.pos1( )
      << " al score: " << ToString( al_score, 3 )
      << " error rate: " << ToString( al_error, 3 )
      << "\n";
  
  return true;
}



/*
 * augment_aligns
 * IsInAllAligns
 *
 * It checks if the alignment between reads of hit1 and hit2 is already
 * in all_aligns_.
 */
bool augment_aligns::IsInAllAligns( int hit1, int hit2 ) const
{
  int read1 = (*hits_)[hit1].query_id;
  int read2 = (*hits_)[hit2].query_id;

  int al_index = all_aligns_->GetAlignsIndex( read1 );
  if ( al_index < 0 )
    return false;

  for (int ii=al_index; ii<(int)all_aligns_->GetNumberAlignments( ); ii++) {
    if ( all_aligns_->GetAlignmentId1( ii ) != read1 )
      break;
    if ( all_aligns_->GetAlignmentId2( ii ) == read2 )
      return true;
  }
  
  return false;
}



/*
 * Overlap
 *
 * Amount of overlap between two reads induced by placing both of them
 * on a master genome. See OffsetOverlap in LookAlign.h for details.
 */
int augment_aligns::Overlap( int hit1, int hit2 ) const
{
  return LookAlignOffsetOverlap( (*hits_)[hit1], (*hits_)[hit2] );
}



/*
 * augment_aligns
 * AlignBandedSW
 *
 * Returns the alignment score or -1 (if an align was not found).
 */
float augment_aligns::AlignBandedSW( int hit1, int hit2, align &al ) const
{
  ForceAssert( hit1 < hit2 );
  
  int overlap = this->Overlap( hit1, hit2 );
  if ( overlap < 1 )
    return -1.0;

  int read1 = (*hits_)[hit1].query_id;
  int read2 = (*hits_)[hit2].query_id;
  int len1 = (*bases_)[read1].size( );
  int len2 = (*bases_)[read2].size( );
  orientation orient1 = (*hits_)[hit1].rc1 ? ReverseOr : ForwardOr;
  orientation orient2 = (*hits_)[hit2].rc1 ? ReverseOr : ForwardOr;
  
  // If overlap is small, then we want to use a smaller bandwidth.
  int default_band = 30;
  int offset = 0;
  if ( !LookAlignOffset( (*hits_)[hit1], (*hits_)[hit2], offset ) )
    return -1.0;
  int band = Min( default_band, overlap / 2 );

  // Align reads.
  int err = 0;
  basevector base1_rc;
  basevector base2_rc;
  bool rc1 = ( orient1 == ReverseOr ) ? true : false;
  bool rc2 = ( orient2 == ReverseOr ) ? true : false;
  if ( rc1 ) {
    base1_rc = (*bases_)[ read1 ];
    base1_rc.ReverseComplement( );
  }
  if ( rc2 ) {
    base2_rc = (*bases_)[ read2 ];
    base2_rc.ReverseComplement( );
  }

  const basevector &bases1 = ( rc1 ) ? base1_rc : (*bases_)[ read1 ];
  const basevector &bases2 = ( rc2 ) ? base2_rc : (*bases_)[ read2 ];
  float sw_score = SmithWatBandedA( bases1, bases2, offset, band, al, err );

  // We may have to reverse alignment.
  if ( orient1 == ReverseOr )
    al.ReverseThis( len1, len2 );
  
  // If align found is a one-base align declare failure.
  int al_len = al.Pos1( ) - al.pos1( );
  if ( al_len < 2 )
    return -1.0;
  
  // Align found, but implied and observed overlaps do not match.
  int observed_overlap = al.Pos1( ) - al.pos1( );
  if ( Abs( observed_overlap - overlap ) > max_discrepancy_ )
    return -1.0;

  // Return alignment score.
  Bool boolRC2 = ( orient1 == orient2 ) ? False : True;
  const basevector &base1 = (*bases_)[read1];
  const basevector &base2 = (*bases_)[read2];
  const qualvector &qual1 = (*quals_)[read1];
  const qualvector &qual2 = (*quals_)[read2];
  return ScoreAlignment( boolRC2, al, base1, qual1, base2, qual2 );

}



/*
 * augment_aligns
 * AlignATB
 *
 * Align the two reads with AlignTwoBasevectors. If an align is found
 * then it returns the score of the alignment, if not it returns -1.
 */
float augment_aligns::AlignATB( int hit1, int hit2, align &al ) const
{
  ForceAssert( hit1 < hit2 );
  
  int read1 = (*hits_)[hit1].query_id;
  int read2 = (*hits_)[hit2].query_id;
  const basevector &base1 = (*bases_)[read1];
  const basevector &base2 = (*bases_)[read2];
  orientation orient1 = (*hits_)[hit1].rc1 ? ReverseOr : ForwardOr;
  orientation orient2 = (*hits_)[hit2].rc1 ? ReverseOr : ForwardOr;
  int RC2 = ( orient1 == orient2 ) ? 0 : 1;
  
  int al_len = AlignTwoBasevectors( base1, base2, al, 0, 1000,
				    1.0, 0, RC2, 0, 8 );
  
  // Alignment not found, return failure.
  if ( al_len < 2 )
    return -1.0;
  
  // We may have to reverse alignment.
  if ( orient1 == ReverseOr ) {
    int len1 = (*bases_)[read1].size( );
    int len2 = (*bases_)[read2].size( );
    al.ReverseThis( len1, len2 );
  }

  // Alignment found, but discrepancy with implied overlap is too large.
  int implied_overlap = this->Overlap( hit1, hit2 );
  int observed_overlap = al.Pos1( ) - al.pos1( );
  if ( Abs( implied_overlap - observed_overlap ) > max_discrepancy_ )
    return -1.0;

  // Alignment found, return score alignment.
  Bool boolRC2 = ( RC2 == 0 ) ? False : True;
  const qualvector &qual1 = (*quals_)[read1];
  const qualvector &qual2 = (*quals_)[read2];
  return ScoreAlignment( boolRC2, al, base1, qual1, base2, qual2 );
  
}



/*
 * augment_aligns
 * ErrorRate
 *
 * Returns alignment's error rate.
 */
float augment_aligns::ErrorRate( int hit1, int hit2, const align &al ) const
{
  int read1 = (*hits_)[hit1].query_id;
  int read2 = (*hits_)[hit2].query_id;
  const basevector &bases1 = (*bases_)[read1];
  const basevector &bases2 = (*bases_)[read2];
  bool is_rc2 = ( (*hits_)[hit1].rc1 == (*hits_)[hit2].rc1 ) ? False : True;
  
  int n_errors = ActualErrors( is_rc2, bases1, bases2, al, 1, 1 );
  int tot_len = al.Pos1( ) - al.pos1( );

  return SafeQuotient( n_errors, tot_len );
}



