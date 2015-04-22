///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "SeqInterval.h"
#include "SupersHandler.h"
#include "Vec.h"
#include "paths/ReadLoc.h"
#include "paths/assisted/CInsertsDB.h"
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

/**
 * CInsertsDB
 * Constructor
 */
CInsertsDB::CInsertsDB( const int MIN_SEP,
			const int MAX_SEP,
			const int *super_id,
			const shandler *supers,
			read_locs_on_disk *locs_parser ) :
  MIN_SEP_ ( MIN_SEP ),
  MAX_SEP_ ( MAX_SEP ),
  supers_ ( supers ),
  locs_parser_ ( locs_parser )
{
  this->SetSuper( super_id );
}
  
/**
 * CInsertsDB
 * SetPointers
 */
void CInsertsDB::SetPointers( const shandler *supers,
			      read_locs_on_disk *locs_parser )
{
  supers_ = supers;
  locs_parser_ = locs_parser;
}

/**
 * CInsertsDB
 * SetSuper
 *
 * Fill locs and loc types (interval_id's set to loc types):
 *   0: valid pair in super (store fw loc only)
 *   1: mate pair is in different super (both reads "close" to super's edge)
 *   2: mate pair unmapped
 */
void CInsertsDB::SetSuper( const int *super_id )
{
  locs_.clear( );
  sis_.clear( );
  if ( ! super_id ) {
    super_id_ = -1;
    return;
  }

  ForceAssert( supers_ && locs_parser_ );

  // Super id.
  super_id_ = *super_id;
  const superb &sup = (*supers_)[super_id_];
  
  // Locs and loc_intervals.
  #pragma omp parallel for
  for (int ii=0; ii<sup.Ntigs( ); ii++) {
    int tig1 = sup.Tig( ii );
    vec<read_loc> tlocs;
    
    #pragma omp critical
    locs_parser_->LoadContig( tig1, tlocs );
    
    vec<read_loc> select;
    vec<seq_interval> seltypes;
    for (int jj=0; jj<tlocs.isize( ); jj++) {
      seq_interval si;
      if ( this->IsValidFw( tlocs[jj], si )   ||
	   this->IsSeparated( tlocs[jj], si ) ||
	   this->IsLoner( tlocs[jj], si ) ) {
	select.push_back( tlocs[jj] );
	seltypes.push_back( si );
      }
    }
    
    #pragma omp critical
    {
      copy( select.begin( ), select.end( ), back_inserter( locs_ ) );
      copy( seltypes.begin( ), seltypes.end( ), back_inserter( sis_ ) );
    }
  }

  // Sort locs.
  SortSync( sis_, locs_ );

}

/**
 * CInsertsDB
 * LocTypeAsString
 */
String CInsertsDB::LocTypeAsString( int ii ) const
{
  if ( sis_[ii].IntervalId( ) == 0 ) return "valid";
  if ( sis_[ii].IntervalId( ) == 1 ) return "separated";
  if ( sis_[ii].IntervalId( ) == 2 ) return "loner";
  return "na";
}

/**
 * CInsertsDB
 * NLocsTotal
 */
int CInsertsDB::NLocsTotal( ) const
{
  return locs_.isize( );
}

/**
 * CInsertsDB
 * NLocsOfType
 *
 * See SetSuper for loc types.
 */
int CInsertsDB::NLocsOfType( int loc_type ) const
{
  ForceAssert( loc_type > -1 && loc_type < 3 );

  int n_valid = 0;
  for (int ii=0; ii<sis_.isize( ); ii++)
    if ( sis_[ii].IntervalId( ) == loc_type ) n_valid++;
  return n_valid;
}

/**
 * CInsertsDB
 * ContigIds
 *
 * It returns the ids of the contigs containing the end reads of the
 * insert ii, with these conventions:
 *
 *   1. a loner insert's mate contig is -2;
 *
 *   2. if insert's mate belongs to another super, then the mate ids
 *      is -1, if the read on this super is rc, or nan, if the read on
 *      this super is fw.
 */
pair<int,int> CInsertsDB::ContigIds( const int ii ) const
{
  pair<int,int> answer = make_pair( -3, -3 );

  const read_loc &loc = locs_[ii];
  if ( this->IsValidType( ii ) ) {
    answer.first = loc.ContigId( );
    answer.second = loc.PartnerContigId( );
    return answer;
  }

  bool is_sep = this->IsSeparatedType( ii );
  
  if ( loc.Fw( ) ) {
    answer.first = loc.ContigId( );
    answer.second = is_sep ? numeric_limits<int>::max( ) : -2;
    return answer;
  }

  if ( loc.Rc( ) ) {
    answer.first = is_sep ? -1 : -2;
    answer.second = loc.ContigId( );
    return answer;
  }
  
  return answer;   // should never get here
}

/**
 * CInsertsDB
 * ContigPos
 *
 * It runs ContigIds and then it converts contig ids into contig pos
 * (position of contig in super). See ContigIds above for the special
 * cases.
 */
pair<int,int> CInsertsDB::ContigPos( const int ii ) const
{
  pair<int,int> answer = this->ContigIds( ii );
  if ( answer.first > -1 && answer.first < numeric_limits<int>::max( ) )
    answer.first = supers_->PosOnSuper( answer.first );
  if ( answer.second > -1 && answer.second < numeric_limits<int>::max( ) )
    answer.second = supers_->PosOnSuper( answer.second );

  return answer;
}

/**
 * CInsertsDB
 * CposToInserts
 *
 * Run this->ContigPos( ) on all inserts of this super, and create a
 * (sorted) map from pairs of contig_pos to the ids of inserts linking
 * those pairs.
 */
vec< triple<int,int,int> > CInsertsDB::CposToInserts( ) const
{
  vec< triple<int,int,int> > result;
  int n_locs = this->Size( );
  result.reserve( n_locs );

  for (int ii=0; ii<n_locs; ii++) {
    pair<int,int> cpos = this->ContigPos( ii );
    triple<int,int,int> result_ii;
    result_ii.first = cpos.first;
    result_ii.second = cpos.second;
    result_ii.third = ii;
    result.push_back( result_ii );
  }

  sort( result.begin( ), result.end( ) );

  return result;
}

/**
 * CInsertsDB
 * LineInfo
 *
 * On each line:
 *   insert_id
 *   window on super
 *   tig_id -> tig_id (or other_super -> tig_id, or tig_id -> other_super)
 *
 * NB: contig ids, not contig pos in supers (note difference with
 * LineInfoAlt below).
 */
vec<String> CInsertsDB::LineInfo( const int ii ) const
{
  vec<String> answer;
  answer.reserve( 3 );

  // Insert id, interval on super, and type.
  answer.push_back( ToString( ii ) );
  
  answer.push_back( "[" + ToString( sis_[ii].Begin( ) ) +
		    ", " + ToString( sis_[ii].End( ) ) + ")" );
  
  // Chain of contigs.
  const superb &sup = (*supers_)[super_id_];
  const read_loc &loc = locs_[ii];
  String chain;
  if ( this->IsValidType( ii ) ) {
    int p1 = supers_->PosOnSuper( loc.ContigId( ) );
    int p2 = supers_->PosOnSuper( loc.PartnerContigId( ) );
    for (int cpos=p1; cpos<=p2; cpos++) {
      chain += "[c" + ToString( sup.Tig( cpos ) ) + "]";
      if ( cpos < p2 ) {
	String str_gap = ToString( sup.Gap( cpos ) );
	String str_dev = ToString( sup.Dev( cpos ) );
	chain += " [" + str_gap + "+/-" + str_dev + "] ";
      }
    }
  }
  else if ( this->IsSeparatedType( ii ) ) {
    String str_c1 = ToString( loc.ContigId( ) );
    String str_s2 = ToString( supers_->ToSuper( loc.PartnerContigId( ) ) );
    if ( loc.Fw( ) ) chain = "[c" + str_c1 + "] -> s" + str_s2;
    else chain = "s" + str_s2 + " -> [c" + str_c1 + "]";
  }
  else if ( this->IsLonerType( ii ) ) {
    String str_c1 = ToString( loc.ContigId( ) );
    if ( loc.Fw( ) ) chain = "[c" + str_c1 + "] -> unmapped";
    else chain = "unmapped -> [c" + str_c1 + "]";
  }
  answer.push_back( chain );
  
  return answer;
}

/**
 * CInsertsDB
 * LineInfoAlt
 *
 * On each line:
 *   super_id
 *   insert_id
 *   tig_pos -> tig_pos (or other_super -> tig_pos, or tig_pos -> other_super)
 *
 * NB: contig pos in supers, not contig ids (note difference with
 * LineInfo above).
 */
vec<String> CInsertsDB::LineInfoAlt( const int ii ) const
{
  vec<String> answer;
  answer.reserve( 3 );
  
  // super_id and insert_id.
  answer.push_back( "s" + ToString( super_id_ ) );
  answer.push_back( ToString( ii ) );
  
  // contig pos.
  const superb &sup = (*supers_)[super_id_];
  const read_loc &loc = locs_[ii];
  String chain;
  if ( this->IsValidType( ii ) ) {
    String str_p1 = ToString( supers_->PosOnSuper( loc.ContigId( ) ) );
    String str_p2 = ToString( supers_->PosOnSuper( loc.PartnerContigId( ) ) );
    chain = "[" + str_p1 + "]" + " -> [" + str_p2 + "]";
  }
  else if ( this->IsSeparatedType( ii ) ) {
    String str_p1 = ToString( supers_->PosOnSuper( loc.ContigId( ) ) );
    String str_s2 = ToString( supers_->ToSuper( loc.PartnerContigId( ) ) );
    if ( loc.Fw( ) ) chain = "[" + str_p1 + "] -> s" + str_s2;
    else chain = "s" + str_s2 + " -> [" + str_p1 + "]";
  }
  else if ( this->IsLonerType( ii ) ) {
    String str_p1 = ToString( supers_->PosOnSuper( loc.ContigId( ) ) );
    if ( loc.Fw( ) ) chain = "[" + str_p1 + "] -> unmapped";
    else chain = "unmapped -> [" + str_p1 + "]";
  }
  answer.push_back( chain );
  
  return answer;
}

/**
 * CInsertsDB
 * BuildCloud
 *
 * si_ids are the ids of the inserts "overlapping" the given gap,
 * r_ids the ids of the actual reads falling in the gap (based on
 * links).
 *
 * NB: r_ids are signed to capture orientation.
 */
void CInsertsDB::BuildCloud( const int gap_id,
			     vec<int> &si_ids,
			     vec<int> *r_ids ) const
{
  si_ids.clear( );
  
  seq_interval win;
  int cid = (*supers_)[super_id_].Tig( gap_id );
  int clen = (*supers_)[super_id_].Len( gap_id );
  int iibeg = supers_->StartOnSuper( cid ) + clen;
  int iiend = supers_->StartOnSuper( cid + 1 );
  win.Set( 0, super_id_, iibeg, iiend );

  // Fill si_ids.
  vec<seq_interval>::const_iterator it;
  int search_beg = Max( 0, iibeg - MAX_SEP_ );
  seq_interval phony_si( 0, super_id_, search_beg, search_beg );
  it = lower_bound( sis_.begin( ), sis_.end( ), phony_si );
  ForceAssert( it != sis_.end( ) );
  int first_id = distance( sis_.begin( ), it );

  for (int jj=first_id; jj<sis_.isize( ); jj++) {
    if ( sis_[jj].Begin( ) > iiend ) break;
    if ( ! this->IsLonerType( jj ) ) continue;
    if ( ! sis_[jj].HasOverlapWith( win ) ) continue;
    si_ids.push_back( jj );
  }
  
  // rids is NULL, return.
  if ( ! r_ids ) return;
  
  // Fill r_ids (end reads of ii too).
  r_ids->clear( );
  r_ids->reserve ( 2 * si_ids.size( ) );
  
  for (int jj=0; jj<si_ids.isize( ); jj++) {
    const read_loc &jj_loc = locs_[ si_ids[jj] ];
    const seq_interval &jj_si = sis_[ si_ids[jj] ];
    const int jj_type = jj_si.IntervalId( );
    
    int cg1 = jj_loc.ContigId( );
    int b1 = supers_->StartOnSuper( cg1 ) + jj_loc.Start( );
    int e1 = supers_->StartOnSuper( cg1 ) + jj_loc.Stop( );
    bool ov1 = ( iibeg <= b1 && b1 < iiend ) || (iibeg < e1 && e1 <= iiend );
    
    // Add only read (partner) if it (its partner) overlap win.
    if ( jj_type == 0 ) {
      int cg2 = jj_loc.PartnerContigId( );
      int b2 = supers_->StartOnSuper( cg2 ) + jj_loc.PartnerStart( );
      int e2 = supers_->StartOnSuper( cg2 ) + jj_loc.PartnerStop( );
      bool ov2 = (iibeg <= b2 && b2 < iiend) || (iibeg < e2 && e2 <= iiend );
      
      bool rc1 = jj_loc.Rc( );
      bool rc2 = jj_loc.PartnerRc( );
      int id1 = rc1 ? - jj_loc.ReadId( ) - 1 : jj_loc.ReadId( );
      int id2 = rc2 ? - jj_loc.PartnerReadId( ) - 1 : jj_loc.PartnerReadId( );

      if ( ov1 ) r_ids->push_back( id1 );
      if ( ov2 ) r_ids->push_back( id2 );
      continue;
    }

    // Add read only if it overlaps win (partner is in different super).
    if ( jj_type == 1 ) {
      int cg2 = jj_loc.PartnerContigId( );
      int s2 = supers_->ToSuper( cg2 );
      ForceAssert( s2 != super_id_ );

      bool rc1 = jj_loc.Rc( );
      int id1 = rc1 ? - jj_loc.ReadId( ) - 1 : jj_loc.ReadId( );

      if ( ov1 ) r_ids->push_back( id1 );
      continue;
    }

    // Missing partner always added, add read only if it overlaps win.
    if ( jj_type == 2 ) {
      ForceAssert( ! jj_loc.PartnerPlaced( ) );
      
      bool rc1 = jj_loc.Rc( );
      bool rc2 = jj_loc.PartnerRc( );
      int id1 = rc1 ? - jj_loc.ReadId( ) - 1 : jj_loc.ReadId( );
      int id2 = rc2 ? - jj_loc.PartnerReadId( ) - 1 : jj_loc.PartnerReadId( );

      if ( ov1 ) r_ids->push_back( id1 );
      r_ids->push_back( id2 );
      continue;
    }
    
    // No other types allowed.
    ForceAssert( 1 == 0 );
  }
  
}

/**
 * CInsertsDB
 * IsValidFw
 * private
 */
bool CInsertsDB::IsValidFw( const read_loc &loc, seq_interval &si ) const
{
  si.Set( -1, super_id_, -1, -1 );

  if ( ! loc.PartnerPlaced( ) ) return false;
  if ( loc.Rc( ) ) return false;
  if ( loc.PartnerFw( ) ) return false;
  
  int t2 = loc.PartnerContigId( );
  if ( supers_->ToSuper( t2 ) != super_id_ ) return false;
  
  int t1 = loc.ContigId( );
  int end1 = loc.Stop( ) + supers_->StartOnSuper( t1 );
  int beg2 = loc.PartnerStart( ) + supers_->StartOnSuper( t2 );
  int sep = beg2 - end1;
  if ( sep < MIN_SEP_ || sep > MAX_SEP_ ) return false;
  
  si.SetIntervalId( 0 );
  si.SetBegin( end1 );
  si.SetEnd( beg2 );
  return true;
}

/**
 * CInsertsDB
 * IsSeparated
 * private
 *
 * Allow for "separated" links only if the total separation is within
 * the MIN_SEP to MAX_SEP range.
 */
bool CInsertsDB::IsSeparated( const read_loc &loc, seq_interval &si ) const
{
  si.Set( -1, super_id_, -1, -1 );

  if ( ! loc.PartnerPlaced( ) ) return false;
  
  int t1 = loc.ContigId( );
  int t2 = loc.PartnerContigId( );
  int s2 = supers_->ToSuper( t2 );
  if ( s2 == super_id_ ) return false;
  
  int sep = 0;
  int end1 = 0;
  int beg1 = 0;
  int s1len = (*supers_)[super_id_].TrueLength( );
  if ( loc.Fw( ) ) {
    end1 = loc.Stop( ) + supers_->StartOnSuper( t1 );
    sep +=  s1len - end1;
  }
  else {
    beg1 = loc.Start( ) + supers_->StartOnSuper( t1 );
    sep += beg1;
  }
  if ( loc.PartnerFw( ) ) {
    int end2 = loc.PartnerStop( ) + supers_->StartOnSuper( t2 );
    sep += (*supers_)[s2].TrueLength( ) - end2;
  }
  else {
    int beg2 = loc.PartnerStart( ) + supers_->StartOnSuper( t2 );
    sep += beg2;
  }
  if ( sep < MIN_SEP_ || sep > MAX_SEP_ ) return false;
  
  si.SetIntervalId( 1 );
  if ( loc.Fw( ) ) {
    si.SetBegin( end1 );
    si.SetEnd( s1len );
  }
  else {
    si.SetBegin( 0 );
    si.SetEnd( beg1 );
  }
  return true;
}

/**
 * CInsertsDB
 * IsLoner
 * private
 */
bool CInsertsDB::IsLoner( const read_loc &loc, seq_interval &si ) const
{
  si.Set( -1, super_id_, -1, -1 );

  if ( loc.PartnerPlaced( ) ) return false;
  
  int t1 = loc.ContigId( );
  int beg1 = loc.Start( ) + supers_->StartOnSuper( t1 );
  int end1 = loc.Stop( ) + supers_->StartOnSuper( t1 );
  int s1len = (*supers_)[super_id_].TrueLength( );
  
  si.SetIntervalId( 2 );
  if ( loc.Fw( ) ) {
    int win_right = Min( end1 + MAX_SEP_, s1len - end1 );
    si.SetBegin( end1 );
    si.SetEnd( win_right );
  }
  else {
    int win_left = Max( 0, beg1 - MIN_SEP_ );
    si.SetBegin( win_left );
    si.SetEnd( beg1 );
  }
  
  return true;
}

