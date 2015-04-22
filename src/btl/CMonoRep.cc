///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "ApplyCorrections.h"
#include "LocsHandler.h"
#include "PrintAlignment.h"
#include "ReadLocation.h"
#include "ScoreAlignment.h"
#include "TokenizeString.h"
#include "String.h"
#include "VecUtilities.h"
#include "btl/MonomerUtilsTALENs.h"
#include "btl/CMonoRep.h"
#include "lookup/LookAlign.h"
#include "pairwise_aligners/AlignTwoBasevectors.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "tiled/AlignsOnConsensus.h"

#include "PrettyPrintTable.h"
#include "tiled/MatchesPerf.h"
#include "tiled/PadsRatio.h"
#include "tiled/TilingPrinter.h"

/**
 * CConsensus
 * Constructor
 */
CConsensus::CConsensus( int n_haplotypes,
			ostream *log ) :
  n_haplotypes_ ( n_haplotypes ),
  log_ ( log ),
  analyzer_ ( 0 ),
  chimeras_ ( 0 )
{
  ForceAssert( n_haplotypes_ == 1 || n_haplotypes_ == 2 );
  this->SetDefaults( );
}

/**
 * CConsensus
 * Destructor
 */
CConsensus::~CConsensus( )
{
  if ( log_ ) *log_ << flush;
}

/**
 * CConsensus
 * SetTilingAnalyzer
 */
void CConsensus::SetTilingAnalyzer( tiling_analyzer *analyzer,
				    const vec<Bool> *chimeras )
{
  analyzer_ = analyzer;
  chimeras_ = chimeras;
  this->ClearAll( );
}

/**
 * CConsensus
 * SetVerbose
 */
void CConsensus::SetVerbose( bool verbose )
{
  verbose_ = verbose;
}

/**
 * CConsensus
 * SetDefaults
 */
void CConsensus::SetDefaults( )
{
  win_size_ = 24;
  max_pads_ratio_ = 0.35;
  good_score_ = 0.01;
  verbose_ = false;
}

/**
 * CConsensus
 * GenerateFixes
 */
void CConsensus::GenerateFixes( vec<corrector> &fixes )
{  
  // Log stream.
  ofstream devnull( "/dev/null" );
  ostream &log = log_ ? *log_ : devnull;

  // Logging start of contig.
  const tiling &tile_set = analyzer_->Tiling( );
  log << "CONTIG_" << tile_set.ContigId( ) << "\n";
  
  // Build and analyze good/bad windows.
  this->ClearAll( );
  this->BuildWindows( );

  // Parse windows.
  for (int ii=0; ii<(int)win_.size( ); ii++) {
    this->SelectLeadReads( ii );
    this->LeadReadsToConsensus( );
  }

  // Ok.
  fixes = fixes_;
  
}

/**
 * CConsensus
 * ClearAll
 */
void CConsensus::ClearAll( )
{
  win_.clear( );
  fixes_.clear( );
  this->ClearTempData( );
}

/**
 * CConsensus
 * ClearTempData
 *
 * Clear all temporary data (temporary rectangles, etc.)
 */
void CConsensus::ClearTempData( )
{
  winsel_ = make_pair( -1, -1 );
  rbases_.clear( );
  rquals_.clear( );
  read_pos_.clear( );
  helpers_.clear( );
  helper_pos_.clear( );
  cscores_.clear( );
}

/**
 * CConsensus
 * BuildWindows
 */
void CConsensus::BuildWindows( )
{
  // Clear windows.
  win_.clear( );
  
  // Tiling.
  const tiling &tile_set = analyzer_->Tiling( );
  const padded_seq &cpads = tile_set.Contig( ).PaddedSeq( );

  // Master window (this may go over the ends of the contig).
  pair<int,int> cwin = analyzer_->CWin( );
  int w_begin = cwin.first;
  int w_end = cwin.second;
  for (int ii=0; ii<tile_set.ReadsCount( ); ii++) {
    pair<int,int> rwin = analyzer_->RWin( ii );
    w_begin = Min( w_begin, rwin.first );
    w_end = Max( w_end, rwin.second );
  }
  int w_len = w_end - w_begin;
  master_win_ = make_pair( w_begin, w_end );
  
  // Define windows.
  int loc_beg = w_begin;
  while ( loc_beg < w_end ) {
    int loc_end = loc_beg + this->NextWinSize( loc_beg );
    win_.push_back( make_pair( loc_beg, loc_end ) );
    loc_beg = loc_end;
  }
}

/**
 * CConsensus
 * NextWinSize
 *
 * Decide how large the next window will be, based on the length of
 * the longest stretch of read with good quality score.
 */
int CConsensus::NextWinSize( int start )
{
  // Defalut min win length (in some cases a shorter length may be returned).
  const int min_win_len = win_size_;

  // Initial cleanup and setup.
  this->ClearTempData( );

  int w_begin = start;
  int w_end = Min( w_begin + 96, master_win_.second );
  int w_len = w_end - w_begin;
  
  base_prob bprob;
  read_score scorer( &bprob );

  // Getting close to the end of contig.
  if ( w_len <= min_win_len )
    return w_len;

  // Generate rectangle.
  analyzer_->Rectangles( w_begin, w_len, read_pos_, rbases_, rquals_ );
  ForceAssert( rbases_.size( ) > 0 );
  
  // Remove from the set of eligible reads those that start after pos.
  vec<Bool> is_valid( rbases_.size( ), False );
  for (int ii=1; ii<(int)rbases_.size( ); ii++) {
    if ( ! IsEmpty( rbases_[ii][0] ) )
      is_valid[ii] = True;
  }
  
  // Remove reads that are too short.
  for (int ii=1; ii<(int)rbases_.size( ); ii++) {
    if ( ! is_valid[ii] )
      continue;
    if ( IsEmpty( rbases_[ii][min_win_len] ) )
      is_valid[ii] = False;
  }
  
  // Local termination codes.
  float term_badscore = 66.6;
  int term_invalid = 0;
  int term_qual = 1;
  int term_alive = 2;
  
  // Score the chunks of reads.
  vec<float> score( rbases_.size( ), 0.0 );
  vec<int> terminated( rbases_.size( ), -1 );
  vec<int> good_len( rbases_.size( ), 0 );
  for (int ii=0; ii<(int)terminated.size( ); ii++) {
    if ( is_valid[ii] ) {
      score[ii] = scorer.Score( rquals_[ii], 0, min_win_len );
      if ( score[ii] > good_score_ )
	terminated[ii] = term_qual;
      else {
	good_len[ii] = min_win_len;
	terminated[ii] = term_alive;
      }
    }
    else {
      score[ii] = term_badscore;
      terminated[ii] = term_invalid;
    }
  }

  // Keep walking along reads.
  for (int cursor=min_win_len+1; cursor<(int)rbases_[0].size( ); cursor++) {
    int n_alive = 0;

    for (int ii=1; ii<(int)terminated.size( ); ii++)
      if ( terminated[ii] == term_alive ) {
	n_alive++;

	// Out of boundary of read.
	if ( IsEmpty( rbases_[ii][cursor] ) ) {
	  terminated[ii] = term_qual;
	  continue;
	}
	
	// A gap.
	if ( ! IsConventionalBase( rbases_[ii][cursor] ) )
	  continue;

	// Eval new score.
	Float old_score = score[ii];
	score[ii] = scorer.Append( rquals_[ii][cursor], old_score );
	if ( score[ii] > good_score_ ) {
	  terminated[ii] = term_qual;
	  continue;
	}

	// Still good.
	good_len[ii] = good_len[ii] + 1;
      }

    if ( n_alive < 1 )
      break;
  }
  
  return Max( min_win_len, Max( good_len ) );
}

/**
 * CConsensus
 * SelectLeadReads
 *
 * Cluster reads into bags (each bag contains reads that confirm each
 * other), and select representatives (lead reads) from each bag.
 *
 * max_discrep: maximum discrepancy between number of gap/empty bases
 *              of a read versus that of consensus for the read to
 *              be good.
 */
void CConsensus::SelectLeadReads( int win_id )
{
  // Heuristics.
  float max_discrep = 0.3;
  
  // Clean up and initial setup.
  this->ClearTempData( );
  winsel_ = win_[win_id];

  base_prob bprob;
  read_score scorer( &bprob );
  
  // Window boundary.
  int wbegin = winsel_.first;
  int wend = winsel_.second;
  int wlength = wend - wbegin;
  
  const tiling &tile_set = analyzer_->Tiling( );
  int contig_id = tile_set.ContigId( );
  
  // Rectangle for whole window.
  analyzer_->Rectangles( wbegin, wlength, read_pos_, rbases_, rquals_ );
  ForceAssert( rbases_.size( ) > 0 );
  
  // Compactify Rectangle.
  vec< vec<char> > rbases_comp( rbases_.size( ) );
  vec< vec<int> > rquals_comp( rquals_.size( ) );

  for (int ii=0; ii<(int)rbases_.size( ); ii++) {
    for (int jj=0; jj<(int)rbases_[ii].size( ); jj++) {
      char theBase = rbases_[ii][jj];
      int theQual = rquals_[ii][jj];
      if ( IsGap( theBase ) || IsEmpty( theBase ) )
	continue;
      
      // Force everything to upper case.
      ToUpperCase( theBase );
      rbases_comp[ii].push_back( theBase );
      rquals_comp[ii].push_back( theQual );
    }
  }

  // Generate scores (helpers[0] is reserved for consensus).
  for (int ii=0; ii<(int)rquals_.size( ); ii++)
    helpers_.push_back( scHelp( rquals_, scorer, ii ) );
  
  // Hack scores for some deviant reads so they sort to the bottom.
  Float min_bracket = 1.0 - max_discrep;
  Float max_bracket = 1.0 + max_discrep;
  int cons_bases = Max( 1, (int)rquals_comp[0].size( ) );

  for (int ii=1; ii<(int)helpers_.size( ); ii++) {

    // Read is chimeric.
    int read_id = analyzer_->Tiling( ).Read( read_pos_[ii] ).Id( );
    if ( chimeras_ && (*chimeras_)[read_id] ) {
      helpers_[ii].Score_ += 1.0;
      continue;
    }
    
    // Read compacted differs too much from consensus compacted length.
    int read_bases = rquals_comp[ii].size( );
    Float ratio = SafeQuotient( read_bases, cons_bases );
    if ( ratio < min_bracket || ratio > max_bracket ) {
      helpers_[ii].Score_ += 1.0;
      continue;
    }
      
    // This is the read beginning or ending.
    int pos = helpers_[ii].Pos_;
    pair<int,int> read_locus = analyzer_->RWin( read_pos_[pos] );
    int read_begin = read_locus.first;
    int read_end = read_locus.second;
    if ( wbegin <= read_begin && read_begin < wend ) {
      helpers_[ii].Score_ += 1.0;
      continue;
    }
    if ( wbegin < read_end && read_end <= wend ) {
      helpers_[ii].Score_ += 1.0;
      continue;
    }
  }

  // Sort helpers.
  sort( helpers_.begin( ) + 1, helpers_.end( ) );
  
  // Partition good reads (two passes, place lead reads first, then the rest).
  {
    vec<Bool> skip( helpers_.size( ), True );
    for (int ii=1; ii<(int)helpers_.size( ); ii++) {
      const vec<char> &consensus = rbases_[0];
      if ( ! this->IsGoodStrong( rbases_[ helpers_[ii].Pos_ ], consensus ) ) {
	skip[ii] = False;
	continue;
      }
      int thePos = -1;
      for (int jj=0; jj<(int)helper_pos_.size( ); jj++) {
	int lead_pos = helpers_[ helper_pos_[jj][0] ].Pos_;
	const vec<char> & lead_read = rbases_comp[lead_pos];
	const vec<char> & this_read = rbases_comp[helpers_[ii].Pos_];
	if ( MatchesPerf( lead_read, this_read ) ) {
	  thePos = jj;
	  break;
	}
      }
      if ( thePos > -1 )
	helper_pos_[thePos].push_back( ii );
      else {
	vec<int> newpos( 1, ii );
	helper_pos_.push_back( newpos );
      }
    }

    for (int ii=1; ii<(int)helpers_.size( ); ii++) {
      if ( skip[ii] )
	continue;
      const vec<char> &consensus = rbases_[0];
      if ( ! this->IsGoodWeak( rbases_[ helpers_[ii].Pos_ ], consensus ) )
	continue;
      int thePos = -1;
      for (int jj=0; jj<(int)helper_pos_.size( ); jj++) {
	int lead_pos = helpers_[ helper_pos_[jj][0] ].Pos_;
	const vec<char> & lead_read = rbases_comp[lead_pos];
	const vec<char> & this_read = rbases_comp[helpers_[ii].Pos_];
	if ( MatchesPerf( lead_read, this_read ) ) {
	  thePos = jj;
	  break;
	}
      }
      if ( thePos > -1 )
	helper_pos_[thePos].push_back( ii );
    }
  }
  
  // Coalesce the bags (treat separately bags with one read only).
  cscores_.clear( );
  cscores_.resize( helper_pos_.size( ), -1.0 );
  for (int ii=0; ii<(int)helper_pos_.size( ); ii++) {
    if ( helper_pos_[ii].size( ) == 1 ) {
      cscores_[ii] = helpers_[ helper_pos_[ii][0] ].Score_;
      continue;
    }
    
    int corepos = helpers_[ helper_pos_[ii][0] ].Pos_;
    vec<int> corequals = rquals_[corepos];

    for (int jj=1; jj<(int)helper_pos_[ii].size( ); jj++) {
      int localpos = helpers_[ helper_pos_[ii][jj] ].Pos_;
      const vec<int> &locquals = rquals_[localpos];

      for (int kk=0; kk<(int)corequals.size( ); kk++) {
	if ( corequals[kk] < 0 ) continue;
	if ( locquals[kk] < 0 ) continue;
	corequals[kk] += locquals[kk];
      }
    }

    for (int kk=0; kk<(int)corequals.size( ); kk++)
      corequals[kk] = Min( fin_qual, corequals[kk] );

    cscores_[ii] = scorer.Score( corequals, 0, corequals.size( ) );
  }

  // Reorganize bags so that those with best coalesced score are at the top.
  {
    vec< pair<float,int> > score2pos;
    for (int ii=0; ii<(int)helper_pos_.size( ); ii++)
      score2pos.push_back( make_pair( cscores_[ii], ii ) );
    sort( score2pos.begin( ), score2pos.end( ) );

    vec<float> newcscores;
    vec< vec<int> > newhelpers;
    for (int ii=0; ii<(int)helper_pos_.size( ); ii++) {
      newcscores.push_back( score2pos[ii].first );
      newhelpers.push_back( helper_pos_[ score2pos[ii].second ] );
    }
    swap( newcscores, cscores_ );
    swap( newhelpers, helper_pos_ );
  }

  // Tag helpers (this is only used in the log).
  for (int ii=0; ii<(int)helper_pos_.size( ); ii++) {
    double cscore = cscores_[ii];
    String tag = "lead_" + ToString( ii );
    tag += ( cscore <= good_score_ ) ? "_ok" : "_bad";
    helpers_[ helper_pos_[ii][0] ].Tag_ = tag;
    for (int jj=1; jj<(int)helper_pos_[ii].size( ); jj++)
      helpers_[ helper_pos_[ii][jj] ].Tag_ = "bag_" + ToString( ii );
  }

  // No effect if verbose_ = false.
  this->PrintCurrentWin( );
}

/**
 * CConsensus
 * LeadReadsToConsensus
 */
void CConsensus::LeadReadsToConsensus( )
{
  // Log stream.
  ofstream devnull( "/dev/null" );
  ostream &log = log_ ? *log_ : devnull;

  // Bhaviour for no coverage: leave it as it is (UntrimToCloseGaps friendly).
  if ( cscores_.size( ) < 1 ) {
    if ( verbose_ )
      log << "Bald spot found (UntrimToCloseGaps?). Skip window.\n\n";
    return;
  }

  // How many haplotypes seen.
  int n_hap_seen = 0;
  for (int ii=0; ii<(int)cscores_.size( ); ii++)
    if ( cscores_[ii] <= good_score_ )
      n_hap_seen++;
  
  // How to run ToConsensusCore.
  int nhap = 1;
  bool is_poor = false;

  // Cluster with lead read is poor quality.
  if ( n_hap_seen < 1 ) {
    if ( verbose_ )
      log << "Poor quality lead read. Using strict mode.\n\n";
    nhap = 1;
    is_poor = true;
  }

  // Number of haplotypes seen matches number of haplotypes expected.
  if ( n_hap_seen == 1 ) {
    if ( verbose_ )
      log << "One haplotype seen, ok.\n\n";
    nhap = 1;
    is_poor = false;
  }
  
  if ( n_hap_seen == 2 && n_haplotypes_ == 2 ) {
    if ( verbose_ )
      log << "Two haplotypes seen, ok.\n\n";
    nhap = 2;
    is_poor = false;
  }
  
  // Too many haplotypes seen: probably an ovecollapsed region.
  if ( n_hap_seen > n_haplotypes_ ) {
    if ( verbose_ )
      log << "BAD: " << n_hap_seen
	  << " haplotypes seen. Using strict mode.\n\n";
    nhap = 1;
    is_poor = false;
  }
  
  // Run consensus core.
  for (int pos=0; pos<(int)rbases_[0].size( ); pos++)
    this->ToConsensusCore( pos, nhap, is_poor );
  
}

/**
 * CConsensus
 * ToConsensusCore
 *
 * look_at: either 1 or 2 (look at the first or the first two lead reads)
 * is_poor: look_at must be 1, and it is not a good quality lead read
 */
void CConsensus::ToConsensusCore( int pos, int look_at, bool is_poor )
{
  // Log stream.
  ofstream devnull( "/dev/null" );
  ostream &log = log_ ? *log_ : devnull;
  
  // Validation.
  ForceAssert( look_at == 1 || look_at == 2 );
  if ( is_poor ) ForceAssert( look_at == 1 );

  // Lead read(s).
  int lead0 = helpers_[ helper_pos_[0][0] ].Pos_;
  int lead1 = look_at == 1 ? lead0 : helpers_[ helper_pos_[1][0] ].Pos_;
  
  // Find replacement for old base/qual.
  int padded_pos = winsel_.first + pos;
  char oldB = rbases_[0][pos];
  int oldQ = rquals_[0][pos];
  
  char base0 = rbases_[lead0][pos];
  char base1 = rbases_[lead1][pos];
  ToUpperCase( base0 );
  ToUpperCase( base1 );
    
  // Bald spot on consensus (do nothing: UntrimToCloseGaps friendly).
  if ( base0 == empty_base ) {
    if ( verbose_ )
      log << " empty lead (UntrimToCloseGaps?). Skip.\n";
    return;
  }

  // Determine consensus qual.
  vec<int> quals0;
  vec<int> quals1;
  for (int jj=1; jj<(int)rbases_.size( ); jj++) {
    if ( IsEmpty( rbases_[jj][pos] ) ) continue;
    char theBase = empty_base;
    int theQual = empty_qual;
    this->ReadBaseQual( read_pos_[jj], padded_pos, theBase, theQual );
    ToUpperCase( theBase );
    
    // This should never happen.
    if ( theQual < 0 )
      continue;

    if ( theBase == base0 ) {
      quals0.push_back( theQual );
      continue;
    }
    if ( theBase == base1 )
      quals1.push_back( theQual );
  }

  // The corrector.
  int c_id = analyzer_->Tiling( ).ContigId( );
  corrector corr( c_id, padded_pos, oldB, empty_base, oldQ, empty_qual );
  
  // Use a stringent criterion for quality score if is_poor = true.
  if ( is_poor ) {
    int qual0 = rquals_[lead0][pos];
    int count = quals0.size( );
    ForceAssert( count > 0 );
    float multiplier = ( 1.0 + ( 0.1 * ( count - 1 ) ) );
    int adjusted_qual = int( float( qual0 ) * multiplier );
    int newQ = ( adjusted_qual > fin_qual ) ? fin_qual : adjusted_qual;
    corr.new_base0_ = base0;
    corr.new_qual0_ = newQ;
  }
  else {
    corr.new_base0_ = base0;
    corr.new_qual0_ = Min( fin_qual, Sum( quals0 ) );
    if ( base1 != base0 ) {
      corr.new_base1_ = base1;
      corr.new_qual1_ = Min( fin_qual, Sum( quals1 ) );
    }
  }

  // No action is taken if corrector new base/qual match old base/qual.
  this->AddFix( corr );
  
}

/**
 * CConsensus
 * ReadBaseQual
 *
 * Set pad qual to the minimum of the flanking base qual. A lower case
 * base (or pad) can be returned. Return false on error.
 */
bool CConsensus::ReadBaseQual( int rPos, int bPos, char &base, int &qual )
{
  // Original base/qual.
  char orig_base = empty_base;
  int orig_qual = empty_qual;
  analyzer_->PrintRead( rPos, bPos, orig_base, orig_qual );
  
  ForceAssert( orig_base != empty_base );
  if ( IsGap( orig_base ) )
    ForceAssert( orig_qual == gap_qual );

  // Base is not a pad.
  if ( ! IsGap( orig_base ) ) {
    base = orig_base;
    qual = orig_qual;
    return true;
  }

  // Base is a pad: look for left and right flanking bases.
  int left_pos = bPos;
  char left_base = orig_base;
  int left_qual = orig_qual;
  while ( left_base != empty_base ) {
    left_pos--;
    analyzer_->PrintRead( rPos, left_pos, left_base, left_qual );
    if ( ! IsGap( left_base ) )
      break;
  }
  
  int right_pos = bPos;
  char right_base = orig_base;
  int right_qual = orig_qual;
  while ( right_base != empty_base ) {
    right_pos++;
    analyzer_->PrintRead( rPos, right_pos, right_base, right_qual );
    if ( ! IsGap( right_base ) )
      break;
  }
  
  // Set base.
  base = orig_base;

  // This should never happen.
  if ( left_base == empty_qual || right_base == empty_qual ) {
    qual = 0;
    return false;
  }
  int min_qual = Min( left_qual, right_qual );
  if ( min_qual < 0 ) {
    qual = 0;
    return false;
  }

  // Set qual and return.
  qual = min_qual;
  return true;
  
}

/**
 * CConsensus
 * AddFix
 *
 * If the new base (or qual) is different from the old one, add the
 * proper fixing event to fixes_. Notice corr will be changed.
 */
void CConsensus::AddFix( corrector &corr )
{
  // Log stream.
  ofstream devnull( "/dev/null" );
  ostream &log = log_ ? *log_ : devnull;

  // Translate to unpadded pos,  and reset gap quals.
  const padded_seq &cpads = analyzer_->Tiling( ).Contig( ).PaddedSeq( );
  int padpos = corr.pos_;
  int unpadpos = cpads.ToUnpaddedPos( padpos );
  
  corr.pos_ = unpadpos;
  if ( IsGap( corr.new_base0_ ) ) corr.new_qual0_ = gap_qual;
  if ( IsGap( corr.new_base1_ ) ) corr.new_qual1_ = gap_qual;

  // Old/new bases and quals.
  char oldB = corr.old_base_;
  int oldQ = corr.old_qual_;
  
  char newB = corr.new_base0_;
  int newQ = corr.new_qual0_;
  
  char newB_2ndhap = corr.new_base1_;
  int newQ_2ndhap = corr.new_qual1_;

  // No change.
  if ( IsEmpty( newB_2ndhap ) && ( newB == oldB && newQ == oldQ ) )
    return;
  
  // Add fix.
  fixes_.push_back( corr );
  
  // Log.
  bool true_snp = ( ! IsEmpty( newB_2ndhap ) );
  String str_snp = true_snp ? " TRUE_SNP" : "";
  log << unpadpos << " (padded_" << padpos << "): "
      << " " << corr.OldBase( )
      << "  ->  " << corr.NewBase( ) << str_snp
      << "\n";
  
}

/**
 * CConsensus
 * PrintCurrentWin
 *
 * Print the given window.
 */
void CConsensus::PrintCurrentWin( ) const
{
  if ( !log_ )
    return;
  if ( !verbose_ )
    return;
  ostream &log = *log_;
  
  log << "\n==================================================\n\n";
  tiling_printer printer ( analyzer_ );
  printer.SetPrintWings( false );
  printer.SetWindowToPrint( winsel_.first, winsel_.second );
  printer.ToStream( log );
  
  vec< vec<String> > table;

  vec<String> line;
  line.push_back( "r_id" );
  line.push_back( "r_score" );
  line.push_back( "bag_score" );
  line.push_back( "n_bases" );
  line.push_back( "r_pos" );
  line.push_back( "tag" );
  table.push_back( line );

  for (int ii=0; ii<(int)helper_pos_.size( ); ii++) {
    const vec<int> &hPos = helper_pos_[ii];
    const double score = cscores_[ii];
    for (int jj=0; jj<(int)hPos.size( ); jj++) {
      int pos = helpers_[ hPos[jj] ].Pos_;
      int read_id = analyzer_->Tiling( ).Read( read_pos_[pos] ).Id( );

      line.clear( );
      line.push_back( ToString( read_id ) );
      line.push_back( ToString( helpers_[ hPos[jj] ].Score_, 6 ) );
      line.push_back( ToString( score, 6 ) );
      line.push_back( ToString( helpers_[ hPos[jj] ].nValid_ ) );
      line.push_back( ToString( helpers_[ hPos[jj] ].Pos_ ) );
      line.push_back( helpers_[ hPos[jj] ].Tag_ );
      table.push_back( line );
    }
  }
  
  BeautifyTable( table );

  for (int ii=0; ii<(int)table.size( ); ii++) {
    for (int jj=0; jj<(int)table[ii].size( ); jj++)
      log << "  " << table[ii][jj];
    log << "\n";
  }
  log << "\n";
  
}

/**
 * CConsensus
 * IsGoodWeak
 *
 * Decide if read is good (does not test for score).
 */
bool CConsensus::IsGoodWeak( const vec<char> &read,
			     const vec<char> &consensus ) const
{
  float pratio = PadsRatio( consensus, read );
  return ( pratio <= max_pads_ratio_ );
}

/**
 * CConsensus
 * IsGoodStrong
 *
 * Decide if read can be trusted as a lead read for the consensus (this 
 * is strictly stronger than IsGoodWeak).
 */
bool CConsensus::IsGoodStrong( const vec<char> &read,
			       const vec<char> &consensus ) const
{
  // Lead read must cover the whole region.
  if ( IsEmpty( read[0] ) || IsEmpty( read.back( ) ) )
    return false;
  
  // Is read good?
  return this->IsGoodWeak( read, consensus );
}


/**
 * class CMonoRep
 * Constructor
 */
CMonoRep::CMonoRep( const String &name, ostream *log ) :
  name_ ( name ),
  log_ ( log )
{
  this->Clear( );
}
  
/**
 * class CMonoRep
 * Clear
 */
void CMonoRep::Clear( )
{
  bases_.clear( );
  quals_.clear( );
  chain_.clear( );
  scores_.clear( );
  gaps_.clear( );
  ids_.clear( );
}

/**
 * class CMonoRep
 * SetFromReference
 */
void CMonoRep::SetFromReference( const String &ref,
				 const vecbvec &primers,
				 const vecbvec &parts,
				 const vecString &parts_ids )
{
  this->Clear( );
  
  vec<String> sids;
  sids.reserve( parts_ids.size( ) );
  for (size_t ii=0; ii<parts_ids.size( ); ii++)
    sids.push_back( ShortName( parts_ids[ii] ) );
  
  gaps_.push_back( 0 );

  vec<String> str_monos;
  Tokenize( ref, str_monos );
  for (size_t ii=0; ii<str_monos.size( ); ii++) {
    String strmon = str_monos[ii];
    vec<String>::iterator it = find( sids.begin( ), sids.end( ), strmon );
    ForceAssert( it != sids.end( ) );
    int pid = distance( sids.begin( ), it );

    int p1 = 0;
    int p2 = PRIMER_LEN + bases_.size( );
    avector<int> gaps( 1, 0 );
    avector<int> lens( 1, (int)parts[pid].size( ) );
    align al( p1, p2, gaps, lens );
    chain_.push_back( look_align( pid, 0, lens( 0 ), 0, false, al, 0, 0, 0 ) );
    scores_.push_back( 1.0 );
    gaps_.push_back( 0 );
    ids_.push_back( pid );

    ForceAssertEq( (int)primers[ii%3].size( ), PRIMER_LEN );
    bases_ = Cat( bases_, primers[ii%3] );
    bases_ = Cat( bases_, parts[pid] );
  }
  
  // Adjust target_length in look_aligns.
  for (int ii=0; ii<chain_.isize( ); ii++)
    chain_[ii].target_length = uint( bases_.size( ) );

  // Finished grade quals.
  quals_.resize( bases_.size( ), FIN_QUAL );
  
}

/**
 * class CMonoRep
 * SetFromHits
 *
 * Store rc of bases (and quals), if rc is true, so all aligns of
 * monomers on sequence will be fw.
 */
void CMonoRep::SetFromHits( const vec<look_align_plus> &hits,
			    const vecbvec &parts,
			    const bvec &bases,
			    const qvec &quals,
			    const int &tid,
			    const bool rc )
{
  // Clear all.
  this->Clear( );
  
  // Bases and quals first.
  bases_ = bases;
  quals_ = quals;
  if ( rc ) {
    bases_.ReverseComplement( );
    quals_.ReverseMe( );
  }
  
  // Pick aligns for selected target, and fill chain_.
  for (int ii=0; ii<hits.isize( ); ii++) {
    if ( hits[ii].target_id != tid ) continue;
    if ( rc != hits[ii].rc1 ) continue;

    look_align_plus hit = hits[ii];
    if ( rc ) {
      hit.a.ReverseThis( hit.query_length, hit.target_length );
      hit.rc1 = false;
    }    
    chain_.push_back( hit );
  }
  order_lookalign_TargetBegin sorter;
  sort( chain_.begin( ), chain_.end( ), sorter );

  // No aligns found: empty all and leave.
  if ( chain_.size( ) < 1 ) {
    this->Clear( );
    if ( log_ ) *log_ << "\nSetFromHits failed on " << name_ << "\n";
    return;
  }

  // Fill scores_, gaps_, and ids_.
  gaps_.push_back( Max( 0, chain_[0].pos2( ) - PRIMER_LEN ) );
  for (int ii=0; ii<chain_.isize( ); ii++) {
    const look_align_plus &hit = chain_[ii];

    // ScoreAlignments returns 1.03077... for perfect aligns.
    const bvec &bpart = parts[hit.query_id];
    qvec qpart( bpart.size( ), FIN_QUAL );
    int err = hit.a.Errors( bpart, bases_ );
    float score = ScoreAlignment( hit.a, bpart, qpart, bases_, quals_ );
    scores_.push_back( err < 1 ? 1.0 : score );

    ids_.push_back( hit.query_id );

    int gap = ( ii == chain_.isize( ) - 1 )
      ? chain_.back( ).target_length - chain_.back( ).Pos2 ( )
      : chain_[ii+1].pos2( ) - chain_[ii].Pos2( ) - PRIMER_LEN;
    gaps_.push_back( gap );
  }
  
}

/**
 * class CMonoRep
 * AddInitialGap
 */
void CMonoRep::AddInitialGap( int gap0 )
{
  this->Clear( );
  gaps_.push_back( gap0 );
}

/**
 * class CMonoRep
 * AddSegment
 *
 * We add a segment containing both primer and monomer, so we need to
 * compute where the new look_align actually starts.
 */
void CMonoRep::AddSegment( const bvec &bases,
			   const qvec &quals,
			   const look_align_plus &hit,
			   const float &score,
			   const int &gap,
			   const int &id )
{
  int len_on_target = hit.a.Pos2( ) - hit.a.pos2( );
  int new_start = bases_.size( ) + bases.size( ) - len_on_target;

  chain_.push_back( hit );
  chain_[ chain_.size( )-1 ].target_id = 0;
  chain_[ chain_.size( )-1 ].SetStartOnTarget( new_start );

  bases_.append( bases.begin( ), bases.end( ) );
  quals_.append( quals.begin( ), quals.end( ) );
  scores_.push_back( score );
  gaps_.push_back( gap );
  ids_.push_back( id );

  for (size_t ii=0; ii<chain_.size( ); ii++)
    chain_[ii].target_length = bases_.size( );
  
}

/**
 * class CMonoRep
 * Merge
 *
 * Set this to the object generated by merging left and right.
 */
bool CMonoRep::Merge( const CMonoRep &left, const CMonoRep &right )
{
  this->Clear( );
  
  // Core sizes, etc.
  const int nleft = left.NMonomers( );
  const int nright = right.NMonomers( );
  
  // Find offset.
  int offset = -666;
  if ( ! Offset( left, right, offset, log_ ) ) {
    if ( log_ ) *log_ << "\nMerge failed (offset not found)\n";
    return false;
  }
  
  // Prepare to build chain walking on the alignment (in monomer space).
  int p1 = ( offset >= 0 ) ? 0 : offset;
  int p2 = ( offset >= 0 ) ? -offset: 0;
  ForceAssert( p1 >= 0 || p2 >= 0 );
  
  const bvec &bleft = left.Bases( );
  const qvec &qleft = left.Quals( );
  const bvec &bright = right.Bases( );
  const qvec &qright = right.Quals( );
  
  // Initial gap.
  this->AddInitialGap( 0 );

  // Keep adding blocks.
  while ( p1 < nleft || p2 < nright ) {
    int sel = 0;
    if ( p2 < 0 || p2 >= nright ) sel = 1;
    else if ( p1 < 0 || p1 >= nleft ) sel = 2;
    else sel = ( left.Score( p1 ) < right.Score( p2 ) ? 1 : 2 );
    
    if ( sel == 1 ) {
      const look_align_plus &al = left.Chain( p1 );
      const int start = Max( 0, al.a.pos2( ) - PRIMER_LEN );
      const int len = al.a.Pos2( ) - start;
      
      bvec bb = bvec( bleft, start , len );
      qvec qq( len, 0 );
      CopyQuals( qleft, start, qq, 0, len );
      const float &score = left.Score( p1 );
      const int &id = left.Id( p1 );

      int gap = -666;
      if ( p1 >= nleft - 2 && p2 >= nright - 2 ) gap = 0;
      else if ( p1 >= nleft - 2 ) gap = right.Gap( p2 + 1 );
      else gap = left.Gap( p1 + 1 );
      
      this->AddSegment( bb, qq, al, score, gap, id ); 
    }
    else {
      const look_align_plus &al = right.Chain( p2 );
      const int start = Max( 0, al.a.pos2( ) - PRIMER_LEN );
      const int len = al.a.Pos2( ) - start;
      
      bvec bb = bvec( bright, start, len );
      qvec qq( len, 0 );
      CopyQuals( qright, start, qq, 0, len );
      const float &score = right.Score( p2 );
      const int &id = right.Id( p2 );

      int gap = -666;
      if ( p1 >= nleft - 2 && p2 >= nright - 2 ) gap = 0;
      else if ( p2 >= nright - 1 ) gap = left.Gap( p1 + 1 );
      else gap = right.Gap( p2 + 1 );
      
      this->AddSegment( bb, qq, al, score, gap, id ); 
    }
    
    p1++;
    p2++;
  }

  // Ok.
  return true;
  
}

/**
 * class CMonoRep
 * Refine
 *
 * Refine consensus. Note that fw2 may be empty. It returns false if
 * not all reads could be realigned to consensus (non a fatal error).
 *
 * NOTE: after you run Refine, you will need to realign the monomers
 * to bases_ and quals_, to update chain_, etc. (ie: run SetFromHits).
 */
bool CMonoRep::Refine( const CMonoRep &fw1,
		       const CMonoRep &fw2,
		       const CMonoRep &rc,
		       const String &out_head )
{
  // HEURISTICS.
  float MAX_ER = 0.35;
  int MIN_AL_LEN = 200;
  int SW_BAND = 48;

  // File names.
  String tiles_file = out_head + ".tiles";
  String corrections_file = out_head + ".fixes";
  String corrections_brief_file = out_head + ".fixes.brief";
  
  // Prepare vecbvec, vecqvec, and read_locs.
  vecbvec rbases;
  vecqvec rquals;
  vec<read_location> rawlocs;

  const bvec &contig = bases_;
  const int lenc = contig.size( );
  
  rbases.push_back( fw1.Bases( ) );
  rquals.push_back( fw1.Quals( ) );
  int len0 = fw1.Bases( ).size( );
  
  int moffset = -1;
  if ( Offset( *this, fw1, moffset, log_ ) ) {
    int start = OffsetBp( *this, fw1, moffset );
    rawlocs.push_back( read_location( 0, len0, 0, start, ForwardOr, lenc ) );
  }

  rbases.push_back( rc.Bases( ) );
  rquals.push_back( rc.Quals( ) );
  int len1 = rc.Bases( ).size( );
  if ( Offset( *this, rc, moffset, log_ ) ) {
    int start = OffsetBp( *this, rc, moffset );
    rawlocs.push_back( read_location( 1, len1, 0, start, ForwardOr, lenc ) );
  }

  if ( ! fw2.IsEmpty( ) ) {
    rbases.push_back( fw2.Bases( ) );
    rquals.push_back( fw2.Quals( ) );
    int len2 = fw2.Bases( ).size( );
    if ( Offset( *this, fw2, moffset, log_ ) ) {
      int start = OffsetBp( *this, fw2, moffset );
      rawlocs.push_back( read_location( 2, len2, 0, start, ForwardOr, lenc ) );
    }
  }
  
  // Turn read_locs into an lhandler (needed by tiling code).
  lhandler locs( rbases.size( ), 1 );
  locs.SetFromVec( rawlocs );

  // Realign reads to consensus.
  vec<t_align> tals;
  for (int loc_id=0; loc_id<locs.Size( ); loc_id++) {
    const read_location &loc = locs[loc_id];
    const int rid = loc.ReadId( );
    const int offset = loc.StartOnContig( );
    const bvec &read = rbases[rid];
    
    align al;
    int err = 0;
    SmithWatBandedA( read, contig, -offset, SW_BAND, al, err );
    
    float f_err = al.Errors( read, contig );
    float f_len = Max( 1.0, double( al.Pos2( ) - al.pos2( ) ) );
    if ( f_err / f_len > MAX_ER ) continue;
    if ( (int)f_len < MIN_AL_LEN ) continue;

    tals.push_back( t_align( rid, loc.Rc( ), al ) );
  }

  tiling tiles( -1, 0 );
  tiles.SetFromAligns( tals );

  // Find fixes and applu corrections.
  vecbvec cbases( 1, bases_ );
  vecqvec cquals( 1, quals_ );

  tiling_analyzer analyzer;
  analyzer.SetPointers( &rbases, &rquals, &cbases, &cquals );
  analyzer.SetTiling( &tiles );
  tiling_printer the_printer( &analyzer );
  ofstream out( tiles_file.c_str( ) );
  the_printer.ToStream( out );
  out.close( );
  
  vec<corrector> fixes;
  CConsensus fixer( 1 );
  fixer.SetTilingAnalyzer( &analyzer );
  fixer.GenerateFixes( fixes );
  WRITE( corrections_file, fixes );
  ofstream brief_out( corrections_brief_file.c_str( ) );
  for (size_t ii=0; ii<fixes.size( ); ii++) {
    if ( fixes[ii].CharOldBase( ) == '.' ) continue;
    if ( fixes[ii].CharOldBase( ) == fixes[ii].CharNewBase( ) ) continue;
    brief_out << fixes[ii];
  }
  brief_out.close( );

  vec<read_location> out_locs = locs.Locs( );
  vec<int> out_flocs = locs.FirstLocs( );
  ApplyCorrections( fixes, out_flocs, out_locs, cbases, cquals );
  
  // Replace bases_ and quals_, and return.
  bases_ = cbases[0];
  quals_ = cquals[0];
  
  if ( log_ && tals.size( ) != rbases.size( ) )
    *log_ << "\nProblem in Refine: " << rbases.size( ) - tals.size( )
	  << " failed to realign\n";
  
  return tals.size( ) == rbases.size( );

}

/**
 * class CMonoRep
 * SaveBasesAndQuals
 */
void CMonoRep::SaveBasesAndQuals( const String &head ) const
{
  const String bases_file = head + ".fastb";
  const String quals_file = head + ".qualb";
  
  vecbvec allbases( 1, bases_ );
  allbases.WriteAll( bases_file );
  
  vecqvec allquals( 1, quals_ );
  allquals.WriteAll( quals_file );
}
  
/**
 * class CMonoRep
 * BlockLength
 */
int CMonoRep::BlockLength( int ii ) const
{
  return PRIMER_LEN + chain_[ii].a.Pos2( ) - chain_[ii].a.pos2( );
}
  
/**
 * class CMonoRep
 * EvalOnRef
 */
void CMonoRep::EvalOnRef( const CMonoRep &ref,
			  const vecString &parts_ids,
			  const String &name,
			  const String &csv_file,
			  ostream &visual ) const
{
  ofstream out( csv_file.c_str( ) );

  // HEURISTICS
  const int SW_BAND = 150;

  // Consts.
  const int ref_size = ref.NMonomers( );
  const int cons_size = this->NMonomers( );

  // Reference line.
  out << ",,,," << name << ",ref,,,,,,";
  for (int ii=0; ii<ref_size; ii++)
    out << ShortName( parts_ids[ref.Id( ii )] ) << ",";
  out << "\n";

  // Find offset (early exit if not found).
  int offset = -1;
  if ( ! Offset( ref, *this, offset, log_ ) ) {
    if ( log_ ) *log_ << "\nFatal error in EvalOnRef: offset failed\n";
    out << name << "0,,,,,," << ref.NMonomers( ) << ",," << "no consensus\n";
    return;
  }
 
  // Offset in base space.
  int offset_bp;
  Offset( ref, *this, offset_bp, log_ );
  
  // Align consensus to reference.
  const bvec &b1 = ref.Bases( );
  const qvec q1( b1.size( ), FIN_QUAL );
  
  align al;
  int err = 0;
  SmithWatBandedA2<uint>( b1, bases_, offset_bp, SW_BAND, al, err, 0, 1, 1, 1 );
  double al_score = ScoreAlignment( al, b1, q1, bases_, quals_ );

  // Print align.
  visual << "ALIGN OF CONSENSUS ON REFERENCE\t"
	 << "score: " << ToString( al_score, 1 ) << "\n"
	 << " reference on top\n"
	 << " on ref:       [" << al.pos1( ) << ", " << al.Pos1( )
	 << ")_" << b1.size( ) << "\n"
	 << " on consensus: [" << al.pos2( ) << ", " << al.Pos2( )
	 << ")_" << bases_.size( ) << "\n";

  PrintVisualAlignment( False, visual, b1, bases_, al, qvec( 0 ) , quals_ );
  visual << "\n";

  // High qual / low qual matches and mismatches.
  int hq_match = 0;
  int lq_match = 0;
  int lq_mis = 0;
  int hq_mis = 0;
  
  // Loop over align.
  int p1 = 0;
  int p2 = 0;
  if ( offset < 0 ) p2 = -offset;
  if ( offset > 0 ) p1 = offset;
  String mono_seq;
  while ( p1 < ref_size && p2 < cons_size ) {
    bool bad_score = this->Score( p2 ) > GOOD_QUAL;
    bool bad_gap = ( p2 < cons_size - 1 && this->Gap( p2+1 ) != 0 );
    bool matches_ref = ( ref.Id( p1 ) == this->Id( p2 ) );
    
    // Try to skip spurious random aligns after the end of the product.
    int matches_sofar = hq_match + lq_match + ( matches_ref ? 1 : 0 );
    bool skip_rest = ( matches_sofar == ref_size );
    
    // Update counters.
    String part_name = ShortName( parts_ids[this->Id( p2 )] );
    if ( matches_ref ) {
      if ( bad_score ) {
	mono_seq += ToLower( part_name );
	lq_match++;
      }
      else {
	mono_seq += part_name;
	hq_match++;
      }
    }
    else {
      if ( bad_score ) {
	mono_seq += ToLower( part_name );
	lq_mis++;
      }
      else {
	mono_seq += part_name;
	hq_mis++;
      }
    }
    
    // Report bad gap.
    if ( ( ! skip_rest ) && bad_gap) {
      mono_seq += ".[" + ToString( this->Gap( p2+1 ) ) + "]";
      hq_mis++;
    }

    // Separator.
    mono_seq += ",";
    
    // Next.
    p1++;
    p2++;
  }


  String str_scorer = "";
  if ( al_score <= PERFECT_QUAL ) str_scorer = "perfect";
  else if ( al_score <= GOOD_QUAL ) str_scorer = "good";
  else if ( al_score <= BAD_QUAL ) str_scorer = "poor";
  else str_scorer = "awful";

  out << hq_match + lq_match - lq_mis - hq_mis << ",,"
      << str_scorer << ",,"
      << name << ",,"
      << hq_match << "," << lq_match << "," << lq_mis << "," << hq_mis << ",,"
      << mono_seq << "\n";
  
  out.close( );

}

/**
 * class CMonoRep
 * Print
 */
void CMonoRep::Print( const vecString &parts_ids, ostream &out ) const
{
  out << "\nRepresentation of " << name_ << " in monomer space:\n\n";
  
  if ( this->IsEmpty( ) ) {
    out << "no monomers found\n";
    return;
  }
  
  out << "(" << this->Gap( 0 ) << ")  ";
  for (int ii=0; ii<(int)this->NMonomers( ); ii++) {
    String name = ShortName( parts_ids[ ids_[ii] ] );
    if ( scores_[ii] > GOOD_QUAL ) name.ToLower( );
    
    // Monomer short name (lower case if not good qual).
    out << name << "  ";
    
    // Print gap only if last, or >0.
    if ( ii == (int)this->NMonomers( )-1 )
      out << "(" << this->Gap( ii+1 ) << ")";
    else
      if ( this->Gap( ii+1 ) > 0 )
	out << "[" << this->Gap( ii+1 ) << "]  ";
  }

  out << "\n";
}

/**
 * class CMonoRep
 * VerbosePrint
 */
void CMonoRep::VerbosePrint( const vecString &parts_ids,
			     const vecbvec &parts,
			     ostream &out ) const
{
  out << "ALIGNS OF MONOMERS ON " << name_ << "\n\n";

  for (int ii=0; ii<chain_.isize( ); ii++) {
    const look_align_plus &hit = chain_[ii];
    const align &al = hit.a;

    out << "Align between [" << hit.pos2( )
	<< ", " << hit.Pos2( )
	<< ") of " << name_
	<< " and " << ShortName( parts_ids[ ids_[ii] ] )
	<< "\tal_score: " << ToString( scores_[ii], 1 )
	<< "\n";

    bvec rc_part;
    if ( hit.Rc1( ) ) {
      rc_part = parts[ ids_[ii] ];
      rc_part.ReverseComplement( );
    }
    const bvec &b1 = hit.Rc1( ) ? rc_part : parts[ ids_[ii] ];

    PrintVisualAlignment( False, out, b1, bases_, al, qvec( 0 ) , quals_ );
    out << "\n";
  }
  
}

/**
 * Offset
 */
bool Offset( const CMonoRep &left,
	     const CMonoRep &right,
	     int &offset,
	     ostream *log,
	     int *matches,
	     int *mismatches )
{
  offset = -666;

  // HEURISTICS.
  const int MIN_LEN = 2;
  
  // Const stuff.
  const int nleft = left.NMonomers( );
  const int nright = right.NMonomers( );
  
  // Special cases.
  bool special_fw1_fw2 = ( left.Name( ) == "fw1" && right.Name( ) == "fw2" );
  bool special_fw1_rc = ( left.Name( ) == "fw1" && right.Name( ) == "rc" );

  // Search space for offsets. HEURISTICS embedded here!
  int range_first = - nright + MIN_LEN;
  int range_last = nleft - MIN_LEN;
  if ( special_fw1_fw2 ) {
    range_first = 5;
    range_last = 7;
  }
  if ( special_fw1_rc ) {
    range_first = 0;
  }

  // Compute score of monomer-alignment.
  vec<SAligner> helper;
  for (int offset=range_first; offset<=range_last; offset++) {
    int p1 = offset;
    int p2 = 0;
    int al_len = 0;
    int al_matches = 0;
    float al_score = .0;
    
    // Go to start of align.
    while ( p1 < 0 || p2 < 0 ) {
      p1++;
      p2++;
    }

    // Slide across blocks.
    while ( p1 < nleft && p2 < nright ) {
      al_len++;
      if ( left.Id( p1 ) != right.Id( p2 ) )
	al_score += 1.0 / Max( left.Score( p1 ), right.Score( p2 ) );
      else
	al_matches++;
      p1++;
      p2++;
    }
    
    // Not enough blocks match.
    if ( al_len < MIN_LEN ) continue;
  
    // Special case (fw1-fw2). HEURISTICS embedded here!
    if ( special_fw1_fw2 && offset == 6 && al_matches >= 2 )
      al_matches += Min( nleft, nright );
    
    // Add to list.
    helper.push_back( SAligner( offset, al_len, al_matches, al_score ) );
  }
  
  // Sort and pick winner (return false if no winner is found).
  sort( helper.begin( ), helper.end( ) );
  if ( helper.size( ) < 1 || helper[0].matches_ < 1 ) {
    if ( log )
      *log << "\nOffset failed - no offset found between " << left.Name( )
	   << " (" << left.NMonomers( ) << " monomers), and " << right.Name( )
	   << " (" << right.NMonomers( ) << " monomers)\n";
    return false;
  }
  
  // Ok, fill matches/mismatches (if given), and return.
  offset = helper[0].offset_;
  if ( matches || mismatches ) {
    if ( matches ) *matches = 0;
    if ( mismatches ) *mismatches = 0;
    const int nleft = left.NMonomers( );
    const int nright = right.NMonomers( );
    
    int p1 = offset >= 0 ? offset : 0;
    int p2 = offset < 0 ? -offset: 0;
    while ( p1 < nleft || p2 < nright ) {
      if ( p1 >= 0 && p1 < nleft && p2 >= 0 && p2 < nright ) {
	if ( matches && left.Id( p1 ) == right.Id( p2 ) ) (*matches)++;
	if ( mismatches && left.Id( p1 ) != right.Id( p2 ) ) (*mismatches)++;
      }
      p1++;
      p2++;    
    }
  }
  
  return true;

}

/**
 * OffsetBp
 */
int OffsetBp( const CMonoRep &left,
	      const CMonoRep &right,
	      const int offset_ms )
{
  if ( offset_ms >= 0 ) {
    const look_align_plus &lhit = left.Chain( offset_ms );
    const look_align_plus &rhit = right.Chain( 0 );
    return lhit.a.pos2( ) - rhit.a.pos2( );
  }

  const look_align_plus &lhit = left.Chain( 0 );
  const look_align_plus &rhit = right.Chain( -offset_ms );
  return lhit.a.pos2( ) - rhit.a.pos2( );
  
}

