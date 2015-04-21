/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Basevector.h"
#include "LocsHandler.h"
#include "ParseSet.h"
#include "PrettyPrintTable.h"
#include "math/Functions.h"
#include "pairwise_aligners/AlignTwoBasevectors.h"
#include "pairwise_aligners/KmerAligner.h"
#include "tiled/CRealign.h"
#include "tiled/TAlign.h"

/**
 * class CRealign
 * Constructor
 *
 * log: optional log stream
 * dbfastalog: optional debug fasta log (save fasta of failed attempts)
 */
CRealign::CRealign( const vecbasevector &cbases,
		    const vecbasevector &rbases,
		    const lhandler &locs,
		    ostream *log,
		    ostream *dbfastalog ) :
  cbases_ ( cbases ),
  rbases_ ( rbases ),
  locs_ ( locs ),
  plog_ ( log ),
  dbfastalog_ ( dbfastalog )
{
  this->SetDefaults( );
}

/**
 * class CRealign
 * SetMinMaxBand
 */
void CRealign::SetMinMaxBand( int min, int max )
{
  ForceAssertLt( 0, min );
  ForceAssertLe( min, max );
  min_band_ = min;
  max_band_ = max;
}

/**
 * class CRealign
 * SetContig
 */
void CRealign::SetContig( int contig_id )
{
  contig_id_ = contig_id;
  set12_ = false;
  aligner24_.SetBases( cbases_[contig_id_] );
}

/**
 * class CRealign
 * PlaceLoc
 *
 * On error it will return false and fill al with a default one-base
 * alignment.
 */
bool CRealign::PlaceLoc( int locpos, alignment &al )
{
  ofstream devnull ( "/dev/null" );
  ostream &log = plog_ ? *plog_ : devnull;
  
  // Generate seeds (this will set locpos_).
  this->FindSeeds( locpos );
  this->ReportSeeds( );
  
  // Determine offset and bandwidth.
  int offset = -1;
  int band = -1;
  if ( ! this->GetOffset( offset, band ) )
    return false;
  
  // Incrementally try with a larger and larger band.
  int n_attempts = 3;
  bool al_found = false;
  for (int attempt=0; attempt<n_attempts; attempt++) {
    int current_band = ( attempt + 1 ) * band;
    
    // Chop contig around the selected target interval.
    int contig_len = cbases_[contig_id_].size( );
    int read_len = locs_[locpos_].LengthOfRead( );
    int beg = Max( 0, offset - current_band );
    int end = Min( contig_len, offset + read_len + current_band );
    int len = end - beg;
    
    cgwin_ = make_pair( beg, end );
    the_contig_.SetToSubOf( cbases_[contig_id_], beg, len );
    
    // Align not found.
    al_found = ( this->RunAlign2Bvecs( al ) );
    if ( ! al_found )
      continue;

    // Adjust alignment.
    int pos1, pos2, errors;
    avector<int> gaps, lens;
    al.Unpack( pos1, pos2, errors, gaps, lens );
    al.Set( pos1 + beg, pos2, errors, gaps, lens );
    break;
  }

  // Align not found, set it to default one-base alignment.
  if ( ! al_found ) {
    this->SetToOneBaseAlign( locpos_, al );
    this->DumpFasta( );
  }
  
  // Return.
  return al_found;
}

/**
 * class CRealign
 * PlaceLoc
 *
 * Wrapper around the real PlaceLoc
 */
bool CRealign::PlaceLoc( int locpos, t_align &tal )
{
  alignment al;
  int read_id = locs_[locpos].ReadId( );
  bool rc = locs_[locpos].OrientationOnContig( ) == ReverseOr;

  bool align_found = this->PlaceLoc( locpos, al );

  int pos1 = 0;
  int pos2 = 0;
  int errors = 0;
  avector<int> gaps;
  avector<int> lens;
  al.Unpack( pos1, pos2, errors, gaps, lens );
  align shortal( pos1, pos2, gaps, lens );
  shortal.Flip( );

  tal.Set( read_id, rc, shortal );
  
  return align_found;
}

/**
 * class CRealign
 * FindSeeds
 */
void CRealign::FindSeeds( int locpos )
{
  // Initial setup.
  locpos_ = locpos;
  seeds_.clear( );

  // Package the_read_.
  int read_id = locs_[locpos].ReadId( );
  bool rc = locs_[locpos].OrientationOnContig( ) == ReverseOr;
  the_read_ = rbases_[read_id];
  if ( rc ) the_read_.ReverseComplement( );
  
  // Try seeding with 24-mers.
  vec<int> rseeds;
  aligner24_.FindPossibleAlignments( the_read_, rseeds );

  // If needed, try again with 12-mers.
  if ( rseeds.size( ) < 1 ) {
    if ( ! set12_ ) {
      aligner12_.SetBases( cbases_[contig_id_] );
      set12_ = true;
    }
    aligner12_.FindPossibleAlignments( the_read_, rseeds );
  }

  // No seeds found.
  if ( rseeds.size( ) < 1 )
    return;
      
  // Package seeds_ and sort.
  sort( rseeds.begin( ), rseeds.end( ) );
  seeds_.push_back( make_pair( rseeds[0], 1 ) );
  for (int ii=1; ii<(int)rseeds.size( ); ii++) {
    if ( rseeds[ii] == rseeds[ii-1] ) seeds_.back( ).second += 1;
    else seeds_.push_back( make_pair( rseeds[ii], 1 ) );
  }
  sort( seeds_.begin( ), seeds_.end( ) );
  
}

/**
 * class CRealign
 * ReportSeeds
 */
void CRealign::ReportSeeds( ) const
{
  if ( ! plog_ ) return;
  ostream &log = *plog_;

  ForceAssert( locpos_ > -1 );

  int floc = locs_.FirstLoc( contig_id_ );
  int nreads = locs_.ReadsInContig( contig_id_ );
  int read_id = locs_[locpos_].ReadId( );
  int read_len = locs_[locpos_].LengthOfRead( );
  String RC = ReverseOr == locs_[locpos_].OrientationOnContig( ) ? "_rc" : "";
  log << "contig_" << contig_id_
      << "." << locpos_ - floc 
      << "/" << nreads - 1
      << " r_" << read_id << RC
      << " len=" << read_len
      << " ";
  
  if ( seeds_.size( ) < 1 ) {
    log << "no_seeds_found\n";
    return;
  }

  int wbeg = seeds_[0].first;
  int wend = seeds_.back( ).first;
  int wlen = wend - wbeg;
  
  int original_start = locs_[locpos_].StartOnContig( );
  int cpos = this->ClosestSeed( );
  int cdist = this->Discrepancy( cpos );
  log << "seeds=[" << wbeg
      << "," << wend
      << ") wedge="
      << wlen << " discrep="
      << cdist << "\n";
  
  if ( verbose_ ) {
    vec< vec<String> > table;
    for (int ii=0; ii<(int)seeds_.size( ); ii++) {
      int dist = Abs( original_start - seeds_[ii].first );
      String str_best = ii == cpos ? "closest" : "";
      
      vec<String> aline;
      aline.push_back( "new=@" + ToString( seeds_[ii].first ) );
      aline.push_back( "orig=@" + ToString( original_start ) );
      aline.push_back( "w=" + ToString( seeds_[ii].second ) );
      aline.push_back( "discrep=" + ToString( dist ) );
      aline.push_back( str_best );
      
      table.push_back( aline );
    }
    
    BeautifyTable( table );
    for(int ii=0; ii<(int)table.size( ); ii++) {
      for (int jj=0; jj<(int)table[ii].size( ); jj++)
	log << "   " << table[ii][jj];
      log << "\n";
    }
  }
  
}

/**
 * class CRealign
 * SetDefaults
 */
void CRealign::SetDefaults( )
{
  min_band_ = 24;
  max_band_ = 512;
  max_discrep_ = 24;
  max_error_ = 0.35;
  verbose_ = false;
  set12_ = false;

  contig_id_ = -1;
  locpos_ = -1;
}

/**
 * class CRealign
 * DumpFasta
 *
 * Dump fasta of read (first) and contig (second) on the dbfastalog stream.
 */
void CRealign::DumpFasta( )
{
  if ( ! dbfastalog_ ) return;
  
  String rc = locs_[locpos_].OrientationOnContig( ) == ReverseOr ? "_rc" : "";
  String read_id = "read_" + ToString( locs_[locpos_].ReadId( ) ) + rc;
  
  String str_beg = ToString( cgwin_.first );
  String str_end = ToString( cgwin_.second );
  String str_win = "_[" + str_beg + "," + str_end + ")";
  String contig_id = "contig_" + ToString( contig_id_ ) + str_win;

  the_read_.Print( *dbfastalog_, read_id );
  the_contig_.Print( *dbfastalog_, contig_id );

  *dbfastalog_ << flush;
}

/**
 * class CRealign
 * GetOffset
 *
 * Select offset and bandwidth (returns false on failure).
 */
bool CRealign::GetOffset( int &offset, int &band ) const
{
  offset = -1;
  band = -1;
  if ( seeds_.size( ) < 1 )
    return false;
  
  // Total weight (how many perfect matches total).
  int tot_weight = 0;
  for (int ii=0; ii<(int)seeds_.size( ); ii++)
    tot_weight += seeds_[ii].second;
  
  // Tag as valid a seed if it is close enough to the original guess.
  vec<Bool> is_valid;
  int tot_valid_weight = 0;
  if ( max_discrep_ < 0 ) {
    tot_valid_weight = tot_weight;
    is_valid.resize( seeds_.size( ), True );
  }
  else {
    is_valid.resize( seeds_.size( ), False );
    for (int ii=0; ii<(int)seeds_.size( ); ii++) {
      if ( Abs( this->Discrepancy( ii ) ) <= max_discrep_ ) {
	tot_valid_weight += seeds_[ii].second;
	is_valid[ii] = True;
      }
    }
  }

  // The best seed.
  int bestpos = -1;
  int minpos = -1;
  int maxpos = -1;
  for (int ii=0; ii<(int)seeds_.size( ); ii++) {
    if ( ! is_valid[ii] )
      continue;
    if ( bestpos < 0 ) {
      bestpos = ii;
      minpos = ii;
      maxpos = ii;
      continue;
    }
    if ( seeds_[ii].second > seeds_[bestpos].second )
      bestpos = ii;
    if ( seeds_[ii].first < seeds_[minpos].first )
      minpos = ii;
    if ( seeds_[ii].first > seeds_[maxpos].first )
      maxpos = ii;
  }
  if ( bestpos < 0 || minpos < 0 || maxpos < 0 )
    return false;
  int wedge = seeds_[maxpos].first - seeds_[minpos].first;
  
  // Fill offset and band, and return.
  offset = seeds_[bestpos].first;
  band = Max( min_band_, Min( wedge, max_band_ ) );
  return true;
}

/**
 * class CRealign
 * ClosestSeed
 */
int CRealign::ClosestSeed( ) const
{
  if ( seeds_.size( ) < 1 ) return -1;
  
  int bpos = 0;
  int dist = Abs( this->Discrepancy( 0 ) );
  for (int ii=1; ii<(int)seeds_.size( ); ii++) {
    int this_dist = Abs( this->Discrepancy( ii ) );
    if ( this_dist > dist ) continue;
    dist = this_dist;
    bpos = ii;
  }
  
  return bpos;
}

/**
 * class CRealign
 * Discrepancy
 */
int CRealign::Discrepancy( int seedpos ) const
{
  int original_start = locs_[locpos_].StartOnContig( );
  int dist = seeds_[seedpos].first - original_start;

  return dist;
}

/**
 * class CRealign
 * SetToOneBaseAlign
 */
void CRealign::SetToOneBaseAlign( int locpos, alignment &al ) const
{
  int c_id = locs_[locpos].Contig( );
  int r_id = locs_[locpos].ReadId( );
  int c_len = cbases_[c_id].size( );

  int pos1 = c_len - 1;
  int pos2 = 0;
  int errors = cbases_[c_id][pos1] == rbases_[r_id][pos2] ? 1 : 0;
  avector<int> gaps( 1, 0 );
  avector<int> lens( 1, 1 );

  al.Set( pos1, pos2, errors, gaps, lens );
}

/**
 * class CRealign
 * RunAlign2Bvecs
 *
 * First try with K=24, then with K=8.
 */
bool CRealign::RunAlign2Bvecs( alignment &al, bool K8 ) const
{
  const bvec &b1 = the_contig_;
  const bvec &b2 = the_read_;
  
  ofstream devnull( "/dev/null" );
  int minOv = 0;
  int maxOv = b1.size( ) + b2.size( );
  float maxErrRate = 1.0;
  ostream &log = devnull;
  int rc = false;
  int mode = 5;
  int K = K8 ? 8 : 24;
  int stretch = 12;
  int nstretch = 6;
  const qvec q1( 0 );
  const qvec q2( 0 );
  float maxScore = 0;
  float maxErr = 0;
  ostream &badlog = devnull;
  bool avoidProm = False;
  int maxCliq8 = 1000;
  int maxAligns8 = 1000;
  int maxErr8 = 1000;
  int locMaxErr = 1000;
  bool altMethod = False;
  int band = 0;
  int minMutmer = 0;
  float ambigThreshold = 3.0;
  double minPolyscoreRatio = 100.0;
  bool swMethod = False;

  if ( plog_ ) {
    int read_id = locs_[locpos_].ReadId( );
    orientation orient = locs_[locpos_].OrientationOnContig( );
    String strK = K8 ? "8" : "24";
    String strClen = ToString( cbases_[contig_id_].size( ) );
    String strWin = ToString( cgwin_.first ) + "," + ToString( cgwin_.second );
    String strOr = ( orient == ForwardOr ) ? "_fw" : "_rc";
    
    *plog_ << "Trying K=" << strK
	   << " c_" << contig_id_
	   << " @[" << strWin
	   << ")_" << strClen
	   << " vs r_" << read_id << strOr
	   << ": ";
  }

  align theAl;
  AlignTwoBasevectors( b1, b2, theAl, minOv, maxOv, maxErrRate, &log, rc,
		       mode, K, stretch, nstretch, q1, q2, maxScore, maxErr,
		       badlog, avoidProm, maxCliq8, maxAligns8, maxErr8,
 		       locMaxErr, altMethod, band, minMutmer, ambigThreshold,
 		       minPolyscoreRatio, swMethod );
  al.Set( theAl, ActualErrors( b1, b2, theAl, 1, 1 ) );

  bool align_found = IsValid( al );
  if ( align_found || K8 )
    return align_found;
  
  return this->RunAlign2Bvecs( al, true );
}

/**
 * class CRealign
 * IsValid
 */
bool CRealign::IsValid( const alignment &al ) const
{
  ofstream devnull ( "/dev/null" );
  ostream &log = plog_ ? *plog_ : devnull;

  const bvec &b1 = the_contig_;
  const bvec &b2 = the_read_;

  // No align found.
  int al_len = al.Pos1( ) - al.pos1( );
  if ( al_len < 1 ) {
    log << "no align found\n";
    return false;
  }

  // Too many errors.
  int al_errors = al.Errors( );
  float error_rate = SafeQuotient( al_errors, al_len );
  if ( max_error_ < 1.0 && error_rate > max_error_ ) {
    log << "high error rate (" << ToString( error_rate, 1 ) << ")\n";
    return false;
  }
  
  // Not a proper align.
  int pos1 = cgwin_.first + al.pos1( );
  int Pos1 = cgwin_.first + al.Pos1( );
  int pos2 = al.pos2( );
  int Pos2 = al.Pos2( );
  int read_len = the_read_.size( );
  int contig_len = cbases_[contig_id_].size( );

  if ( pos1 > 0 && pos2 > 0 ) {
    log << "improper (pos1=" << pos1 << ", pos2=" << pos2 << ")\n";
    return false;
  }
  if ( Pos1 < contig_len && Pos2 < read_len ) {
    String str1 = "contig_len-Pos1=" + ToString( contig_len - Pos1 );
    String str2 = "read_len-Pos2=" + ToString( read_len - Pos2 );
    log << "improper (" << str1 << ", " << str2 << ")\n";
    
    return false;
  }
  
  // Ok.
  log << "proper!\n";
  return true;
}

