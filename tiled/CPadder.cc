/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "Qualvector.h"
#include "tiled/CharBaseUtil.h"
#include "tiled/CPadder.h"
#include "tiled/PaddedSeq.h"
#include "tiled/Tiling.h"
#include "tiled/TilingAnalyzer.h"
#include "tiled/TilingPrinter.h"

/**
 * CPadder
 * Constructor
 */
CPadder::CPadder( const vecbasevector &cbases,
		  const vecbasevector &rbases,
		  ostream *log ) :
  cbases_ ( cbases ),
  rbases_ ( rbases ),
  log_ ( log ),
  before_out_ ( 0 ),
  after_out_ ( 0 )
{ }

/**
 * CPadder
 * SetOutStreams
 */
void CPadder::SetOutStreams( ostream *before_out, ostream *after_out )
{
  before_out_ = before_out;
  after_out_ = after_out;
}

/**
 * CPadder
 * BlockJustifyLeft
 *
 * Set tiling and run BlockJustifyLeft( ) on all reads.
 */
void CPadder::BlockJustifyLeft( tiling *theTiling )
{
  // Print consensus before shifting.
  if ( before_out_ ) this->PrintTiling( *theTiling, before_out_ );

  // Run BlockJustifyLeft.
  this->SetTiling( theTiling );
  for (int rpos=0; rpos<theTiling_->ReadsCount( ); rpos++)
    while ( this->BlockJustifyLeft( rpos ) > 0 ) ;

  // Print consensus after shifting.
  if ( after_out_ ) this->PrintTiling( *theTiling_, after_out_ );
}

/**
 * CPadder
 * BlockJustifyLeft
 *
 * All pads are treated in blocks (if a pad belongs to a group of
 * three, then all three pads will be shifted at the same time), in
 * other words, no group of pads will be broken. The read is
 * identified by its position in the tiling. It returns the number of
 * shifting events performed.
 */
int CPadder::BlockJustifyLeft( int pos )
{
  ofstream devnull ( "/dev/null" );
  ostream &log = log_ ? *log_ : devnull;

  // Keep track of events and shifts.
  int n_events = 0;
  padded_seq &rpads = theTiling_->Read( pos ).PaddedSeq( );

  // Loop over pads' blocks.
  int padblock_begin = 0;
  while ( padblock_begin < rpads.PadsCount( ) ) {
    int begin = padblock_begin;
    int padblock_end = begin + 1;
    for (int ii=begin+1; ii<rpads.PadsCount( ); ii++) {
      if ( rpads.Pad( ii-1 ) + 1 < rpads.Pad( ii ) )
	break;
      padblock_end = ii + 1;
    }
    int npads = padblock_end - padblock_begin;

    // It may happen that the pads have been pushed all the way on the left.
    if ( rpads.Pad( padblock_begin ) == 0 ) {
      padblock_begin = padblock_end;
      continue;
    }

    // Cannot push the block too much left (only >= wall is allowed).
    int wall = 0;
    if ( padblock_begin > 0 )
      wall = rpads.Pad( padblock_begin - 1 ) + 1;
    ForceAssert( wall >= 0 && wall < rpads.Pad( padblock_begin ) );
    
    // Get contig and read as vec<char>'s.
    tiling_analyzer theAnalyzer;
    theAnalyzer.SetPointers( &rbases_, 0, &cbases_, 0 );
    theAnalyzer.SetTiling( theTiling_ );

    pair<int,int> rwin = theAnalyzer.RWin( pos );
    vec<int> rpos;
    rpos.push_back( -1 );
    rpos.push_back( pos );
    int post = rwin.first;
    int width = rwin.second - rwin.first;

    vec< vec<char> > rb;
    vec< vec<int> > rq;
    theAnalyzer.TamedRectangles( post, width, rpos, rb, rq );
    ForceAssert( rb.size( ) == 2 );

    // Find shift amount.
    int shift_amt = 0;
    for (int ii=rpads.Pad(padblock_begin)-1; ii>=wall; ii--) {
      char rbase = rb[1][ii];
      char cbase = rb[0][ii+npads];
      
      if ( IsEmpty( rbase ) || IsEmpty( cbase ) )
	continue;
      ToUpperCase( rbase );
      ToUpperCase( cbase );
      if ( rbase != cbase )
	break;
      
      shift_amt++;
    }

    if ( shift_amt > 0 ) {
      n_events++;
      for (int ii=padblock_begin; ii<padblock_end; ii++)
	rpads[ii] += -shift_amt;
      log << theTiling_->Contig( ).PrettyName( 0 ) << "   "
	  << theTiling_->Read( pos ).PrettyName( 0 ) << "   <"
	  << post + rpads.Pad( padblock_begin ) << ","
	  << post + rpads.Pad( padblock_end - 1 ) << ">   shift="
	  << shift_amt << "\n";
    }

    padblock_begin = padblock_end;
  }
  
  // Ok, return.
  return n_events;
}

/**
 * CPadder
 * PrintAlign
 */
void CPadder::PrintAlign( int pos, ostream &out ) const
{
  tiling_analyzer theAnalyzer;
  theAnalyzer.SetPointers( &rbases_, 0, &cbases_, 0 );
  theAnalyzer.SetTiling( theTiling_ );

  pair<int,int> rwin = theAnalyzer.RWin( pos );
  
  int width = 79;
  vec<int> posts;
  int pos_on_master = rwin.first;
  while ( pos_on_master < rwin.second ) {
    posts.push_back( pos_on_master );
    pos_on_master += width;
  }

  const padded_seq &cpads = theTiling_->Contig( ).PaddedSeq( );
  out << "Align of " << theTiling_->Read( pos ).PrettyName( 0 )
      << " on " << theTiling_->Contig( ).PrettyName( 0 )
      << "   padded_[" << rwin.first
      << ", " << rwin.second
      << ")   unpadded_[" << cpads.ToUnpaddedPos( rwin.first )
      << ", " << cpads.ToUnpaddedPos( rwin.second )
      << ")\n\n";
  
  vec<int> rpos;
  rpos.push_back( -1 );
  rpos.push_back( pos );
  vec< vec<char> > rb;
  vec< vec<int> > rq;
  for(int ii=0; ii<(int)posts.size( ); ii++) {
    int post = posts[ii];
    theAnalyzer.TamedRectangles( post, width, rpos, rb, rq );

    vec<char> tags( rb[0].size( ), ' ' );
    if ( rb.size( ) == 2 ) {
      for (int kk=0; kk<(int)rb[0].size( ); kk++) {
	char cbase = rb[0][kk];
	char rbase = rb[1][kk];
	if ( ! IsEmpty( cbase ) ) ToUpperCase( cbase );
	if ( ! IsEmpty( rbase ) ) ToUpperCase( rbase );
	if ( cbase == rbase ) continue;
	if ( IsEmpty( cbase ) || IsEmpty( rbase ) ) continue;
	if ( IsGap( cbase ) && IsGap( rbase ) ) continue;
	if ( IsGap( cbase ) || IsGap( rbase ) ) { tags[kk] = '|'; continue; }
	tags[kk] = '*';
      }
    }

    for (int kk=0; kk<(int)tags.size( ); kk++)
      out << tags[kk];
    out << "\n";
    for (int jj=0; jj<(int)rb.size( ); jj++) {
      for (int kk=0; kk<(int)rb[jj].size( ); kk++)
	out << rb[jj][kk];
      out << "\n";
    }
    out << "\n";
  }
  out << "end_of_align\n\n";
}

/**
 * CPadder
 * PrintTiling
 */
void CPadder::PrintTiling( const tiling &tiles, ostream *out ) const
{
  tiling_analyzer analyzer;
  analyzer.SetPointers( &rbases_, 0, &cbases_, 0 );
  analyzer.SetTiling( &tiles );
  tiling_printer printer( &analyzer );
  printer.SetPrintQuals( false );
  printer.ToStream( *out );
}

