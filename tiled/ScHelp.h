/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef SC_HELP_H
#define SC_HELP_H

#include "ReadScore.h"
#include "Vec.h"

/**
 * class scHelp
 *
 * Score Helper: used to keep track of chunks of reads over a given
 * consensus window.
 */
class scHelp {

public:
  
  scHelp( ) :
    Score_ ( -1.0 ), nValid_ ( -1 ), Pos_ ( -1 ), Tag_ ( "" ) { }
  
  scHelp( Float score, int n_valid, int pos, const String &tag ) :
    Score_ ( score ), nValid_ ( n_valid ), Pos_ ( pos ), Tag_ ( tag ) { }
  
  scHelp( const vec< vec<int> > &quals, const read_score &scorer, int pos ) {
    this->Set( quals, scorer, pos );
  }
  
  // Set from a vec<int>.
  void Set( const vec<int> &quals, const read_score &scorer, int pos ) {
    Score_= scorer.Score( quals, 0, (int)quals.size( ) );
    nValid_ = 0;
    for (int ii=0; ii<(int)quals.size( ); ii++)
      if ( quals[ii] != empty_qual && quals[ii] != gap_qual )
	nValid_++;
    Pos_ = pos;
    Tag_ = "";
  }

  // Set from given Rectangle.
  void Set( const vec< vec<int> > &quals, const read_score &scorer, int pos ) {
    this->Set( quals[pos], scorer, pos );
  }
  
  // Sort by:  Score / nValid (high to low) / pos.
  friend bool operator< ( const scHelp &left, const scHelp &right ) {
    if ( left.Score_ == right.Score_ ) {
      if ( left.nValid_ == right.nValid_ ) {
	return ( left.Pos_ < right.Pos_ );
      }
      return ( left.nValid_ > right.nValid_ );
    }
    return ( left.Score_ < right.Score_ );
  }

  // Simple operator<<.
  friend ostream& operator<< ( ostream &out, const scHelp &helper ) {
    out << ToString( helper.Score_, 6 ) << "\t"
	<< helper.nValid_ << "\t"
	<< helper.Pos_ << "\t"
	<< helper.Tag_;
    return out;
  }


public:
  
  float Score_;  // score of the read subset
  int nValid_;   // count of real bases (i.e. not empty or gap)
  int Pos_;      // position in a given rectangle
  String Tag_;   // all-purpose tag

};

#endif
