/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef BUILD_SCAFFOLD_LOCS_H
#define BUILD_SCAFFOLD_LOCS_H

#include "Superb.h"
#include "graph/Digraph.h"
#include "math/HoInterval.h"
#include "paths/reporting/CLinkBundle.h"

/**
 * class CSloc (super location)
 *
 * Place a super onto a local coordinate system.
 */
class CSloc {

public:

  CSloc( ) :
    id_ ( -1 ), rc_ ( false ),
    begin_ ( 0 ), end_ ( 0 ),
    dev_ ( 0 ), weight_ ( 0 ), score_ ( 0.0 ),
    spread_ ( ho_interval( 0, 0 ) ),
    spread2_ ( ho_interval( 0, 0 ) ) { }

  CSloc( int id, bool rc, int begin, int end,
	 int dev, int weight, float score,
	 ho_interval spread,
	 ho_interval spread2 ) :
    id_( id ), rc_( rc ),
    begin_ ( begin ), end_ ( end ),
    dev_ ( dev ), weight_ ( weight ), score_ ( score ),
    spread_ ( spread ), spread2_ ( spread2 ) { }

  void Set( int id, bool rc, int begin, int end,
	    int dev, int weight, int score,
	    ho_interval spread, ho_interval spread2 ) {
    id_ = id;
    rc_ = rc;
    begin_ = begin;
    end_ = end; 
    dev_ = dev;
    weight_ = weight;
    score_ = score;
    spread_ = spread;
    spread2_ = spread2;
  }

  int Id( ) const { return id_; }
  bool Rc( ) const { return rc_; }
  int Begin( ) const { return begin_; }
  int End( ) const { return end_; }
  int Dev( ) const { return dev_; }
  int Weight( ) const { return weight_; }
  float Score( ) const { return score_; }
  ho_interval Spread( ) const { return spread_; }
  ho_interval Spread2( ) const { return spread2_; }

  int OverlapWith( const CSloc &other, const double &MAX_STRETCH = 3.0 ) {
    int sdl = (int)( MAX_STRETCH * (double)this->Dev( ) );
    int sdr = (int)( MAX_STRETCH * (double)other.Dev( ) );
    int left = Max( this->Begin( ) - sdl, other.Begin( ) + sdr );
    int right = Min( this->End( ) - sdl , other.End( ) + sdr );
    return Max( 0, right - left );
  }

  vec<String> TableLine( ) const {
    vec<String> result( 6, "." );

    String strw1 = ToString( spread_.Start( ) );
    String strw2 = ToString( spread_.Stop( ) );
    result[0] = ToString( id_ ) + ( rc_ ? "rc" : "fw" );
    result[1] = "[" + ToString( begin_ ) + "," + ToString( end_ ) + ")";
    if ( weight_ == 0 && score_ == 0.0 ) return result;

    result[2] = " +/- " + ToString( dev_ );
    result[3] = "w=" + ToString( weight_ );
    result[4] = "s=" + ToString( score_, 2 );
    result[5] = "[" + strw1 + ", " + strw2 + ")";
    return result;
  }

  friend bool operator< ( const CSloc &left, const CSloc &right ) {
    if ( left.Begin( ) != right.Begin( ) )
      return ( left.Begin( ) < right.Begin( ) );
    if ( left.End( ) != right.End( ) )
      return ( left.End( ) < right.End( ) );
    if ( left.Weight( ) != right.Weight( ) )
      return ( left.Weight( ) > right.Weight( ) );
    if ( left.Score( ) != right.Score( ) )
      return ( left.Score( ) < right.Score( ) );
    if ( left.Rc( ) != right.Rc( ) )
      return ( left.Rc( ) < right.Rc( ) );
    return ( left.Id( ) < right.Id( ) );
    // Forget about spread_ (should do it for completeness).
  }


private:

  int id_;
  bool rc_;
  int begin_;
  int end_;
  int dev_;
  int weight_;
  float score_;
  ho_interval spread_;
  ho_interval spread2_;

};

/**
 * BuildScaffoldLocs
 *
 * The set of all and only the edges connected (directly) to the given
 * super (the "center").
 */
void BuildScaffoldLocs( vec<CSloc> &locs,
			const int s1,			
			const digraphE<CLinkBundle> &graph,
			const vec<superb> &supers,
			const int MIN_LINKS );

#endif
