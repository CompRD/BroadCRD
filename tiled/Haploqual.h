/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef HAPLOQUAL_H
#define HAPLOQUAL_H

#include <set>
#include "math/Functions.h"
#include "Vec.h"
#include "tiled/PaddedSeq.h"

/**
 * class haploqual
 *
 * Say that reads pile up to form consensus. An haploqual is a way
 * to look at a part of a slice-column: each instance of an haploqual
 * stores one base and the quality scores for that base.
 *
 * Warning: quality scores are scored as int.
 */
class haploqual {
  
public:
  
  haploqual( ) : base_ ( 'N' ) { }

  haploqual( char base ) : base_ ( base ) { }

  haploqual( char base, const vec<int> &quals ) :
    base_ ( base ), quals_ ( quals ) { }
  
  void Set( char base, const vec<int> &quals ) {
    base_ = base;
    quals_ = quals;
  }

  void AddQual( int qual ) { quals_.push_back( qual ); }

  char Base( ) const { return base_; }

  int SumQual( ) const { return ( quals_.size( ) > 0 ) ? Sum( quals_ ) : 0; }

  int MaxQual( ) const { return ( quals_.size( ) > 0 ) ? Max( quals_ ) : 0; }

  int MinQual( ) const { return ( quals_.size( ) > 0 ) ? Min( quals_ ) : 0; }

  int Coverage( ) const { return (int)quals_.size( ); }
  
  int ConsensusQual( ) const { 
    if ( quals_.size( ) == 0 )
      return 0;
    float multiplier = ( 1.0 + ( 0.1 * ( quals_.size( ) - 1 ) ) );
    int adjusted_qual = int( float( this->MaxQual( ) ) * multiplier );
    return ( adjusted_qual > fin_qual ) ? fin_qual : adjusted_qual;
  }
  
  
private:
  
  char base_;
  vec<int> quals_;

};

/**
 * order_haploqual_Consensus
 * ordering functor
 */
struct order_haploqual_Consensus :
  public binary_function<const haploqual&, const haploqual&, bool>
{
public:
  bool operator() ( const haploqual *left, const haploqual *right ) {
    return ( left->ConsensusQual( ) < right->ConsensusQual( ) );
  }
};

/**
 * order_haploqual_Coverage
 * ordering functor
 */
struct order_haploqual_Coverage :
  public binary_function<const haploqual&, const haploqual&, bool>
{
public:
  bool operator() ( const haploqual *left, const haploqual *right ) {
    return ( left->Coverage( ) < right->Coverage( ) );
  }
};

/**
 * order_haploqual_SumQual
 * ordering functor
 */
struct order_haploqual_SumQual :
  public binary_function<const haploqual&, const haploqual&, bool>
{
public:
  bool operator() ( const haploqual *left, const haploqual *right ) {
    return ( left->SumQual( ) < right->SumQual( ) );
  }
};

#endif
