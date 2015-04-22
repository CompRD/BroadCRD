///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "paths/reporting/CLinkBundle.h"

/**
 * CLinkBundle
 * Constructor
 */
CLinkBundle::CLinkBundle( ) :
  sep_ ( 0 ),
  dev_ ( 0 ),
  weight_ ( 0 ),
  score_ ( 0.0 ),
  win1_ ( make_pair( 0.0, 0.0 ) ),
  win2_ ( make_pair( 0.0, 0.0 ) )
{ }

/**
 * CLinkBundle
 * Constructor
 */
CLinkBundle::CLinkBundle( int sep, int dev, int weight, float score,
 			  pair<double,double> win1,
			  pair<double,double> win2 ) :
  sep_ ( sep ),
  dev_ ( dev ),
  weight_ ( weight ),
  score_ ( score ),
  win1_ ( win1 ),
  win2_ ( win2 )
{ }
  
/**
 * CLinkBundle
 * Constructor
 */
void CLinkBundle::Set( int sep, int dev, int weight, float score,
		       pair<double,double> win1,
		       pair<double,double> win2 )
{
  sep_ = sep;
  dev_ = dev;
  weight_ = weight;
  score_ = score;
  win1_ = win1;
  win2_ = win2;
}

/**
 * CLinkBundle
 * ToString
 */
String CLinkBundle::AsString( bool brief ) const 
{
  if ( brief )
    return
      ToString( sep_ ) + "+/-" + ToString( dev_ )
      + " w=" + ToString( weight_ ) + " s=" + ToString( score_, 2 );

  return
    ToString( sep_ ) + "+/-" + ToString( dev_ )
    + " w=" + ToString( weight_ ) + " s=" + ToString( score_, 2 )
    + " win1=[" + ToString( win1_.first, 2 ) + ","
    + ToString( win1_.second, 2 ) + ")"
    + " win2=[" + ToString( win2_.first, 2 ) + ","
    + ToString( win2_.second, 2 ) + ")";
}

/**
 * CLinkBundle
 * CombinedScore
 *
 * It returns an overall combined score of Score( ) and Weight( ). The
 * returned value is in the range (0, sqrt( Weight( ) ], where higher
 * is better.
 *
 * HEURISTICS: use CENTER to decide where the real distribution is
 * centered (in theory it should be 1.0).
 *
 * HEURISTICS: use DAMPING to damp the effect of the exponential
 * (at 1.0 it would decay too fast).
 */
double CLinkBundle::CombinedScore( ) const 
{
  const double CENTER = 1.0;
  const double DAMPING = 100.0;
  double N = (double)weight_;
  double sig = ( score_ - CENTER ) * ( score_ - CENTER );
  return sqrt( N ) * exp( - N * sig / DAMPING );
}

/**
 * Template initializations
 */
#include "graph/DigraphTemplate.h"
template digraphE<CLinkBundle>::digraphE();
template void digraphE<CLinkBundle>::AddVertices(int);
template void digraphE<CLinkBundle>::DeleteEdges(const vec<int>&);
template void digraphE<CLinkBundle>::DeleteEdges(const vec<int>&, const vec<int>&);
template void digraphE<CLinkBundle>::DeleteEdgesAtVertex(int);
template void digraphE<CLinkBundle>::DeleteEdgeFrom(int, int);
template const CLinkBundle& digraphE<CLinkBundle>::EdgeObject(int) const;
template CLinkBundle const& digraphE<CLinkBundle>::EdgeObjectByIndexFrom(int, int) const;
template int digraphE<CLinkBundle>::EdgeObjectCount() const;
template int digraphE<CLinkBundle>::EdgeObjectIndexByIndexFrom(int, int) const;
template int digraphE<CLinkBundle>::EdgeObjectIndexByIndexTo(int, int) const;
template vec<CLinkBundle> const& digraphE<CLinkBundle>::Edges() const;
template vec<int> digraphE<CLinkBundle>::EdgesBetween(int, int) const;
template vec<CLinkBundle> digraphE<CLinkBundle>::EdgeObjectsBetween(const int, const int) const;
template vec<vec<int> > const& digraphE<CLinkBundle>::FromEdgeObj() const;
template vec<int> const& digraphE<CLinkBundle>::FromEdgeObj(int) const;
template void digraphE<CLinkBundle>::GiveEdgeNewFromVx(int, int, int);
template void digraphE<CLinkBundle>::GiveEdgeNewToVx(int, int, int);
template void digraphE<CLinkBundle>::Initialize(const vec< vec<int> >&, const vec< vec<int> >&, const vec<CLinkBundle>&, const vec< vec<int> >&, const vec< vec<int> >&, const Bool);
template int digraphE<CLinkBundle>::InputFromOutputTo(int, int) const;
template int digraphE<CLinkBundle>::InputToOutputFrom(int, int) const;
template void digraphE<CLinkBundle>::readBinary(BinaryReader&);
template vec<int> digraphE<CLinkBundle>::RemoveDeadEdgeObjects();
template vec<vec<int> > const& digraphE<CLinkBundle>::ToEdgeObj() const;
template vec<int> const& digraphE<CLinkBundle>::ToEdgeObj(int) const;
template void digraphE<CLinkBundle>::ToLeft(vec<int>&) const;
template void digraphE<CLinkBundle>::ToRight(vec<int>&) const;
template void digraphE<CLinkBundle>::Used(vec<unsigned char>&) const;
template void digraphE<CLinkBundle>::writeBinary(BinaryWriter&) const;
