/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Superb.h"
#include "math/HoInterval.h"
#include "paths/BuildScaffoldLocs.h"
#include "paths/reporting/CLinkBundle.h"

/**
 * BuildScaffoldLocs
 */
void BuildScaffoldLocs( vec<CSloc> &locs,
			const int s1,			
			const digraphE<CLinkBundle> &graph,
			const vec<superb> &supers,
			const int MIN_LINKS )
{
  locs.clear( );

  // Vertex id and length of super.
  int v1 = 2 * s1;
  int len1 = supers[s1].TrueLength( );

  // The center.
  ho_interval bogus( 0, 0 );
  locs.resize( 1, CSloc( s1, false, 0, len1, 0, 0, 0.0, bogus, bogus ) );

  // Add "to" edges.
  vec<int> tos = graph.To( v1 );
  for (int ii=0; ii<tos.isize( ); ii++) {
    int v2 = tos[ii];
    int s2 = v2 / 2;
    int len2 = supers[s2].TrueLength( );
    bool rc2 = ( 1 == v2 % 2 );
    vec<int> bundle_ids = graph.EdgesBetween( v2, v1 );
    ForceAssertEq( bundle_ids.isize( ), 1 );
    const CLinkBundle &bundle = graph.EdgeObject( bundle_ids[0] );
    int weight = bundle.Weight( );
    if ( weight < MIN_LINKS ) continue;
    int sep = bundle.Sep( );
    int beg = - sep - len2;
    int end = beg + len2;
    int dev = bundle.Dev( );
    float score = bundle.Score( );
    pair<double,double> spread = bundle.Win2( );
    pair<double,double> spread2 = bundle.Win1( );
    int w1_beg = spread.first * len1 / 100.0;
    int w1_end = spread.second * len1 / 100.0;
    int w2_beg = spread2.first * len2 / 100.0;
    int w2_end = spread2.second * len2 / 100.0;
    ho_interval win1( w1_beg, w1_end );
    ho_interval win2( w2_beg, w2_end );
    CSloc newloc( s2, rc2, beg, end, dev, weight, score, win1, win2 );
    locs.push_back( newloc );
  }

  // Add "from" edges.
  vec<int> froms = graph.From( v1 );
  for (int ii=0; ii<froms.isize( ); ii++) {
    int v2 = froms[ii];
    int s2 = v2 / 2;
    int len2 = supers[s2].TrueLength( );
    bool rc2 = ( 1 == v2 % 2 );
    vec<int> bundle_ids = graph.EdgesBetween( v1, v2 );
    ForceAssertEq( bundle_ids.isize( ), 1 );
    const CLinkBundle &bundle = graph.EdgeObject( bundle_ids[0] );
    int weight = bundle.Weight( );
    if ( weight < MIN_LINKS ) continue;
    int sep = bundle.Sep( );
    int beg = sep + len1;
    int end = beg + len2;
    int dev = bundle.Dev( );
    float score = bundle.Score( );
    pair<double,double> spread = bundle.Win1( );
    pair<double,double> spread2 = bundle.Win2( );
    int w1_beg = spread.first * len1 / 100.0;
    int w1_end = spread.second * len1 / 100.0;
    int w2_beg = spread2.first * len2 / 100.0;
    int w2_end = spread2.second * len2 / 100.0;
    ho_interval win1( w1_beg, w1_end );
    ho_interval win2( w2_beg, w2_end );
    CSloc newloc( s2, rc2, beg, end, dev, weight, score, win1, win2 );
    locs.push_back( newloc );
  }

  // Sort.
  sort( locs.begin( ), locs.end( ) );

}
