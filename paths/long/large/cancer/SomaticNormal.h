///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Return hg38 coordinates for regions of the genome that are often rearranged
// in normal cells.  These would usually involve antibodies.

#ifndef SOMATIC_NORMAL_H
#define SOMATIC_NORMAL_H

#include "CoreTools.h"
#include "math/HoInterval.h"

// Report intervals on grch38, or hg19 if alt = True.

inline void SomaticNormal( vec< pair<int,ho_interval> >& sn, const Bool alt = False )
{    sn.clear( );

     if ( !alt )
     {    sn.push( 1,  ho_interval(  88859000,  89313000 ) ); // 2:88.859M-89.313M
          sn.push( 13, ho_interval( 105745000, 106351000 ) ); // 14:105.745M-106.351M
               }
     else 
     {    
          // via https://genome.ucsc.edu/cgi-bin/hgLiftOver
          // could not translate chr14 interval, however the interval given below
          // maps to inside the above interval on chr14; it could probably be
          // extended with a little more work

          sn.push( 1,  ho_interval(  89158513,  89612757 ) );
          sn.push( 13, ho_interval( 106215480, 106506183 ) );    }    }

void FindBadEdges( const HyperBasevectorX& hb, const vec<int>& inv,
     const vec<vec<vec<vec<int>>>>& lines, const vec< vec< pair<int,int> > >& hits,
     vec<Bool>& bad, const Bool alt = False );

#endif
