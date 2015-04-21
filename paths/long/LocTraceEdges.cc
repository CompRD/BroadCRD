///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * AlignToSHBV.cc
 *
 *  Created on: Sep 3, 2013
 *      Author: blau
 */
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
#include "paths/long/LocTraceEdges.h"
#include <omp.h>
#include <iomanip>
#include "MainTools.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "paths/long/EvalByReads.h"
#include <queue>

vec<vec<read_place > > LocTraceEdges( const HyperBasevector& hb_fw, const vecbasevector& bases, const vecqualvector& quals, const double delta_qsum, const int64_t L, bool debug /*=false*/) {
    const int64_t minQual=3;
    const int infinity = 1000000000;
    const double prox=delta_qsum*1000;
    HyperBasevector hb_rc(hb_fw);
    hb_rc.Reverse( );
    vec<int> to_left_fw, to_left_rc;
    hb_fw.ToLeft(to_left_fw), hb_rc.ToLeft(to_left_rc);
    vec<int> to_right_fw, to_right_rc;
    hb_fw.ToRight(to_right_fw), hb_rc.ToRight(to_right_rc);

    const int64_t nReads=bases.size();

    vecbasevector x_fw, x_rc;
    for ( int i = 0; i < hb_fw.EdgeObjectCount( ); i++ ) x_fw.push_back( hb_fw.EdgeObject(i) );
    for ( int i = 0; i < hb_rc.EdgeObjectCount( ); i++ ) x_rc.push_back( hb_rc.EdgeObject(i) );

    vec<vec<read_place>> places_log(nReads);
    VecIntPairVec locs_fw, locs_rc;
    CreateGlocs( x_fw, L, locs_fw );
    CreateGlocs( x_rc, L, locs_rc );
    #pragma omp parallel for schedule(dynamic)
    for ( size_t id = 0; id < bases.size( ) ; id++ ){
        int n = KmerId( bases[id], L, 0 );
        int qual_sum = infinity;
        FindPlaces( bases[id], quals[id], n, hb_fw, hb_rc, to_right_fw, to_right_rc, locs_fw, locs_rc, places_log[id], qual_sum, minQual, prox);
        Sort(places_log[id],[](const read_place&L,const read_place&R){return L.Qsum() < R.Qsum();});
    }

    size_t r_digits=1;
    for(size_t tmp=nReads/10; tmp>0 ; ++r_digits,tmp/=10){};
    r_digits=max(size_t(6),r_digits);
    if ( debug ) {
	std::cout << "# SHBV.K()   L delta_qsum nQueries" << "\n# "
		<< std::setw(8) << hb_fw.K()<< " "
		<< std::setw(3) << L << " "
		<< std::setw(10) << delta_qsum << " "
		<< nReads << "\n"
		<< "#query place/nPlaces fw/rc-edges 1st_edge.offset,2nd_edge,3rd_edge... [q-sum of mismatches]" << std::endl;
	for ( size_t id = 0; id < bases.size( ); id++ ){
	    const auto& places=places_log[id];
	    const size_t np = places.size( );
	    if( np==0)                       { std::cout << std::setw(r_digits) << id+1 << " " << 0 << "/" << np <<"\n"; }
	    for ( size_t i = 0; i < np; i++ ){ std::cout << std::setw(r_digits) << id+1 << " " << i+1 << "/" << np << "\t" << places[i] << "\n";}
	}
    }

    return places_log;
}
