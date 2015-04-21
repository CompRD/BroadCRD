///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file TestMultiVec.cc
 * \author neilw
 * \date October 10, 2013
 *
 * \brief
 */
#include "MainTools.h"
#include "MultiVec.h"
#include <iostream>
#include <fstream>
#include "Vec.h"


int main( int argc, char* argv[] )
{
	vec<size_t> a{1,2,3};
	vec<size_t> b{8,9,10,11};
	vec<size_t> c{25,26,27,28,29};

        std::vector< vec<size_t> > collection = {a,c,b};

	MultiVec< vec<size_t> > mv(collection);

	cout << "size=" << mv.size() << endl;
	for ( size_t i = 0; i < mv.size(); ++i )
		cout << " mv[" << i << "]=" << mv[i] << endl;

        mv[6]=1;
        mv[7]=2;

        cout << "---- changed 6->1 and 7->2" << endl;
	for ( size_t i = 0; i < mv.size(); ++i )
		cout << " mv[" << i << "]=" << mv[i] << endl;

        cout << "front is " << mv.front() << endl;
        cout << "back is " << mv.back() << endl;

//        cout << "about to do something bad..." << endl;
//        cout << mv[mv.size()] << endl;

        return 0;
}
