///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file EditDistance.h
 * \author tsharpe
 * \date Oct 19, 2012
 *
 * \brief
 */
#ifndef EDITDISTANCE_H_
#define EDITDISTANCE_H_

#include <algorithm>
#include <cstddef>

 // iterators must have comparable value_types
template<class Itr, class Itr2>
unsigned editDistance( Itr const& beg, Itr const& end,
                        Itr2 itr2, Itr2 const& end2 )
{
    using std::distance;
    size_t width = distance(beg, end);
    unsigned* rows[2];
    rows[0] = new unsigned[width * 2];
    rows[1] = rows[0] + width;
    unsigned* pRow = rows[0];
    for ( size_t idx = 1; idx <= width; ++idx )
        *pRow++ = idx;
    size_t whichRow = 0;
    unsigned firstVal = 0;
    unsigned prevVal = 0;
    for ( ; itr2 != end2; ++itr2 )
    {
        typename Itr2::value_type test = *itr2;
        unsigned* pLastRow = rows[whichRow];
        pRow = rows[whichRow ^= 1];
        unsigned diagVal = firstVal;
        prevVal = ++firstVal;
        for ( Itr itr(beg); itr != end; ++itr )
        {
            unsigned s = diagVal + (*itr != test);
            unsigned s2 = prevVal + 1;
            unsigned s3 = (diagVal = *pLastRow++) + 1;
            using std::min;
            prevVal = min(s,min(s2,s3));
            *pRow++ = prevVal;
        }
    }
    delete[] rows[0];
    return prevVal;
}

#endif /* EDITDISTANCE_H_ */
