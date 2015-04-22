///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PAIR_H
#define PAIR_H

/**

\file Pair.h

Functions to help with managing std::pair class

*/


#include <utility>
#include <iostream>

template <class F, class S>
inline std::ostream& operator<<( std::ostream& os, std::pair<F,S> const& pr )
{ return os << '(' << pr.first << ',' << pr.second << ')'; }

#endif //PAIR_H
