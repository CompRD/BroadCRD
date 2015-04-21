///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef DATA_ORIENTATION
#define DATA_ORIENTATION

#include <iostream>

using std::ostream;

enum Orientation {
  orient_FW,
  orient_RC
};

inline
Orientation Flip( const Orientation anOrientation )
{ return ( anOrientation == orient_FW ? orient_RC : orient_FW ); }

inline
Orientation operator+ ( const Orientation lhs, const Orientation rhs )
{
  if ( lhs == rhs ) 
    return orient_FW;
  else
    return orient_RC;
}

inline 
ostream & operator<< ( ostream & out, const Orientation anOrientation )
{
    return out << ( anOrientation == orient_FW ? "+" : "-" );
}
    
#endif
