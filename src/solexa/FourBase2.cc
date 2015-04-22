/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "solexa/FourBase2.h"

int four_base2::Call( ) const
{    if ( A( ) >= C( ) && A( ) >= G( ) && A( ) >= T( ) ) return 0;
     if ( C( ) >= G( ) && C( ) >= T( ) ) return 1;
     if ( G( ) >= T( ) ) return 2;
     return 3;    }

float four_base2::CallQuality( ) const
{    int call = Call( );
     float call_value = base_[call], next_base = 0.0;
     if ( call_value == 0 ) return 1;
     for ( int i = 0; i < 4; i++ )
     {    if ( i == call ) continue;
          next_base = Max( next_base, base_[i] );    }
     if ( next_base == 0 ) return 1000000000;
     return call_value / next_base;    }

ostream& operator<<( ostream& out, const four_base2& b )
{    return out << "[A=" << b.A( ) << ",C=" << b.C( ) << ",G=" << b.G( )
          << ",T=" << b.T( ) << "]";    }

void Call( const VecFourBase2Vec& I, vecbasevector& bases )
{    longlong nbases = 0;
     int nseqs = I.size( );
     for ( size_t i = 0; i < I.size( ); i++ )
          nbases += I[i].size( );
     bases.clear( );
     bases.Reserve( nbases/16 + nseqs, nseqs );
     for ( size_t i = 0; i < I.size( ); i++ )
     {    static basevector b;
          b.resize( I[i].size( ) );
          for ( int j = 0; j < b.isize( ); j++ )
               b.Set( j, I[i][j].Call( ) );
          bases.push_back(b);    }    }


#include "feudal/SmallVecDefs.h"
#include "feudal/OuterVecDefs.h"
template class SmallVec< four_base2, MempoolAllocator<four_base2> >;
template class OuterVec<FourBase2Vec>;
