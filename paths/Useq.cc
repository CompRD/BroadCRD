///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "paths/Useq.h"

int useq::Kmers( ) const
{    Assert( unibases_ != 0 );
     int kmers = 0;
     for ( int j = 0; j < N( ); j++ )
          kmers += Unibase( U(j) ).size( ) - ( K( ) - 1 );
     return kmers;    }

int useq::Len( ) const
{    return Kmers( ) + K( ) - 1;    }

void useq::ReverseMe( )
{    u_.ReverseMe( );
     for ( int j = 0; j < N( ); j++ )
          u_[j] = (*to_rc_)[ u_[j] ];    }

const vecbasevector* useq::unibases_(0);
const vec<int>* useq::to_rc_(0);
int useq::K_(0);

void ExtendAlignmentRight( const useq& x1, const useq& x2, const int max_delta, 
     vec< pair<int,int> >& a )
{    int p1 = a.back( ).first, p2 = a.back( ).second;
     if ( p1 == x1.N( ) - 1 || p2 == x2.N( ) - 1 ) return;
     const vecbasevector& unibases = x1.Unibases( );
     int K = x1.K( ), n1 = 0, n2 = 0;
     for ( int j2 = p2 + 1; j2 < x2.N( ); j2++ )
     {    for ( int j1 = p1 + 1; j1 < x1.N( ); j1++ )
          {    if ( x1.U(j1) == x2.U(j2) )
               {    if ( Abs(n1-n2) <= max_delta )
                    {    a.push( j1, j2 );
                         ExtendAlignmentRight( x1, x2, max_delta, a );    
                         return;    }    }
               n1 += unibases[ x1.U(j1) ].isize( ) - (K-1);    }
          n2 += unibases[ x2.U(j2) ].isize( ) - (K-1);    }    }

void ExtendAlignmentLeft( const useq& x1, const useq& x2, const int max_delta, 
     vec< pair<int,int> >& a )
{    int p1 = a.back( ).first, p2 = a.back( ).second;
     if ( p1 == 0 || p2 == 0 ) return;
     const vecbasevector& unibases = x1.Unibases( );
     int K = x1.K( ), n1 = 0, n2 = 0;
     for ( int j2 = p2 - 1; j2 >= 0; j2-- )
     {    for ( int j1 = p1 - 1; j1 >= 0; j1-- )
          {    if ( x1.U(j1) == x2.U(j2) )
               {    if ( Abs(n1-n2) <= max_delta )
                    {    a.push_back( make_pair( j1, j2 ) );
                         ExtendAlignmentLeft( x1, x2, max_delta, a );    
                         return;    }    }
               n1 += unibases[ x1.U(j1) ].isize( ) - (K-1);    }
          n2 += unibases[ x2.U(j2) ].isize( ) - (K-1);    }    }

ualign::ualign( const useq& x1, const useq& x2, const int p1, const int p2, 
     const int max_delta )
{    ties_.push( p1, p2 );
     ExtendAlignmentRight( x1, x2, max_delta, ties_ );
     ties_.ReverseMe( );
     ExtendAlignmentLeft( x1, x2, max_delta, ties_ );
     ties_.ReverseMe( );    }

void ualign::Print( ostream& out, const useq& x1, const useq& x2 ) const
{    int w = 0, M = x1.Unibases( ).size( ) - 1;
     while( M > 0 )
     {    M /= 10;
          w++;    }
     const int s = 2;
     vec< vec<int> > fields(2);
     for ( int j = 0; j < ties_.front( ).first; j++ )
          fields[0].push_back( x1.U(j) );
     for ( int j = 0; j < ties_.front( ).second; j++ )
          fields[1].push_back( x2.U(j) );
     while( fields[0].size( ) < fields[1].size( ) )
          fields[0].push_front(-1);
     while( fields[1].size( ) < fields[0].size( ) )
          fields[1].push_front(-1);
     for ( int r = 0; r < ties_.isize( ); r++ )
     {    fields[0].push_back( x1.U( ties_[r].first ) );
          fields[1].push_back( x2.U( ties_[r].second ) );
          if ( r < ties_.isize( ) - 1 )
          {    for ( int u = ties_[r].first + 1; u < ties_[r+1].first; u++ )
                    fields[0].push_back( x1.U(u) );
               for ( int u = ties_[r].second + 1; u < ties_[r+1].second; u++ )
                    fields[1].push_back( x2.U(u) );    }
          while( fields[0].size( ) < fields[1].size( ) )
               fields[0].push_back(-1);
          while( fields[1].size( ) < fields[0].size( ) )
               fields[1].push_back(-1);    }
     for ( int u = ties_.back( ).first + 1; u < x1.N( ); u++ )
     {    fields[0].push_back( x1.U(u) );
          fields[1].push_back(-1);    }
     for ( int u = ties_.back( ).second + 1; u < x2.N( ); u++ )
     {    fields[1].push_back( x2.U(u) );
          fields[0].push_back(-1);    }
     for ( int l = 0; l < 2; l++ )
     {    if ( l == 1 )
          {    for ( int j = 0; j < fields[l].isize( ); j++ )
               {    if ( j > 0 ) out << String( s, ' ' );
                    if ( fields[0][j] >= 0 
                         && fields[0][j] == fields[1][j] )
                    {    out << String( w, '|' );    }
                    else out << String( w, ' ' );    }
               out << "\n";    }
          for ( int j = 0; j < fields[l].isize( ); j++ )
          {    if ( j > 0 ) out << String( s, ' ' );
               if ( fields[l][j] == -1 ) out << String( w, ' ' );
               else
               {    String y = ToString( fields[l][j] );
                    while ( y.isize( ) < w ) y = "0" + y;
                    out << y;    }    }
          out << "\n";    }    }

