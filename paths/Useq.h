///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef USEQ_H
#define USEQ_H

#include "CoreTools.h"
#include "Basevector.h"

class useq {

     public:

     useq( ) { }
     useq( const vec<int>& u ) : u_(u) { }

     const vecbasevector& Unibases( ) const
     {    Assert( unibases_ != 0 );
          return *unibases_;    }
     const basevector& Unibase( int k ) const
     {    Assert( unibases_ != 0 );
          return Unibases( )[k];    }
     int ToRc( int u ) const
     {    Assert( to_rc_ != 0 );
          return (*to_rc_)[u];    }
     void SetUnibases( const int K, const vecbasevector& unibases, 
          const vec<int>& to_rc ) 
     {    K_ = K;
          unibases_ = &unibases;
          to_rc_ = &to_rc;    }

     int N( ) const { return u_.size( ); }
     const vec<int>& U( ) const { return u_; }
     int U( int k ) const { return u_[k]; }
     int K( ) const { return K_; }

     int Len( ) const; // return length in bases
     int Kmers( ) const; // return length in kmers

     void ReverseMe( );

     friend ostream& operator<<( ostream& out, const useq& x )
     {    for ( int j = 0; j < x.N( ); j++ )
          {    if ( j > 0 ) out << " ";
               out << x.U(j);    }
          return out;    }

     private:

     vec<int> u_;       // unibase ids
     static int K_;        
     const static vecbasevector* unibases_;
     const static vec<int>* to_rc_;

};

class ualign {

     public:

     ualign( ) { }

     ualign( const useq& x1, const useq& x2, const int p1, const int p2, 
          const int max_delta );

     pair<int,int> Tie( int j ) const { return ties_[j]; }
     const vec< pair<int,int> >& Ties( ) const { return ties_; }

     void Print( ostream& out, const useq& x1, const useq& x2 ) const;

     friend Bool operator==( const ualign& a1, const ualign& a2 )
     {    return a1.ties_ == a2.ties_;    }

     friend Bool operator<( const ualign& a1, const ualign& a2 )
     {    return a1.ties_ < a2.ties_;    }

     friend Bool operator>( const ualign& a1, const ualign& a2 )
     {    return a1.ties_ > a2.ties_;    }

     private:

     vec< pair<int,int> > ties_;

};

#endif
