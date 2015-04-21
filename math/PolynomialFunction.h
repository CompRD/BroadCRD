/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// This file defines the polynomial_function class and the 
// piecewise_polynomial_function class.
// 
// A polynomial_function represents a function R --> R given by a polynomial.
//
// A piecewise_polynomial_function represents a function [a,b] --> R that is 
// piecewise polynomial.  It is defined by numbers a = p1 < p2 < ... < pn = b and 
// by n-1 polynomials pf.

#ifndef POLYNOMIAL_FUNCTION_H
#define POLYNOMIAL_FUNCTION_H

#include "CoreTools.h"
#include "math/Functions.h"

class polynomial_function {

     public:

     polynomial_function( ) { }

     void Initialize( const vec<double>& coeff )
     {    ForceAssert( coeff.nonempty( ) );
          coeff_ = coeff;    }

     polynomial_function( const vec<double>& coeff )
     {    Initialize(coeff);    }

     Bool Defined( ) const { return coeff_.nonempty( ); }

     Bool Defined( const double x ) const { return Defined( ); }

     double operator( ) ( const double x ) const;

     int Degree( ) const
     {    ForceAssert( coeff_.nonempty( ) );
          return coeff_.size( ) - 1;    }

     // Shift(a) changes f(x) into the new polynomial f(x+a).

     void Shift( const double a );

     // Scale(c) changes f(x) into the new polynomial f(cx).

     void Scale( const double c );

     // Definite and indefinite integrals.

     void Integrate( polynomial_function& I ) const;
     double Integrate( const double a, const double b ) const;

     friend ostream& operator<<( ostream& o, const polynomial_function& f );

     private:

     vec<double> coeff_;

};

class piecewise_polynomial_function {

     public:

     piecewise_polynomial_function( ) { }

     void Initialize( const vec<double>& p, const vec<polynomial_function>& pf );

     piecewise_polynomial_function( 
          const vec<double>& p, const vec<polynomial_function>& pf )
     {    Initialize( p, pf );    }

     Bool Defined( ) const { return p_.nonempty( ); }

     Bool Defined( const double x ) const 
     {    return Defined( ) && x >= p_.front( ) && x <= p_.back( );    }

     void RequireDefined( const double x ) const;

     // Loc(x): for x lying within the domain of the function, return the largest k
     // such that x >= p_k.

     int Loc( const double x ) const;

     double operator( ) ( const double x ) const;

     const polynomial_function& Poly( int k ) const { return pf_[k]; }

     double Juncture( int k ) const { return p_[k]; }

     // Inverse: return x such that (*this)(x) = y, provided that (*this) is a
     // strictly increasing function and y lies in its range.  This could be
     // implemented to work for strictly decreasing functions (i.e. for any
     // one-to-one function).  The implementation is likely very inefficient, as
     // compared to what is possible.  The increasing function requirement is not
     // checked.

     double Inverse( const double y );

     Bool Valid( const double eps );

     // Definite integral.

     double Integrate( double a, double b ) const;

     friend ostream& operator<<(ostream& o, const piecewise_polynomial_function& f);

     private:

     vec<double> p_; 
     vec<polynomial_function> pf_;

};

#endif

// PiecewiseLinear Given points (x1,y1),...,(xn,yn), where x1 < x2 < ... < xn, 
// return the piecewise linear function f defined on [x1,xn] such that yi = f(xi) 
// for all i.

void PiecewiseLinear( const vec<double>& x, const vec<double>& y,
     piecewise_polynomial_function& f );

// NaturalCubicSpline.  Given points (x1,y1),...,(xn,yn), where x1 < x2 < ... < xn, 
// return the piecewise cubic function f defined on [x1,xn] such that yi = f(xi) for 
// all i, such that the pieces agree up to order two at the junctures, and such that 
// the second derivatives are zero at the endpoints.

void NaturalCubicSpline( const vec<double>& x, const vec<double>& y,
     piecewise_polynomial_function& f );

// CubicSolution.  Let v = (a,b,c,d).  Return the real root of the equation
// a + bx + cx^2 + dx^3 = 0, provided that the equation has exactly one real root.
// Otherwise, assert.

double CubicSolution( const vec<double>& v );
