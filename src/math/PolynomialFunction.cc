/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "math/Matrix.h"
#include "math/PolynomialFunction.h"

void polynomial_function::Shift( const double a )
{    ForceAssert( Defined( ) );
     vec<double>& c = coeff_;
     int n = c.size( );
     vec<double> cnew( n, 0 ), pow( n, 0 );
     pow[0] = 1;
     for ( int i = 0; i < n; i++ )
     {    for ( int j = 0; j <= i; j++ )
               cnew[j] += c[i] * pow[j];
          if ( i < n-1 )
          {    pow[i+1] = 1;
               for ( int j = i; j >= 1; j-- )
                    pow[j] = a * pow[j] + pow[j-1];
               pow[0] *= a;    }    }
     c = cnew;    }

void polynomial_function::Scale( const double c )
{    ForceAssert( Defined( ) );
     int n = coeff_.size( );
     for ( int i = 0; i < n - 1; i++ )
     {    for ( int j = i + 1; j < n; j++ )
               coeff_[j] *= c;    }    }

// The implementation of NaturalCubicSpline is derived from the math given in
// J. Stoer and R. Bulirsch, Introduction to Numerical Analysis, Second Edition,
// Springer-Verlag, 1993.

void NaturalCubicSpline( const vec<double>& x, const vec<double>& y,
     piecewise_polynomial_function& f )
{    ForceAssertGe( x.size( ), 2u );
     ForceAssertEq( x.size( ), y.size( ) );
     ForceAssert( x.UniqueOrdered( ) );
     int n = x.size( );
     vec<double> h(n), z(n), l(n-1), u(n), d(n), M(n);
     for ( int i = 1; i < n; i++ )
     {    h[i] = x[i] - x[i-1];
          z[i] = y[i] - y[i-1];    }
     for ( int i = 1; i < n-1; i++ )
     {    l[i] = h[i+1] / ( h[i] + h[i+1] );
          u[i] = 1.0 - l[i];    }
     l[0] = u[n-1] = 0.0;
     for ( int i = 1; i < n-1; i++ )
          d[i] = ( 6.0 / (h[i] + h[i+1]) ) * ( z[i+1] / h[i+1] - z[i] / h[i] );
     d[0] = d[n-1] = 0.0;
     matrix<double> T( n, n, 0 );
     for ( int i = 0; i < n; i++ )
          T(i,i) = 2.0;
     for ( int i = 0; i < n-1; i++ )
     {    T(i,i+1) = l[i];
          T(i+1,i) = u[i+1];    }
     Solve( T, d, M );
     vec<polynomial_function> pf(n-1);
     for ( int i = 0; i < n-1; i++ )
     {    vec<double> coeff(4);
          coeff[0] = y[i];
          coeff[1] = -M[i] * h[i+1] / 2.0
               + z[i+1] / h[i+1] - h[i+1]/6.0 * ( M[i+1] - M[i] );
          coeff[2] = M[i]/2.0;
          coeff[3] = ( M[i+1] - M[i] ) / ( 6.0 * h[i+1] );
          pf[i].Initialize(coeff);
          pf[i].Shift( -x[i] );    }
     f.Initialize( x, pf );    }

void PiecewiseLinear( const vec<double>& x, const vec<double>& y,
     piecewise_polynomial_function& f )
{    ForceAssertGe( x.size( ), 2u );
     ForceAssertEq( x.size( ), y.size( ) );
     ForceAssert( x.UniqueOrdered( ) );
     int n = x.size( );
     vec<polynomial_function> pf(n-1);
     for ( int i = 0; i < n-1; i++ )
     {    vec<double> coeff(2);
          coeff[1] = (y[i+1] - y[i]) / (x[i+1] - x[i]);
          coeff[0] = y[i] - coeff[1] * x[i];
          pf[i].Initialize(coeff);    }
     f.Initialize( x, pf );    }

double CubicSolution( const vec<double>& v )
{    double a = v[0], b = v[1], c = v[2], d = v[3];
     ForceAssert( d != 0.0 );
     double u = -2.0*c*c*c + 9.0*b*c*d - 27.0*a*d*d;
     double q = -c*c + 3.0*b*d;
     double q3 = q*q*q;
     double D = 4.0*q3 + u*u;
     ForceAssert( D > 0.0 );
     double cr2 = pow( 2.0, 1.0/3.0 );
     double x = pow( u + sqrt(D), 1.0/3.0 );
     return -c/(3.0*d) - cr2 * q / ( 3.0 * d * x ) + x / ( 3.0 * cr2 * d );    }

double polynomial_function::operator( ) ( const double x ) const
{    ForceAssert( Defined( ) );
     double y = 0.0;
     for ( int i = coeff_.isize( ) - 1; i >= 0; i-- )
     {    y *= x;
          y += coeff_[i];    }
     return y;    }

ostream& operator<<( ostream& o, const polynomial_function& f )
{    Bool first_print = True;
     for ( int i = f.coeff_.isize( ) - 1; i >= 0; i-- )
     {    double c = f.coeff_[i];
          if ( c == 0.0 ) continue;
          if ( !first_print )
          {    if ( c > 0.0 ) o << " + ";
               else 
               {    o << " - ";
                    c = -c;    }    }
          first_print = False;
          if ( i == 0 ) o << c;
          else 
          {    if ( c != 1.0 ) o << c << " ";
               o << "x";
               if ( i != 1 ) o << "^" << i;    }    }
     if (first_print) o << 0;
     return o;    }

void piecewise_polynomial_function::Initialize( 
     const vec<double>& p, const vec<polynomial_function>& pf )
{    p_ = p;
     pf_ = pf;
     ForceAssertGe( p.size( ), 2u );
     ForceAssertEq( p.size( ) - 1, pf.size( ) );
     for ( int i = 1; i < p.isize( ); i++ )
          ForceAssert( p[i] > p[i-1] );    }

void piecewise_polynomial_function::RequireDefined( const double x ) const
{    if ( !Defined(x) )
     {    cout << "\npiecewise_polynomial_function called with arg "
               << x << ",\nbut arg must lie between " << p_.front( )
               << " and " << p_.back( ) << "\n";    }
     ForceAssert( Defined(x) );    }

int piecewise_polynomial_function::Loc( const double x ) const
{    RequireDefined(x);
     int k = upper_bound( p_.begin( ), p_.end( ), x ) - p_.begin( ) - 1;
     if ( k == p_.isize( ) - 1 ) --k;
     return k;    }

double piecewise_polynomial_function::operator( ) ( const double x ) const
{    return pf_[ Loc(x) ](x);    }

double piecewise_polynomial_function::Integrate( double a, double b ) const
{    Bool negate = False;
     if ( a == b ) return 0.0;
     if ( b < a )
     {    swap( a, b );
          negate = True;    }
     int ia = Loc(a), ib = Loc(b);
     double sum = 0.0;
     if ( ia == ib ) sum = Poly(ia).Integrate( a, b );
     else
     {    sum += Poly(ia).Integrate( a, Juncture(ia+1) );
          ++ia;
          if ( ib != p_.isize( ) - 1 && b != Juncture(ib+1) )
          {    sum -= Poly(ib).Integrate( b, Juncture(ib+1) );
               ++ib;    }
          for ( int j = ia; j < ib; j++ )
               sum += Poly(j).Integrate( Juncture(j), Juncture(j+1) );    }
     if (negate) sum = -sum;
     return sum;    }

double piecewise_polynomial_function::Inverse( const double y )
{    int n = p_.size( );
     const piecewise_polynomial_function& f = *this;
     double L = p_[0], R = p_[n-1];
     ForceAssert( f(L) <= y && y <= f(R) );
     double last_diff = R - L + 1.0;
     while(1)
     {    double diff = R - L;
          if ( diff < 0 || diff >= last_diff ) return L;
          last_diff = diff;
          double M = ( L + R ) / 2.0;
          if ( y == f(M) ) return M;
          if ( y < f(M) ) R = M;
          else L = M;    }    }

Bool piecewise_polynomial_function::Valid( const double eps )
{    ForceAssert( Defined( ) );
     for ( int i = 1; i < p_.isize( ) - 1; i++ )
          if ( Abs( pf_[i-1]( p_[i] ) - pf_[i]( p_[i] ) ) > eps ) return False;
     return True;    }

ostream& operator<<( ostream& o, const piecewise_polynomial_function& f )
{    for ( int i = 0; i < f.pf_.isize( ); i++ )
          o << "[" << f.p_[i] << "," << f.p_[i+1] << "]: " << f.pf_[i] << "\n";
     return o;    }

void polynomial_function::Integrate( polynomial_function& I ) const
{    I.coeff_.resize( coeff_.size( ) + 1 ); 
     I.coeff_[0] = 0.0; // arbitrary choice
     for ( int n = 0; n < coeff_.isize( ); n++ )
          I.coeff_[n+1] = coeff_[n] / double(n+1);    }    

double polynomial_function::Integrate( const double a, const double b ) const
{    polynomial_function I;
     Integrate(I);
     return I(b) - I(a);    }

