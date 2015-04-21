/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// ComputeUniqueMolecules.  Estimate number of unique molecules in a library.

#include "MainTools.h"

double f( double x, double c, double n )
{    return c/x - 1 + exp(-n/x);    }

int main( int argc, char *argv[] )  
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_Double_Doc(n, "number of aligned pairs");
     CommandArgument_Double_Doc(c, "number of distinct aligned pairs");
     CommandArgument_Bool_OrDefault(QUIET, False);
     EndCommandArguments;

     if ( c >= n )
     {    cout << "You must have n > c." << endl;
          exit(1);     }

     if ( !QUIET )
     {    cout << (longlong)roundl(c) << " distinct pairs" << endl;
          cout << "fraction of unique pairs = " << PERCENT_RATIO( 3, c, n ) 
               << endl;    }

     double m = 1.0, M = 100.0;

     ForceAssert( f(m*c,c,n) > 0 );
     while( f(M*c,c,n) >= 0 ) M *= 10.0;

     for ( int i = 0; i < 40; i++ )
     {    double r = (m+M)/2.0;
          double u = f( r * c, c, n );
          if ( u == 0 ) break;
          else if ( u > 0 ) m = r;
          else if ( u < 0 ) M = r;    }
     longlong x = (longlong)( roundl( c * (m+M)/2.0 ) );
     if ( !QUIET ) cout << "estimated library size = ";
     cout << x << endl;
     if (QUIET) exit(0);

     cout << "\nSuppose that we increased the amount of sequence, by some amount,\n"
          << "say z-fold.  (So z = 1 means we don't add any, z = 2 means we\n" 
          << "double, etc.)  Here we predict the number of distinct pairs that\n"
          << "would be obtained, as compared the number that we currently have, c.";
     cout << "\n\n" << "estimate of (distinct pairs in z times given reads) / "
          << "(distinct pairs in given reads)\n";
     for ( int z = 1; z <= 10; z++ )
     {    cout << "z = " << z << ": " 
               << setprecision(3) << x * ( 1 - exp(-(z*n)/x) ) / c << endl;    }    }
