/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"

#include "math/Functions.h"
#include "graphics/BarGraph.h"
#include "graphics/HistogramsForDummies.h"

#include <functional>

// Generates a histogram from the given data (a file of integers, one
// per line).
//
// SHIFT: shift data so there are no negative bins

int main( int argc, char *argv[] )
{
  RunTime();

     BeginCommandArguments;
     CommandArgument_String(DATA);
     CommandArgument_Int_OrDefault(DIGITS_TO_RIGHT, 4);
     CommandArgument_Bool_OrDefault( SHIFT, True );
     CommandArgument_String(OUTHEAD);
     EndCommandArguments;

  Ifstream( in, DATA );
  vec<int> x;
  while(1)
  {
    String s;
    in >> s;
    if ( !in ) break;
    x.push_back( int(round( s.Double( ) )) );
  }
  sort( x.begin(), x.end() );

  if ( SHIFT && x.front() < 0 )
  {
    cout << "Had to shift everything by " << -x.front() << "." << endl;
    transform( x.begin(), x.end(),
  	       x.begin(),
  	       bind2nd( plus<int>(), -x.front() ) );
  }

  float cutoff = 0;
  float mean = 0;
  float stdev = 0;

  String cutoff_str( ( argc > 2 ? argv[2] : "" ) );
  if ( cutoff_str.IsInt() )
  {
    cutoff = cutoff_str.Int();
  }

  if ( cutoff > 0 )
  {
    String mean_str( ( argc > 3 ? argv[3] : "" ) );
    String stdev_str( ( argc > 4 ? argv[4] : "" ) );

    if ( mean_str.IsInt() && 
	 stdev_str.IsInt() )
    {
      mean = mean_str.Int();
      stdev = stdev_str.Int();

      cout << "Pretending mean is " << mean << " and stdev is " << stdev << "." << endl;
    }
  }

  vec<ns_psplot::freetext> labels;

  GenerateHistogram( labels, x, OUTHEAD, "", cutoff, mean, stdev, DIGITS_TO_RIGHT );

  // System( "gv " + tempfile + ".eps" );
  // Remove( tempfile + ".eps" );
}
