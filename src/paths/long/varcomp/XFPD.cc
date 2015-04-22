///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Get error bars for DISCOVAR het FP rate.

#include "FastIfstream.h"
#include "MainTools.h"
#include "random/Random.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(DATE);
     EndCommandArguments;

     fast_pipe_ifstream in( "cat /wga/scr4/vartests/" 
          + DATE + "/test_xhet*.report | grep false" );
     String line, loc, junk;
     int count;
     vec<int> counts( 100, 0 );
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          istringstream iline( line.c_str( ) );
          iline >> loc >> junk >> count;
          int pos = loc.Before( "-" ).Int( );
          counts[ ( pos - 10000000 ) / 1000000 ] += count;    }
     const int passes = 1000;
     const double base = 553;
     vec<double> fps;
     for ( int pass = 0; pass < passes; pass++ )
     {    int X = 0;
          for ( int j = 0; j < 100; j++ )
               X += counts[ randomx( ) % 100 ];
          double FP = 100.0 * ( X / 100000000.0 ) / ( 1.8 * base / 1000000.0 );
          fps.push_back(FP);    }
     cout << setprecision(2) << StdDev( fps, Mean(fps) ) << endl;    } 
