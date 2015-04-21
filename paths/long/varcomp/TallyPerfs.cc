///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Used by N50_perfect_stretch.

#include "MainTools.h"
#include "math/Functions.h"

int main( )
{
     RunTime( );

     String line;
     vec<int> perfs;
     Bool plus = False;
     while(1)
     {    getline( cin, line );
          if ( cin.fail( ) ) break;
          if ( line == "+" ) plus = True;
          else
          {    int n = line.Int( );
               if (plus)
               {    if ( perfs.nonempty( ) ) perfs.back( ) += n;
                    else perfs.push_back(n);
                    plus = False;    }    
               else perfs.push_back(n);    }    }
     Sort(perfs);

     cout << "N50 perfect stretch = " << N50(perfs) << endl;
     cout << "Q" << ToString( log10(N50(perfs))*10, 1 ) << endl;    }
