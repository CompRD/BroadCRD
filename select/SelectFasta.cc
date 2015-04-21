///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// SelectFasta: select some of the entries from a fasta file, of bases or
// quality scores.

#include <algorithm>
#include <fstream>
#include "FastIfstream.h"
#include "MainTools.h"
#include "ParseSet.h"

int main( int argc, char *argv[] )
{
     RunTime();

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(INPUT);
     CommandArgument_String(OUTPUT);
     CommandArgument_String(IDS);
     CommandArgument_Bool_OrDefault_Doc(UPCASE,False,"Optionally set bases to upper case.");
     EndCommandArguments;

     ofstream os(OUTPUT.c_str());
     vec<int> ids;
     ParseIntSet( IDS, ids );
     vec<Bool> used( ids.size( ), False );

     fast_ifstream in( PRE + "/" + INPUT );
     int count = 0;
     String line;
     while( Sum(used) < used.isize( ) )
     {    if ( in.fail( ) ) break;
          getline( in, line );
          if ( in.fail( ) ) break;
          if ( line.size( ) == 0 || line[0] != '>' ) continue;
          int p = BinPosition( ids, count );
          if ( p >= 0 )
          {    used[p] = True;
               os << line << "\n";
               char c;
               while(1)
               {    in.peek(c);
                    if ( c == '>' ) break;
                    if ( in.fail( ) ) break;
                    getline( in, line );
                    if ( in.fail( ) ) break;
                    if ( UPCASE ) line.ToUpper();
                    os << line << "\n";    }    }
          ++count;    }
    os.close(); }
