/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "PrettyPrintTable.h"
#include "String.h"
#include "Vec.h"
#include "math/Functions.h"
#include "reporting/PrintLengthsNstats.h"

/**
 * PrintLengthsNstats
 */
void PrintLengthsNstats( const vec<int> &lens,
			 const String &name,
			 ostream &out )
{
  // Make a local reverse-sorted copy.
  vec<int> rslens = lens;
  ReverseSort( rslens );

  // No object found.
  int n_lens = rslens.size( );
  if ( n_lens == 0 ) {
    out << name << "   NONE\n"<< endl;
    return;
  }
  
  // Indexes.
  longlong total_len = 0;
  for (int ii=0; ii<n_lens; ii++)
    total_len += rslens[ii];

  longlong sum = 0;
  int N25index = 0;
  while ( sum < total_len / 4 ) {
    sum += rslens[N25index];
    N25index++;
  }
  int N50index = N25index;
  while ( sum < total_len / 2 ) {
    sum += rslens[N50index];
    N50index++;
  }
  int N75index = N50index;
  while ( sum < 3 * total_len / 4 ) {
    sum += rslens[N75index];
    N75index++;
  }
  int N98index = N75index;
  while ( sum < ( 98 * total_len ) / 100 ) {
    sum += rslens[N98index];
    N98index++;
  }

  // Prepare table.
  vec< vec<String> > table;
  
  vec<String> line;
  for (int ii=0; ii<5; ii++) {
    String str_N = "";
    String str_largeness = "";
    int index = -1;
    if ( ii < 4 ) {
      switch ( ii ) {
      case 0 :
	str_N = "N98";
	str_largeness = "large";
	index = N98index;
	break;
      case 1:
	str_N = "N75";
	str_largeness = "Large";
	index = N75index;
	break;
      case 2:
	str_N = "N50";
	str_largeness = "LARGE";
	index = N50index;
	break;
      case 3:
	str_N = "N25";
	str_largeness = "LARGEST";
	index = N25index;
	break;
      }
    }

    line.clear( );

    if ( name != "" )
      line.push_back( name );

    if ( ii < 4 ) {
      line.push_back( str_N );
      line.push_back( ToString( rslens[index-1] ) );
      line.push_back( ToString( index ) );
      line.push_back( str_largeness );
    }
    else {
      line.push_back( "total" );
      line.push_back( ToString( total_len ) );
      line.push_back( ToString( n_lens ) );
      line.push_back( "TOTAL" );
    }
    table.push_back( line );
  }

  vec<Bool> rjustify( table[0].size( ), False );
  int beg = name == "" ? 0 : 1;
  int end = name == "" ? 3 : 4;
  for (int ii=beg; ii<end; ii++)
    rjustify[ii] = True;
  
  BeautifyTable( table, &rjustify );

  // Send to stream.
  for (int ii=0; ii<(int)table.size( ); ii++) {
    for (int jj=0; jj<(int)table[ii].size( ); jj++) {
      String space = "";
      if ( jj < 3 ) space =  "     ";
      else if ( jj == 3 ) space = " ";
      
      out << table[ii][jj] << space;
    }
    out << "\n";
  }
  
}
