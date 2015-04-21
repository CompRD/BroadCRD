/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


/** Make sure that a tab (or other) separated table looks pretty.

\file SameSizeColumns.cc

This can function as a filter (from stdin to stdout)
if no IN or OUT parameters are given.

By default, it will work in place and modify the IN file.

SEP (default tab) is the input table separator, 
and OUTPUT_SEP (default two spaces)
is the output table separator.
*/

#include "MainTools.h"
#include <sstream>

void ProcessRows(ostream & os, const vec<vec<String> > & rows,
		 const String & OUTPUT_SEP, bool leftJustify);

int main( int argc, char *argv[] ) {

  RunTime( );

  BeginCommandArgumentsAcceptEmptyArgList;
  CommandArgument_String_OrDefault(SEP, "\t");
  CommandArgument_String_OrDefault(IN, "");
  CommandArgument_String_OrDefault(JUSTIFY, "L");
  CommandArgument_String_OrDefault(OUT, IN);
  CommandArgument_String_OrDefault(OUTPUT_SEP, "  ");
  EndCommandArguments;

  ostringstream os;

  istream * is = 0;
  if (IN.empty()) { is = &cin; }
  else is = new ifstream(IN.c_str());

  bool leftJustify = ('l' == tolower(JUSTIFY[0])); 

  vec< vec<String> > rows;
  String line;
  while(1) {    
    getline( *is, line );
    if ( !(*is) ) break;
    if (line.empty()) { 
      ProcessRows(os, rows, OUTPUT_SEP, leftJustify);
      rows.clear();
      os << "\n";
    }
    else {
      rows.push_back(vec<String>());
      Tokenize(line, rows.back(), SEP);
    }
  }
  if (!IN.empty()) delete is; //cleanup and close it.
  ProcessRows(os, rows, OUTPUT_SEP, leftJustify);

  ostream * out = 0;
  if (OUT.empty()) { out = &cout; }
  else out = new ofstream(OUT.c_str());
     
  (*out) << os.str();
  if (!OUT.empty()) delete out; //cleanup and close it.

  return 0;
}

void ProcessRows(ostream & os, const vec<vec<String> > & rows,
		 const String & OUTPUT_SEP, bool leftJustify) {
  vec<int> sizes;
  int csize=0;
  for (int col=0; true; ++col) {
    csize = -1;
    for (int i=0; i != rows.isize(); ++i) {
      if (rows[i].isize() <= col) continue;
      csize = max(csize, rows[i][col].isize());
    }
    if (-1 == csize) break;
    else sizes.push_back(csize);
  }

  for (int i=0; i != rows.isize(); ++i) {
    for (int col = 0; col != rows[i].isize(); ++col) {
      vec<char> spaces(sizes[col] - rows[i][col].size(), ' ');
      if (leftJustify) {
	os << rows[i][col] << String(spaces.begin(),spaces.end()) << OUTPUT_SEP;
      }
      else {
	os << String(spaces.begin(),spaces.end()) << rows[i][col] << OUTPUT_SEP;
      }
    }
    os << "\n";
  }
}
    
  


