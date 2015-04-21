// Copyright (c) 2003 Broad Institute/Massachusetts Institute of Technology
//

#include "String.h"
#include "STLExtensions.h"
#include "system/System.h"
#include "select/ParseFasta.h"



/*
 * ParseFasta
 */
void ParseFasta( const String &in_file,
		 const vec<String> &names,
		 int &n_found,
		 int &n_parsed,
		 ostream &out )
{
  ForceAssert( is_sorted( names.begin( ), names.end( ) ) );

  n_found = 0;
  n_parsed = 0;

  istream *in = 0;
  procbuf *pipe = 0;
  if ( in_file.Contains( ".gz" ) ) {
    string command = "gzip -dc " + in_file;
    pipe = new procbuf( command.c_str(), ios::in );
    in = new istream( pipe );
  }
  else
    in = new ifstream( in_file.c_str() );
  
  String line;
  while( 1 ) {
    getline( *in, line );
    if ( in->fail( ) ) break;
    if ( line.size( ) == 0 || line[0] != '>' ) continue;
    n_parsed++;
    String read_name = line.After( ">" );
    while ( read_name.Contains( " " ) )
      read_name = read_name.After( " " );

    if ( binary_search( names.begin( ), names.end( ), read_name ) ) {
      char c;
      n_found++;
      out << line << "\n";
      while( 1 ) {
	if ( in->peek( ) == '>' ) break;
	if ( in->fail( ) ) break;
	getline( *in, line );
	if ( in->fail( ) ) break;
	out << line << "\n";
      }
    }
  }

  if ( in )
    delete ( in );
  if ( pipe )
    delete ( pipe );
}
