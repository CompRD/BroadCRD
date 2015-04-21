///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// ValidateEfasta: check to see if a given file is in efasta format.

#include "FastIfstream.h"
#include "MainTools.h"
#include "TokenizeString.h"
#include "efasta/EfastaTools.h"

#define Err(message)                                      \
{    cout << message << endl << "\nInvalid.\n" << endl;   \
     return 1;    }

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(F);
     EndCommandArguments;

     fast_ifstream in(F);
     String line;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) 
          {    if ( line.size( ) == 0 )
                    Err( "Illegal empty file." );
               break;    }
          if ( line.size( ) == 0 ) Err( "Encountered line of length 0." );
          if ( line[0] != '>' )
          {    Err( "See line = '" << line << "', which was expected "
                    << "to start with >." );    }
          for ( int i = 0; i < line.isize( ); i++ )
          {    char x = line[i];
               if ( x < 33 && x > 126 && x != ' ' && x != '\t' )
               {    Err( "See illegal character '" << x << "' in line = '"
                         << line << "'." );    }    }
          if ( line.size( ) < 2 || line[1] == ' ' || line[1] == '\t' )
               Err( "Illegal header line '" << line << "'." );
          vec<String> lines;
          Bool eof = False;
          while(1)
          {    char c;
               in.peek(c);
               if ( in.fail( ) ) { eof = True; break; }
               if ( c == '>' ) break;
               getline( in, line );
               lines.push_back(line);    }
          String all;
          int64_t all_size = 0;
          for ( size_t i = 0; i < lines.size( ); i++ )
               all_size += lines[i].size( );
          all.reserve(all_size);
          for ( size_t i = 0; i < lines.size( ); i++ )
               all.append( lines[i] );
          ValidateEfastaRecord(lines);
          vec<char> seps;
          seps.push_back( ',' );
          for ( size_t i = 0; i < all.size( ); i++ )
          {    while( i < all.size( ) && all[i] != '{' ) i++;
               if ( i == all.size( ) ) break;
               size_t j;
               for ( j = i+1; j < all.size( ); j++ )
                    if ( all[j] == '}' ) break;
               String stuff;
               for ( size_t l = i + 1; l < j; l++ )
                    stuff.push_back( all[l] );
               vec<String> parts;
               TokenizeStrictly( stuff, seps, parts );
               int n = parts.size( );
               UniqueSort(parts);
               if ( parts.isize( ) != n )
                    Err( "Duplicated item in {" << stuff << "}" );
               i = j - 1;    }
          if (eof) break;    }
     cout << "Valid.\n" << endl;    }
