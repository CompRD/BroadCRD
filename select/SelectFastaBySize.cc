// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

// SelectFastaBySize: select records from a fasta file having size between
// MIN and MAX bases.

#include "MainTools.h"
#include "FastIfstream.h"

const unsigned int UNDEFINED = 4000000000u;

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(FILE);
     CommandArgument_UnsignedInt_OrDefault(MIN, 0);
     CommandArgument_UnsignedInt_OrDefault(MAX, UNDEFINED);
     EndCommandArguments;

     fast_ifstream in(FILE);
     String line;
     vec<String> lines(1000);
     while(1)
     {    int linecount = 0;
          if ( in.fail( ) ) break;
          getline( in , line );
          if ( in.fail( ) ) break;
          if ( linecount == (int) lines.size( ) ) lines.resize( linecount * 2 );
          lines[linecount++] = line;
          while(1)
          {    char c;
               in.peek(c);
               if ( in.fail( ) || c == '>' ) break;
               getline( in, line );
               if ( linecount == (int) lines.size( ) ) lines.resize( linecount * 2 );
               lines[linecount++] = line;    }
          unsigned int basecount = 0;
          for ( int i = 1; i < linecount; i++ )
               basecount += lines[i].size( );
          if ( basecount < MIN ) continue;
          if ( MAX != UNDEFINED && basecount > MAX ) continue;
          for ( int i = 0; i < linecount; i++ )
               cout << lines[i] << "\n";    }    }
