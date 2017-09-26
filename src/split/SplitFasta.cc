// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

// Split a fasta file into records, each of size n, discarding remnants of 
// size less than n.

#include "MainTools.h"
#include "FastIfstream.h"

#include "random/Random.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(IN);
     CommandArgument_UnsignedInt(n);
     CommandArgument_Int_OrDefault_Doc(RDIV, 1, 
          "randomly keep only 1/RDIV records");
     EndCommandArguments;

     String line, header;
     fast_ifstream in(IN);
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          ForceAssert( line.Contains( ">", 0 ) );
          header = line;
          static vec<char> bases;
          bases.clear( );
          while(1)
          {    char c;
               in.peek(c);
               if ( c == '>' || in.fail( ) ) break;
               getline( in, line );
               if ( in.fail( ) ) break;
               for ( int j = 0; j < (int) line.size( ); j++ )
                    bases.push_back( line[j] );    }
          int chunk = 1;
          for ( int u = 0; u < bases.isize( ); u += n )
          {    if ( u + (int) n > bases.isize( ) ) break;
               if ( RDIV == 1 || randomx( ) % RDIV == 0 )
               {    cout << header << "." << chunk++ << "\n";
                    for ( int j = 0; j < (int) n; j++ )
                    {    if ( j > 0 && j % 80 == 0 ) cout << "\n";
                         cout << bases[u+j];    }
                    cout << "\n";    }    }
          if ( in.fail( ) ) break;    }    }
