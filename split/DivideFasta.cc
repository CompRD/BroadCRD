// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

// DivideFasta: divide a fasta file into roughly N pieces of roughly equal size.
//
// If SPLIT_OUT=True, new files have names OUT.1, OUT.2, ....
// Otherwise, new file is OUT.
//
// If BREAK_ANYWHERE=False, only break at fasta record boundaries.

#include "MainTools.h"
#include "FastIfstream.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(IN);
     CommandArgument_String(OUT);
     CommandArgument_UnsignedInt(N);
     CommandArgument_Bool_OrDefault(BREAK_ANYWHERE, True);
     CommandArgument_Bool_OrDefault(SPLIT_OUT, False);
     EndCommandArguments;

     longlong n = FileSize(IN)/N;
     fast_ifstream in(IN);
     String line;
     int filecount = 1, chunk_count = 1;
     String last_header;
     if ( !SPLIT_OUT ) Remove(OUT);
     while(1)
     {    longlong charcount = 0;
          String outfile = OUT;
          if (SPLIT_OUT) 
          {    outfile += "." + ToString(filecount++);
               Remove(outfile);    }
          ofstream out(outfile.c_str( ), ios::app);
          while(1)
          {    getline( in , line );
               if ( in.fail( ) ) break;
               if ( line.Contains( ">", 0 ) )
               {    out << line << " (CHUNK 1)\n";
                    last_header = line;
                    chunk_count = 1;    }
               else out << line << "\n";
               charcount += line.size( ) + 1;
               if ( charcount > n )
               {    char c;
                    in.peek(c);
                    if ( in.fail( ) || c == '>' ) break;    
                    if (BREAK_ANYWHERE)
                    {    ++chunk_count;
                         out << last_header << " (CHUNK " << chunk_count << ")\n";
                         break;    }    }    }
          if ( in.fail( ) ) break;    }    }
