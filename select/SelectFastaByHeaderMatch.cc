// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

const char * DOC = "Select the fasta records from a file IN whose header "
"line contains a given string M.  The file IN may be gzipped.";

#include "MainTools.h"
#include "FastIfstream.h"

int main(int argc, char ** argv) 
{
     BeginCommandArguments;
     CommandDoc(DOC);
     CommandArgument_String_Doc(IN, "Input fasta file.");
     CommandArgument_String_Doc(M, "String to test for (NOT regular expression).");
     CommandArgument_Bool_Abbr_OrDefault_Doc(COMPLEMENT,C,False,"If true, select lines *not* matching M");
     
     EndCommandArguments;

     String line;
     String com;
     if ( IN.Contains( ".gz", -1 ) ) com = "zcat " + IN;
     else com = "cat " + IN;
     fast_pipe_ifstream in(com);
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          ForceAssert( line.Contains( ">", 0 ) );
          Bool keep = line.Contains(M);
	  if (COMPLEMENT) keep = !keep;
          if (keep) cout << line << "\n";
          while(1)
          {    char c;
               in.peek(c);
               if ( in.fail( ) || c == '>' ) break;
               getline( in, line );
               if ( in.fail( ) ) break;
               if (keep) cout << line << "\n";    }
          if ( in.fail( ) ) break;    }    }
