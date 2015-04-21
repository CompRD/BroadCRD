// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

// Qualb: convert a fasta quality score file into a qualb file.

#include <ctype.h>

#include <strstream>

#include "MainTools.h"
#include "FastIfstream.h"
#include "Qualvector.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(FILE);
     EndCommandArguments;

     String filename = PRE + "/" + FILE, new_filename;
     if ( filename.Contains( ".fasta", -1 ) )
          new_filename = filename.Before ( ".fasta" ) + ".qualb";
     else if ( filename.Contains( ".qual", -1 ) )
          new_filename = filename.Before ( ".qual" ) + ".qualb";
     else FatalErr( "FILE has unrecognized suffix." );

     vecqualvector Qs;

     int total_seqs = 0;
     longlong total_bases = 0;
     for ( int pass = 1; pass <= 2; pass++ )
     {    if ( pass == 2 ) Qs.Reserve( total_bases, total_seqs );
          fast_ifstream quals(filename);
          String line;
          qualvector q;
          Bool first = True;
          while(1)
          {    getline( quals, line );
               if ( quals.fail( ) )
               {    if ( pass == 1 ) 
                    {    total_bases += q.size( );
                         ++total_seqs;    }
                    if ( pass == 2 ) Qs.push_back(q);
                    q.clear( );
                    break;    }
               ForceAssert( line.size( ) > 0 );
               if ( line[0] == '>' )
               {    if ( !first )
                    {    if ( pass == 1 ) 
                         {    total_bases += q.size( );
                              ++total_seqs;    }
                         if ( pass == 2 ) Qs.push_back(q);    }
                    first = False;
                    q.clear( );
                    while(1)
                    {    char c;
                         quals.peek(c);
                         if ( quals.fail( ) || c == '>' ) break;
                         getline( quals, line );
                         for ( int j = 0; j < (int) line.size( ); j++ )
                         {    if ( !isspace(line[j]) && !isdigit(line[j]) )
                                   cout << "illegal line: \"" << line
                                        << "\"" << endl;
                              ForceAssert( isspace(line[j]) || isdigit(line[j]) );
                                   }
                         istrstream i( line.c_str( ) );
                         while(1)
                         {    int n;
                              i >> n;
                              if ( i.fail( ) ) break;
                              q.push_back(n);    }    }    }
               if ( quals.fail( ) )
               {    if ( pass == 1 ) 
                    {    total_bases += q.size( );
                         ++total_seqs;    }
                    if ( pass == 2 ) Qs.push_back(q);
                    q.clear( );
                    break;    }    }    }
          
     Qs.WriteAll( new_filename );    }
