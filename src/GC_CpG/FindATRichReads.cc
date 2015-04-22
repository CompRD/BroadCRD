// Copyright (c) 2000, 2001 Whitehead Institute for Biomedical Research

// FindATRichReads: find trimmed reads having >= 70% AT.

#include "Basevector.h"
#include "MainTools.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_Double_OrDefault(MIN_AT, 0.7);
     EndCommandArguments;

     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String source_dir = run_dir +
          ( IsDirectory( run_dir + "/preFindSeeds" ) ? "/preFindSeeds" : "" );

     vecbasevector reads;
     reads.ReadAll( source_dir + "/reads.fastb" );
     for ( int i = 0; i < (int) reads.size( ); i++ )
     {    basevector& b = reads[i];
          if ( b.size( ) == 0 ) continue;
          int AT = 0;
          for ( int j = 0; j < (int) b.size( ); j++ )
          {    unsigned char c = b[j];
               if ( as_base(c) == 'A' || as_base(c) == 'T' ) ++AT;    }
          if ( float(AT) / float( b.size( ) ) >= MIN_AT ) 
               cout << i << "\n";    }    }
