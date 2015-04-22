// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

// CatVecKmerPath: given a list of files IN containing vecKmerPaths, concatenate
// them to create a single file OUT.

#include "MainTools.h"
#include "ParseSet.h"
#include "paths/KmerPath.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(IN);
     CommandArgument_String(OUT);
     EndCommandArguments;

     vec<String> in;
     ParseStringSet( IN, in );
     vecKmerPath p;
     for ( int i = 0; i < in.isize( ); i++ )
     {    vecKmerPath q( in[i] );
          p.Append( q, 0, q.size( ) );    }
     p.WriteAll(OUT);    }
