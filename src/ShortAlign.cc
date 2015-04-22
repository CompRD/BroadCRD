// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#include "CoreTools.h"
#include "ShortAlign.h"

void BuildIndex( const vec<short_align>& aligns, vec<int>& index, int nseq )
{    index.resize_and_set( nseq + 1, -1 );
     index[0] = 0;
     for ( int i = 0; i < aligns.isize( ); i++ )
          index[ aligns[i].Id1( ) + 1 ] = i + 1;
     for ( int i = 1; i <= nseq; i++ )
          if ( index[i] < 0 ) index[i] = index[i-1];    }
