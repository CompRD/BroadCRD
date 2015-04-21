// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#include "CoreTools.h"
#include "math/Functions.h"
#include "lookup/LookAlignIndex.h"
#include "Vec.h"
#include "system/file/FileReader.h"
#include "feudal/BinaryStream.h"

void PullAligns( const String& dir, int target, int start, int stop, vec<char>& buf,
     vec<unsigned char>& mult, vec<QLT_index>* indexp )
{    vec<QLT_index> index0;
     if ( indexp == 0 ) 
     {    BinaryReader::readFile( dir + "/all_hits.index", &index0 );
          indexp = &index0;    }
     const vec<QLT_index>& index = *indexp;
     vec< pair<int, int> > target_start( index.size( ) );
     for ( int i = 0; i < index.isize( ); i++ )
          target_start[i] 
               = make_pair(index[i].target_id, index[i].start_on_target);
     pair<int, int> request1( target, start );
     int low = lower_bound( target_start.begin( ), target_start.end( ), request1 )
          - target_start.begin( );
     pair<int, int> request2( target, stop );
     int high = upper_bound( target_start.begin( ), target_start.end( ), request2 )
          - target_start.begin( );
     int ilow = Max( 0, low-1 ), ihigh = Min( high, index.isize( ) - 1 );
     while ( ilow >= 1 )
     {    if ( index[ilow-1].target_id == index[ilow].target_id )
          {    if ( index[ilow-1].stop_on_target <= start ) break;
               --ilow;    }
          else
          {    --ilow;
               break;    }    }
     longlong fbegin = index[ilow].start_on_hits_file, fend;
     if ( ihigh != index.isize( ) - 1 ) fend = index[ihigh].start_on_hits_file;
     else fend = FileSize( dir + "/all_hits" );
     buf.resize( fend - fbegin );
     if ( true )
     {
         FileReader fr( (dir + "/all_hits").c_str() );
         fr.seek( fbegin );
         fr.read( &buf[0], fend - fbegin );
     }

     int mbegin = index[ilow].alignment_id, mend;
     if ( ihigh != index.isize( ) - 1 ) mend = index[ihigh].alignment_id;
     else mend = BinaryVecNumElements( dir + "/all_hits.mult" );
     mult.ReadRange( dir + "/all_hits.mult", mbegin, mend );   }
