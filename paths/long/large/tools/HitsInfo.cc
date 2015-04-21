///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// HitsInfo. Find edges that align to a particular part of the genome.

#include "Basevector.h"
#include "Intvector.h"
#include "MainTools.h"
#include "Qualvector.h"
#include "VecUtilities.h"
#include "paths/long/ReadPath.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault_Doc(DIR, ".", 
          "looks for DIR/a.{fastb,inv,paths,paths.inv}");
     CommandArgument_String_Doc(X, 
          "region, in form a:b-c, where a is a chromosome number");
     EndCommandArguments;

     int chr = X.Before( ":" ).Int( ) - 1;
     int start = X.Between( ":", "-" ).Int( );
     int stop = X.After( "-" ).Int( );
         
     vec< vec< pair<int,int> > > hits;
     BinaryReader::readFile( DIR + "/a.aligns", &hits );

     vecbasevector edges( DIR + "/a.fastb" );

     vec< triple<int,int,int> > places;
     for ( int i = 0; i < hits.isize( ); i++ )
     for ( int j = 0; j < hits[i].isize( ); j++ )
     {    if ( hits[i][j].first != chr ) continue;
          places.push( hits[i][j].second, i, j );    }
     Sort(places);
     vec<int> es;
     for ( int i = 0; i < places.isize( ); i++ )
     {    int pos = places[i].first;
          if ( pos >= start && pos <= stop )
          {    int e = places[i].second;
               // index of placement of edge e, usually 0
               int s = places[i].third; 
               cout << pos << " - " << pos + edges[e].isize( )
                    << "   " << e << "[" << s+1 << "]\n";    
               es.push_back(e);    }    }    
     int ip = 0;
     for ( int i = 0; i < es.isize( ); i++ )
          if ( i == 0 || ( i > 0 && es[i] != es[i-1] ) ) es[ip++] = es[i];
     es.resize(ip);
     cout << "edges = " << printSeq(es) << endl;    }
