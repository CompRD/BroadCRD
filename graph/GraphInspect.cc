///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// GraphInspect.  Rudimentary interactive tool to inspect a digraph object.

#include "MainTools.h"
#include "Set.h"
#include "graph/Digraph.h"
#include "feudal/BinaryStream.h"

int main( int argc, char *argv[] )
{
     BeginCommandArguments;
     CommandArgument_String(G);
     EndCommandArguments;

     digraph X;
     BinaryReader::readFile( G, &X );
     const uint max_edges = 1000000;
     const uint max_edges_to_show = 100;
     while(1)
     {    cout << "\nu=";
          int u;
          cin >> u;
          for ( int i = 0; i < X.To(u).isize( ); i++ )
               cout << X.To(u)[i] << " --> u" << endl;
          for ( int i = 0; i < X.From(u).isize( ); i++ )
               cout << "u --> " << X.From(u)[i] << endl;    
          set<int> C0, C;
          C0.insert(u);
          while( C0.size( ) > 0 && C0.size( ) + C.size( ) <= max_edges )
          {    int x = *( C0.begin( ) );
               C0.erase( C0.begin( ) );
               C.insert(x);
               for ( int i = 0; i < X.To(x).isize( ); i++ )
               {    int w = X.To(x)[i];
                    if ( !Member( C, w ) ) C0.insert(w);    }
               for ( int i = 0; i < X.From(x).isize( ); i++ )
               {    int w = X.From(x)[i];
                    if ( !Member( C, w ) ) C0.insert(w);    }    }
          vec<int> V;
          for ( set<int>::iterator i = C0.begin( ); i != C0.end( ); i++ )
               V.push_back( *i );
          for ( set<int>::iterator i = C.begin( ); i != C.end( ); i++ )
               V.push_back( *i );
          UniqueSort(V);
          if ( V.size( ) >= max_edges ) 
               cout << "component has " << max_edges << "+ edges" << endl;
          else if ( V.size( ) > max_edges_to_show )
          {    cout << "component has " << V.size( ) << " edges" << endl;    }
          else
          {    cout << "component = {" << V[0];
               for ( int j = 1; j < V.isize( ); j++ )
                    cout << "," << V[j];
               cout << "}" << endl;    }    }    }
