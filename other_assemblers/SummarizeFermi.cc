///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Stupid program to print a couple of fermi assembly stats.

#include "MainTools.h"
#include "math/HoInterval.h"
#include "other_assemblers/ConnectionAssembly.h"
#include "other_assemblers/ExtractFermi.h"

int main(int argc, char *argv[])
{
     RunTime( );
     double clock = WallClockTime( );

     BeginCommandArguments;
     CommandArgument_String(IN_DIR);
     EndCommandArguments;

     connection_assembly A;
     ExtractFermi( IN_DIR, A );
     cout << A.bases.size( ) << " fermi edges" << endl;

     vecbasevector G( IN_DIR + "/genome.fastb" );

     vec< triple< int, ho_interval, vec<int> > > matches;
     MatchRef( A, G, matches );    
     for ( int j = 0; j < matches.isize( ); j++ )
     {    cout << matches[j].first << "." << matches[j].second << " :: ";
          const vec<int>& v = matches[j].third;
          for ( int l = 0; l < v.isize( ); l++ )
          {    if ( l > 0 ) cout << ",";
               if ( v[l] >= 0 ) cout << v[l] << "fw";
               else cout << -v[l]-1 << "rc";    }
          cout << "\n";    }

     const int max_gap_bases = 1000;
     Bool perf = False;
     for ( int j = 0; j < matches.isize( ); j++ )
     {    ho_interval h = matches[j].second;
          int gap_bases = h.Start( ) + 1000000 - h.Stop( );
          if ( gap_bases < max_gap_bases ) perf = True;    }
     cout << "fermi is " << ( perf ? "perfect" : "imperfect" ) << endl;

     const digraphE<connection> & g = A.G;
     const vecbasevector & bases = A.bases;

     int v_count = g.N();
     cout << "Vertices: " << v_count << endl;
     cout << "Edges: " << g.EdgeObjectCount() << endl;
     
     vec<String> v_labels(v_count);
     vec<int> v_list(v_count);
     for (int i = 0; i < v_count; ++i) {
       cout << "-------- " << i << "--------" << endl;
       A.PrintConnections(cout, i);
       v_list[i] = i;
       v_labels[i] = ToString(i) + " [" + ToString(bases[i].size()) + "]";
     }

     ofstream out;
     String filename = "fermi.dot";
     OpenOfstream( out, filename );

     g.DOT_vl(out, v_labels, v_list);
     out.close();

     HyperBasevector hb;
     vec<int> inv;
     ToHbv(A, 12, hb, inv, true); 

    }
