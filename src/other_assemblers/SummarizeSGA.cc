///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Stupid program to print a couple of SGA assembly stats.

#include "MainTools.h"
#include "math/HoInterval.h"
#include "other_assemblers/ConnectionAssembly.h"
#include "other_assemblers/ExtractSGA.h"

int main(int argc, char *argv[])
{
     RunTime( );
     double clock = WallClockTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(IN_DIR,"location of SGA Assembly");
     CommandArgument_Bool_OrDefault_Doc(VARIANTS, True, "Include variants.");
     EndCommandArguments;

     connection_assembly A;
     ExtractSGA( IN_DIR, A, VARIANTS );
     cout << A.bases.size( ) << " SGA vertices" << endl;

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
     cout << "SGA is " << ( perf ? "perfect" : "imperfect" ) << endl;

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
     String filename = "SGA.dot";
     OpenOfstream( out, filename );

     g.DOT_vl(out, v_labels, v_list);
     out.close();

     std::cout << "after dot" <<std::endl;

     digraphE<BaseVec2KmerPath> sequence_graph;

     MakeSequenceGraph(sequence_graph,A);
     std::cout << "after MakeSequenceGraph"<<std::endl;

     SimplifySequenceGraph(sequence_graph);
     HyperKmerize(sequence_graph, 10);
//     std::cout << "after Simplification"<<std::endl;

     vec< vec<String>> edge_labels(sequence_graph.N());
     vec< String> vertex_labels(sequence_graph.N());

     String fa_filename = "sequence_graph.fasta";
     Ofstream( fa_out, fa_filename );

     int nVertices=0;

     for(size_t ii=0;ii<edge_labels.size();++ii){
         vertex_labels[ii]=ToString(ii);
         int nn = sequence_graph.FromSize(ii);
         edge_labels[ii].resize(nn);
         for(int jj = 0 ; jj < nn ; ++jj){
             auto& bases =  sequence_graph.EdgeObjectByIndexFrom(ii,jj);
             int length = bases.isize();
             if( length==0){
                 edge_labels[ii][jj] = 'X';
             }
             else{
                 edge_labels[ii][jj] = ToString(nVertices)+" ("+ToString(length)+")";
                 bases.Print(fa_out,nVertices);
                 ++nVertices;
//                     (length)?(ToString(length)):(sequence_graph.EdgeObjectByIndexFrom(ii,jj).ToString());
             }
         }
     }
     fa_out.close();






     String sg_filename = "sequence_graph.dot";
     Ofstream( sg_out, sg_filename );
//     sequence_graph.DOT(sg_out,edge_labels);
     sequence_graph.DOT_vl(sg_out,vertex_labels,"",vec<vec<String>>(), vec<String>(), edge_labels,vec<String>());
     sg_out.close();


}
