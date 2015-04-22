///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/// GenomeUnipathStats.  Generate unipaths from a given genome, for a given K,
/// and display some statistics for them.

#include "Basevector.h"
#include "FeudalMimic.h"
#include "MainTools.h"
#include "graph/Digraph.h"
#include "math/Functions.h"
#include "paths/HyperKmerPath.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/Unipath.h"
#include <map>

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(GENOME, "program looks for GENOME.fastb");
     CommandArgument_String_OrDefault_Doc(DOT, "", 
          "generate DOT-style graph in this file");
     CommandArgument_Int_OrDefault(DOT_MIN_COMPONENT, 0);
     CommandArgument_Bool_OrDefault(DOT_REMOVE_HANGING_ENDS, False);
     CommandArgument_Bool_OrDefault(DOT_SIMPLE, False);
     CommandArgument_Bool_OrDefault(DOT_LABEL_CONTIGS, True);
     CommandArgument_Bool_OrDefault(DOT_LABEL_VERTICES, False);
     CommandArgument_Bool_OrDefault(DOT_LABEL_EDGES, False);
     CommandArgument_String_OrDefault(DOT_FASTA, "");
     CommandArgument_Bool_OrDefault(DOT_COV, False);
     CommandArgument_Int(K);
     CommandArgument_Bool_OrDefault_Doc(DUMP_SIZES, False, 
          "dump unipath sizes in bp");
     CommandArgument_Bool_OrDefault(FW_ONLY, False);
     EndCommandArguments;

     // Generate unipaths.

     vecbasevector genome( GENOME + ".fastb" );
     vecKmerPath paths, paths_rc, unipaths;
     vec<big_tagged_rpint> pathsdb, unipathsdb;
     if ( !FW_ONLY ) ReadsToPathsCoreY( genome, K, paths, paths_rc, pathsdb );
     else
     {    ReadsToPathsCoreY( genome, K, paths );
          CreateDatabase( paths, paths_rc, pathsdb );    }
     Unipath( paths, paths_rc, pathsdb, unipaths, unipathsdb );
     vec<int> ulen( unipaths.size( ) );
     for ( size_t id = 0; id < unipaths.size( ); id++ )
          ulen[id] = K - 1 + unipaths[id].KmerCount( );

     // Map them back to the genome.

     vec< vec<int> > unipath_tigs( unipaths.size( ) );
     vec< vec<int> > unipath_starts( unipaths.size( ) );
     vec< vec<int> > genome_starts( genome.size( ) );
     vec< vec<int> > genome_starts_rc( genome.size( ) );
     cout << Date( ) << ": building starts" << endl; // XXX
     for ( size_t i = 0; i < genome.size( ); i++ )
     {    int x = 0, xrc = 0;
          for ( int j = 0; j < paths[i].NSegments( ); j++ )
          {    genome_starts[i].push_back(x);
               x += paths[i].Length(j);    }
          if ( !FW_ONLY )
          {    for ( int j = 0; j < paths_rc[i].NSegments( ); j++ )
               {    genome_starts_rc[i].push_back(xrc);
                    xrc += paths_rc[i].Length(j);    }    }    }
     cout << Date( ) << ": mapping back" << endl; // XXX
     int64_t total_kmers = 0, sub_100_kmers = 0;
     for ( size_t id = 0; id < unipaths.size( ); id++ )
     {    const KmerPath& u = unipaths[id];
          longlong l = u.Start(0);
          vec<longlong> con;
          Contains( pathsdb, l, con );
          for ( int j = 0; j < con.isize( ); j++ )
          {    const big_tagged_rpint& t = pathsdb[ con[j] ];
               int tig = t.ReadId( );
               longlong start;
               if ( t.Fw( ) ) 
               {    start = l - t.Start( ) + genome_starts[tig][ t.PathPos( ) ];    }
               else 
               {    start = l - t.Start( ) + genome_starts_rc[tig][ t.PathPos( ) ];
                    start = genome[tig].isize( ) - ( start + ulen[id] );    }
               unipath_tigs[id].push_back(tig);
               unipath_starts[id].push_back(start);    
               total_kmers += u.KmerCount( );
               if ( u.KmerCount( ) < 100 ) sub_100_kmers += u.KmerCount( );
               }    }
     cout << "\nfraction of placed kmers in unipaths shorter than "
          << "100 kmers: " << PERCENT_RATIO( 3, sub_100_kmers, total_kmers )
          << "\n\n";

     // Generate statistics.

     cout << Date( ) << ": generating stats" << endl << endl; // XXX
     vec<int> sizes;
     for ( size_t i = 0; i < unipaths.size( ); i++ )
          sizes.push_back( unipaths[i].KmerCount( ) + K - 1 );
     Sort(sizes);
     cout << "There are " << unipaths.size( ) << " unipaths, having N50 size "
          << N50(sizes) << " bp.\n";
     cout << "Coverage statistics:\n";
     cout << "(Note that an n-kmer unipath covers n+K-1 bases, which can "
          << "be misleading.)\n";
     int top = 1;
     for ( int d = 0; d <= 5; d++ )
     {    vecbitvector cov;
          Mimic( genome, cov );
          int count = 0;
          for ( size_t id = 0; id < unipaths.size( ); id++ )
          {    if ( ulen[id] >= top )
               {    ++count;
                    for ( int j = 0; j < unipath_starts[id].isize( ); j++ )
                    {    int tig = unipath_tigs[id][j];
                         int start = unipath_starts[id][j];
                         for ( int k = start; k < start + ulen[id]; k++ )
                         {    if ( k < 0 || static_cast<unsigned>(k) >= cov[tig].size( ) ) continue;
                              cov[tig].Set( k, True );    }    }    }    }
          cout << "coverage by " << count << " unipaths >= " << top << "bp: "
               << ToString( 100.0 * Coverage(cov), 2 ) << "%" << endl;
          top *= 10;    }


     // Generate copy-number statistics
     {
       cout << "\nCoverage statistics by unipaths with copy numbers bigger than 1:\n";
       map<int,int> cn_freqs;
       vecbitvector cov;
       cov.resize( genome.size( ) );
       for (size_t ii=0; ii<genome.size( ); ii++)
	 cov[ii].resize( Max( 0, genome[ii].isize( ) - K + 1 ) );
       int count = 0;
       for ( size_t id = 0; id < unipaths.size( ); id++ ){
	 int cn = unipath_starts[id].isize( );
	 cn_freqs[ cn ]++;
	 if ( unipath_starts[id].isize( ) <= 1 )
	   continue;
	 count++;
	 for ( int j = 0; j < unipath_starts[id].isize( ); j++ ){    
	   int tig = unipath_tigs[id][j];
	   int start = unipath_starts[id][j];
	   for ( int k = start; k < start + ulen[id] - K + 1; k++ ){    
	     if ( k < 0 || static_cast<unsigned>(k) >= cov[tig].size( ) ) continue;
	     cov[tig].Set( k, True );    
	   }    
	 }    
       }
       cout << "coverage by " << count << " repetitive unipaths "
	    << ToString( 100.0 * Coverage(cov), 4 ) << "%" << endl;
       cout << "\ncopy-number unipaths count:\n";
       for ( map<int,int>::iterator it = cn_freqs.begin(); it != cn_freqs.end(); it++ ){
	 cout.width(9);
	 cout .fill(' ');
	 cout << it -> first << " " ;
	 cout.width(10);
	 cout .fill(' ');
	 cout << it -> second << "\n";
       }
       cout << endl;
     }
    

     // Generate graph.

     if ( DOT != "" )
     {    cout << Date( ) << ": Writing to DOT file" << endl;
          Ofstream( out, DOT );
          digraph A;
          BuildUnipathAdjacencyGraph( paths, paths_rc, pathsdb, unipaths,
               unipathsdb, A );
          HyperKmerPath h;
          BuildUnipathAdjacencyHyperKmerPath( K, A, unipaths, h );
          if (DOT_REMOVE_HANGING_ENDS)
               RemoveHangingEnds( h, &KmerPath::KmerCount, 250, 5.0 );
          if ( DOT_MIN_COMPONENT > 0 ) h.RemoveSmallComponents(DOT_MIN_COMPONENT);
          h.RemoveDeadEdgeObjects( );
          h.RemoveEdgelessVertices( );
          h.RemoveUnneededVertices( );

          vec<String> edge_labels_extra( h.EdgeObjectCount( ) );
          if (DOT_COV)
          {    for ( int e = 0; e < h.EdgeObjectCount( ); e++ )
               {    const KmerPath& p = h.EdgeObject(e);
                    longlong cov = 0;
                    for ( int j = 0; j < p.NSegments( ); j++ )
                    {    for ( longlong x = p.Start(j); x <= p.Stop(j); x++ )
                         {    vec<longlong> con;
                              Contains( pathsdb, x, con );
                              cov += con.size( );    }    }
                    ostringstream out;
                    if ( p.KmerCount( ) < 1000 )
                         out << "(" << p.KmerCount( ) << " b)";
                    out << "[" << setiosflags(ios::fixed) << setprecision(1) 
                         << double(cov) / double( p.KmerCount( ) ) << "]";
                    edge_labels_extra[e] = out.str( );    }    }

          if ( DOT_FASTA != "" )
          {    KmerBaseBrokerBig kbb( K, paths, paths_rc, pathsdb, genome );
               Ofstream( out, DOT_FASTA )
               for ( int i = 0; i < h.EdgeObjectCount( ); i++ )
               {    kbb.Seq( h.EdgeObject(i) ).Print( 
                         out, "edge_" + ToString(i) );    }    }

          if (DOT_SIMPLE) h.PrintSummaryDOT0(out);
          else
          {    h.PrintSummaryDOT0w( out, DOT_LABEL_CONTIGS, DOT_LABEL_VERTICES,
                    DOT_LABEL_EDGES, NULL, False, 
                    DOT_COV ? &edge_labels_extra : NULL );    }    }

     // Dump unipath sizes.

     if (DUMP_SIZES)
     {    cout << "Unipath sizes:\n";
          for ( int i = 0; i < sizes.isize( ); i++ )
               cout << sizes[i] << "\n";    }
     
     cout << Date( ) << ": Done!" << endl;
}
