///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// UnipathDot.  Generate DOT file from unipaths.

#include "Equiv.h"
#include "MainTools.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerPath.h"
#include "paths/Unipath.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(HEAD);
     CommandArgument_Int(K);
     CommandArgument_Bool_OrDefault(LABEL_CONTIGS, True);
     CommandArgument_Bool_OrDefault(LABEL_EDGES, True);
     EndCommandArguments;

     // Load data.

     String KS = ToString(K);
     vecKmerPath paths( HEAD + ".paths.k" + KS );
     vecKmerPath pathsrc( HEAD + ".paths_rc.k" + KS );
     vecKmerPath unipaths( HEAD + ".unipaths.k" + KS );
     BREAD2( HEAD + ".pathsdb.k" + KS, vec<tagged_rpint>, pathsdb );
     BREAD2( HEAD + ".unipathsdb.k" + KS, vec<tagged_rpint>, unipathsdb );

     // Build hyperkmerpath.

     digraph A;
     BuildUnipathAdjacencyGraph( paths, pathsrc, pathsdb, unipaths, unipathsdb, A );
     HyperKmerPath h;
     BuildUnipathAdjacencyHyperKmerPath( K, A, unipaths, h );

     // Generate dot file for it.

     Ofstream( out, HEAD + ".unipaths.k" + KS + ".dot" );
     Ofstream( lout, HEAD + ".unipaths.k" + KS + ".dot.log" );
     h.PrintSummaryDOT0w( out, LABEL_CONTIGS, False, LABEL_EDGES, 0, False );    

     // Print some summary statistics about it.

     vec< vec<int> > comp;
     h.Components(comp);
     lout << comp.size( ) << " components of unipath graph:\n";
     for ( int i = 0; i < comp.isize( ); i++ )
     {    int edges = 0;
          for ( int j = 0; j < comp[i].isize( ); j++ )
               edges += h.From( comp[i][j] ).size( );
          lout << "\n[" << i+1 << "] has " << edges << " edges\n";    }    }
