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
#include "paths/long/LongReadsToPaths.h"
#include <map>

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(GENOME, "program looks for GENOME");
     CommandArgument_String_Doc(HB, "outputs HyperBasevector");
     CommandArgument_Bool_OrDefault_Doc(FW_ONLY, false, "only output the FW strand");
     CommandArgument_String_OrDefault_Doc(DOT,"","outputs DOT");
     CommandArgument_String_OrDefault_Doc(FASTB,"","outputs FASTB");
     CommandArgument_String_OrDefault_Doc(FASTA,"","outputs FASTA");
     CommandArgument_Int(K);
     EndCommandArguments;

     // Generate unipaths.

     vecbasevector genome( GENOME );
     vecKmerPath paths, paths_rc;
     HyperBasevector hbv;
     HyperBasevectorX hb;
     digraph A;
     HyperKmerPath h;
     unsigned const COVERAGE = 5u;
     uint loglevel = 1;
     LongReadsToPaths( genome, K, COVERAGE, loglevel, false,
           &hbv, &h, &paths, &paths_rc );
     cout << Date( ) << ": constructing hb" << endl;
     hb = HyperBasevectorX(hbv);
     cout << Date() << ": Writing HBX to " << HB << endl;
     BinaryWriter::writeFile( HB, hb );

     if ( DOT != "" ) 
     {    cout << Date() << ": Writing DOT to " << DOT << endl;
	  vec<vec<String> > edge_labels;
	  for ( int v = 0; v < h.N(); v++ ) 
          {    edge_labels.push_back(vec<String>() );
	       for ( size_t j = 0; j < hbv.From(v).size(); j++ ) {
		    int e = hbv.EdgeObjectIndexByIndexFrom(v, j);
		    std::ostringstream s;
		    s << BaseAlpha(e) << " (" << hbv.EdgeObject( e ).size() <<")";
		    edge_labels.back().push_back( s.str() );    }    }
	  ofstream dot(DOT);
	  hbv.DOT(dot, edge_labels);
	  dot.close( );    }
     if ( FASTB != "" ) 
     {    cout << Date() << ": Writing FASTB to " << FASTB << endl;
	  hbv.DumpFastb(FASTB);    }

     if ( FASTA != "" ) 
     {    cout << Date() << ": Writing FASTA to " << FASTA << endl;
	  hbv.DumpFasta(FASTA);    }
     
     cout << Date( ) << ": Done!" << endl;
     Scram(0);
}
