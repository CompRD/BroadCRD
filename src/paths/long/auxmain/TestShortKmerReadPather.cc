///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

const char *DOC =
  "Code to test finding read paths using short kmers.";


#include "MainTools.h"
#include "Basevector.h"
#include "feudal/PQVec.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/ShortKmerReadPather.h"
#include "paths/long/ExtendReadPath.h"


int main(int argc, char *argv[])
{
    RunTime( );
    
    BeginCommandArguments;
    CommandDoc(DOC);
    CommandArgument_String_Doc(READ_HEAD,
      "Base and quality score files: READ_HEAD.{fastb,qualb}");
    CommandArgument_String_Doc(PATHS, 
      "ReadPath file containing paths on the HBV");
    CommandArgument_String_OrDefault_Doc(HBV, "",
      "HyperBasevector (HBV) file");
    CommandArgument_String_OrDefault_Doc(READ_IDS, "",
      "Optional list of read IDs to display, or a file containing the list.");
    CommandArgument_Int_OrDefault_Doc(K, 24,
      "Short kmer size");
    EndCommandArguments;

    int num_threads = configNumThreads(0);
   
    String hbv_file = HBV;
    String hbv_paths_file = PATHS;

    // Import HBV if available, then build to and from index
    // Alternatively, use the edges vecbasevector
    HyperBasevector hbv;
    vec<int> to_left, to_right;
    BinaryReader::readFile( hbv_file, &hbv );
    hbv.ToLeft(to_left);
    hbv.ToRight(to_right);
    int edge_count = hbv.EdgeObjectCount();

    // Import read paths (selected)
    ReadPathVec hbv_paths;
    vec<int> ids;
    if (READ_IDS != "") {
	if (IsRegularFile(READ_IDS)) {
	    Ifstream( in, READ_IDS);
	    ids.ReadFromTextStream(in);
	} else {
	    ParseIntSet( READ_IDS, ids, /*sort*/ false );
	}
	hbv_paths.SparseRead(hbv_paths_file, ids);
    } else {
	hbv_paths.ReadAll(hbv_paths_file);
	ids = vec<int>(hbv_paths.size(), vec<int>::IDENTITY);
    }

    // Find empty paths
    vec<size_t> ids2;
    for ( auto read_id : ids )
	if (hbv_paths[read_id].empty())
	    ids2.push_back(read_id);

    // Import reads and quals, just the ones we are interested in
    vecbasevector bases;
    bases.SparseRead( READ_HEAD + ".fastb", ids2);
    VecPQVec quals;
    quals.SparseRead( READ_HEAD + ".qualp", ids2);

    ShortKmerReadPather::FindPaths(bases, quals, hbv, hbv_paths, 24, num_threads, ids2 , true);

}
