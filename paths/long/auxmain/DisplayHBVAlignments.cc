///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

const char *DOC =
  "Visualize readpaths on a hyperbasevector.";

// Visualization of a readpath.
// Requires an HBV (or edges fastb), read fastb and qualb, and a paths file
// Will display all readpaths, or just a single one as specified by READ
// For a single read an alternative offset and edge path can be given as:
//   offset:edge1,edge2,...
//
// The visualization lines are:
//   Quality scores
//   Mismatches
//   Read sequence
//   Path sequence
//   Edge overlap
//   <optional edges of the path>
//
// If the path starts or stops within the read, the overlap line uses
// [ or ] to indicate the edge is a source or a sink, or a + if it is
// connected to additional edges not in the the path.
// 
// If your only want to examine a couple of paths, the code runs faster if 
// supply an EDGE fastb rather than an HBV. 


#include "Basevector.h"
#include "Intvector.h"
#include "MainTools.h"
#include "Qualvector.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/ReadPathTools.h"
#include "feudal/VirtualMasterVec.h"
#include "paths/long/ExtendReadPath.h"


bool ValidatePathFromEdgesOnly(vec<int> edge_list, const int edge_count, 
			       String& message) {
    for (auto edge_id : edge_list)
	if (edge_id >= edge_count) {
	    message = "ERROR - Invalid edge ID: " + ToString(edge_id);
	    return false;
	}
    message = "";
    return true;
}


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
      "HyperBasevector (HBV) file (or supply EDGES instead - faster for subset of reads)");
    CommandArgument_String_OrDefault_Doc(EDGES, "", 
      "Fastb file containing edges in the HyperBasevector (or supply HBV instead)");
    CommandArgument_Int_OrDefault_Doc(K, -1, 
      "HyperBasevector K - only required if HBV is NOT supplied");
    CommandArgument_String_OrDefault_Doc(READ_IDS, "",
      "Optional list of read IDs to display.");
    CommandArgument_String_OrDefault_Doc(ALT_PATH, "", 
      "Alternative path to display for a single READ_IDS in the form:  offset:edge1,edge2,..." );
    CommandArgument_Bool_OrDefault_Doc(EXTEND_PATH, False, 
      "For a single read, extend the supplied ALT_PATH rightwards first");
    CommandArgument_Bool_OrDefault_Doc(DEBUG_EXTENSION, False, 
      "When using ALT_PATH, display extension debug info");
    CommandArgument_Bool_OrDefault_Doc(SHOW_EDGES, False, 
      "Show each overlapping edge in a path");
    CommandArgument_Bool_OrDefault_Doc(SHOW_RULER, False, 
      "Show a ruler counting bases");
    CommandArgument_Bool_OrDefault_Doc(USE_HBV, True, 
      "Use the HyperBasevector if both the HBV and EDGES are available.");
    CommandArgument_Bool_OrDefault_Doc(FILTER_SINGLES, False, 
      "Do not display paths that contain only a single edge");
    CommandArgument_String_OrDefault_Valid_Doc(DISPLAY_MODE, "paths", "{paths,validate}",
      "Display paths visualization or just validate them");
    EndCommandArguments;

    if (HBV == "" && EDGES == "") 
	InputErr("You must supply either and HBV or a EDGES file");

    
    bool use_hbv = (USE_HBV && (HBV != "")) || (EDGES == "");

    if (!use_hbv && K == -1) 
	InputErr("You must supply K if you don't use an HBV file");

    String hbv_file = HBV;
    String edges_file = EDGES;
    String hbv_paths_file = PATHS;

    // Parse display mode
    bool view_path = (DISPLAY_MODE == "paths");
    bool validate = (DISPLAY_MODE == "validate");

    // Import HBV if available, then build to and from index
    // Alternatively, use the edges vecbasevector
    HyperBasevector hbv;
    vec<int> to_left, to_right;
    VirtualMasterVec<basevector>* edges = nullptr;
    int edge_count;
    if (use_hbv) {
	BinaryReader::readFile( hbv_file, &hbv );
	hbv.ToLeft(to_left);
	hbv.ToRight(to_right);
	edge_count = hbv.EdgeObjectCount();
    } else {
	edges = new VirtualMasterVec<basevector>(edges_file);
	edge_count = edges->size(); 
    } 
    // Import readpath
    VirtualMasterVec<ReadPath> hbv_paths(hbv_paths_file);
   
    // Load reads.
    VirtualMasterVec<basevector> bases( READ_HEAD + ".fastb" );
    VirtualMasterVec<qualvector> quals( READ_HEAD + ".qualb" );

    vec<int> ids(bases.size(), vec<int>::IDENTITY);
    if (READ_IDS != "") {
	ids.clear();
	if (IsRegularFile(READ_IDS)) {
	    Ifstream( in, READ_IDS);
	    ids.ReadFromTextStream(in);
	} else {
	    ParseIntSet( READ_IDS, ids, /*sort*/ false );
	}
    }
     
    // Some stats
    int empty_count = 0;
    int invalid_count = 0;

    // just make sure that we have only one id if we use ALT_PATH
    if (ALT_PATH != "" && ids.size() != 1)
        FatalErr("you may only specify one read if you specify ALT_PATH");

    // Examine ReadPaths
    for ( size_t i = 0; i < ids.size( ); i++ ) {
	size_t read_id = ids[i];

	bool found_path = hbv_paths[read_id].size();
	ReadPath read_path = hbv_paths[read_id];

	if (!found_path)
	    empty_count++;

        // Examine a single read with an ALT_PATH
        if (ALT_PATH != "" && ids.size() == 1) {
            int offset = ALT_PATH.SafeBefore(":").Int();
            String path_str = ALT_PATH.SafeAfter(":");
            vec<int> edges;
            ParseIntSet(path_str, edges, /*sort*/false);
            IntVec x;
            for (auto i : edges)
                x.push_back(i);
            (IntVec&) read_path = x;
            read_path.setOffset(offset);
	    found_path = true;
	    if (EXTEND_PATH) {
		ExtendReadPath::attemptLeftRightExtension( read_path, bases[read_id], quals[read_id],
							   hbv, to_left, to_right, DEBUG_EXTENSION);
	    }
        }

	// apply filtereing
	if (FILTER_SINGLES && read_path.size() < 2 ) continue;

	int offset = read_path.getOffset();
	vec<int> edge_list(read_path.begin(), read_path.end());

	// validate the path
	bool valid = true;
	String message;
	if (found_path)
	    if (use_hbv)
		valid =  ValidateReadPath(hbv, to_left, to_right, offset, edge_list,
					  message, bases[read_id].size());
	    else
		valid = ValidatePathFromEdgesOnly(edge_list, edge_count, message);

	if (!valid)
	    invalid_count++;

	if (validate && !valid) { 
	    // validation path only
	    cout << "Path " << i << " = " << offset << ":" 
		 << printSeq( read_path ) << "  " << message <<endl;

	} else if (view_path) {  
	    // display the path
	    cout << "[" << read_id << "]";
	    cout << "  hbv_path= " << offset << ":"<< printSeq( edge_list );
	    cout << endl << endl;;
	    
	    if (found_path && valid) {

		if (use_hbv) {  // use the HBV as the source of edges
		    
		    DisplayReadPath( cout, hbv, to_left, to_right,
				     offset, edge_list,
				     bases[read_id], quals[read_id], SHOW_EDGES, SHOW_RULER);
		
		} else {  // use the edge fastb as the source of edges (limits output)
		    vecbvec path_edges;
		    for (auto edge_id : read_path) {
			path_edges.push_back((*edges)[edge_id]);
		    }
		    DisplayReadPath( cout, path_edges, K,
				     offset, vec<int>(path_edges.size(), vec<int>::IDENTITY),
				     bases[read_id], quals[read_id], SHOW_EDGES, SHOW_RULER);
		}
		
	    } else if (!valid)
		cout << message << endl;
	    else 
		cout << "No hbv path." << endl;
	
	    cout << endl << endl << endl;
	}
    }
    
    if (validate) {
	cout << "Graph edges    : " << edge_count << endl;
	cout << "Reads          : " << bases.size() << endl;
	cout << "Examined paths : " << ids.size() << endl;
	cout << "Empty paths    : " << empty_count << endl;
	cout << "Invalid paths  : " << invalid_count << endl;
    }
}
