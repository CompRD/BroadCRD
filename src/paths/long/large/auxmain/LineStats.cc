///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Takes a DISCOVAR assembly and converts the lines it contains into contigs
// by flattening ambiguities, picking the first alternative. Will ignore RC
// lines, even in the case when the graph is broken and not symmetric.


#include "MainTools.h"
#include "VecUtilities.h"
#include "STLExtensions.h"
#include "graph/DigraphTemplate.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/Lines.h"


int main( int argc, char *argv[] ) {
    RunTime( );

    BeginCommandArguments;
    CommandArgument_String_Doc(HBV_HEAD, "HyperBasevector (HBV) file");
    CommandArgument_Bool_OrDefault_Doc(LINE_IDS, False, "Write out Line ID");
    EndCommandArguments;

    String hbv_file = HBV_HEAD + ".hbv";
    String lines_file = HBV_HEAD + ".lines";

    // Load the hyperbasevector (hbv)
    cout << "Loading HBV file" << endl;
    HyperBasevector hbv;
    BinaryReader::readFile( hbv_file, &hbv );
    int K = hbv.K();

    // Load the lines file
    cout << "Loading lines file" << endl;
    vec<vec<vec<vec<int>>>> lines;
    BinaryReader::readFile( lines_file, &lines );
    int line_count = lines.size();  // number of lines
    
    // Display the lines for diagnositic purposes.
    vec<int> empty_lines;
    vec<int> partially_empty_lines;
    vec<int> gapped_lines;

    int line_id = 0;
    for (const auto& line : lines) {
	for (const auto& paths : line) {
	    if (paths.size() == 1 && paths.front().empty()) 
		gapped_lines.push_back(line_id);
	    else {
		for (const auto& path : paths) {
		    for (int edge: path) { 
			if (hbv.EdgeLength(edge) == 0) {
			    if (line.size() == 1)
				empty_lines.push_back(line_id);
			    else 
				partially_empty_lines.push_back(line_id);
			}
		    }
		}
	    }
	}
	line_id++;
    }
    
    UniqueSort(partially_empty_lines);
    int gap_count = gapped_lines.size();
    UniqueSort(gapped_lines);

    cout << "Lines: " << line_count << endl;
    cout << "Edges: " << hbv.E() << endl;
    cout << "Gaps : " << gap_count << endl;
    cout << "Gap lines  : " << gapped_lines.size() << endl;
    cout << "Empty lines : " << empty_lines.size() << endl;
    cout << "Invalid lines : " << partially_empty_lines.size() << endl;

    if (LINE_IDS) {
	cout << "Invalid lines   " << printSeq(partially_empty_lines) << endl;
	cout << "Empty lines " << printSeq(empty_lines) << endl;
	cout << "Gapped lines " << printSeq(gapped_lines) << endl;
    }
}

    
