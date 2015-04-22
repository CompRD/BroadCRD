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
#include "paths/UnibaseUtils.h"
#include "paths/long/large/Lines.h"

// Write the line RC map to standard out
void display_map( const map<int, int>& inv_map, int line_size ) {
    for (auto itr = inv_map.begin(); itr != inv_map.end(); ++itr) {
	double frac = static_cast<double>(itr->second) / line_size;
	cout << itr->first << " (" << frac << "), ";
    }
    cout << endl;
}

// Count the number of edges in a line
vec<int> edges_in_line(const vec<vec<vec<int>>>& line) {
    vec<int> edges;
    for (const auto& paths : line) 
	for (const auto& path : paths) 
	    edges.append(path);
    UniqueSort(edges);
    return edges;
}

// Write the line to standard out
void display_line(const vec<vec<vec<int>>>& line) {
    vec<int> edges;
    for (const auto& paths : line) {
	cout << " {";
	for (const auto& path : paths) {
	    cout << " (";
	    for (int edge: path) 
		cout << edge << " ";
	    cout << ") ";
	}
	cout << "} ";
    }
}


// Flatten a line to a linear path by picking the first edge at cells and SNPs
// Note that gap edges are represented by -1
vec<int> flatten_line(const vec<vec<vec<int>>>& line) {
    vec<int> edges;
    for (const auto& paths : line) {
	if (paths.size() == 1 && paths.front().empty()) 
	    edges.push_back(-1); // Gap edge
	else
	    edges.append(paths.front()); // Alt path or SNP (general)
    }
    return edges;
}

int path_length(const vec<int>& path, const vecbasevector& edges, const int K, const int gap_size = 20) {
    int length = 0;
    bool gap = true;   // treat first edge as following a gap
    for (int edge_id : path) {
	if (gap) {  // last edge was a gap
	    gap = false;
	    length += edges[edge_id].size();
	} else if (edge_id == -1) {  // this edge is a gap 
	    gap = true;
	    length += gap_size;
	} else
	    length += edges[edge_id].size() - (K-1);
    }
    return length;
}

// Convert a list of edges to fasta sequence, replacing zero length edges with Ns
String gapped_fasta(const vec<int>& path, const vecbasevector& edges, const int K, const String gap_str) {
    int length = path_length(path, edges, K, gap_str.size());
    String bases;
    if (length <= 0)
	return bases;
    bases.reserve(length);
    bool gap = true;   // treat first edge as following a gap
    for (int edge_id : path) {
	if (gap) {  // last edge was a gap
	    gap = false;
	    bases.append(edges[edge_id].ToString());
	} else if (edge_id == -1) {  // this edge is a gap 
	    gap = true;
	    bases.append(gap_str);
	} else
	    bases.append(edges[edge_id].ToString(), K - 1, edges[edge_id].size() - (K -1)) ;
    }
    return bases;
}


int main( int argc, char *argv[] )
{
    RunTime( );

    BeginCommandArguments;
    CommandArgument_String_Doc(HBV_HEAD, "HyperBasevector (HBV) file");
    CommandArgument_Int_OrDefault_Doc(GAP_SIZE, 20, "Number of Ns to use to represent a gap");
    CommandArgument_Int_OrDefault_Doc(MIN_LENGTH, 1000, "Ignore lines smaller than this (including Ns)");
    CommandArgument_Bool_OrDefault_Doc(DEBUG, False, "Display debug information");
    CommandArgument_Bool_OrDefault_Doc(WARN, False, "Display warnings");
    EndCommandArguments;

    String hbv_file = HBV_HEAD + ".hbv";
    String inv_file = HBV_HEAD + ".inv";
    String lines_file = HBV_HEAD + ".lines";
    String fastb_file = HBV_HEAD + ".contigs.fastb";
    String fasta_file = HBV_HEAD + ".contigs.fasta";

    // Load the hyperbasevector (hbv)
    // This is the graphical respresentation of the assembly, where
    // edges represent sequence, and nodes just provide adjacency info.
    // It is a subclass of DigraphE
    cout << "Loading HBV file" << endl;
    HyperBasevector hbv;
    BinaryReader::readFile( hbv_file, &hbv );
    int edge_count = hbv.EdgeObjectCount();  // number of edges in the graph
    int K = hbv.K();

    // Extract edge sequences from the graph (hbv) and store them in a
    // vecbasevector, with one entry per edge. The vecbasevector is
    // a vector of basevectors (sequences of ATC and G).
    // It uses a special data structure which is optimized for nested vectors,
    // and handles its own memory allocations.
    vecbvec edges( hbv.Edges().begin(),hbv.Edges().end() );

    // Load the lines file
    // Lines are linear-ish paths through the graph.
    // They contain bubbles and other small features but do not branch globally.
    // Uses a quad vector nested as follows:
    // lines = vec<line>, and each
    // line = vec<paths>, and each
    // paths = vec<path>, and each
    // path = vec<edge Ids>
    cout << "Loading lines file" << endl;
    vec<vec<vec<vec<int>>>> lines;
    BinaryReader::readFile( lines_file, &lines );
    int line_count = lines.size();  // number of lines
    cout << "Found " << line_count << " lines consisting of " << edges.size() << " edges" << endl;
    
    // Compute the edge involution
    // The assembly graph is symmetric and each edge is represented twice,
    // Corresponding to the forward and reverse complement base sequence.
    // The involution provides a mapping of an edge ID onto its complement.
    cout << "Computing involution" << endl;
    vec<int> inv;
    UnibaseInvolution( edges, inv );

    // Display the lines for diagnositic purposes.
    if (DEBUG) {
	for (size_t line_id = 0; line_id < lines.size(); ++line_id) {
	    cout << line_id << ": ";
	    display_line(lines[line_id]);
	    cout << endl;
	}
    }

    // build edge/line lookup tables
    cout << "Building lookup tables" << endl;
    vec<vec<int>> line_2_edges(line_count);
    for (size_t line_id = 0; line_id < lines.size(); ++line_id) 
	line_2_edges[line_id] = edges_in_line(lines[line_id]);
    vec<int> edge_2_line(edge_count, -1);    
    for (size_t line_id = 0; line_id < lines.size(); ++line_id) 
        for (auto edge : line_2_edges[line_id]) 
	    edge_2_line[edge] = line_id;

    // Compute potential line inversions. This should be a one to one relationship
    // but due to problems in the assembly graph it might not be.
    // Vec < Map < RC line, RC line size> >
    cout << "Building line involution table" << endl;
    int empty_edge_count = 0;
    vec<map<int, int>> inv_map(lines.size());
    vec<int> line_size(lines.size(), 0);
    for (size_t line_id = 0; line_id < lines.size(); ++line_id) {
	map<int, int>& this_inv_map = inv_map[line_id];
	for (auto edge : line_2_edges[line_id]) {
	    if ( edges[edge].size() == 0) { // Some edges are zero length - ignore them
		empty_edge_count++;
		continue;
	    }
	    this_inv_map[edge_2_line[inv[edge]]] += edges[edge].size();
	    line_size[line_id] += edges[edge].size();
	}
    }
    cout << "WARNING: Found " << empty_edge_count << " empty edges" << endl;

    // Create an index into the lines by order of size (largest to smallest)
    vec<int> size_index(line_size.size(),vec<int>::IDENTITY);
    vec<int> size_copy(line_size);
    ReverseSortSync(size_copy, size_index);
    

    // Bookkeeping
    vec<int> chosen; // Lines to use when building contigs 
    set<int> used;   // Lines used (either chosen or their RC)

    // Pick lines to turn into contigs, discarding their RC.
    // More complicated than it should be due to lack of symmetry in the graph
    // Start with the largest line and work downwards
    cout << "Eliminating RC lines." << endl;
    for (int line_id : size_index) {
	
	if (used.find(line_id) != used.end())  // Skip if we have already used this line's RC
	    continue;

	// Find the RC lines of this line, skip if there isn't one but log the error
	const auto& this_line = inv_map[line_id];
	if (this_line.empty()) { 
	    if (WARN)
		cout << "Missing: " << line_id << " (" << line_size[line_id] << ") has no RC)" << endl;
	    continue;
	}

	// Find the largest RC line (usually only one unless something is wrong)
	const auto& max_value = max_element(this_line.begin(), this_line.end(), [] (const pair<int,int>& one, const pair<int,int>& two) {return one.second < two.second;});

	// Check that the we haven't used any of the RC lines before
	int inv_line_id = max_value->first;
	auto& inv_line = inv_map[inv_line_id];
	const auto& inv_max_value = max_element(inv_line.begin(), inv_line.end(), [] (const pair<int,int>& one, const pair<int,int>& two) {return one.second < two.second;});
	for (const auto& value : inv_line)
	if (used.find(inv_line_id) != used.end())
	    continue;

	// Log potential errors, but don't stop
	if ( (inv_max_value->first != line_id || this_line.size() != 1 || inv_line.size() != 1) && WARN ){
	    cout << "Mismatch: " << line_id << " (" << line_size[line_id] << ") v " << inv_line_id << " (" << line_size[inv_line_id] << ")" << endl;
	    cout << line_id << "-> ";
	    display_map(this_line, line_size[line_id]);
	    cout << inv_line_id << "-> ";
	    display_map(inv_line, line_size[inv_line_id]);
	}

	// If we have reached this point this is a good line to use to build a contig
	chosen.push_back(line_id);
	// Update the list of used lines
	used.insert(line_id);
	std::transform(this_line.begin(), this_line.end(), std::inserter(used, used.begin()), [] (const pair<int,int>& value) { return value.first;});
	
    }
    
    // Write a list of chosen line IDs for diagnostic purposes.
    if (DEBUG)  {
	cout << "Remaining lines: ";
	for (int line_id : chosen) {
	    cout << line_id << ", ";
	}
	cout << endl;
    }

    // Build the contigs from the chosen lines and write to a fasta file
    cout << "Writing lines as contigs" << endl;
    Ofstream( out, fasta_file );
    String gap_str(GAP_SIZE, 'N');
    int empty_lines = 0;
    int contig_id = 0;
    for (int line_id : chosen) {
	// flatten the line into a single contig
	vec<int> path = flatten_line(lines[line_id]);
	String fasta = gapped_fasta(path, edges, K, gap_str);
	// Ignore empty lines (and warn)
	if (fasta.size() == 0) 
	    empty_lines++;
	// Only keep those contigs greater than 1000 bases long
	else if (fasta.isize() >= MIN_LENGTH) {
	    if (DEBUG)
		cout << contig_id++ << " (" << line_id << "): " << printSeq(path) << endl;
	    // Write the contig to fasta file
	    out << ">" << ToString(contig_id) << ":" << line_id << "\n";
	    out << fasta << "\n";
	    contig_id++;
	    
	}
    }

    if (empty_lines != 0)
	cout << "WARNING: found " << empty_lines << " empty lines" << endl;
    cout << "Flattened " << contig_id << " lines of size " << MIN_LENGTH << " or greater" << endl;
	    
    cout << Date( ) << ": Done." << endl;
}

    
