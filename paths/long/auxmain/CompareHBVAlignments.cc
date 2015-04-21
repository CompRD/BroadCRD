///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

const char *DOC =
  "Compare and visualize readpaths on a hyperbasevector.";

// Requires an HBV, fastb, qualb and one or two paths files.
// Will compare two readpath files, or one readpath file against internally
// computed ReadPlaces0.
// Will display all readpaths, or just a single one as specified by READ.
// You can filter the paths to display using FILTER_MATCHING and FILTER_SINGLE.
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


#include "Basevector.h"
#include "Intvector.h"
#include "MainTools.h"
#include "Qualvector.h"
#include "paths/HyperBasevector.h"
#include "paths/long/EvalByReads.h"
#include "paths/long/PlaceReads0.h"
#include "paths/long/ReadPath.h"
#include "paths/long/ReadPathTools.h"


bool ValidatePath(const HyperBasevector& hbv, const vec<int>& to_left,
		  const vec<int>& to_right, const int offset,
		  const vec<int>& edge_list) {

    String message;
    bool valid = ValidateReadPath(hbv, to_left, to_right, offset, edge_list,
				  message);
    cout << message << endl;
    return valid;
}


int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandDoc(DOC);
     CommandArgument_String_Doc(READ_HEAD,
      "Base and quality score files: READ_HEAD.{fastb,qualb}");
     CommandArgument_String_Doc(PATHS1, 
      "ReadPath file containing paths on the HBV");
     CommandArgument_String_Doc(PATHS2, 
      "ReadPath file containing paths on the HBV");
     CommandArgument_String_Doc(HBV,
      "HyperBasevector (HBV) file");
    CommandArgument_String_OrDefault_Doc(READ_IDS, "",
      "Optional list of read IDs to display.");
     CommandArgument_Bool_OrDefault_Doc(SHOW_EDGES,False, 
      "Show each overlapping edge in a path");
     CommandArgument_Bool_OrDefault_Doc(SHOW_RULER,False, 
      "Show a ruler counting bases");
     CommandArgument_Bool_OrDefault_Doc(FILTER_SINGLES,False, 
      "Do not display paths that contain only a single edge");
     CommandArgument_Bool_OrDefault_Doc(FILTER_MATCHING, True,
      "Do not display paths that are identical");
     CommandArgument_String_OrDefault_Valid_Doc(DISPLAY_MODE, "both", "{first,second,both}",
      "Display visualizations of first, second or both sets of paths");
     EndCommandArguments;

     // Parse display mode
     bool view_first = (DISPLAY_MODE == "first" || DISPLAY_MODE == "both");
     bool view_second = (DISPLAY_MODE == "second" || DISPLAY_MODE == "both");
     bool use_places = (PATHS2 == "");

     // Import hbv and read paths from ReadQGrapher

     String hbv_file = HBV;
     String paths1_file =  PATHS1;
     String paths2_file =  PATHS2;

     HyperBasevector hbv;
     BinaryReader::readFile( hbv_file, &hbv );
     ReadPathVec paths1, paths2;
     paths1.ReadAll(paths1_file);
     if (paths2_file != "")
       paths2.ReadAll(paths2_file);

     vec<int> to_left, to_right;
     hbv.ToLeft(to_left);
     hbv.ToRight(to_right);
     
     // Write hbv as dot file

     Ofstream( out, hbv_file + ".dot" );
     hbv.PrintSummaryDOT0w( out, True, False, True );

     // Load reads.

     vecbasevector bases;
     vecqualvector quals;
     bases.ReadAll( READ_HEAD + ".fastb", True );
     quals.ReadAll( READ_HEAD + ".qualb", True );

     // Select subset of paths to display
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

     // Find alignments using PlaceReads0
     vec< vec<read_place> > read_places( bases.size( ) );
     if (use_places)
       PlaceReads0( hbv, bases, quals, read_places );

     // Some stats
     int empty_match = 0;
     int non_empty_match = 0;
     int paths1_found = 0;
     int paths2_found = 0;

     // Compare ReadPaths
     for ( size_t i = 0; i < ids.size( ); i++ ) {
	size_t read_id = ids[i];

       ReadPath tmp;
       if (use_places) {      
	 read_place place = (read_places[read_id].size() != 0 ?  read_places[read_id][0] : read_place());
	 tmp.setOffset(	place.P());
	 for (auto& j : place.E()) 
	   tmp.push_back(j);
       }
       
       ReadPath& read_path1 = paths1[read_id];
       ReadPath& read_path2 = (use_places ? tmp : paths2[read_id]);

       bool found_path1 = read_path1.size();
       bool found_path2 = read_path2.size();

       // Compute stats
       if (!found_path2 && !found_path1)
	 empty_match++;
       if (found_path2 && found_path1 && read_path1 == read_path2 )
	 non_empty_match++;
       if (found_path2)
	 paths2_found++;
       if (found_path1)
	 paths1_found++;
       
       // apply filtereing
       if (FILTER_SINGLES && paths1[read_id].size() < 2 && paths2[read_id].size() < 2) continue;
       if (FILTER_MATCHING && !found_path2 && !found_path1) continue;
       if (FILTER_MATCHING && found_path2 && read_path1 == read_path2 ) continue;

	 int offset1 = read_path1.getOffset();
	 int offset2 = read_path2.getOffset();
	 vec<int> edge_list1(read_path1.begin(), read_path1.end());
	 vec<int> edge_list2(read_path2.begin(), read_path2.end());


	 cout << "[" << read_id << "]";
	 cout << "  hbv_path1= " << offset1 << ":"<< printSeq( edge_list1 );
	 cout << "  hbv_path2= " << offset2 << ":"<< printSeq( edge_list2 );
	 cout << endl << endl;;

	 if (view_first && found_path1) {
	     if (ValidatePath(hbv, to_left, to_right, offset1, edge_list1 ))
		 DisplayReadPath( cout, hbv, to_left, to_right, offset1, edge_list1,
				  bases[read_id], quals[read_id], SHOW_EDGES, SHOW_RULER);
	 }
	 else if (view_first)
	   cout << "No valid first hbv path found." << endl;
	 
	 cout << endl;

	 if (view_second && found_path2) {
	     if (ValidatePath(hbv, to_left, to_right, offset2, edge_list2)) 
		 DisplayReadPath( cout, hbv, to_left, to_right, offset2, edge_list2,
				  bases[read_id], quals[read_id], SHOW_EDGES, SHOW_RULER);
	 } else if (view_second)
	   cout << "No valid second hbv path found." << endl;
	 
	 cout << endl << endl;
     }
     
     cout << "Graph edges : " << hbv.EdgeObjectCount() << endl;
     cout << "Reads       : " << bases.size() << endl;
     cout << "Matching    : " << non_empty_match << "  (non_empty)" << endl;
     cout << "Matching    : " << empty_match << "  (empty)" << endl;
     cout << "Conflicting : " << bases.size() - empty_match - non_empty_match << "" << endl;
     cout << "HBV paths   : " << paths1_found << endl;
     cout << "Places      : " << paths2_found << endl;
}
