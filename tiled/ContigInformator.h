// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

#include "CoreTools.h"
#include "ReadLocation.h"
#include "ReadPairing.h"
#include "lookup/LookAlign.h"

#include <functional>

//  Calculate the read coverage of a contig in consecutive, non-overlapping 
//  'window'-sized regions.  Display window (in bases along the contig) and coverage.
//  coverage calculated by #bases-in-window/window size
//
void ContigCoverage( vec<read_location> &vecLocs,
		     int contig,
		     int window );


// For UNIQUELY PLACED reads only,
// print info about a read if its partner should be in the same contig but isn't
// info includes (for both the read and its partner): 
//          placement on human
//          read location from assembly
//          separation, stdev from pairs file
//
// a read partner should be in the same contig as a read if the distance from the
// read to the end of the contig (orientation dependent) is > sep + 4*std for the 
// insert.
void ReadAndPartnerInfo( vec<read_pairing> &pairs,
			 vec<int> &pairs_index,
			 vec<read_location> &locs,
			 vec<int> &locs_by_id,
			 vec<genome_pos> &uniq_human,
			 vec<int> &individual,
			 vecString &read_names,
			 vec<int> &read_lens,
			 int contig );



// return -1 if partner not placed (error)
// return 0 if partner not placed in same contig but should be
// return 1 if partner placed in same contig
//			 
int IsPartnerPlaced( read_location &loc,
		     vec<int> &locs_by_id,
		     vec<read_pairing> &pairs,
		     vec<int> &pairs_index );
