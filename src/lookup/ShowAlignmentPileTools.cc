/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "Bitvector.h"
#include "CoreTools.h"
#include "lookup/LookAlign.h"
#include "lookup/AlignCollector.h"
#include "lookup/ShowAlignmentPileTools.h"

#include <map>

bool EarlierAlignment(const pair< look_align, basevector > &a, const pair< look_align, basevector > &b) {
	return a.first.pos2() < b.first.pos2();
}

class GridModifier {
	public:
		GridModifier(int gapsize, unsigned int begin, unsigned int end) : _gapsize(gapsize), _begin(begin), _end(end) {}
		int _gapsize;
		unsigned int _begin;
		unsigned int _end;
};

void PileOfAlignments::AddAlignment(look_align &newla, basevector &newread) {
			// If we're given a read with an insertion, expand the grid so that the inserted base appears in the space
			if (newla.a.Nblocks() > 1) {
				for (int block = 0; block < newla.a.Nblocks(); block++) {
					if (newla.a.Gaps(block) < 0) {
						int elementsize = abs(newla.a.Gaps(block)) + 1;
						int gapplacement = newla.pos2() - g._start - 1;

						// Figure out exactly where the space needs to go
						for (int subblock = 0; subblock < block; subblock++) {
						    gapplacement += newla.a.Lengths(subblock);
						    if  ( newla.a.Gaps(subblock) > 0 ) gapplacement += newla.a.Gaps(subblock);
						}

						// If we've already expanded this region and the old expansion is smaller than the one we're about to perform, overwrite it
						if (grid[gapplacement] < elementsize) {
							grid[gapplacement] = elementsize;
						}
					}
				}
			}

			// Store our info (in an easily sortable form for later)
			pair< look_align, basevector > a_and_r(newla, newread);
			alignments_and_reads.push_back(a_and_r);
		}

void PileOfAlignments::Print(ofstream &stream, bool as_html, vecbitvector basesUsed, String display_unused) {
			// Sort by start position
			sort(alignments_and_reads.begin(), alignments_and_reads.end(), EarlierAlignment);

			// First, print out the reference in a nice grid that makes room for insertions in the read (looks like gap in reference)
			int refgridsize = 0;

			vec<String> refgrid(grid.size());
			for (unsigned int gridpos = 0; gridpos < grid.size(); gridpos++) {
				refgridsize += grid[gridpos];

				refgrid[gridpos] = as_base(ref[gridpos]);

				for (int element = 1; element < grid[gridpos]; element++) {
					refgrid[gridpos] += "+";
				}
			}

			if (as_html) { stream << "<tr>"; }
			for (unsigned int gridpos = 0; gridpos < grid.size(); gridpos++) {
				stream << (as_html ? "<td>" : "") << refgrid[gridpos] << (as_html ? "</td>" : "");
			}
			stream << (as_html ? "<td>" : "\t") << "(" << (g._contig) << ":" << g._start << "-" << (g._start + refgrid.size()) << " reference)" << (as_html ? "</td>" : "") << endl;
			if (as_html) { stream << "</tr>"; }

			// Now print the reads into the grid
			for (unsigned int readnum = 0; readnum < alignments_and_reads.size(); readnum++) {

			        // these are the original read and its alignment in thier entiriety
				basevector bv = alignments_and_reads[readnum].second;
				look_align alignment = alignments_and_reads[readnum].first; 
				const look_align & orig_alignment = alignments_and_reads[readnum].first; 

				if (alignment.rc1) { bv.ReverseComplement(); }

				// let's take care of the alignments that extend beyond the reference stretch shown.

				pair<unsigned int, unsigned int> deleted_bases(0,0); // how many bases on each end of the basevector were overhanging from the grid and got chopped off?

				if ( static_cast<unsigned int>(alignment.StartOnTarget()) < g._start || 
				     static_cast<unsigned int>(alignment.EndOnTarget()) > g._end ) {
				  deleted_bases=BuildEffectiveAlignment(alignment,bv);
				}

				// This is how much left-sided buffer will we need
				int refoffset = alignment.StartOnTarget() - g._start;

				// Start by putting everything in the grid
				vec<String> bvgrid(grid.size(), " ");
				for (unsigned int cycle = 0, gridpos = refoffset; gridpos < bvgrid.size() && cycle < bv.size(); cycle++, gridpos++) {
					bvgrid[gridpos] = as_base(bv[cycle]);

					// if this base was not used in coverage... lowercase it
					if (basesUsed.size() > 0) {
					  // we want read offset on the original read to be able to retrieve the used status,
					  // while our vector bv may have been chopped if it was overhanging from the grid!
					  // here we take care of the possibility that some bases were dropped:
					  int readOffset = !alignment.rc1 ? (cycle+deleted_bases.first) : (bv.size() + deleted_bases.second  - cycle - 1);
					  if ( basesUsed[alignment.QueryId()][readOffset] == 0) {
					    if (display_unused == "dot") {
					      bvgrid[gridpos] = ".";
					    } else if (display_unused == "lower") {
					      bvgrid[gridpos] = ToLower(bvgrid[gridpos]);
					    }
					  }
					} 
				}

				vec<GridModifier> gms;

				for (int block = 0, offset = refoffset; block < alignment.a.Nblocks(); 
				     offset += alignment.a.Lengths(block), block++) {

					if (alignment.a.Gaps(block) != 0) {
						gms.push_back(GridModifier(alignment.a.Gaps(block), 
									   offset, 
									   offset + abs(alignment.a.Gaps(block))
									   )
							      );
						// additional shift: before we merged read's
						// bases appropriately, each insertion gap introduces
						// additional offset in bvgrid that we need to take care of
						if ( alignment.a.Gaps(block) < 0 ) offset -= alignment.a.Gaps(block);
					}
				}

				for (int gmindex = gms.size() - 1; gmindex >= 0; gmindex--) {
					GridModifier gm = gms[gmindex];

					if (gm._gapsize > 0) {
					  bvgrid.insert(bvgrid.begin()+gm._begin, gm._gapsize, "-");
					} else if (gm._gapsize < 0) {
					  vec<String>::iterator position = bvgrid.begin()+gm._begin - 1;
					  for (vec<String>::iterator it = bvgrid.begin() + gm._begin; 
					       it != bvgrid.begin()+gm._end; it++) {
							(*position) += (*it);
					  }
					}
				}

				for (int gmindex = gms.size() - 1; gmindex >= 0; gmindex--) {
					GridModifier gm = gms[gmindex];

					if (gm._gapsize < 0) {
					  bvgrid.erase(bvgrid.begin()+gm._begin, bvgrid.begin()+gm._end);
					}
				}

				unsigned int oldsize = bvgrid.size();
				bvgrid.resize(grid.size());
				for (unsigned int pos = oldsize; pos < bvgrid.size(); pos++) {
					bvgrid[pos] = " ";
				}

				// Print out the finalized read grid
				if (as_html) { stream << "<tr>"; }
				for (int gridpos = 0; gridpos < grid.isize(); gridpos++) {
					     String estr(grid[gridpos], (gridpos < refoffset + bv.isize() && gridpos >= refoffset ? '+' : ' '));
                    String element(estr.c_str());
					element.replace(0, bvgrid[gridpos].size(), bvgrid[gridpos]);
					stream << (as_html ? "<td>" : "") << element << (as_html ? "</td>" : "") << flush;
				}

				// Print a little meta information (on the right-hand side so that we don't have to do extra work to get everything to line up nicely)
				stream << (as_html ? "<td>" : "\t") << "(" << (orig_alignment.target_id) << ":" << orig_alignment.pos2() << "-" << orig_alignment.Pos2() 
				       << " " << orig_alignment.query_id << " " << (orig_alignment.StartOnTarget() - (int)g._start) << " " << (alignment.rc1 ? "RC" : "FW");
				for (int block = 0; block < orig_alignment.a.Nblocks(); block++) {
					stream << " [" << block << " " << orig_alignment.a.Gaps(block) << " " << orig_alignment.a.Lengths(block) << "]";
				}
				stream << ")" << (as_html ? "</td>" : "") << endl;
				if (as_html) { stream << "</tr>"; }
			}
		}


pair<unsigned int, unsigned int> PileOfAlignments::BuildEffectiveAlignment(look_align & alignment, basevector & bv) {

		 unsigned int bases_to_chop_left = 0;
		 unsigned int bases_to_chop_right = 0;
		 int block = 0;

		 if ( static_cast<unsigned int>(alignment.StartOnTarget()) < g._start ) {
		    look_align effective_alignment = alignment;
                    unsigned int block_start = effective_alignment.StartOnTarget(); // start of the alignment block on target reference
		    
		    // lets find where in the alignment (in which block) we have the start of the displayed genomic stretch:
		    for ( ; block < effective_alignment.a.Nblocks() ; block++ ) {
		        if ( effective_alignment.a.Gaps(block) < 0 ) {
			    cout << "Sorry, reads that span outside the displayed reference stretch and have insertions can not be dealt with yet" << endl;
			    exit(1);
			}

			block_start += effective_alignment.a.Gaps(block); // jump over the deletion on the read
		      
			// now block_start points to the correct position on the whole genome reference, where the current block starts
			if ( block_start + effective_alignment.a.Lengths(block) > g._start ) break; // we located the block spanning over the leftmost end of the displayed reference stretch

			block_start += effective_alignment.a.Lengths(block);
			bases_to_chop_left += effective_alignment.a.Lengths(block); // current block in its entiriety is outside of the displayed reference chunk and have to be chopped off
		    }

		    unsigned int block_offset = 0; // by default, we are going to copy all blocks block+1, ... into blocks 1,...

		    // first block (0) requires special treatment
		    if ( block_start <= g._start ) {
		        unsigned int overhang = g._start - block_start; // block overhangs the reference by this many bases
			// block starts to the left of the displayed reference's start; 
			// just update the first alignment block (being first, it's gap is already 0)

			effective_alignment.a.SetLength(0, effective_alignment.a.Lengths(block) - overhang);
			bases_to_chop_left += overhang; // we also will need to chop those overhanging bases

		    } else {
		        // that's a trickier case: displayed reference stretch starts in the middle of the gap
		        // on the original alignment. Hack: use first block of length 0 to enforce drawing the gap

		        effective_alignment.a.SetLength(0,0);
			effective_alignment.a.SetGap(1,block_start-g._start); // if g._start falls into a gap, we are guaranteed that there are at least 2 blocks!
			effective_alignment.a.SetLength(1,effective_alignment.a.Lengths(block));
			block_offset = 1; // we used blocks 0 and 1 already, got to copy remaining blocks into blocks 2,...
		    }

		    // we took care of the block that got us to g._start and copied it manually; now let's copy remaining blocks unchanged:
		    for ( int i = block+1 ;  i < effective_alignment.a.Nblocks() ; i++ ) {

		        // if block_offset=0 (i.e. we had block_start < g._start, we just copy back block+1 --> 1, block+2 --> 2, etc (if block==0, we waste time copying 1-->1, 2-->2 etc, but no harm!)
		        // if block offset=1 (i.e. g._start was in the gap, we already set the fake block of length 0, then we used block 1 to
		        // copy data from block=block, now we copy block+1 --> 2, block+2-->3, etc. Note that since g._start fell into the gap, we are
		        // guaranteed that block>=1, so in the worst case block=1 we copy values onto themselves (2-->2, 3-->3, etc); if block > 1 we do copy back.
		      effective_alignment.a.SetGap(i-block+block_offset, effective_alignment.a.Gaps(i));
		      effective_alignment.a.SetLength(i-block+block_offset, effective_alignment.a.Lengths(i));
		    }

		    effective_alignment.query_length -= bases_to_chop_left;
		    effective_alignment.a.SetNblocks( effective_alignment.a.Nblocks() - block+block_offset );
		    effective_alignment.SetStartOnTarget(g._start);
		    
		    basevector b;
		    b.SetToSubOf(bv,bases_to_chop_left,bv.size()-bases_to_chop_left);
		    // done with left side
		    alignment = effective_alignment;
		    bv = b;
		 }
		  // **********************************
		  // now let's take care of the right end of the alignment possibly protruding out of the reference grid:
		  // **********************************
		 
		 if ( static_cast<unsigned int>(alignment.EndOnTarget()) > g._end ) {
		 
		    look_align effective_alignment = alignment;

		    unsigned int block_end = effective_alignment.EndOnTarget(); // one past end of the alignment block on target reference
		    block = effective_alignment.a.Nblocks()-1;

		    // lets find where in the alignment (in which block) we have the end of the displayed genomic stretch:
		    for ( ; block >= 0  ; block-- ) {
		        if ( effective_alignment.a.Gaps(block) < 0 ) {
			    cout << "Sorry, reads that span outside the displayed reference stretch and have insertions can not be dealt with yet" << endl;
			    exit(1);
			}

			if ( block_end - effective_alignment.a.Lengths(block) < g._end ) break; // we located the block spanning over the rightmost end of the displayed reference stretch

			block_end -= effective_alignment.a.Gaps(block); // jump over the deletion on the read
			block_end -= effective_alignment.a.Lengths(block); // jump over the deletion on the read

			bases_to_chop_right += effective_alignment.a.Lengths(block); // current block in its entiriety is outside of the displayed reference chunk and have to be chopped off
		    }

		    effective_alignment.query_length -= bases_to_chop_right;

		    // first block (0) requires special treatment
		    if ( block_end >=  g._end ) {
		        unsigned int overhang = block_end - g._end; // block overhangs the reference by this many bases

			// block ends to the left of the displayed reference's end; 
			// just update the block's length:

			effective_alignment.a.SetLength(block, effective_alignment.a.Lengths(block) - overhang);
			bases_to_chop_right += overhang; // we also will need to chop those overhanging bases
			effective_alignment.a.SetNblocks(block+1); // delete all blocks after the one we found to cover reference grid end
				    
		    } else {
		        // that's a trickier case: displayed reference stretch ends in the middle of the gap
		        // on the original alignment. Hack: use last block of length 0 to enforce drawing the gap

		        effective_alignment.a.SetNblocks(block+2); // we keep the block we found (last before grind end) unchanged and add additional block
			effective_alignment.a.SetLength(block+1,0); // that extra block has zero length
			effective_alignment.a.SetGap(block+1,g._end - block_end); // draw a gap till the end of the grid
		    }

		    basevector b;
		    b.SetToSubOf(bv,0,bv.size()-bases_to_chop_right);
		    
		    // done with left side
		    alignment = effective_alignment;
		    bv = b;
		 }
		 return make_pair(bases_to_chop_left, bases_to_chop_right);
	       }

MaxErrDiffAlignCollector LoadAlignments(String qltoutfile, unsigned int size, int max_errors) {
	MaxErrDiffAlignCollector aligns(0, max_errors, size);
	LoadLookAligns(qltoutfile, aligns);
	aligns.Consolidate();

	return aligns;
}

bool IsArray(String str) {
	return (str.Contains("{") && str.Contains("}")) ? 1 : 0;
}

GenomeCoordinateRange ParseCoordinate(String coordinate, unsigned int WINDOW_SIZE, pair<bool, unsigned int> & spec) {
	unsigned int contig = 0, start = 0, end = 0;

	contig = coordinate.Before(":").Int();

	if (coordinate.Contains("-")) {
		start = coordinate.After(":").Before("-").Int();
		end   = coordinate.After("-").Int();
		spec.first=false; // no, the complete interval was specified
	} else {
		start = coordinate.After(":").Int() ;
		end   = coordinate.After(":").Int() + WINDOW_SIZE;
		spec.first=true; // yes, a single coordinate was specified
		spec.second = start;
		start -= WINDOW_SIZE;
	}

	GenomeCoordinateRange g(contig, start, end);

	return g;
}

vec<GenomeCoordinateRange> LoadCoordinates(String COORDINATES, unsigned int WINDOW_SIZE, vec< pair<bool, unsigned int> > & specs ) {
	vec<GenomeCoordinateRange> coordinates;
 
	pair<bool, unsigned int> spec;

	if (IsRegularFile(COORDINATES) || IsArray(COORDINATES)) {
		String coordStr;
		vec<String> coordStrs;

		if (IsRegularFile(COORDINATES)) {
			string coordstr;
			ifstream coordfile(COORDINATES.c_str());

			while (getline(coordfile, coordstr)) {
				coordStr = coordstr;
				if (!coordStr.Contains("#")) {
					coordStrs.push_back(coordStr);
				}
			}
			coordfile.close();
		} else if (IsArray(COORDINATES)) {
			coordStr = COORDINATES.Between("{", "}");
			Tokenize(coordStr, ',', coordStrs);
		}

		for (unsigned int i = 0; i < coordStrs.size(); i++) {
		  coordinates.push_back(ParseCoordinate(coordStrs[i], WINDOW_SIZE, spec));
		  specs.push_back(spec);
		}
	} else {
	  coordinates.push_back(ParseCoordinate(COORDINATES, WINDOW_SIZE,spec));
	  specs.push_back(spec);
	}

	return coordinates;
}
