/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef SHOW_ALIGNMENT_PILE_TOOLS_H
#define SHOW_ALIGNMENT_PILE_TOOLS_H

#include "Basevector.h"
#include "Bitvector.h"
#include "MainTools.h"
#include "lookup/LookAlign.h"
#include "lookup/AlignCollector.h"

class GenomeCoordinateRange {
	public:
		GenomeCoordinateRange(unsigned int contig = 0, unsigned int start = 0, unsigned int end = 0) : _contig(contig), _start(start), _end(end) {}
		bool operator<(const GenomeCoordinateRange &other) const {
			return (_contig < other._contig || _start < other._start || _end < other._end);
		}

		unsigned int _contig;
		unsigned int _start;
		unsigned int _end;
};

class PileOfAlignments {
	public:
		PileOfAlignments(GenomeCoordinateRange newg, basevector &newref) : g(newg), ref(newref) {
			// The grid for this set of coordinates has as many elements as the section of the reference.
			grid.resize(ref.size());
			for (unsigned int i = 0; i < ref.size(); i++) {
				grid[i] = 1;
			}
		}

        void AddAlignment(look_align &newla, basevector &newread);

        void Print(ofstream &stream, bool as_html = 0, vecbitvector basesUsed = vecbitvector(0), String display_unused = "dot");

	private:
		GenomeCoordinateRange g;
		basevector &ref;
		vec< pair<look_align, basevector> > alignments_and_reads;
		vec<int> grid;

        pair<unsigned int, unsigned int> BuildEffectiveAlignment(look_align & alignment, basevector & bv);

};

MaxErrDiffAlignCollector LoadAlignments(String qltoutfile, unsigned int size, int max_errors = 4);

vec<GenomeCoordinateRange> LoadCoordinates(String COORDINATES, unsigned int WINDOW_SIZE, vec< pair<bool, unsigned int> > & specs );

#endif
