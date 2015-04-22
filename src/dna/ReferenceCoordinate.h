#ifndef _REFERENCE_COORDINATE_H
#define _REFERENCE_COORDINATE_H

#include "Basevector.h"
#include "feudal/FeudalFileReader.h"

#include <cstring>
#include <fstream>
#include <sstream>

// Forward declaration for the mapper between linear and genomic coordinates
class ReferenceCoordinateMapper;

/**
    Class: ReferenceCoordinate

    Class to represent a genomic coordinate (contig:offset) or linear coordinate, and easily convert between the two.  This can be used as a loop variable, making it very easy to access a reference position and a LargeSparseArray position.
*/

struct ReferenceCoordinate {

    ReferenceCoordinate() : linear_offset(0), contig(0), offset(0), rcm(0) {}

    void Shift(long shift);

    void operator++(int) { Shift(1); }
    void operator++() { Shift(1); }

    void operator--(int) { Shift(-1); }
    void operator--() { Shift(-1); }

    void operator+=(long shift) { Shift(shift); }
    void operator-=(long shift) { Shift(-shift); }

    bool operator<(ReferenceCoordinate g) const { return (linear_offset < g.linear_offset); }
    bool operator<=(ReferenceCoordinate g) const { return (linear_offset <= g.linear_offset); }

    bool operator>(ReferenceCoordinate g) const { return (linear_offset > g.linear_offset); }
    bool operator>=(ReferenceCoordinate g) const { return (linear_offset >= g.linear_offset); }

    bool operator==(ReferenceCoordinate g) const { return (linear_offset == g.linear_offset); }
    bool operator!=(ReferenceCoordinate g) const { return (linear_offset != g.linear_offset); }

    String ToString()
    {
        return ::ToString(this->contig) + ":" + ::ToString(this->offset);
    }

    friend ostream &operator<<(ostream &stream, const ReferenceCoordinate &gc) {
        stream << gc.contig << ":" << gc.offset;
        return stream;
    }

    unsigned int linear_offset;
    int contig;
    unsigned int offset;
    const ReferenceCoordinateMapper *rcm;
};

inline bool EarlierCoordinate(const ReferenceCoordinate &a, const ReferenceCoordinate &b) {
    return a.linear_offset < b.linear_offset;
}

/**
    Class: ReferenceRange

    Represents a reference range.  Can be interpreted and displayed as a linear range (A-B) or a genomic range (C1:O1-C2:O2).
*/

struct ReferenceRange {

    String ToString()
    {
        return this->start.ToString() + "-" + this->stop.ToString();
    }

    friend ostream &operator<<(ostream &stream, const ReferenceRange &rr) {
        stream << rr.start << "-" << rr.stop;
        return stream;
    }

    String AsLinearRange(void) const {
        ostringstream ss;
        ss << start.linear_offset << "-" << stop.linear_offset;

        return ss.str();
    }

    String AsGenomicRange(void) const {
        return ::ToString(start.contig) + ":" + ::ToString(start.offset) + "-" + ::ToString(stop.contig) + ":" + ::ToString(stop.offset);
    }

    /// Extracts bases specified by this range from the reference <ref> into the basevector <b>.
    /// Previous content of b, if any, will be destroyed
    void ExtractBases(const vecbasevector & ref, basevector &b) {
      b.resize(0);
      for ( int i = start.contig; i <= stop.contig ; i++ ) {
        basevector t ;
        unsigned int start_on_contig, stop_on_contig;
        start_on_contig =  ( ( i == start.contig ) ? start.offset : 0 );
        stop_on_contig =  ( ( i == stop.contig ) ? stop.offset : ref[i].size() );
        t.SetToSubOf( ref[i], start_on_contig, stop_on_contig - start_on_contig );
        b = Cat(b,t);
      }
    }

    /// Computes count of GC bases inside this range.
    ///
    /// @par ref the actual reference sequence on which this range is specified
    unsigned int GcBases(const vecbasevector & ref) {

      unsigned int gc = 0;

      for ( int i = start.contig; i <= stop.contig ; i++ ) {
        unsigned int start_on_contig, stop_on_contig;
        start_on_contig =  ( ( i == start.contig ) ? start.offset : 0 );
        stop_on_contig =  ( ( i == stop.contig ) ? stop.offset : ref[i].size() );
        gc += ref[i].GcBases(start_on_contig, stop_on_contig );
      }
      return gc;
    }

    ReferenceCoordinate start;
    ReferenceCoordinate stop;
};


inline String ToString(const vec<ReferenceRange> &vrr) {
    String ranges = "{";
    for (unsigned int rangeindex = 0; rangeindex < vrr.size(); rangeindex++) {
        ranges += vrr[rangeindex].AsLinearRange();

        if (rangeindex < vrr.size() - 1) {
            ranges += ",";
        }
    }

    ranges += "}";

    return ranges;
}

inline bool EarlierRange(const ReferenceRange &a, const ReferenceRange &b) {
    return a.start.linear_offset < b.start.linear_offset;
}

/**
    Class: ReferenceCoordinateMapper

    Takes a reference and allows the easy conversion between linear and genomic coordinates.  Also parses coordinate ranges.

    If the reference is specified to the class as just a filename, the full reference will not be loaded into memory.  Rather, the fastb control block and contig sizes will be read directly.

    Be aware: the conversion of contig:offset to linear_offset is more efficient as it is just an array lookup and an addition.  In contrast, the conversion of a linear_offset to a contig:offset requires a binary search over an array to find the right offset range and then a subsequent subtraction.  Keep this in mind when writing a program that loops over a range that needs to be converted to the alternative format in the inner loop.
*/

class ReferenceCoordinateMapper {
    public:
        /// Loads the information from a reference file required to efficiently map from one coordinate system to another.  This method is memory-efficient: the reference is *not* loaded.
        ReferenceCoordinateMapper(String referencefile) {
            unsigned int size;
            FeudalFileReader rdr(referencefile.c_str());

            unsigned long nnn = rdr.getNElements();
            contigoffsets.resize(nnn);
            linear_size = 0;
            for ( unsigned long iii = 0; iii < nnn; ++iii )
            {
                memcpy(&size,rdr.getFixedData(iii,sizeof(size)),sizeof(size));
                contigoffsets[iii] = linear_size;
                linear_size += size;
            }
        }

	/// Default ctor
	ReferenceCoordinateMapper() {}

        /// Loads the information from a reference file required to efficiently map from one coordinate system to another.  This method requires the reference to have already been loaded.
        ReferenceCoordinateMapper(const vecbasevector &reference) {
	  LoadReference(reference);
        }

	/// Loads the information from the specified pre-loaded reference; the mapper will be reinitialized, old reference information, if any,
        /// will be lost.
	void LoadReference(const vecbasevector &reference) {
            contigoffsets.resize(reference.size());

            linear_size = 0;
            for (size_t contig = 0; contig < reference.size(); contig++) {
                contigoffsets[contig] = linear_size;
                linear_size += reference[contig].size();
            }
        }

        /// Converts a contig and offset to a linear_offset and returns it via the reference variable.
        void LinearOffset(unsigned int contig, unsigned int offset, unsigned int &linear_offset) const {
            linear_offset = contigoffsets[contig] + offset;

            // We're lax in the actual parsing of ranges where we rationalize user input to reasonable numbers.  However, in the program itself, we're strict on the allowable coordinates we allow.
            Assert(linear_offset <= linear_size);
        }

        /// Converts a linear_offset to a contig and offset and returns them via the reference variables.
        void ContigAndOffset(unsigned int linear_offset, int &contig, unsigned int &offset) const {
            bool found = 0;
            unsigned int start = 0;
            unsigned int stop = contigoffsets.size();
            unsigned int probe;
            unsigned int contig_offset;

            do {
                probe = start + (stop - start)/2;
                contig_offset = contigoffsets[probe];

                if      (contig_offset == linear_offset || probe == start || probe == stop) { found = 1; }
                else if (contig_offset <  linear_offset) { start = probe; }
                else if (contig_offset >  linear_offset) { stop  = probe; }
            } while (!found);

            // We're lax in the actual parsing of ranges where we rationalize user input to reasonable numbers.  However, in the program itself, we're strict on the allowable coordinates we allow.
            contig = probe;
            Assert(contig < contigoffsets.isize());

            offset = (linear_offset - contigoffsets[probe]);
            Assert(offset <= NumBasesInContig(contig));
        }

        /// Returns a ReferenceCoordinate object given a linear_offset.  This is less efficient than using the operator() method.
        ReferenceCoordinate operator[](unsigned int linear_offset) const {
            ReferenceCoordinate gc;

            gc.linear_offset = linear_offset;
            ContigAndOffset(linear_offset, gc.contig, gc.offset);
            gc.rcm = this;

            return gc;
        }

        /// Returns a ReferenceCoordinate object given a contig and offset.  This is more efficient than using the operator[] method.
        ReferenceCoordinate operator()(int contig, unsigned int offset) const {
            ReferenceCoordinate gc;

            gc.contig = contig;
            gc.offset = offset;
            LinearOffset(contig, offset, gc.linear_offset);
            gc.rcm = this;

            return gc;
        }

        /// Parses a single range.  Allowable forms are L1-L2, C: (for full chromosome), C1:O1-C2:O2, C1:O1-O2, L1, L1+w, C1:O1, C1:O1+w,
        /// or any string not containing a dash (interpreted as the entire reference).
        ReferenceRange ParseRangeString(const String & range, unsigned int window_size = 50) const {
            ReferenceRange gr;

            if (range.Contains("-")) {
                String start = range.Before("-");
                String stop  = range.After("-");

                if (start.Contains(":")) {
                    gr.start.contig = start.Before(":").Int();
                    gr.start.offset = start.After(":").Int();

                    if (stop.Contains(":")) {
                        gr.stop.contig = stop.Before(":").Int();
                        gr.stop.offset = stop.After(":").Int();
                    } else {
                        gr.stop.contig = gr.start.contig;
                        gr.stop.offset = stop.Int();
                    }

                    // Do a bunch of range checking.  If the specified range is out of bounds, try to rationalize it
                    // to the end of the contig or the end of the reference.  The only exception is the start contig,
                    // which *must* be an accessible contig.
                    Assert(gr.start.contig < contigoffsets.isize());
                    if (gr.start.offset > NumBasesInContig(gr.start.contig)) {
                        gr.start.offset = NumBasesInContig(gr.start.contig);
                    }

                    if (gr.stop.contig > contigoffsets.isize()) {
                        gr.stop.contig = contigoffsets.isize() - 1;
                        gr.stop.offset = NumBasesInContig(gr.stop.contig);
                    } else if (gr.stop.offset > NumBasesInContig(gr.stop.contig)) {
                        gr.stop.offset = NumBasesInContig(gr.stop.contig);
                    }

                    LinearOffset(gr.start.contig, gr.start.offset, gr.start.linear_offset);
                    LinearOffset(gr.stop.contig, gr.stop.offset, gr.stop.linear_offset);
                } else {
                    gr.start.linear_offset = start.Int();

                    // Allow an offset that's out of range to be rationalized back to the end of the reference
                    gr.stop.linear_offset = (stop.Int() < NumBases()) ? stop.Int() : NumBases();

                    ContigAndOffset(gr.start.linear_offset, gr.start.contig, gr.start.offset);
                    ContigAndOffset(gr.stop.linear_offset, gr.stop.contig, gr.stop.offset);
                }
            } else if (range.Contains(":")) {
                unsigned int window = window_size;
                unsigned int startoffset;
                bool full_chr = false;
                if (range.Contains("+")) {
                    window = range.After("+").Int();
                    startoffset = range.After(":").Before("+").Int();
                } else {
                    String offset_str = range.After(":");
                    if ( offset_str.empty() ) {
                      startoffset = 0;
                      full_chr = true; // full chromosome requested using "C:" format
                    }
                    else startoffset = offset_str.Int();
                }

                gr.start.contig = range.Before(":").Int();
                gr.stop.contig  = gr.start.contig;

                if ( full_chr ) {
                  // we could set 'window' to something ridiculously large to enforce the assignments below to result
                  // from the same two lines of code in the 'else' clause, but this way it is cleaner:
                  gr.start.offset = 0;
                  gr.stop.offset = NumBasesInContig( gr.start.contig );
                } else {
                  gr.start.offset = (startoffset >= window) ? startoffset - window : 0;
                  gr.stop.offset  = (startoffset + window <= NumBasesInContig(gr.start.contig)) ? startoffset + window : NumBasesInContig(gr.start.contig);
                }

                LinearOffset(gr.start.contig, gr.start.offset, gr.start.linear_offset);
                LinearOffset(gr.stop.contig, gr.stop.offset, gr.stop.linear_offset);
            } else if ((range.Contains("+") && range.Before("+").IsInt()) || range.IsInt()) {
                unsigned int window = window_size;
                unsigned int startlinearoffset;
                if (range.Contains("+")) {
                    window = range.After("+").Int();
                    startlinearoffset = range.Before("+").Int();
                } else {
                    startlinearoffset = range.Int();
                }

                gr.start.linear_offset = (startlinearoffset >= window) ? startlinearoffset - window : 0;
                gr.stop.linear_offset  = (startlinearoffset + window <= NumBases()) ? startlinearoffset + window : NumBases();

                ContigAndOffset(gr.start.linear_offset, gr.start.contig, gr.start.offset);
                ContigAndOffset(gr.stop.linear_offset, gr.stop.contig, gr.stop.offset);
            } else {
                gr.start.contig = 0;
                gr.start.offset = 0;
                gr.start.linear_offset = 0;

                gr.stop.contig = contigoffsets.size() - 1;
                gr.stop.offset = NumBasesInContig(gr.stop.contig);
                gr.stop.linear_offset = NumBases();
            }

            gr.start.rcm = this;
            gr.stop.rcm = this;

            // Barf on a bogus range.
            Assert(gr.start < gr.stop);

            return gr;
        }

        /// Parses a vector of ranges.  Ranges may be specified on the command-line or in a file.
        vec<ReferenceRange> ParseRanges(const vec<String> & sranges, unsigned int window_size = 50) const {
            vec<ReferenceRange> vranges;

            if (sranges.size() > 0) {
                for (unsigned int rangeindex = 0; rangeindex < sranges.size(); rangeindex++) {
                    if (IsRegularFile(sranges[rangeindex])) {
                        ifstream rangestream(sranges[rangeindex].c_str());
                        String rangestring;

                        while (getline(rangestream, rangestring)) {
                            if (!rangestring.Contains("#")) {
                                vranges.push_back(ParseRangeString(rangestring, window_size));
                            }
                        }

                        rangestream.close();
                    } else {
                        vranges.push_back(ParseRangeString(sranges[rangeindex], window_size));
                    }
                }

                // There's a lot of code that's much more efficient if it doesn't need to backtrack to process a range.  Here, we sort the ranges and make explicit that condition.  This may not be a Good Thing(tm) overall, but I'm goin' with it for now.
                sort(vranges.begin(), vranges.end(), EarlierRange);
            } else {
                vranges.push_back(ParseRangeString("all"));
            }

            return vranges;
        }

        /// Parses a vector of ranges.
        vec<ReferenceRange> ParseRanges(const String & file_name) const
        {
            vec<String> ranges;
            Ifstream(in, file_name);
            while (! in.eof())
            {
                string s;
                getline(in, s);
                if (s == "") { continue; }
                ranges.push_back(String(s));
            }
            return this->ParseRanges(ranges);
        }


        /// Parses a single coordinate string.
        ReferenceCoordinate ParseCoordinateString(const String & scoord) const {
            if (scoord.Contains(":")) {
                int contig = scoord.Before(":").Int();
                unsigned int offset = scoord.After(":").Int();

                return (*this)(contig, offset);
            }

            unsigned int linear_offset = scoord.Int();
            return (*this)[linear_offset];
        }

        /// Parses a vector of coordinates.  Coordinates may be specifiedon the command-line or in a file.
        vec<ReferenceCoordinate> ParseCoordinates(const vec<String> & scoords) const {
            vec<ReferenceCoordinate> vcoords;

            for (unsigned int coordindex = 0; coordindex < scoords.size(); coordindex++) {
                if (IsRegularFile(scoords[coordindex])) {
                    ifstream coordstream(scoords[coordindex].c_str());
                    String coordstring;

                    while (getline(coordstream, coordstring)) {
                        if (!coordstring.Contains("#")) {
                            vcoords.push_back(ParseCoordinateString(coordstring));
                        }
                    }

                    coordstream.close();
                } else {
                    vcoords.push_back(ParseCoordinateString(scoords[coordindex]));
                }
            }

            sort(vcoords.begin(), vcoords.end(), EarlierCoordinate);

            return vcoords;
        }

        /// Number of bases in reference
        unsigned int NumBases(void) const { return linear_size; }

        /// Number of contigs in reference
        int NumContigs(void) const { return contigoffsets.size(); }

        /// Number of bases in contig
        unsigned int NumBasesInContig(int contig) const { return ((contig == contigoffsets.isize() - 1) ? linear_size : contigoffsets[contig+1]) - contigoffsets[contig]; }

    private:
        vec<unsigned int> contigoffsets;
        unsigned int linear_size;
};

// This is declared here rather than with the object defintion because we needed the definition of ReferenceCoordinateMapper before this method would work.
inline void ReferenceCoordinate::Shift(long shift) {
    linear_offset += shift;

    if (offset + shift < rcm->NumBasesInContig(contig)) {
        offset += shift;
    } else {
        rcm->ContigAndOffset(linear_offset, contig, offset);
    }
}

#endif
