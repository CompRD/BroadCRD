/////////////////////////////////////////////////////////////////////////////
////                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
////       This software and its documentation are copyright (2007) by the   //
////   Broad Institute/Massachusetts Institute of Technology.  All rights    //
////   are reserved.  This software is supplied without any warranty or      //
////   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
////   can be responsible for its use, misuse, or functionality.             //
///////////////////////////////////////////////////////////////////////////////

#ifndef BAITMAP_H
#define BAITMAP_H

#include <string.h>

#include "Basevector.h"
#include "Boolvector.h"
#include "CoreTools.h"
#include "math/Functions.h"
#include "math/HoInterval.h"

#include <map>

class BaitMap
{
public:

    BaitMap() { regions.resize(0); regions.reserve(100000); }
    BaitMap(const String& file_name) { regions.resize(0); regions.reserve(100000); this->read(file_name); }
    BaitMap(vec<HoIntervalWithId>& regions) { this->regions = regions; }

    void init_from_reference(vecbasevector& ref) // make a baitmap that just covers the reference
    {
        regions.resize(0);
        regions.reserve(ref.size());
        
        for (size_t i = 0; i < ref.size(); i++)
        {
            this->add(i, 0, ref[i].size()-1);
        }
    }

    /// Tries to add the sequence interval [<start>, <end>] from contig <contig> to this map;
    /// if an interval with exactly the same parameters is already stored, the method will
    /// return and new redundant copy will not be added. If the interval is added,
    /// its name will be auto-generated.
    void add(int contig, longlong start, longlong end);

    /// Tries to add the named sequence interval [<start>, <end>] from contig <contig> to this map;
    /// if an interval with exactly the same parameters is already stored, the method will
    /// return and new redundant copy will not be added. Otherwise the interval is added
    /// and assigned the specified <name>
    void add(int contig, longlong start, longlong end, const String & name);

    /// clear this map
    void clear() { regions.resize(0); regions.reserve(100000); }

    /// reserve space for <n> intervals
    void reserve(int n) { regions.reserve(n); }

    /// Translates an interval on the original genome g_(contig, start, end) into an interval
    /// on the "reference" represented by the collectoin of intervals stored in this map: r_contig
    /// is the index of the stored interval and r_start, r_end are the positions within that interval
    /// onto which g_star and g_end map. If positions g_start, g_end from contig g_contig
    /// do not map inside the same stored interval ("reference contig"), this method returns
    /// false and the content of output variables r_* is not chanhed.
    bool lookup_from_genome(int  g_contig, longlong  g_start, longlong  g_end,
                            int* r_contig, longlong* r_start, longlong* r_end);

    /// Translates a pair of positions r_(start, end) on a <r_contig>-th stored interval
    /// into pair of offsets g_(start, end) on contig g_contig of the original genome, to
    /// which the stored intervals refer
    bool lookup_from_reference(int  r_contig, longlong  r_start, longlong  r_end,
                               int* g_contig, longlong* g_start, longlong* g_end);

    /// Translates position on the original genome (to which the stored
    /// intervals refer) to the position on the stored intervals (i.e.
    /// the index of the stored interval r_contig, and offset within that
    /// interval, r)
    bool lookup_from_genome(int  g_contig, longlong  g,
                            int* r_contig, longlong* r);


    /// Translates position on the on the stored intervals (i.e.
    /// offset r within r_contigh-th stored interval) into the
    /// location on th eoriginal genom, to which the stored intervals refer
    bool lookup_from_reference(int  r_contig, longlong  r,
                               int* g_contig, longlong* g);

    bool lookup_from_genome(int  g_contig, ho_interval  g,
                            int* r_contig, ho_interval& r);
    bool lookup_from_reference(int  r_contig, ho_interval  r,
                               int* g_contig, ho_interval& g);

    /// Writes intervals stored in this map into uncompressed text file
    /// with specified name using the format "contig_id start end", one line
    /// per interval.
    bool write(String file_name);

    /// Reads intervals into this map from a file; the format used
    /// is text file with one "contig_id start end" line per interval.
    /// This method \em adds intervals from the file to the currently
    /// stored collection of intervals, i.e. the old content is preserved
    bool read(String file_name, bool skip_reverse_index=false);

    /// writes names of the stored intervals into a file with
    /// specified name in uncompressed text format.
    bool write_names(String file_name);

    /// writes names of the stored intervals into a file with
    /// specified name in binary format (fastb.names).
    bool write_names_bin(String bin_file_name);

    /// Returns i-th interval stored in this map
    HoIntervalWithId& operator[] (unsigned int i) { return regions[i]; }

    /// Returns i-th interval stored in this map
    const HoIntervalWithId& operator[] (unsigned int i) const { return regions[i]; }

    /// Returns the number of intervals stored in this map
    size_t size() const { return regions.size(); }

    /// Merge each group of overlapping stored intervals into one interval,
    /// and store the resulting set of disjoint intervals into <out>. Old content
    /// of <out> will be lost, and the content of this map will stay unchanged.
    void Merge(BaitMap& out, bool skip_reverse_index=false);

    /// Fills out the passed map <out> with intervals stored in this map
    /// widened by <size> bases on both sides (if an interval can not be
    /// widened to the left by <size> bases because contig start is reached,
    /// then it will be widened till contig start). An old content of <out>
    /// is destroyed. This map is left unchanged.
    void Widen(int size, BaitMap& out, bool skip_reverse_index=false);

    /// Fills <contig_length> vector with lengths of each intrerval
    /// stored in this map. <contig_length> will be resized to match
    /// the number of intervals stored in this map, and the order of
    /// the returned length vaules matches the order of intervals.
    void Dimensions(vec<longlong>& contig_lengths);

    /// Removes from this map all the intervals that overlap by <overlap>
    /// or more bases with any of the intervals stored in <other_map>.
    vec<size_t> remove_intersection(const BaitMap& other_map, int overlap);

    /// Makes a vec< vec<char> > mask indicating which positions are under the mask.
    void MakeMask(vec< vec<char> >& mask);

    /// Makes a vec<boolvector> mask indicating which positions are under the mask.
    /// The ends of each contig are extended by padding.
    void MakeBoolVectorMask(vec<boolvector> *mask, int padding=0);

    /// Extracts actual nucleotide sequences as specified by stored intervals form genome
    /// reference <reference> and stores them into collection of sequences <output>.
    /// Sequences in <output> are syncronized with intervals stored by this class (i.e. the order
    /// is preserved)
    /// @param[in] reference genome sequence on which the stored intervals are defined
    /// @param[out] output subsequences from the reference genome defined by the stored intervals
      void extract_sequences(const vecbasevector& reference, vecbasevector& output, bool trim = false);

    /// Returns name previously assigned to the i-th interval
    String name(int i) { return this->names[i]; }

    /// Returns <true> if the sequence interval <interval> is fully contained
    /// within one of the intervals already stored in this map
    bool Contains(HoIntervalWithId& interval)
    {
        for (int i = 0; i < regions.isize(); i++)
        {
            if (regions[i].Contains(interval))
            {
                return true;
            }
        }
        return false;
    }

    template<class IDX>
    void permute(const vec<IDX> & permutation) {
    	PermuteVec(this->regions, permutation);
    	PermuteVec(this->names, permutation);
    }

private:
    BaitMap(BaitMap const& other); // unimplemented -- no copying
    BaitMap& operator=(BaitMap const& other); // unimplemented -- no copying

    vec<HoIntervalWithId> regions;  ///< keeps sequence intervals stored in this map as (start, stop+1, contig)
    map<HoIntervalWithId, bool> regions_unique;    ///< every interval is also saved here (for quicker search)
    vec<String> names;    ///< names of each sequence interval; synchronized with <regions>
    map< pair<int, longlong>, pair<int, longlong> > reverse_index;
              ///< map (original genome location) -> (location on the set of stored sequence intervals)

   /// Fills the map (genome location) -> (sequence interval location), where
   /// genome location is given as (contig, offset on the contig) in terms of the
   /// original genome sequence the intervals were built from, and sequence interval location
   /// is (interval id, position within the interval) with interval id being the index
   /// into the <regions> vector. Each genome location can map onto multiple intervals
   /// if the intervals overlap; only one mapping will be stored by this method
    void compute_reverse_index();
};

#endif
