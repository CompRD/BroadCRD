#include <string.h>
#include <cassert>

#include "Basevector.h"
#include "MainTools.h"
#include "CoreTools.h"
#include "math/Functions.h"
#include "util/BaitMap.h"
#include <fstream>

/// Tries to add the sequence interval [<start>, <end>] from contig <contig> to this map;
/// if an interval with exactly the same parameters is already stored, the method will
/// return and new redundant copy will not be added. If the interval is added,
/// its name will be auto-generated.
void BaitMap::add(int contig, longlong start, longlong end)
{
    char buf[4096];
    sprintf(buf, "bait-%lu", this->size());
    String name(buf);
    this->add(contig, start, end, name);
}

/// Tries to add the named sequence interval [<start>, <end>] from contig <contig> to this map;
/// if an interval with exactly the same parameters is already stored, the method will
/// return and new redundant copy will not be added. Otherwise the interval is added
/// and assigned the specified <name>
void BaitMap::add(int contig, longlong start, longlong end, const String & name)
{
    HoIntervalWithId r(start, end+1, contig);

#if 0
#if 0
    // This linear search could be improved, of course.
    for (int i = 0; i < regions.isize(); i++)
    {
        if (regions[i].id == r.id)
        {
            if (regions[i] == r)
            {
                // This is a duplicate.
                return;
            }
        }
    }
#else
    if (regions_unique.find(r) != regions_unique.end())
    {
        // this is a duplicate.
        return;
    }
#endif
#endif

    regions.push_back(r);
    regions_unique[r] = true;
    names.push_back(name);
}

/// Translates an interval on the original genome g_(contig, start, end) into an interval
/// on the "reference" represented by the collectoin of intervals stored in this map: r_contig
/// is the index of the stored interval and r_start, r_end are the positions within that interval
/// onto which g_star and g_end map. If positions g_start, g_end from contig g_contig
/// do not map inside the same stored interval ("reference contig"), this method returns
/// false and the content of output variables r_* is not chanhed.
bool BaitMap::lookup_from_genome(int  g_contig, longlong  g_start, longlong  g_end,
                                 int* r_contig, longlong* r_start, longlong* r_end)
{

    if (this->reverse_index.size() <= 0)
    {
        printf("BaitMap ERROR: trying to lookup from genome but reverse index hasn't been computed.\n"); fflush(stdout);
    }
    ForceAssert(this->reverse_index.size() > 0);

    if (this->size() == 0)
    {
        *r_contig = g_contig;
        *r_start  = g_start;
        *r_end    = g_end;
        return true;
    }

    pair<int, longlong> g_1(g_contig, g_start);
    pair<int, longlong> g_2(g_contig, g_end);

    if ((reverse_index.count(g_1) != 0) &&
        (reverse_index.count(g_2) != 0))
    {
        pair<int, longlong> r_1 = reverse_index[g_1];
        pair<int, longlong> r_2 = reverse_index[g_2];

        if (r_1.first != r_2.first)
        {
            printf("ERROR A: %d != %d\n", r_1.first, r_2.first);
            fflush(stdout);
            return false;
        }

        if((r_2.second - r_1.second) != (g_2.second - g_1.second))
        {
            printf("ERROR B: %ld != %ld\n", (r_2.second - r_1.second), (g_2.second - g_1.second));
            fflush(stdout);
            return false;
        }

        //printf("OK: %d\n", r_1.first); fflush(stdout);

        //assert(r_1.first == r_2.first);
        //assert((r_2.second - r_1.second) == (g_2.second - g_1.second));

        *r_contig = r_1.first;
        *r_start  = r_1.second;
        *r_end    = r_2.second;

        return true;
    }

    return false;
}

/// Translates a pair of positions r_(start, end) on a <r_contig>-th stored interval
/// into pair of offsets g_(start, end) on contig g_contig of the original genome, to
/// which the stored intervals refer
bool BaitMap::lookup_from_reference(int  r_contig, longlong  r_start, longlong  r_end,
                                    int* g_contig, longlong* g_start, longlong* g_end)
{
    if (this->size() == 0)
    {
        *g_contig = r_contig;
        *g_start  = r_start;
        *g_end    = r_end;
        return true;
    }

    *g_contig = regions[r_contig].id;
    *g_start  = r_start + regions[r_contig].Start();
    *g_end    = r_end   + regions[r_contig].Start();

    return true;
}

/// Translates position on the original genome (to which the stored
/// intervals refer) to the position on the stored intervals (i.e.
/// the index of the stored interval r_contig, and offset within that
/// interval, r)
bool BaitMap::lookup_from_genome(int  g_contig, longlong  g,
                                 int* r_contig, longlong* r)
#if 1
{
    if (this->size() == 0)
    {
        *r_contig = g_contig;
        *r  = g;
        return true;
    }

    pair<int, longlong> g_1(g_contig, g);

    if (reverse_index.count(g_1) != 0)
    {
        pair<int, longlong> r_1 = reverse_index[g_1];

        //printf("OK: %d\n", r_1.first); fflush(stdout);

        //assert(r_1.first == r_2.first);
        //assert((r_2.second - r_1.second) == (g_2.second - g_1.second));

        *r_contig = r_1.first;
        *r        = r_1.second;

        return true;
    }

    return false;
}
#else
{
    longlong g2 = g + 1;
    longlong r2;
    return lookup_from_genome(g_contig, g, g2, r_contig, r, &r2);
}
#endif


/// Translates position on the on the stored intervals (i.e.
/// offset r within r_contigh-th stored interval) into the
/// location on th eoriginal genom, to which the stored intervals refer
bool BaitMap::lookup_from_reference(int  r_contig, longlong  r,
                                    int* g_contig, longlong* g)
{
    longlong r2 = r + 1;
    longlong g2;
    return lookup_from_reference(r_contig, r, r2, g_contig, g, &g2);
}


/// Reads intervals into this map from a file; the format used
/// is text file with one "contig_id start end" line per interval.
/// This method \em adds intervals from the file to the currently
/// stored collection of intervals, i.e. the old content is preserved
bool BaitMap::read(String file_name, bool skip_reverse_index)
{
  std::ifstream is(file_name.c_str());

  String line;
  while ( getline(is,line) )
  {
      DeleteLeadingWhiteSpace(line);
      DeleteTrailingWhiteSpace(line);
      if (line.size() == 0) { /* printf("LINE: \"%s\"\n", line.c_str()); fflush(stdout); */ break; }
      vec<String> tokens;
      Tokenize(line, tokens);

      longlong contig = tokens[0].Int();
      longlong start  = tokens[1].Int();
      longlong end    = tokens[2].Int();

      if ( tokens.size() > 3 ) {
          this->add(contig,start,end,tokens[3]); // add interval with name
      } else {
          this->add(contig, start, end); // add interval and generate name automatically
      }
  }
  ForceAssert(is.eof());

  if (skip_reverse_index == false) { this->compute_reverse_index(); }

  //printf("BaitMap::read() : read (%u,%u) intervals from %s.\n", (unsigned int)(this->size()), (unsigned int)(this->regions.size()), file_name.c_str()); fflush(stdout);

  return true;
}

/// Writes intervals stored in this map into uncompressed text file
/// with specified name using the format "contig_id start end", one line
/// per interval.
bool BaitMap::write(String file_name)
{
    std::ofstream out(file_name.c_str());
    ForceAssert(out);

    for (unsigned int i = 0; i < this->size(); i++)
    {
      HoIntervalWithId region = (*this)[i];

      longlong contig = region.id;
      longlong start  = region.Start();
      longlong end    = (region.Stop() - 1);

      out << contig << ' ' << start << ' ' << end << '\n';
    }

    out.close();
    ForceAssert(out);

    return true;
}

/// writes names of the stored intervals into a file with
/// specified name in uncompressed text format.
bool BaitMap::write_names(String file_name)
{
    std::ofstream out(file_name.c_str());
    ForceAssert(out);

    for (unsigned int i = 0; i < this->size(); i++)
    {
      out << names[i] << '\n';
    }

    out.close();
    ForceAssert(out);

    return true;
}


/// writes names of the stored intervals into a file with
/// specified name in binary format (fastb.names).
bool BaitMap::write_names_bin(String bin_file_name)
{

    vecString vnames;
    vnames.resize(this->size());
    for (unsigned int i = 0; i < this->size(); i++)
    {
      vnames[i] = names[i];
    }
    vnames.WriteAll(bin_file_name);

    return true;
}

/// Merge each group of overlapping stored intervals into one interval,
/// and store the resulting set of disjoint intervals into <out>. Old content
/// of <out> will be lost, and the content of this map will stay unchanged.
void BaitMap::Merge(BaitMap& out, bool skip_reverse_index)
{
    vec<HoIntervalWithId> x = regions;
    Sort(x, LessById);
    ::Merge(x);
    out.regions = x;
    out.regions_unique.clear();
    out.names.resize(out.size());
    char buf[4096];
    for ( long unsigned int i = 0 ; i < out.size() ; i++ ) {
      out.regions_unique[x[i]]=true;
      sprintf(buf, "Merged-contig-%lu", i);
      out.names[i] = buf;
    }
    if (!skip_reverse_index) { out.compute_reverse_index(); }
}

/// Fills out the passed map <out> with intervals stored in this map
/// widened by <size> bases on both sides (if an interval can not be
/// widened to the left by <size> bases because contig start is reached,
/// then it will be widened till contig start). If size==0, this map will
/// be copied to the target <out> map, but trimming on the left side will still
/// occur for the intervals that extend beyond the contig start.
/// The names of new, widened  intervals will be set to "<old_name>_widened_<size>"
/// if size!=0 or old names will be used if size==0.
/// An old content of <out> is destroyed and this map is left unchanged.
void BaitMap::Widen(int size, BaitMap& out, bool skip_reverse_index)
{

  out.regions.resize(0);

  for (unsigned int i = 0; i < this->size(); i++)
  {
      HoIntervalWithId region = (*this)[i];

      longlong contig = region.id;
      longlong start  = region.Start() - size;
      longlong end    = (region.Stop() - 1) + size;

      if (start < 0) { start = 0; }

      if ( size != 0 ) {
          out.add(contig, start, end, names[i]+"_widened_"+ToString(size));
      } else {
          out.add(contig, start, end, names[i]);
      }
  }
  if (!skip_reverse_index) { out.compute_reverse_index(); }
}

/// Fills <contig_length> vector with lengths of each intrerval
/// stored in this map. <contig_length> will be resized to match
/// the number of intervals stored in this map, and the order of
/// the returned length vaules matches the order of intervals.
void BaitMap::Dimensions(vec<longlong>& contig_lengths)
{
    contig_lengths.resize(this->size());

    for (unsigned int i = 0; i < this->size(); i++)
    {
        HoIntervalWithId region = (*this)[i];
        contig_lengths[i] = region.Length();
    }
}

/// Removes from this map all the intervals that overlap by <overlap>
/// or more bases with any of the intervals stored in <other_map>.
/// Returns the indexes of the removed baits in the original collection
vec<size_t> BaitMap::remove_intersection(const BaitMap& other_map, int overlap)
{
	vec<size_t> removed;
    vec<HoIntervalWithId> new_regions;
    vec<String> new_names;
    new_regions.reserve(regions.size());
    new_names.reserve(regions.size());

    unsigned int my_size = this->size();
    unsigned int other_size = other_map.size();

    for (unsigned int i = 0; i < my_size; i++)
    {
        unsigned int j;
        for (j = 0; j < other_size; j++)
        {
            HoIntervalWithId& A = (*this)[i];
            const HoIntervalWithId& B = other_map[j];

            if ((A.id == B.id) && (Overlap(A, B) >= overlap))
            {
                break;
            }
        }

        if (j == other_map.size()) {
            new_regions.push_back((*this)[i]);
            new_names.push_back(this->name(i));
        }
        else {
        	removed.push_back(i);
        }
    }

    regions = new_regions;
    names = new_names;
    compute_reverse_index();

    return removed;
}

/// Extracts actual nucleotide sequences as specified by stored intervals form genome
/// reference <reference> and stores them into collection of sequences <output>.
/// Sequences in <output> are syncronized with intervals stored by this class (i.e. the order
/// is preserved). If an interval extends beyond the reference contig it is supposed to come from,
/// and <trim> is true, then the interval will be trimmed to the contig boundaries on one or both
/// sides as needed
///
/// @param[in] reference genome sequence on which the stored intervals are defined
/// @param[out] output subsequences from the reference genome defined by the stored intervals
/// @param trim if true, trim the intervals that run beyond the boundaries of their source contigs
/// to those boundaries; if trim is false (default) the runaway intervals currently result in
/// undefined behavior and possibly segmentation fault
void BaitMap::extract_sequences(const vecbasevector& reference, vecbasevector& output, bool trim)
{
    output.resize(0);
    output.reserve(this->size());

    for (unsigned int i = 0; i < this->size(); i++)
    {
        HoIntervalWithId& interval = (*this)[i];
        int start = interval.Start();
        int stop = interval.Stop();

        if (stop > reference[interval.id].isize())
        {
            printf("ERROR: interval extends to the right beyond reference.\n");
            printf("BAIT %d named %s contig %d from %d to %d\n",
                        i, this->names[i].c_str(), interval.id, interval.Start(), interval.Stop());
            if ( trim ) stop = reference[interval.id].isize();
        }
        if (start < 0 )
        {
            printf("ERROR: interval extends to the left beyond reference.\n");
            printf("BAIT %d named %s contig %d from %d to %d\n",
                        i, this->names[i].c_str(), interval.id, interval.Start(), interval.Stop());
            if ( trim ) start = 0;
        }

        if (start < stop && stop <= reference[interval.id].isize()) {
            basevector v;
            v.SetToSubOf(reference[interval.id], start, stop-start);
            output.push_back(v);
        }
    }
}

/// Fills the map (genome location) -> (sequence interval location), where
/// genome location is given as (contig, offset on the contig) in terms of the
/// original genome sequence the intervals were built from, and sequence interval location
/// is (interval id, position within the interval) with interval id being the index
/// into the <regions> vector. Each genome location can map onto multiple intervals
/// if the intervals overlap; only one mapping will be stored by this method [POSSIBLE BUG?]
void BaitMap::compute_reverse_index()
{
    this->reverse_index.clear();

    for (int i = 0; i < regions.isize(); i++)
    {
        for (longlong j = 0; j < regions[i].Length(); j++)
        {
            pair<int, longlong> genomic_location(regions[i].id, regions[i].Start() + j);
            pair<int, longlong> map_location(i, j);

            this->reverse_index[genomic_location] = map_location;

            pair<int, longlong> test_key(regions[i].id, regions[i].Start() + j);
            pair<int, longlong> test_value = this->reverse_index[test_key];
            assert(test_key == genomic_location);
            assert(test_value == map_location);
        }
    }
}

void BaitMap::MakeMask(vec< vec<char> >& mask)
{
    for (size_t i = 0; i < this->size(); i++)
    {
        longlong contig = (*this)[i].id;
        longlong start  = (*this)[i].Start();
        longlong stop   = (*this)[i].Stop();

        if (mask.isize() <= contig) { mask.resize(contig+1); }
        if (mask[contig].isize() <= stop) { mask[contig].resize(stop+1); }

        for (longlong j = Min((longlong)0,start); j <= stop; j++)
        {
            mask[contig][j] = 1;
        }
    }
}

void BaitMap::MakeBoolVectorMask(vec<boolvector> *mask, int padding)
{
    mask->resize(0);

    for (size_t i = 0; i < this->size(); i++)
    {
        longlong contig = (*this)[i].id;
        longlong start  = (*this)[i].Start() - padding;
        longlong stop   = (*this)[i].Stop() + padding;

        if (mask->isize() <= contig) { mask->resize(contig+1); }
        if ((*mask)[contig].isize() <= stop) { (*mask)[contig].resize(stop, false); }

        for (longlong j = Min((longlong)0,start); j < stop; j++)
        {
            (*mask)[contig][j] = 1;
        }
    }
}

