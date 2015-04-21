///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "PrintAlignment.h"
#include "VecUtilities.h"
#include "lookup/LookAlign.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "pairwise_aligners/LongReadsAligner.h"
#include "reporting/SnapToGrid.h"
#include "util/CSmallKmers.h"

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

/**
 * SanitizeAlign
 */
void SanitizeAlign(const bvec &target, const bvec &query, align &al)
{
  // HEURISTICS.  
  const int WIN_SIZE = 64;      // sliding window size
  const double MAX_ER = 0.45;   // max error rate

  // Align length.
  int align_len = al.Pos1() - al.pos1();
  for (int ii=0; ii<al.Nblocks(); ii++)
    align_len += Max(0, al.Gaps(ii));

  // Align shorter than window size.
  if (align_len < WIN_SIZE) {
    align empty_align;
    al = empty_align;
    return;
  }

  // Errors counter ([0]=undef, [-1]=match, [1]=ins, [2]=del, [3]=mis).
  vec<int> errors(align_len, 0);
  vec<bool> is_gap(align_len, false);
  {
    int pos = 0;
    int pos1 = al.pos1();
    int pos2 = al.pos2();
    for (int ii=0; ii<al.Nblocks(); ii++) {
      int gap = al.Gaps(ii);
      int len = al.Lengths(ii);

      // Insertion.
      if (gap < 0) {
	for (int jj=0; jj<-gap; jj++) {
	  is_gap[pos] = true;
	  errors[pos] = 1;
	  pos++;
	  pos1++;
	}
      }
      
      // Deletion.
      if (gap > 0) {
	for (int jj=0; jj<gap; jj++) {
	  is_gap[pos] = true;
	  errors[pos] = 2;
	  pos++;
	  pos2++;
	}
      }

      // Match or mismatch.
      for (int jj=0; jj<len; jj++) {
	errors[pos] = (query[pos1] == target[pos2]) ? -1 : 3;
	pos++;
	pos1++;
	pos2++;
      }
    }
  }

  // Error counters for sliding windows.
  vec<int> win_errors(align_len, 0);
  win_errors[0] = errors[0] > 0 ? 1 : 0;
  for (int pos=1; pos<WIN_SIZE; pos++)
    win_errors[pos] = win_errors[pos-1] + (errors[pos] > 0 ? 1 : 0);

  for (int pos=WIN_SIZE; pos<align_len; pos++) {
    win_errors[pos]
      = win_errors[pos-1]
      + (errors[pos] > 0 ? 1 : 0)
      - (errors[pos-WIN_SIZE] > 0 ? 1 : 0);
    if (win_errors[pos] < 0) win_errors[pos] = 0;
  }
  
  // Convert each sliding window error rate to a bool flag.
  vec<bool> win_good(align_len, false);
  for (int pos=WIN_SIZE-1; pos<align_len; pos++)
    win_good[pos] = SafeQuotient(win_errors[pos], WIN_SIZE) <= MAX_ER;

  // Find good windows.
  vec< pair<int,int> > good_segments;
  int pos = WIN_SIZE - 1;
  while (pos < align_len) {
    while (pos < align_len && ! win_good[pos]) pos++;
    int end = pos + 1;
    while (end < align_len && win_good[end]) end++;
    good_segments.push_back(make_pair(pos, end));
    pos = end;
  }

  // Find longest segment.
  pair<int,int> winner;
  {
    vec< pair<int,int> > len2sid;
    for (int ii=0; ii<good_segments.isize(); ii++) {
      int len = good_segments[ii].second - good_segments[ii].first;
      len2sid.push_back(make_pair(len, ii));
    }
    sort(len2sid.rbegin(), len2sid.rend());
    winner = good_segments[ len2sid[0].second ];
  }

  // Build new align (slide off gaps first).
  int alp1 = winner.first - WIN_SIZE + 1;
  int alp2 = winner.second;
  while (alp1 < align_len && is_gap[alp1]) alp1++;
  while (alp2 > 0 && is_gap[alp2-1]) alp2--;
  
  int final_pos1 = -1;
  int final_pos2 = -1;
  vec<bool> used_block(al.Nblocks(), false);
  {
    int pos = 0;
    int pos1 = al.pos1();
    int pos2 = al.pos2();
    for (int ii=0; ii<al.Nblocks(); ii++) {
      int gap = al.Gaps(ii);
      int len = al.Lengths(ii);
      if (gap < 0) {
	pos += -gap;
	pos1 += -gap;
      }
      if (gap > 0) {
	pos += gap;
	pos2 += gap;
      }
      if (alp1 <= pos && pos < alp2) {
	used_block[ii] = true;
	if (final_pos1 < 0) {
	  final_pos1 = pos1;
	  final_pos2 = pos2;
	}
      }
      pos += len;
      pos1 += len;
      pos2 += len;
    }
  }

  avector<int> final_gaps;
  avector<int> final_lens;

  int n_final = 0;
  for (int ii=0; ii<al.Nblocks(); ii++) {
    if (used_block[ii]) {
      n_final++;
      final_gaps.Append(al.Gaps(ii));
      final_lens.Append(al.Lengths(ii));
    }
  }
  if (n_final > 0) final_gaps(0) = 0;
  
  // Return (default empty align).
  align final;
  if (n_final > 0)
    final.Set(final_pos1, final_pos2, final_gaps, final_lens);
  
  swap(final, al);
  
}

/**
 * ClusterPlacements
 */
void ClusterPlacements(const vec< pair<int,int> > &placements,
                       vec<SRange> &ranges)
{
  ForceAssert(is_sorted(placements.begin(), placements.end()));
  
  // HEURISTICS.
  const int SEED_RADIUS = 24;
  const int MIN_REQUIRED = 18;

  // Clean up and sort seeds. Early exit of there are no placements.
  ranges.clear();
  if (placements.size() < 1) return;

  // Starts of clusters.
  int n_max_clusters = 0;
  vec<bool> starter(placements.size(), false);
  starter[0] = true;
  for (int ii=1; ii<placements.isize(); ii++) {
    int this_tid = placements[ii].first;
    int prev_tid = placements[ii-1].first;
    int this_start = placements[ii].second;
    int prev_start = placements[ii-1].second;
    if (this_tid != prev_tid || this_start - prev_start > SEED_RADIUS) {
      starter[ii] = true;
      n_max_clusters++;
    }
  }

  // Create ranges.
  ranges.reserve(n_max_clusters);
  int p1 = 0;
  while (p1 < starter.isize()) {
    while (p1 < starter.isize() && ! starter[p1]) p1++;
    int p2 = p1 + 1;
    while (p2 < starter.isize() && ! starter[p2]) p2++;
    if (p2 - p1 >= MIN_REQUIRED) {
      int tid = placements[p1].first;
      int start = placements[p1].second;
      int stop = placements[p2-1].second;
      int weight = p2 - p1;
      ForceAssertEq(tid, placements[p2-1].first);
      ForceAssertLe(start, stop);
      ranges.push_back(SRange(tid, start, stop, weight));      
    }
    p1 = p2;
  }

}

/**
 * FindRanges
 */
void FindRanges(const bvec &query,
                const CSmallKmers &finder,
                vec<SRange> &ranges)
{
  ranges.clear();
  
  // Collect seeds.
  int K_SM = finder.K();
  size_t n_seeds = 0;
  for (int pos=0; pos<(int)query.size()-K_SM; pos++)
    n_seeds += finder.NPlacements(finder.ToInt(query, pos));
  
  // Implied placement of query on target.
  vec< pair<int,int> > placements;
  placements.reserve(n_seeds);
  for (int pos=0; pos<(int)query.size()-K_SM; pos++) {
    int id = finder.ToInt(query, pos);
    for (int ii=0; ii<(int)finder.NPlacements(id); ii++) {    
      int tid = finder.TargetId(id, ii);
      int match = finder.TargetPos(id, ii);
      placements.push_back(make_pair(tid, match - pos));
    }
  }
  sort(placements.begin(), placements.end());
  
  // Cluster implied placements.
  ClusterPlacements(placements, ranges);
  
}




/**
 * LongReadsCoverages
 */
void LongReadsCoverages(const vecbvec &targets,
                        const vecbvec &queries,
                        const vec<look_align> &aligns,
                        const int sample_size,
                        ostream &out,
                        const bool csv,
                        const String nickname)
{

  // Extrapolate coverage based on sample_size.
  ForceAssertGt(sample_size, 0);
  float scale = float(queries.size()) / float(sample_size);


  // Bins.
  const vec<double> bins = MkVec(0., 1., 2., 3., 4., 5.);
  const int nbins = bins.size();

  // The windows (aligned portions).
  vec<int> al_lens(queries.size(), 0);
  for (int ii = 0; ii < aligns.isize(); ii++) {
    const look_align &al = aligns[ii];
    int qid = al.query_id;
    int win = al.a.Pos2() - al.a.pos2();
    if (win > al_lens[qid])
      al_lens[qid] = win;
  }
  
  // Coverages on each bin (not cumulative).
  vec<int>      n_reads(nbins, 0);
  vec<longlong> n_bases(nbins, 0);
  for (int qid = 0; qid < al_lens.isize(); qid++) {
    if (al_lens[qid] > 0) {
      int ibin = PlaceInBin(double(al_lens[qid]) / 1000.0, bins);
      if (ibin >= 0 && ibin < nbins) {
        ForceAssertGe(al_lens[qid], 1000*bins[ibin]);
        n_bases[ibin] += longlong(al_lens[qid]);
        n_reads[ibin] += 1;
      }
    }
  }
  

  // Coverages on each bin (cumulative).
  for (int ibin = nbins - 2; ibin >= 0; ibin--) {
    n_bases[ibin] += n_bases[ibin + 1];
    n_reads[ibin] += n_reads[ibin + 1];
  }


  // Total target length.
  longlong tot_tlen = 0;
  for (size_t ii = 0; ii < targets.size(); ii++)
    tot_tlen += targets[ii].size();

  cout << "tot_tlen= " << tot_tlen << endl;
  cout << "scale= " << scale << endl;
  
  // output results: 
  
 
  cout << "-----------------------------------------------------" << endl;



  cout << setw(14) << " bins";
  for (int i = 0; i < nbins; i++)
    cout << " " << setw(10) << bins[i];
  cout << "  " << nickname << endl;

  cout << setw(14) << " scale n_reads";
  for (int i = 0; i < nbins; i++)
    cout << " " << setw(10) << longlong(scale * n_reads[i]);
  cout << "  " << nickname << endl;
  
  cout << setw(14) << " L0 x n_reads";
  for (int i = 0; i < nbins; i++) {
    const int len_bin = 1000 * bins[i];
    cout << " " << setw(10) << longlong(scale * (len_bin * n_reads[i]));
  }
  cout << "  " << nickname << endl;

  cout << setw(14) << " n_bases";
  for (int i = 0; i < nbins; i++)
    cout << " " << setw(10) << longlong(scale * n_bases[i]);
  cout << "  " << nickname << endl;

  cout << setw(14) << " coverage";
  for (int i = 0; i < nbins; i++)
    cout << " " << setw(10) << longlong(scale * n_bases[i] / float(tot_tlen));
  cout << "  " << nickname << endl;



  cout << setw(14) << " n links";
  for (int i = 0; i < nbins; i++) {
    const int len_bin = 1000 * bins[i];
    cout << " " << setw(10) << longlong(scale * (n_bases[i] - len_bin * n_reads[i]) / float(tot_tlen));
  }
  cout << "  " << nickname << endl;

  cout << "-----------------------------------------------------" << endl;


  
  // Generate table.
  vec< vec<String> > table;
  vec<String> line;
  line.push_back(">= of (Kb)");
  for (int ii=0; ii<bins.isize(); ii++) 
    line.push_back(ToString(bins[ii]));
  table.push_back(line);
  line.clear();
  line.push_back("by count");
  for (int ii=0; ii<bins.isize(); ii++)
    line.push_back(ToString(n_reads[ii]));
  table.push_back(line);
  line.clear();
  line.push_back("cov (x)");
  for (int ii=0; ii<bins.isize(); ii++) {
    double cov = scale * SafeQuotient(n_bases[ii], tot_tlen);
    line.push_back(ToString(cov, 1));
  }
  table.push_back(line);

  // Print table.
  out << "LONG READ COVERAGE AS FUNCTION OF READ LENGTH (ALIGNED PORTIONS ONLY)"
      << "\n\n";
  if (scale != 1.0)
    out << "NB: coverage is extrapolated from a set of "
	<< sample_size << " reads, out of a\n"
	<< "    set of " << queries.size() << " reads (corresponding to "
	<< ToString(100.0 / scale, 1) << "% of the total)\n"
	<< "\n";
  if (csv)
    for (int ii=0; ii<table.isize(); ii++) {
      for (int jj=0; jj<table[ii].isize(); jj++) {
	out << table[ii][jj] << ",";
	if (jj == table[ii].isize() - 1)
	  out << (ii == table.isize() - 1 ? nickname + "\n" : "\n");
      }
    }
  else
    PrintTabular(out, table, 3, "lrrrrrrrr");
  out << endl;
  
}







/**
 * LongReadsAligner
 */
void LongReadsAligner(const int K_SM,
                      const vecbvec &targets,
                      const vecbvec &queries,
                      vec<look_align> * look_aligns_p,
                      const vec<int> *select,
                      const int VERBOSE,
                      ostream *log)
{
  // HEURISTICS.
  const double TO_NEXT_RATIO = 0.75;  // filter out reads with too many
  const int TO_ALIGN_MAX = 4;         //   locations with similar seed count

  const double MAX_BAND_R = .2;  // max band for SW: query length * MAX_BAND_R
  const int MIN_BAND = 12;       // min band for SW (overrides MAX_BAND_R)
  const int MIS = 3;             // cost for (affine) SW: mimsatch,
  const int GAP_OPEN = 3;        //   gap opening,
  const int GAP_EXT = 1;         //   and gap extension

  const int RADIUS = 36;         // used to decide if two aligns overlap
  

  // Create maps of K_SM-mers for contigs.
  CSmallKmers finder(K_SM, &targets);

  // Select ids to align.
  const size_t n_selected = (select != 0) ? select->size() : queries.size();
  vec<bool> is_selected(queries.size(), (select != 0) ? false : true);
  if (select != 0)
    for (size_t ii = 0; ii < select->size(); ii++)
      is_selected[(*select)[ii]] = true;
  
  // Log start date.
  if (log)  *log << Date() << ": aligning " << n_selected << " reads." << endl;
  

  // Loop over all queries.
  #pragma omp parallel for
  for (int query_id = 0; query_id < is_selected.isize(); query_id++) {
    if (is_selected[query_id]) {
      
      const bvec &query_fw = queries[query_id];
      
      if ((int)query_fw.size() >= K_SM) {

        bvec query_rc = query_fw;
        query_rc.ReverseComplement();

        // Find ranges for both fw and rc copies, and sort them.
        vec<SRange> ranges;
        vec<SRange> ranges_rc;
        FindRanges(query_fw, finder, ranges);
        FindRanges(query_rc, finder, ranges_rc);
        for (int ii = 0; ii < ranges_rc.isize(); ii++)
          ranges_rc[ii].weight_ = - ranges_rc[ii].weight_;
        copy(ranges_rc.begin(), ranges_rc.end(), back_inserter(ranges));
        sort(ranges.begin(), ranges.end());
    
        int break_point = -1;
        if (ranges.isize() == 1)
          break_point = 1;
        else {
          for (int ii = 1; ii < ranges.isize(); ii++) {
            if (ii > TO_ALIGN_MAX) break;
            double ratio = TO_NEXT_RATIO * double(Abs(ranges[ii-1].weight_));
            if (double(Abs(ranges[ii].weight_)) < ratio) {
              break_point = ii;
              break;
            }
          }
        }
        if (ranges.size() > 0)
          if (break_point < 0 && ranges.isize() <= TO_ALIGN_MAX)
            break_point = ranges.isize();
    
        // Extra log on ranges.
        if (VERBOSE > 0) {
#pragma omp critical
          for (int ii=0; ii<ranges.isize(); ii++)
            ranges[ii].PrintInfo(query_id, cout);
        }

        // No valid locations for this read (or its rc).
        if (break_point < 0) continue;

        // Align read.

        for (int loc_id = 0; loc_id < break_point; loc_id++) {
          const SRange &range = ranges[loc_id];
          const bvec &target = targets[range.target_];
          const bvec &query = range.weight_ < 0 ? query_rc : query_fw;
      
          int offset = (range.start_ + range.stop_) / 2;
          int band = range.stop_ - range.start_;   // generous band
      
          // Refuse to align if band is too big.
          int max_band = MAX_BAND_R * double(query.size());
          if (band > max_band) continue;

          // Ovveride band, if too small.
          band = Max(MIN_BAND, band);

          // Chop target.
          int start = Max(0, offset - band);
          int len = Min((int)query.size() + 2 * band,
                         (int)target.size() - start);
          bvec chop(target, start, len);

          // Do align.
          alignment ali;
          SmithWatAffine(query, chop, ali, true, true, MIS, GAP_OPEN, GAP_EXT);
      
          // Convert alignment.
          int alip1, alip2, errors;
          avector<int> aligaps, alilenghts;
          ali.Unpack(alip1, alip2, errors, aligaps, alilenghts);
          align al(alip1, start + alip2, aligaps, alilenghts);
      
          // Sanitize alignment (check align is not empty).
          SanitizeAlign(target, query, al);

          if (al.Nblocks() >= 1) {
        
            // Turn into look_align.
            int tid = range.target_;
            int qlen = query.size();
            int tlen = target.size();      
            bool rc1 = (range.weight_ < 0);
            look_align theAlign(query_id, tid, qlen, tlen, rc1, al, 0, 0, 0);
            theAlign.ResetFromAlign(al, query, target);

            #pragma omp critical
            look_aligns_p->push_back(theAlign);
        
            // Log event.
            if (VERBOSE >= 2) {
              #pragma omp critical
              {
                cout << "Align of r" << query_id
                     << (range.weight_ < 0 ? "[-]" : "[+]")
                     << " (" << qlen << " bases)  vs  "
                     << "c" << range.target_ << "\n"
                     << " alignment:   "
                     << "[" << theAlign.pos1()
                     << ", " << theAlign.Pos1()
                     << ")_" << qlen
                     << "   to   "
                     << "[" << theAlign.pos2()
                     << ", " << theAlign.Pos2()
                     << ")_" << tlen
                     << "\n"
                     << " align_len=" << theAlign.Pos1() - theAlign.pos1()
                     << " (" << theAlign.Errors() << " errors,"
                     << " err_rate=" << ToString(theAlign.ErrorRate(), 2) << ")\n"
                     << " based on " << Abs(range.weight_) << " seeds, with "
                     << "offset=" << offset << " and band=" << band << "\n";
	
                PrintVisualAlignment(true, cout, query, target, al);
      
              }
            }
          }
        }
      }
    }
  }
  if (log) *log << endl;

  // Sort aligns.
  if (log) *log << Date() << ": reloading and sorting aligns" << endl;
  sort(look_aligns_p->begin(), look_aligns_p->end());

  // Remove duplicates.
  if (log) *log << Date() << ": removing duplicates" << endl;
  {
    size_t na = look_aligns_p->size();
    vec<look_align> clean_aligns;
    clean_aligns.reserve(na);
    for (size_t ii = 1; ii < na; ii++) {
      const look_align & prev = (*look_aligns_p)[ii-1];
      const look_align & curr = (*look_aligns_p)[ii];
      
      if (prev.target_id != curr.target_id ||
          Abs(prev.a.pos1() - curr.a.pos1()) >= RADIUS ||
          Abs(prev.a.Pos1() - curr.a.Pos1()) >= RADIUS ||
          Abs(prev.a.pos2() - curr.a.pos2()) >= RADIUS ||
          Abs(prev.a.Pos2() - curr.a.Pos2()) >= RADIUS)
        
        clean_aligns.push_back(curr);
    }
    *look_aligns_p = clean_aligns;
  }

}
