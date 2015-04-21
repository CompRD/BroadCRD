///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Equiv.h"
#include "Intvector.h"
#include "PairsManager.h"
#include "PrintAlignment.h"
#include "Qualvector.h"
#include "ScoreAlignment.h"
#include "TokenizeString.h"
#include "btl/CBarcodes.h"
#include "btl/ClusterOffsets.h"
#include "btl/DedupWinner.h"
#include "pairwise_aligners/KmerAligner.h"
#include "pairwise_aligners/SmithWatBandedA.h"

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

/**
 * CBarcodes
 * Clear
 */
void
CBarcodes::Clear() {
    indexes_.clear();
    indexes_.resize(4, -1);
    pids_.clear();
}

/**
 * CBarcodes
 * SetOne
 */
void
CBarcodes::SetOne(const int ii, const int index) {
    indexes_[ii] = index;
}

/**
 * CBarcodes
 * SetIndexes
 */
void
CBarcodes::SetIndexes(const IntVec &indexes) {
    for (size_t ii = 0; ii < indexes.size(); ii++)
        indexes_[ii] = indexes[ii];
}

/**
 * CBarcodes
 * SetIndexes
 *
 * It requires as_string to be like this: aa.bb.cc.dd (where aa,... dd
 * are intergers, >= -1).
 */
void
CBarcodes::SetIndexes(const String &as_string) {
    vec<String> tokens;
    Tokenize(as_string, '.', tokens);
    ForceAssertEq( (int)tokens.size( ), 4);
    IntVec indexes(4, -1);
    for (int ii = 0; ii < 4; ii++)
        indexes[ii] = tokens[ii].Int();
    this->SetIndexes(indexes);
}

/**
 * CBarcodes
 * SetAll
 */
void
CBarcodes::SetAll(const int64_t pid, const IntVec &indexes) {
    pids_.clear();
    pids_.resize(1, pid);
    this->SetIndexes(indexes);
}

/**
 * CBarcodes
 * Reserve
 */
void
CBarcodes::Reserve(const size_t n_pids) {
    pids_.reserve(n_pids);
}

/**
 * CBarcodes
 * AddPairId
 */
void
CBarcodes::AddPairId(const int64_t pid) {
    pids_.push_back(pid);
}

/**
 * CBarcodes
 * Rank
 */
int
CBarcodes::Rank() const {
    return IntVecRank(indexes_);
}

/**
 * CBarcodes
 * CompactName
 */
String
CBarcodes::CompactName() const {
    String name = "";
    for (size_t ii = 0; ii < indexes_.size(); ii++)
        name += ToString(indexes_[ii]) + (ii < indexes_.size() - 1 ? "." : "");
    return name;
}

/**
 * CBarcodes
 * Write
 */
void
CBarcodes::Write(ostream &out) const {
    for (size_t ii = 0; ii < indexes_.size(); ii++)
        out << indexes_[ii] << (ii < indexes_.size() - 1 ? "\t" : "\n");
    out << pids_.size() << "\n";
    for (size_t ii = 0; ii < pids_.size(); ii++)
        out << pids_[ii] << "\n";
}

/**
 * CBarcodes
 * PrintOneLine
 */
void
CBarcodes::PrintOneLine(ostream &out) const {
    out << "r" << this->Rank() << "  [";
    for (size_t ii = 0; ii < indexes_.size(); ii++)
        out << indexes_[ii] << (ii < indexes_.size() - 1 ? "  " : "");
    out << "]  " << pids_.size() << " pair(s): ";
    if (pids_.size() <= 6)
        for (size_t ii = 0; ii < pids_.size(); ii++)
            out << pids_[ii] << (ii < pids_.size() - 1 ? ", " : "\n");
    else
        out << pids_[0] << ", " << pids_[1] << ", (" << pids_.size() - 4
                << " more pairs), " << pids_[pids_.size() - 2] << ", "
                << pids_[pids_.size() - 1] << "\n";
}

/**
 * CBarcodes
 * operator==
 */
bool
operator==(const CBarcodes &left, const CBarcodes &right) {
    for (size_t ii = 0; ii < left.indexes_.size(); ii++)
        if (left.indexes_[ii] != right.indexes_[ii])
            return false;
    return true;
}

/**
 * CBarcodes
 * operator!=
 */
bool
operator!=(const CBarcodes &left, const CBarcodes &right) {
    return (!(left == right));
}

/**
 * CBarcodes
 * operator<
 */
bool
operator<(const CBarcodes &left, const CBarcodes &right) {
    for (size_t ii = 0; ii < left.indexes_.size() - 1; ii++) {
        if (left.indexes_[ii] < right.indexes_[ii])
            return true;
        if (left.indexes_[ii] > right.indexes_[ii])
            return false;
    }
    return (left.indexes_.back() < right.indexes_.back());
}

/**
 * IntVecRank
 */
int
IntVecRank(const IntVec &indexes) {
    int indexes_size = indexes.size();
    for (int ii = 0; ii < indexes_size; ii++)
        if (indexes[ii] < 0)
            return ii;
    return indexes_size;
}

/**
 * IntVecsConsistent
 */
bool
IntVecsConsistent(const IntVec &left, const IntVec &right) {
    if (left.size() != right.size())
        return false;
    for (size_t ii = -0; ii < left.size(); ii++) {
        if (left[ii] < 0 || right[ii] < 0)
            continue;
        if (left[ii] != right[ii])
            return false;
    }
    return true;
}

/**
 * WriteVecBarcodes
 */
void
WriteVecBarcodes(const String &outfile, const vec<CBarcodes> &barcodes) {
    ofstream out(outfile.c_str());
    for (size_t ii = 0; ii < barcodes.size(); ii++)
        barcodes[ii].Write(out);
    out.close();
}

/**
 * LoadVecBarcodes
 */
void
LoadVecBarcodes(const String &infile, vec<CBarcodes> &barcodes) {
    barcodes.clear();

    ifstream in(infile.c_str());

    // Load one set at a time.
    while (in) {

        // Go to head of block (leave, if there are no blocks).
        IntVec indexes(4, -2);
        while (in) {
            String line;
            vec<String> tokens;
            getline(in, line);
            if (!in)
                break;
            Tokenize(line, tokens);
            if (tokens.size() != 4)
                continue;
            for (int ii = 0; ii < 4; ii++)
                indexes[ii] = tokens[ii].Int();
            break;
        }

        bool head_found = true;
        for (int ii = 0; ii < 4; ii++) {
            if (indexes[ii] == -2) {
                head_found = false;
                break;
            }
        }
        if (!head_found)
            break;

        // Number of pairs.
        size_t n_pairs;
        in >> n_pairs;

        // Pair ids.
        CBarcodes newbar;
        newbar.SetIndexes(indexes);
        newbar.Reserve(n_pairs);
        for (size_t ii = 0; ii < n_pairs; ii++) {
            size_t pair_id;
            in >> pair_id;
            newbar.AddPairId(pair_id);
        }

        // Add newbar to barcodes.
        barcodes.push_back(newbar);

    }

    in.close();
}

/**
 * BuildVecBarcodes
 */
void
BuildVecBarcodes(const VecIntVec &raw_barcodes, const IntVec &exp_ranks,
        vec<CBarcodes> &barcodes) {
    barcodes.clear();

    // Consts.
    const size_t n_pairs = raw_barcodes.size();

    // Build initial barcodes (one pid each), and sort.
    barcodes.resize(n_pairs);
    for (size_t ii = 0; ii < n_pairs; ii++)
        barcodes[ii].SetAll(ii, raw_barcodes[ii]);
    sort(barcodes.begin(), barcodes.end());

    // Find "firsts".
    vec<size_t> firsts;
    firsts.reserve(barcodes.size());

    firsts.push_back(0);
    for (size_t ii = 1; ii < n_pairs; ii++)
        if (barcodes[ii] != barcodes[ii - 1])
            firsts.push_back(ii);

    // Build final barcodes.
    vec<CBarcodes> result;
    result.reserve(firsts.size());

    for (size_t fid = 0; fid < firsts.size(); fid++) {
        size_t bid = firsts[fid];
        size_t n_bars = 1;

        // Initial loop (so we can reserve).
        for (size_t jj = bid + 1; jj < n_pairs; jj++) {
            if (barcodes[jj] != barcodes[bid])
                break;
            else
                n_bars++;
        }

        // New CBarcodes (reserve memory for all pair ids).
        CBarcodes newbar = barcodes[bid];
        newbar.Reserve(n_bars);

        // Add pair ids to new CBarcodes.
        for (size_t jj = bid + 1; jj < n_pairs; jj++) {
            if (barcodes[jj] != barcodes[bid])
                break;
            else
                newbar.AddPairId(barcodes[jj].PairId(0));
        }
        result.push_back(newbar);
    }

    // Swap and go.
    swap(barcodes, result);

}

/**
 * DiscardNonFullRank
 */
void
DiscardNonFullRank(const int FULL_RANK, const VecIntVec &raw_barcodes,
        const IntVec &exp_ranks, vec<int> &deleted) {
    const size_t n_pairs = raw_barcodes.size();

    for (size_t pid = 0; pid < n_pairs; pid++) {
        if (deleted[pid] != NOT_DELETED)
            continue;

        const int exprank = exp_ranks[pid];
        const int rank = IntVecRank(raw_barcodes[pid]);

        if (rank < 1) {
            deleted[pid] = (exprank < 1) ? IDX_A_CUT : IDX_A_MISS;
            continue;
        }
        if (FULL_RANK < 2)
            continue;

        if (rank < 2) {
            deleted[pid] = (exprank < 2) ? IDX_B_CUT : IDX_B_MISS;
            continue;
        }
        if (FULL_RANK < 3)
            continue;

        if (rank < 3) {
            deleted[pid] = (exprank < 3) ? IDX_C_CUT : IDX_C_MISS;
            continue;
        }
        if (FULL_RANK < 4)
            continue;

        if (rank < 4) {
            deleted[pid] = (exprank < 4) ? IDX_D_CUT : IDX_D_MISS;
            continue;
        }
    }

}

/**
 * SetCapEnds
 */
void
SetCapEnds(const vecbvec&reads
         ,const bvec& cap_bases
         ,const vec<minial>& caps
         ,vec<short>& CapEnds
         ){
    CapEnds.assign(reads.size(),0);
    const int cap_length=cap_bases.isize();
    ForceAssert(cap_length>0);
    #pragma omp parallel for
    for( size_t cc=0 ; cc<caps.size(); ++cc){
        const minial& loc_min = caps[cc];
        if(loc_min.Rc()) continue;
        const auto rid = loc_min.ReadId();
        const auto offset=loc_min.Offset();
        if(   offset >=0 && offset+cap_length <= reads[rid].isize()
           && std::equal(cap_bases.begin(),cap_bases.end()
                        ,reads[rid].begin()+offset)
          ){
            CapEnds[rid] = offset+cap_length;
        }
    }
}

/**
 * DiscardShortGenomic
 */
void
DiscardShortGenomic(const int MIN_GENLEN, const PairsManager &pairs,
        const vecbvec &reads, const vec<minial> &caps,
        const vec<CBarcodes> &barcodes, vec<short>&vCapEnds, vec<int> &klens, vec<int> &deleted) {
    const size_t n_bars = barcodes.size();
    const size_t n_reads = pairs.nReads();

    klens.clear();
    klens.resize(n_reads, 0);

    // Build map read_id to cap_align_id.
    vec<int64_t> caps_idx(n_reads, -1);
    for (int64_t ii = 0; ii < (int64_t) caps.size(); ii++) {
        int64_t rid = caps[ii].ReadId();
        if (caps_idx[rid] == -1)
            caps_idx[rid] = ii;
        else
            caps_idx[rid] = -2;
    }

    // Loop over all barcodes.
    for (size_t bar_id = 0; bar_id < n_bars; bar_id++) {
        const CBarcodes &bar = barcodes[bar_id];
        const size_t n_pairs = bar.NPairs();

        // Loop over all pairs in barcode.
        for (size_t ii = 0; ii < n_pairs; ii++) {
            const size_t pid = bar.PairId(ii);
            if (deleted[pid] != NOT_DELETED)
                continue;

            size_t id1 = pairs.ID1(pid);
            size_t id2 = pairs.ID2(pid);

            if( vCapEnds[id1]!=0 ){
                if (reads[id1].isize() >= vCapEnds[id1] + MIN_GENLEN){
                    vCapEnds[id2]=0;
                    continue;
                }
                else
                    vCapEnds[id1]=0;
            }
            if( vCapEnds[id2]!=0 ){
                if (reads[id2].isize() >= vCapEnds[id2] + MIN_GENLEN)
                    continue;
                else
                    vCapEnds[id2]=0;
            }

            const minial *al1 = (caps_idx[id1] > -1) ? &caps[caps_idx[id1]] : 0;
            const minial *al2 = (caps_idx[id2] > -1) ? &caps[caps_idx[id2]] : 0;
            const bool fwcap_1 = (al1 && !al1->Rc());
            const bool fwcap_2 = (al2 && !al2->Rc());

            // Check input.
            ForceAssert( al1 || al2);
            if (al1 && al2)
                ForceAssert( al1->Rc( ) != al2->Rc( ));

            // Cap (and indexes) found on both reads.
            if (al1 && al2) {
                const minial *sel_al = (fwcap_1) ? al2 : al1;
                const int sel_id = (fwcap_1) ? id2 : id1;
                const int genlen = sel_al->offset_;

                if (genlen < MIN_GENLEN)
                    deleted[pid] = TOO_SHORT;
                else
                    klens[sel_id] = genlen;

                continue;
            }

            // Cap (and indexes) found on one read only.
            if (al1) {
                ForceAssertGe( (int)reads[id2].size( ), MIN_GENLEN);
                klens[id2] = reads[id2].size();
            }
            else {
                ForceAssertGe( (int)reads[id1].size( ), MIN_GENLEN);
                klens[id1] = reads[id1].size();
            }

        } // loop over all pairs

    } // loop over all barcodes

}

/**
 * DiscardDuplicates
 */
void
DiscardDuplicates(const PairsManager &pairs, const vecbvec &reads,
        const vecqvec &quals, const vec<CBarcodes> &barcodes, const vec<short> &vCapEnds,
        const vec<int> &klens, vec<int> &deleted) {
    const size_t n_reads = pairs.nReads();
    const size_t n_bars = barcodes.size();

    // HEURISTICS.
    const int max_mis = 12;
    const int max_ind = 2;
    const double max_er = 0.45;
    const int band = max_ind + 2;
    const float max_al_score = 100.0;

    // SANTEMP
#ifdef SKIP
    CBarcodes keyb;
    keyb.SetOne( 0, 97 );
    keyb.SetOne( 1, 72 );
    keyb.SetOne( 2, 89 );
    keyb.SetOne( 3, 165 );
#endif

    // Loop over all barcodes.
    for (size_t bar_id = 0; bar_id < n_bars; bar_id++) {
        const CBarcodes &bar = barcodes[bar_id];
        const size_t n_pairs = bar.NPairs();

        // SANTEMP
#ifdef SKIP
        if ( bar == keyb ) {
            bar.PrintOneLine( cout );
            cout << "\n" << endl;
            for (int ii=0; ii<(int)bar.NPairs( ); ii++)
            cout << pairs.ID1( bar.PairId( ii ) ) << ",";
            cout << "\n" << endl;
        }
#endif

        // Build local chunks.
        vecbvec chunks;
        vecqvec qchunks;
        vec<size_t> to_pid;
        chunks.reserve(n_pairs);
        qchunks.reserve(n_pairs);
        to_pid.reserve(n_pairs);
        
        bool bHasLocalSequence=false;

        for (size_t ii = 0; ii < n_pairs; ii++) {
            const size_t pid = bar.PairId(ii);
            if (deleted[pid] != NOT_DELETED)
                continue;

            size_t id1 = pairs.ID1(pid);
            size_t id2 = pairs.ID2(pid);
            ForceAssert( vCapEnds[id1] != 0 || vCapEnds[id2] != 0 || klens[id1] > 0 || klens[id2] > 0 );

            qvec qual;
            if(vCapEnds[id1]!=0){
                const size_t klen = reads[id1].size()-vCapEnds[id1];
                qual.resize(klen);
                CopyQuals(quals[id1], 0, qual, 0, klen, true);
                chunks.push_back(bvec(reads[id1], vCapEnds[id1], klen));
                chunks.back().ReverseComplement();
                bHasLocalSequence=true;
            }
            else if(vCapEnds[id2]!=0){
                const size_t klen = reads[id2].size()-vCapEnds[id2];
                qual.resize(klen);
                CopyQuals(quals[id2], 0, qual, 0, klen, true);
                chunks.push_back(bvec(reads[id2], vCapEnds[id2], klen));
                chunks.back().ReverseComplement();
                bHasLocalSequence=true;
            }
            else{
                int id = klens[id1] > 0 ? id1 : id2;
                int klen = klens[id1] > 0 ? klens[id1] : klens[id2];
                qual.resize(klen);

                CopyQuals(quals[id], 0, qual, 0, klen);
                chunks.push_back(bvec(reads[id], 0, klen));
            }

            qchunks.push_back(qual);
            to_pid.push_back(pid);
        }

        // Nothing left for this barcode.
        const int nchunks = chunks.size();
        if (nchunks < 1)
            continue;

        // Length of shortest chunk (needed to decide which K to use).
        int shorter_clen = chunks[0].size();
        for (size_t ii = 1; ii < chunks.size(); ii++)
            shorter_clen = Min((int) chunks[ii].size(), shorter_clen);

        // Create a KmerAligner. HEURISTICS embedded here.
        const int K = (shorter_clen < 50 || bHasLocalSequence) ? 12 : 24;
        KmerAligner<12> aligner12;
        KmerAligner<24> aligner24;
        aligner12.SetKmerStep(2);
        aligner24.SetKmerStep(4);

        if (K == 12)
            aligner12.SetBases(chunks);
        else if (K == 24)
            aligner24.SetBases(chunks);
        else
            ForceAssert( 1 == 0);

        // Create an equiv rel object, by connecting duplicated reads.
        equiv_rel equiv(nchunks);

        // Loop over all chunks.
        for (int ii = 0; ii < nchunks; ii++) {
            vec<int> kapos;
            vec<int> kaids;
            if (K == 12)
                aligner12.FindPossibleAlignments(chunks[ii], kapos, kaids);
            else if (K == 24)
                aligner24.FindPossibleAlignments(chunks[ii], kapos, kaids);

            // Separate seeds (only test ii against jj, for jj>ii).
            vec<vec<int> > all_offsets(nchunks);
            for (size_t kid = 0; kid < kaids.size(); kid++) {
                int jj = kaids[kid];
                if (jj > ii)
                    all_offsets[jj].push_back(-kapos[kid]);
            }

            // Cluster offsets between ii and jj (jj>ii).
            vec<vec<int> > clusters(nchunks);
            for (int jj = ii + 1; jj < nchunks; jj++) {
                if (all_offsets[jj].size() > 0) {
                    sort(all_offsets[jj].begin(), all_offsets[jj].end());
                    ClusterOffsets(1, all_offsets[jj], clusters[jj]);
                }
            }

            // Decide (by alignment) if ii and jj are connected.
            for (int jj = ii + 1; jj < nchunks; jj++) {

                // Chunks ii and jj may already belong to the same orbit.
                if (equiv.Equiv(ii, jj))
                    continue;

                // There may be several clusters between ii and jj.
                bool ii_jj_align = false;
                for (int kk = 0; kk < clusters[jj].isize(); kk++) {
                    const int off = clusters[jj][kk];

                    // Align ii and jj, according to offset kk.
                    align al;
                    int err = 0;
                    SmithWatBandedA(chunks[ii], chunks[jj], off, band, al, err);

                    if (al.Pos1() - al.pos1() < K)
                        continue;

                    unsigned int iilen = chunks[ii].size();
                    unsigned int itlen = chunks[jj].size();
                    bool embedded1 = (al.pos1() == 0
                            && (uint) al.Pos1() == iilen);
                    bool embedded2 = (al.pos2() == 0
                            && (uint) al.Pos2() == itlen);
                    if (!(embedded1 || embedded2))
                        continue;

                    int n_err = al.Errors(chunks[ii], chunks[jj]);
                    look_align theLA(ii, jj, iilen, itlen, false, al, 0, n_err,
                            0);
                    if (!theLA.IsProper())
                        continue;

                    int n_muts = theLA.Mutations();
                    int n_ind = theLA.Indels();
                    float score = ScoreAlignment(al, chunks[ii], qchunks[ii],
                            chunks[jj], qchunks[jj]);
                    bool test1 = (n_muts <= max_mis && n_ind <= max_ind);
                    bool test2 = (score <= max_al_score);
                    bool test3 = (theLA.ErrorRate() <= max_er);

                    if (!(test1 || test2 || test3))
                        continue;

                    // Found valid align between ii and jj.
                    equiv.Join(ii, jj);
                    ii_jj_align = true;
                    break;
                }

                // One align found, skip eventual remaining (ii-jj) offsets.
                if (ii_jj_align)
                    break;
            }

        }

        // Select the best representative from each orbit, delete the rest.
        vec<int> reps;
        equiv.OrbitReps(reps);

        for (int ii = 0; ii < reps.isize(); ii++) {
            vec<int> orbit;
            equiv.Orbit(reps[ii], orbit);
            if (orbit.size() < 2)
                continue;

            int winner = DedupWinner(qchunks, orbit);
            for (int jj = 0; jj < orbit.isize(); jj++)
                if (orbit[jj] != winner)
                    deleted[to_pid[orbit[jj]]] = IS_DUPLICATE;
        }

    }  // loop over all barcodes

}

/**
 * SaveIndexedfasta
 */
void
SaveIndexedFasta(const PairsManager &pairs, const vecbvec &reads,
        const vecqvec &quals, const String &out_head,
        const vec<CBarcodes> &barcodes, const vec<short>& vCapEnds, const vec<int> &klens,
        const vec<int> &deleted) {
    // File names.
    String out_fasta_file = out_head + ".fasta";
    String out_qual_file = out_head + ".qual";
    String out_counts_file = out_head + ".barcounts";

    // Open output streams, and dump bases and quals.
    ofstream bout(out_fasta_file.c_str());
    ofstream qout(out_qual_file.c_str());
    ofstream countsout(out_counts_file.c_str());

    size_t current = 0;
    for (size_t barid = 0; barid < barcodes.size(); barid++) {
        const CBarcodes &bar = barcodes[barid];
        const int n_pairs = bar.NPairs();

        int count = 0;
        for (int ii = 0; ii < n_pairs; ii++) {
            const size_t pid = bar.PairId(ii);
            if (deleted[pid])
                continue;
            count++;

            const size_t id1 = pairs.ID1(pid);
            const size_t id2 = pairs.ID2(pid);

            ForceAssert( vCapEnds[id1] != 0 || vCapEnds[id2] != 0 || klens[id1] > 0 || klens[id2] > 0 );

            size_t id, id_bc;
            bvec chunk_bases;
            qvec chunk_quals;
            if(vCapEnds[id1]!=0){
                FatalErr("not fully implemented");
                id = id1;
                chunk_bases.assign(reads[id1].begin()+vCapEnds[id1], reads[id1].end());
                chunk_bases.ReverseComplement();
                const size_t klen = reads[id1].size()-vCapEnds[id1];
                chunk_quals.resize(klen);
                CopyQuals(quals[id1], 0, chunk_quals, 0, klen, true);
            }
            else if(vCapEnds[id2]!=0){
                FatalErr("not fully implemented");
                id = id2;
                chunk_bases.assign(reads[id2].begin()+vCapEnds[id2], reads[id2].end());
                chunk_bases.ReverseComplement();
                const size_t klen = reads[id2].size()-vCapEnds[id2];
                chunk_quals.resize(klen);
                CopyQuals(quals[id2], 0, chunk_quals, 0, klen, true);
            }
            else{
                id = klens[id1] > 0 ? id1 : id2;
                if( klens[id1] > 0){
                    id=id1;
                    id_bc=id2;
                }
                else{
                    id=id2;
                    id_bc=id1;

                }
                chunk_quals.resize(klens[id]);
                CopyQuals(quals[id], 0, chunk_quals, 0, klens[id]);
                chunk_bases.assign(reads[id].begin(), reads[id].begin()+klens[id]);
            }

//            const size_t id = klens[id1] > 0 ? id1 : id2;

            String name = ToString(current) + "_" + bar.CompactName() + "_"
                                                  + ToString(id_bc) + "_"
                                                  + ToString(id) + "_"
                                                  + ToString(0) + "_"
                                                  + ToString(klens[id]) + "_"
                                                  ;
            current++;

//            bvec chunk_bases(reads[id], 0, klens[id]);
//            qvec chunk_quals(klens[id], 0);
//            CopyQuals(quals[id], 0, chunk_quals, 0, klens[id]);

            chunk_bases.Print(bout, name);
            Print(qout, chunk_quals, name);
        }

        if (count > 0) {
            countsout << bar.CompactName() << "\t" << count << "\n";
            bout << "\n";
            qout << "\n";
        }
    }

    // Close streams.
    bout.close();
    qout.close();

}

/**
 * ReportDeleted
 */
void
ReportDeleted(const vec<int> &deleted, const PairsManager &pairs,
        ostream &out) {
    const int n_deletions = 16;
    const int n_pairs = pairs.nPairs();

    vec<size_t> totals(n_deletions, 0);
    for (size_t ii = 0; ii < deleted.size(); ii++)
        totals[deleted[ii]]++;
    for (int ii = 2; ii < n_deletions; ii++)
        totals[ii] += totals[ii - 1];

    vec<vec<String> > table;
    table.push_back(
            MkVec(String(""), ToStringAddCommas(n_pairs), String("100.0%"),
                    String("total pairs in input")));

    for (int ii = 1; ii < n_deletions; ii++) {
        size_t n_remaining = n_pairs - totals[ii];
        double ratio = SafeQuotient(n_remaining, n_pairs);

        String legend = "";
        if (ii == 1)
            legend = "low quality pairs";
        else if (ii == 2)
            legend = "unplaced cap";
        else if (ii == 3)
            legend = "multiple placement of cap on the same read";
        else if (ii == 4)
            legend = "cap placed illogically on both reads";
        else if (ii == 5)
            legend = "cap on both reads, but inconsistent indexes";
        else if (ii == 6)
            legend = "not enough space for A, before cap";
        else if (ii == 7)
            legend = "there is space for A, but A is not found";
        else if (ii == 8)
            legend = "not enough space for B, before A";
        else if (ii == 9)
            legend = "there is space for B, but B is not found";
        else if (ii == 10)
            legend = "not enough space for C, before B";
        else if (ii == 11)
            legend = "there is space for C, but C is not found";
        else if (ii == 12)
            legend = "not enough space for D, before C";
        else if (ii == 13)
            legend = "there is space for D, but D is not found";
        else if (ii == 14)
            legend = "genomic portion too short";
        else if (ii == 15)
            legend = "after removal of duplicated pairs";

        table.push_back(
                MkVec("[" + ToString(ii) + "]", ToStringAddCommas(n_remaining),
                        ToString(100. * ratio, 1) + "%", legend));
    }

    out << "\nCAP AND INDEXES PLACEMENT STATISTICS\n\n"
            << "Each line contains the number of pairs remaining after\n"
            << "a test has been performed (tests are run sequentially).\n"
            << "The last line, therefore, contains the number of usable pairs.\n\n";

    PrintTabular(out, table, 3, "rrrl");
    out << endl;

}

