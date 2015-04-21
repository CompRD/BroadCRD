///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
#ifndef READ_ERROR_ANALYZER_H
#define READ_ERROR_ANALYZER_H

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "ReadError.h"
#include "RefLocus.h"
#include "ParallelVecUtilities.h"

// Find the canonical form of read errors (edits) that will generate the
// observed read from the original read that contains errors. The canonical
// form contains minimal number of edits(ins,sub,del) for  generating the observed
// read. In case there are several equally good minima, the canonical form is
// chosen by the following rules:
// 1. As much head and tail bases matched as possible
// 2. Preference match > sub > del > ins, counting from the back of the sequence.
//
// Example: 
//
// If your original and observed reads are:
// a0 = BaseVec("AACCCGTG");
// b0 = BaseVec("ATCCGAG");
// The returned vector shall contains edits 
//    ReadError(1 DEL)
//    ReadError(2 SUB T)
//    ReadError(6 SUB A)
// which correspond to the alignment
//   AACCCGTG
//   A-TCCGAG 
// instead of the alternative alignments that has the same edit distance:
//   AACCCGTG
//   AT-CCGAG  (edits: 2 SUB T, 3 DEL, 6 SUB A, not chosen because positon 3 prefers matching according to rule 2)
// or 
//   AACCCGTG
//   -ATCCGAG  (edits: 0 DEL, 2 SUB T, 6 SUB A, not chosen because rule 1 )

void FindCanonicalReadError( const basevector& original, const basevector& observed, 
        vec<ReadError> *p_edits,
        bool MaxHeadMatching = true,
        int D_Score = 1 , int I_Score = 1, int S_Score = 1
        );


// =========================================================================
// A friend teller can tell or predict if two reads are friend and 
// return their relative positions.
// =========================================================================

class FriendTeller {
public:
    FriendTeller( const RefLocusVec& ref_loc_vec, 
            const vecbasevector& genome,
            int MinOverlap=100 ) 
        :  ref_loc_vec_( ref_loc_vec ) , genome_(genome)
    {   this->InitFriends(MinOverlap); }

    //// interface function needed to be provided by implementation
    size_t Gid (size_t rid) const { return ref_loc_vec_[rid].mRefID; }
    size_t Start (size_t rid) const { return ref_loc_vec_[rid].mOffset; }
    size_t Stop (size_t rid) const { return ref_loc_vec_[rid].mOffset + ref_loc_vec_[rid].mLength; }
    Bool   Rc (size_t rid) const { return ref_loc_vec_[rid].mRC; }
    size_t NReads() const { return ref_loc_vec_.size(); }
    const vecbvec Genome() const { return genome_; }

    int Overlap ( size_t i, size_t j ){
        if ( Gid( i ) != Gid( j ) ) return 0;
        return IntervalOverlap( (signed)Start( i ), (signed)Stop( i ), 
                (signed)Start( j ), (signed)Stop( j ) );
    }

    // relative shift of the read i aligned to read j ( read j fw )
    int Offset( int i, int j ) const { 
        ForceAssertEq( Gid(i), Gid(j) );
        if ( ! Rc( j ) )
            return Start( i ) - Start( j );
        else
            return Stop( j ) - Stop( i );
    }

    void InitFriends (int MinOverlap) {
        vec<RefLocus> locs( ref_loc_vec_.begin(), ref_loc_vec_.end() );
        vec<size_t> rids( NReads(), vec<size_t>::IDENTITY );
        ParallelSortSync( locs, rids );
        flist_.resize( NReads() ) ;
        for ( size_t i = 0; i < locs.size(); ++i ) {
            size_t gid = locs[i].mRefID;
            size_t start = locs[i].mOffset;
            size_t stop = locs[i].mOffset + locs[i].mLength;
            for ( size_t j = i+1; j < locs.size(); j++ ) {
                if ( locs[j].mOffset < stop - MinOverlap && locs[i].mRefID == gid ) {
                    flist_[ rids[i] ].push_back( rids[j] );
                    flist_[ rids[j] ].push_back( rids[i] );
                } else {
                    break;
                }
            }
        }
        //for ( size_t i = 0; i < NReads(); ++i ) 
        //    for ( size_t j = i+1; j < NReads(); ++j ) {
        //        if ( Overlap( i, j ) >= MinOverlap ) {
        //            flist_[i].push_back( j );
        //            flist_[j].push_back( i );
        //        }
        //    }
        for ( size_t i = 0; i < flist_.size(); ++i ) 
            Sort( flist_[i] );
    }

    const vec<int>& FriendList( int i ) const { return flist_[i]; }

    bool IsFriend(size_t i, size_t j) const { return BinMember( flist_[i], j ); }
    Bool IsRc(size_t i, size_t j) const { return Rc(i) ^ Rc( j ); }
    Bool IsRc(size_t i) const { return Rc( i ); }

    void GetExtendedRead( size_t read_id, int left_ext, int right_ext, 
            basevector* p_read ) const {
        bool rc = Rc( read_id );
        int start = Start( read_id ) - ( rc ? right_ext : left_ext ); 
        int stop = Stop( read_id )   + ( rc ? left_ext  : right_ext );
        p_read->SetToSubOf( Genome()[Gid( read_id )], start, stop - start );
        if ( rc ) p_read->ReverseComplement();
    }

    void GetActualRead( size_t read_id, basevector* p_read ) const 
    {   GetExtendedRead (read_id, 0, 0, p_read);   }

    void GetActualReads( vecbasevector* p_reads ) const {
        p_reads->clear();
        p_reads->reserve( NReads() );
        for ( size_t read_id = 0; read_id < NReads(); ++read_id ) {
            bvec read;
            GetActualRead( read_id, &read );
            p_reads->push_back( read );
        }
    }

private:
    vec< vec<int> > flist_;
    const RefLocusVec& ref_loc_vec_;
    const vecbasevector& genome_;
    // disable copy and assign constructors
    FriendTeller( const FriendTeller& );
    FriendTeller& operator= ( const FriendTeller& );
};

// Friend teller that use RefLocusVec generated from ReadSimulatorSimpleGenerator
class FriendTellerFromRefLocusVec : public FriendTeller {
public:

    // interface function needed to be provided by implementation

private:
};



// ====================================
// An analyzer for ReadErrorCorrectors
// ====================================

class ReadErrorAnalyzer {
public:
    struct ResultT {
        size_t read_id;
        int n_errors;
        int n_missed;
        int n_added;
        String log;
        ResultT( ) : read_id(-1), n_errors(-1), n_missed(-1), n_added(-1), log() {};
        ResultT( int read_id1, int n_errors1, int n_missed1, int n_added1, const String& log1 )
            : read_id(read_id1), n_errors(n_errors1), n_missed(n_missed1),
            n_added(n_added1), log(log1) {};
    };

    ReadErrorAnalyzer( const vecbasevector& truth, const vecbasevector& reads, 
                            const vecbasevector& corrected ) 
        : truth_(truth), reads_(reads), corrected_(corrected), 
          p_friend_teller_(NULL), verbosity_(0) {};
    ~ReadErrorAnalyzer( ) {}

    // Analyze all reads and print summary
    void Analyze();
    // Analyze only one read
    ResultT AnalyzeRead( size_t rid ) const;
    // Analyze all and print errors
    void AnalyzeAll( ) const;

    // == Visualize the read errors by pileup of the read with its friends ==

    // Need a friend teller in order to perform read pileup
    void AddFriendTeller ( FriendTeller* p_teller) {   p_friend_teller_ = p_teller;   }
    // Create a pile of the friend reads over the founder read
    vec<String> Pileup( int read_id ) const;
    vec<String> PileupPerfectReads( int read_id ) const;
    vec<String> PileupActualReads( int read_id ) const;

    // == Logging options
    void SetVerbosity( int v ) { verbosity_ = v; };

    // == Static utility functions ==

    // Makeup the read erros in the original sequence, for building pile of the reads.
    static vec<String> PileupReadErrors( const basevector& read, const vec<vec<ReadError> >& edits_vec,
                                  const vec<int>& starts, const vec<int>& stops );
    // Return the sequence string where the read errors are marked up.
    static vec<String> MarkupReadErrors( const basevector& read, const vec<ReadError>& edits );

private:
    size_t NReads() const { return reads_.size(); }
    // disable copy and assign constructors
    ReadErrorAnalyzer(const ReadErrorAnalyzer&);
    ReadErrorAnalyzer& operator= (const ReadErrorAnalyzer&);

    const vecbasevector& truth_;
    const vecbasevector& reads_;
    const vecbasevector& corrected_;
    const FriendTeller* p_friend_teller_;
    int verbosity_;
};

#endif
