///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
#include "simulation/ReadErrorAnalyzer.h"
#include "VecUtilities.h"
#include "math/Array.h"
#include "util/TextTable.h"
#include <iterator>

void FindCanonicalReadError( const basevector& original, const basevector& observed, 
                              vec<ReadError> *p_edits,
                              bool MaxHeadMatching,
                              int D_Score , int I_Score, int S_Score )
{
    const int U_Score = 0; 
    //ForceAssertGt( original.size(), 0u );
    //ForceAssertGt( observed.size(), 0u );
    int M0 = original.size(), N0 = observed.size();
    // trim off same head and tail only  if matching score is zereo
    int len_head = 0;
    int len_tail = 0;
    if ( MaxHeadMatching ) 
        while ( len_head < M0 && len_head < N0 
                && original[ len_head ] == observed[ len_head ] ) 
            len_head++;
    while ( M0 - len_tail > len_head && N0 - len_tail > len_head
            && original[ M0- len_tail -1 ] == observed[ N0 - len_tail  -1 ] ) 
        len_tail++; 
    //LOG1 << "Trimming " << len_head << " and " << len_tail << " from both end. " << endl;
    // match the middle parts of the two sequences
    BaseVec a( original, len_head, M0 - len_head - len_tail );
    BaseVec b( observed, len_head, N0 - len_head - len_tail );
    int M = a.size(), N = b.size();
    RecArray<int> s( M+1, N+1 );                  // scoring matrix, uninitialized
    for ( int i = 0; i <= M; ++i ) s[i][0] = D_Score * i;
    for ( int j = 0; j <= N; ++j ) s[0][j] = I_Score * j;
    for ( int i = 1; i <= M; i++ )
        for ( int j = 1; j <= N; j++ ) {
            int sub_score = s[i-1][j-1];
            if ( a[i-1] != b[j-1] ) sub_score += S_Score;
            else sub_score += U_Score;
            int del_score = s[i-1][j] + D_Score;
            int ins_score = s[i][j-1] + I_Score;
            s[i][j] = std::min( {sub_score, del_score, ins_score} );
        }

    // Now find the edits that generate seq b from seq a
    // Remember to add back the trimming from front
    vec<ReadError> edits;
    int i = M , j = N ;
    //For debug purpose
    //String aa = a.ToString();
    //String bb = b.ToString();
    //string align_a, align_b;
    while ( i > 0 && j > 0 ){
        int score = s[i][j];
        // substitution
        if ( score == s[i-1][j-1] + ( a[i-1] == b[j-1] ? U_Score : S_Score ) ) {
            //align_a = aa[i-1] + align_a;
            //align_b = bb[j-1] + align_b;
            if ( a[i-1] != b[j-1] )
                edits.push_back( ReadError( len_head + i-1, ReadError::SUBSTITUTION, b[j-1], a[i-1] ) );
            i--, j--;
        }
        else if ( score == s[i-1][j] + D_Score ) {
            //align_a = aa[i-1] + align_a;
            //align_b = '-'  + align_b;
            edits.push_back( ReadError( len_head + i-1, ReadError::DELETION, ReadError::GAP_CODE, a[i-1] ) );
            i--;
        }
        else if ( score == s[i][j-1] + I_Score ) {
            //align_a = '-'  + align_a;
            //align_b = bb[j-1] + align_b;
            edits.push_back( ReadError( len_head + i, ReadError::INSERTION, b[j-1], ReadError::GAP_CODE ) );
            j--;
        }
        else { cout << "Error at position " << i << " " << j << endl; CRD::exit(1); }
    }

    while ( i > 0 ) {
        //align_a = aa[i-1] + align_a;
        //align_b = '-'  + align_b;
        edits.push_back( ReadError( len_head + i-1, ReadError::DELETION, ReadError::GAP_CODE, a[i-1] ) );
        i--;
    }

    while ( j > 0 ) {
        //align_a = '-'  + align_a;
        //align_b = bb[j-1] + align_b;
        edits.push_back( ReadError( len_head + i, ReadError::INSERTION, b[j-1], ReadError::GAP_CODE ) );
        j--;
    }
    p_edits->assign( edits.rbegin(), edits.rend() );
}

namespace {

String MarkupReadErrorsWithMapping( const basevector& read, const vec<ReadError>& edits, 
        const vec<int> pos_mapping, unsigned start= 0, unsigned stop=-1, bool ShowUnchanged = true ) {
    // initialize empty string spanning the whole region
    String str( pos_mapping.back() + 1, ' ' );
    start = min( start, read.size() );
    stop = min( stop, read.size() );
    if ( stop <= start ) return str;
    // from start to stop
    for( size_t pos = pos_mapping[start]; pos < (size_t)pos_mapping[stop-1]+1; ++pos ) 
        str[pos] = '|';
    for( size_t pos = start; pos < stop; ++pos ) 
        str[ pos_mapping[pos] ] = ( ShowUnchanged ? as_base( read[pos] ) : '.' ) ;
    for ( size_t i = 0; i < edits.size(); ++i ) {
        const ReadError& ed = edits[i];
        int pos = ed.getLocation();
        switch ( ed.getType() ) {
            case ReadError::SUBSTITUTION :
                str[ pos_mapping[pos] ] = tolower( as_base( ed.getReadBase() ) );
                break;
            case ReadError::DELETION :
                str[ pos_mapping[pos] ] = '-';
                break;
            case ReadError::INSERTION :
                int ins_pos = pos_mapping[pos] - 1;
                while( ins_pos >= 0 && str[ins_pos] == '|' ) --ins_pos;
                ++ins_pos;
                ForceAssertEq( str[ins_pos], '|' );
                str[ ins_pos ] = tolower( as_base( ed.getReadBase() ) );
                break;
        }
    }
    return str;
}

}

vec<String> ReadErrorAnalyzer::MarkupReadErrors( const basevector& read, const vec<ReadError>& edits ) {
    vec< vec<ReadError> > edit_vec( 1, edits );
    vec<int> starts(1, 0), stops(1, read.size() );
    return PileupReadErrors( read, edit_vec, starts, stops );
}

vec<String> ReadErrorAnalyzer::PileupReadErrors( const basevector& read, const vec<vec<ReadError> >& edits_vec,
                              const vec<int>& starts, const vec<int>& stops ) {
    vec<int> pos_mapping( read.size(), vec<int>::IDENTITY );
    for ( size_t i = 0; i < edits_vec.size(); ++i ) {
        for ( size_t j = 0; j < edits_vec[i].size(); ++j ) {
            if ( edits_vec[i][j].getType() == ReadError::INSERTION ) {
                for( size_t x = edits_vec[i][j].getLocation(); x < read.size(); x++ )
                    pos_mapping[x]++;
            }
        }
    }
    vec<String> pileup;
    vec<ReadError> empty_edits;
    pileup.push_back( MarkupReadErrorsWithMapping( read, empty_edits, pos_mapping ) );
    for ( size_t i = 0; i < edits_vec.size(); ++i ) {
        String str = MarkupReadErrorsWithMapping( read,  edits_vec[i], pos_mapping, starts[i], stops[i], false );
        pileup.push_back( str );
    }
    return pileup;
}

void ReadErrorAnalyzer::Analyze () {
    const int VERBOSITY = verbosity_;
    int n_perfect_reads = 0;
    int n_total_errors = 0;
    int n_missed_errors = 0;
    int n_added_errors = 0;
    #pragma omp parallel for
    for ( size_t i = 0; i < corrected_.size(); ++i ) {
        // No error in the read
        vec<ReadError> edit1, edit2;
        if ( reads_[i] != truth_[i] )  
            FindCanonicalReadError( truth_[i], reads_[i], &edit1 );
        if ( corrected_[i] != truth_[i] )  
            FindCanonicalReadError( truth_[i], corrected_[i], &edit2 );
        vec<ReadError> fixed, added;
        // errors that are fixed during error correction
        set_difference( edit1.begin(), edit1.end(), edit2.begin(), edit2.end(), 
                        back_inserter( fixed ) );
        // errors that are added during error correction
        set_difference( edit2.begin(), edit2.end(), edit1.begin(), edit1.end(), 
                        back_inserter( added ) );
        #pragma omp critical 
        {
            n_total_errors += edit1.size();
            n_missed_errors += edit1.size() - fixed.size();
            n_added_errors += added.size();
            if ( edit1.empty() ) n_perfect_reads++;
        }
        if ( VERBOSITY <= 0 )
            continue;
        if ( VERBOSITY < 3 && !added.size() && fixed.size() == edit1.size() )
            continue;
        cout << "========== read " << i << " =============\n";
        if ( VERBOSITY >= 2 ) {
            // Print detailed information of all the
            cout << "actual seq: " << truth_[i].ToString() << endl;
            cout << "observed  : " << reads_[i].ToString() << endl;
            cout << "corrected : " << corrected_[i].ToString() << endl;
            cout << "The actual read errors are          : " ;
            for ( size_t i = 0; i < edit1.size(); ++i ) 
                cout << edit1[i] << " ";
            cout << endl;
            cout << "The read errors after correction are: " ;
            for ( size_t i = 0; i < edit2.size(); ++i ) 
                cout << edit2[i] << " ";
            cout << endl;
        }
        cout << "    Found " << edit1.size() << " errors. " << endl;
        if ( added.size() > 0 )
            cout << "    Additional " << added.size() << " error introduced at read " << i << endl;
        if ( fixed.size() < edit1.size() )
            cout << "    Missed " << edit1.size() - fixed.size() << " error at read " << i << endl;
    }
    // Generating Report
    TextTable table;
    table << DoubleLine;
    table << "Total number of reads" << Tab << reads_.size() << EndRow; 
    table << "Number of perfect reads" << Tab << n_perfect_reads << EndRow;
    table << "Total number of errors" << Tab << n_total_errors << EndRow;
    table << "Number of missed errors" << Tab << n_missed_errors << EndRow;
    table << "Number of added errors" << Tab << n_added_errors << EndRow;
    table << DoubleLine;
    table.Print(cout, 5, "lr" );
}

ReadErrorAnalyzer::ResultT ReadErrorAnalyzer::AnalyzeRead (size_t rid ) const {
    ostringstream out;
    out << "========== read " << rid << " =============\n";
    vec<ReadError> edit1, edit2;
    if ( reads_[rid] != truth_[rid] )  
        FindCanonicalReadError( truth_[rid], reads_[rid], &edit1 );
    if ( corrected_[rid] != truth_[rid] )  
        FindCanonicalReadError( truth_[rid], corrected_[rid], &edit2 );

    vec<ReadError> fixed, added;
    // errors that are fixed during error correction
    set_difference( edit1.begin(), edit1.end(), edit2.begin(), edit2.end(), 
                    back_inserter( fixed ) );
    // errors that are added during error correction
    set_difference( edit2.begin(), edit2.end(), edit1.begin(), edit1.end(), 
                    back_inserter( added ) );

    // Print detailed information of all the
    // out << "actual seq: " << truth_[rid].ToString() << endl;
    // out << "observed  : " << reads_[rid].ToString() << endl;
    // out << "corrected : " << corrected_[rid].ToString() << endl;
    // out << "The actual read errors are          : " ;
    // for ( size_t i = 0; i < edit1.size(); ++i ) 
    //     out << edit1[i] << " ";
    // out << endl;
    // out << "The read errors after correction are: " ;
    // for ( size_t i = 0; i < edit2.size(); ++i ) 
    //     out << edit2[i] << " ";
    // out << endl;

    vec< vec<ReadError> > edits;
    edits.push_back( edit1 );
    edits.push_back( edit2 );
    basevector read_truth = truth_[rid];
    vec<String> lines = PileupReadErrors( read_truth, edits, vec<int>(2,0), vec<int>(2,-1) );
    for ( size_t i = 0; i < lines.size(); ++i ) 
        out << lines[i] << endl;

    out << "    Found " << edit1.size() << " errors. " << endl;
    if ( added.size() > 0 )
        out << "    Additional " << added.size() << " error introduced at read " << rid << endl;
    if ( fixed.size() < edit1.size() )
        out << "    Missed " << edit1.size() - fixed.size() << " error at read " << rid << endl;

    return ResultT( rid, edit1.size(), edit1.size() - fixed.size(), added.size(), out.str() );
}

void ReadErrorAnalyzer::AnalyzeAll () const {
    vec< ResultT > results( NReads() );
    #pragma omp parallel for
    for ( size_t i = 0; i < results.size(); ++i ) 
        results[i] = AnalyzeRead( i );
    for ( size_t i = 0; i < results.size(); ++i ) {
        const ResultT& res = results[i];
        if ( res.n_missed != 0 || res.n_added != 0 )
            cout << res.log;
    }
}

vec<String> ReadErrorAnalyzer::PileupActualReads( int read_id ) const {
    ForceAssert( p_friend_teller_ != NULL );
    // All the friend reads ( including self)
    vec<int> friend_ids = p_friend_teller_->FriendList( read_id );
    friend_ids.push_back( read_id );
    bool rc1 =  p_friend_teller_->IsRc(read_id);       // is the founder rc aligned?

    vec<int> r_starts; // offset of the friend read start in the founder forward direction
    for ( size_t i = 0; i < friend_ids.size(); ++i ) 
        r_starts.push_back( p_friend_teller_->Offset( friend_ids[i], read_id ) );
    vec<int> r_stops; // offset of the friend read end in the founder forward direction
    for ( size_t i = 0; i < friend_ids.size(); ++i ) 
        r_stops.push_back( r_starts[i] + truth_[ friend_ids[i] ].size() );
    SortSync( r_starts, r_stops,friend_ids );

    int start_ext = -r_starts.front();  // left_ext >= 0
    int stop_ext = Max( r_stops ) - truth_[read_id].size();
    basevector extend_read; 
    p_friend_teller_->GetExtendedRead( read_id, start_ext, stop_ext, &extend_read );
    // in the new extended read
    for ( size_t i = 0; i < r_starts.size(); ++i ) 
        r_starts[i] += start_ext;
    for ( size_t i = 0; i < r_stops.size(); ++i ) 
        r_stops[i] += start_ext;

    // pile with markuped reads
    vec<vec<ReadError> > read_error_vec;
    for ( size_t i = 0; i < friend_ids.size(); ++i ) {
        basevector fread = reads_[ friend_ids[i] ];
        basevector fread_truth = truth_[ friend_ids[i] ];
        if ( p_friend_teller_->IsRc( read_id, friend_ids[i] ) ) {
                fread.ReverseComplement();
                fread_truth.ReverseComplement();
        }
        vec<ReadError> errors;
        FindCanonicalReadError( fread_truth, fread, &errors );
        for ( size_t j = 0; j < errors.size(); ++j ) 
            errors[j].setLocation( errors[j].getLocation() + r_starts[i] );
        read_error_vec.push_back( errors );
    }
    vec<String> pile = PileupReadErrors( extend_read, read_error_vec, r_starts, r_stops );
    return pile;
}

