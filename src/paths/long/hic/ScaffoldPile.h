///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Sep 2, 2014 - <crdhelp@broadinstitute.org>
//

// efficient structure for holding a list of scaffolds (chains of line ids)
// from which we intend to modify, merge, insert, and delete.

#ifndef _SCAFFOLDPILE_H
#define _SCAFFOLDPILE_H

#include <unordered_map>

class ScaffoldPile {
public:

    // -1 is marker for the end
    // adjacency is a pair of left, right "pointers" to the previous and next
    // line in the scaffold
    struct adjacency : public pair<int,int> {
        adjacency() { this->first = -1; this->second = -1; }
        int left() const { return this->first; }
        int right() const { return this->second; }
        int& left() { return this->first; }
        int& right() { return this->second; }
        bool single() const { return this->left() == -1 && this->right() == -1; }
        bool start() const { return this->left() == -1; }
        bool end() const { return this->right() == -1; }
        adjacency(int left, int right) { this->first = left; this->second = right; }
    };

    // ** the fundamental, underlying datastructure
    typedef std::unordered_map<int, adjacency > PileType;

    ScaffoldPile() = delete;
    explicit ScaffoldPile(vec<int> lineInv) : mLineInv(lineInv) {};

    // return a copy of the scaffold containing lineId
    vec<int> operator()( const int lineId ) const { return scaffold(lineId); }

    // check the structure for consistency
    void check(bool fatal=false) const {
        bool bad = false;
        for ( auto itr = mPile.begin(); !bad && itr != mPile.end(); ++itr ) {
            auto const& key = itr->first;
            auto const& adj = itr->second;

            if ( adj.first != -1 ) {
                // check left
                if ( !this->has_key(adj.first) ) {
                    bad = true;
                    cout << "key " << key << " has unknown left adjacency "
                            << adj.first << endl;
                } else if ( mPile.at(adj.first).second != key ) {
                    bad = true;
                    cout << "key " << key << " has left adjacency " << adj.first <<
                            ", but that node has right adjancency of " <<
                            mPile.at(adj.first).second << endl;
                }
            }

            if ( adj.second != -1 ) {
                // check right
                if ( !this->has_key( adj.second) ) {
                    bad = true;
                    cout << "key " << key << " has unknown right adjacency "
                            << adj.second << endl;
                } else if ( mPile.at(adj.second).first != key ) {
                    bad = true;
                    cout << "key " << key << " has right adjacency " << adj.second <<
                            ", but that node has left adjacency of " <<
                            mPile.at(adj.second).first << endl;
                }
            }

            if ( bad ) {
                cout << "inverse of key " << key << " is "
                        << mLineInv[key] << endl;
            }
        }

        if ( bad && fatal ) {
            cout << Date() << ": terminating" << endl;
            TracebackThisProcess();
        }
    }

    vec<int> singletons() const {
        // return a list of singletons
        vec<int> lineIds;
        for ( auto itr = mPile.begin(); itr != mPile.end(); ++itr )
            if ( itr->second.single() ) lineIds.push_back(itr->first);
        return lineIds;
    }

    vec<int> starts() const {
        vec<int> starts;
        for ( auto itr = mPile.begin(); itr != mPile.end(); ++itr )
            if ( itr->second.start() ) starts.push_back(itr->first);
        return starts;
    }

    int start(const int lineId) const {
//        cout << "call find in start" << endl;
        auto itr = mPile.find(key(lineId));
        while ( itr->second.left() != -1 ) itr = mPile.find(itr->second.left());
        return itr->first;
    }

    vec<int> scaffold( const int lineId ) const {
        vec<int> scaffold;
//        cout << "call find in scaffold" << endl;
        auto itr = mPile.find(key(lineId));
        ForceAssert(itr != mPile.end());
        while ( itr->second.left() != -1 ) {
            itr = mPile.find(itr->second.left());
            if ( itr == mPile.end() ) this->check(true);
        }
        scaffold.push_back( itr->first );
        while ( itr->second.right() != -1 ) {
            itr = mPile.find(itr->second.right());
            if ( itr == mPile.end() ) this->check(true);
            scaffold.push_back(itr->first);
        }

        return scaffold;
    }

    vec<vec<int>> scaffolds() const {
        vec<vec<int>> scaffolds;
        for ( auto itr = mPile.begin(); itr != mPile.end(); ++itr )
            if ( itr->second.start() )
                scaffolds.push_back( scaffold(itr->first) );
        return scaffolds;
    }

    int Max() const {
        return mLineInv.size();
    }

    // merge the scaffolds contining these lines
    void merge(vec<int> const& lineIds ) {
        if ( lineIds.size() < 2 ) return;       // nothing to do

        vec<int> starts;
        for ( auto lineId : lineIds ) {
            int lineId_start = start(lineId);
            if ( Member(starts, lineId_start) ) FatalErr("a loop in merging");
            starts.push_back(lineId_start);
        }

        cout << "building scaffs" << endl;
        vec<vec<int>> scaffs( lineIds.size() );
        for ( size_t i = 0; i < lineIds.size(); ++i )
            scaffs[i] = this->scaffold(lineIds[i]);
        cout << "checking ends" << endl;
        for ( size_t i = 0; i < scaffs.size(); ++i ) {
            ForceAssert( mPile.at(scaffs[i].front()).start() );
            ForceAssert( mPile.at(scaffs[i].back()).end() );
        }
        cout << "checking seen" << endl;
        vec<bool> seen( this->Max(), false );
        for ( size_t i = 0; i < scaffs.size(); ++i )
            for ( size_t j = 0; j < scaffs[i].size(); ++j ) {
                if ( seen[ scaffs[i][j] ] )
                    FatalErr("attempt to merge scaffolds that intersect");
                seen[scaffs[i][j]] = true;
            }
        // check that the beginning of each scaffold is -1
        // as is the end
        cout << "merging..." << endl;
        for ( size_t i = 1; i < scaffs.size(); ++i )
            this->join(scaffs[i-1].back(), scaffs[i].front());
    }

    // join two ends
    void join( int left, int right ) {
        ForceAssertEq(mPile.at(left).right(),-1);
        ForceAssertEq(mPile.at(right).left(),-1);
        mPile.at(left).right() = right;
        mPile.at(right).left() = left;
    }

    // RC the scaffold containing lineId
    vec<int> rc(int lineId) {
        auto sc = del(lineId);   // delete the scaffold, but remember it in s
        vec<int> rc_sc(sc.rbegin(), sc.rend());
        for ( auto& line : rc_sc ) line = mLineInv[line];
        add(rc_sc);
        return rc_sc;
    }

    vec<int> del(int lineId) {
        auto s = scaffold(lineId);
        ForceAssert(mPile.at(s.front()).start());
        ForceAssert(mPile.at(s.back()).end());
        for ( auto const line : s ) {
            ForceAssert(!this->has_key(mLineInv[line]));
            mPile.erase(line);
        }
        return s;       // return the deleted scaffold, as a convenience
    }

    void add(vec<int> const& scaff) {
        for ( auto const line : scaff ) {
            ForceAssert(!this->has_key(line));
            ForceAssert(!this->has_key(mLineInv[line]));
        }
        size_t n = scaff.size();
        for ( size_t i = 0; i < n; i++ ) {
            int left = (i > 0) ? scaff[i-1] : -1;
            int right = ( i+1 < n ) ? scaff[i+1] : -1;
            mPile[scaff[i]] = adjacency(left,right);
        }
    }

    bool has_key( int lineId ) const {
//        cout << "START call find in has_key for lineId " << lineId << endl;
        auto ret = mPile.find(lineId) != mPile.end();
//        cout << "FINIS call find in has_key for lineId " << lineId << endl;
        return ret;
    }

    bool has_key_or_inv( int lineId ) const {
//        cout << "call find in has_key_or_inv" << endl;
        return ( mPile.find(lineId) != mPile.end() ||
                 mPile.find(inv(lineId)) != mPile.end() );
    }

    int inv(int lineId ) const { return mLineInv[lineId]; }

    // key - normalize the lineId by returning the corresponding
    // element of { lineId, inv[lineId] } that is present in the
    // current pile.  Otherwise abort.
    int key(int lineId) const {
        if ( has_key(lineId) ) return lineId;
        else if ( has_key(inv(lineId)) ) return inv(lineId);
        else FatalErr("BUG: tried to pull a scaffold line that doesn't exist");
    }

    void clear() { mPile.clear(); }

    PileType::iterator begin() { return mPile.begin(); }
    PileType::iterator end() { return mPile.end(); }

private:
    PileType mPile;
    vec<int> mLineInv;
};

#endif  //_SCAFFOLDPILE_H
