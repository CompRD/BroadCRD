/*
 * AnalyzeBridges.cc
 *
 *  Created on: Dec 9, 2013
 *      Author: blau
 */

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
#include <limits>
#include <fstream>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <unordered_map>
#include <functional>
#include <iterator>
#include <omp.h>

#include "Basevector.h"
#include "FetchReads.h"
#include "PrintAlignment.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "paths/LongReadPatchOptimizer.h"


namespace {
    // a much faster hack of SmithWatAffine in SmithWatAffine.cc
    class SmithWatAffineResource_t{
    public:
        SmithWatAffineResource_t(size_t n1,size_t n2):mN1(n1),mN2(n2){this->reset(n1,n2);};
        void reset(size_t n1,size_t n2){
            mN1=n1;
            mN2=n2;
            sx.assign(mN1*mN2,0);
            sy.assign(mN1*mN2,0);
            sz.assign(mN1*mN2,0);
            xf.assign(mN1*mN2,0);
            yf.assign(mN1*mN2,0);
            zf.assign(mN1*mN2,0);
        }

        unsigned int& score_x(size_t row,size_t col){return sx[row*mN2+col];};
        unsigned int& score_y(size_t row,size_t col){return sy[row*mN2+col];};
        unsigned int& score_z(size_t row,size_t col){return sz[row*mN2+col];};

        unsigned char& x_from(size_t row,size_t col){return xf[row*mN2+col];};
        unsigned char& y_from(size_t row,size_t col){return yf[row*mN2+col];};
        unsigned char& z_from(size_t row,size_t col){return zf[row*mN2+col];};

        vec<unsigned char>& x_from(){return xf;};
        vec<unsigned char>& y_from(){return yf;};
        vec<unsigned char>& z_from(){return zf;};
        size_t n1()const{return mN1;};
        size_t n2()const{return mN2;};

    private:
        vec<unsigned int> sx;
        vec<unsigned int> sy;
        vec<unsigned int> sz;

        vec<unsigned char> xf;
        vec<unsigned char> yf;
        vec<unsigned char> zf;
        size_t mN1,mN2;
    };
unsigned int SmithWatAffine_loc( const basevector& S, const basevector& T, alignment& a
                               , SmithWatAffineResource_t& resource
                               , bool penalize_left_gap, bool penalize_right_gap, const int mismatch_penalty, const int gap_open_penalty, const int gap_extend_penalty )
{
     const int Infinity = 100000000;
     ForceAssertGt( S.size(), 0u );
     ForceAssertGt( T.size(), 0u );

     unsigned int n = S.size( ), N = T.size( );

     //     ForceAssertLe( n, N );

     avector<char> s, t;
     s.resize(n);
     for ( unsigned int i = 0; i < n; i++ )
          s(i) = S[i];
     t.resize(N);
     for ( unsigned int i = 0; i < N; i++ )
          t(i) = T[i];

     int best_score = Infinity;
     resource.reset(n+1,N+1);

     resource.score_x(0,0) = 0;
     resource.score_y(0,0) = Infinity;
     resource.score_z(0,0) = Infinity;
     resource.x_from(0,0) = 's';
     resource.y_from(0,0) = 's';
     resource.z_from(0,0) = 's';

     for ( unsigned int i = 1; i <= n; i++ )
     {    resource.score_x(i,0) = Infinity;
          resource.score_y(i,0) = Infinity;
          resource.score_z(i,0) = gap_open_penalty + gap_extend_penalty * i;
          resource.x_from(i,0) = 's';
          resource.y_from(i,0) = 's';
          resource.z_from(i,0) = 's';   }

     for ( unsigned int j = 1; j <= N; j++)
       {  resource.score_x(0,j) = Infinity;
          resource.score_y(0,j) = (penalize_left_gap ? gap_open_penalty + gap_extend_penalty * j : 0);
          resource.score_z(0,j) = Infinity;
          resource.x_from(0,j) = 's';
          resource.y_from(0,j) = 's';
          resource.z_from(0,j) = 's';   }

     for ( unsigned int i = 1; i <= n; i++ )
     {   for ( unsigned int j = 1; j <= N; j++ )
          {    unsigned int x_x = resource.score_x(i-1,j-1) + mismatch_penalty * ( s(i-1) != t(j-1) );
               unsigned int x_y = resource.score_y(i-1,j-1) + mismatch_penalty * ( s(i-1) != t(j-1) );
               unsigned int x_z = resource.score_z(i-1,j-1) + mismatch_penalty * ( s(i-1) != t(j-1) );
               unsigned int y_x = resource.score_x(i,j-1) + (i != n || penalize_right_gap ? gap_open_penalty : 0);
               unsigned int y_y = resource.score_y(i,j-1) + (i != n || penalize_right_gap ? gap_extend_penalty : 0);
               unsigned int y_z = Infinity; //score_z[i][j-1] + gap_open_penalty;
               unsigned int z_x = resource.score_x(i-1,j) + gap_open_penalty;
               unsigned int z_y = Infinity; //score_y[i-1][j] + gap_open_penalty;
               unsigned int z_z = resource.score_z(i-1,j) + gap_extend_penalty;

               resource.score_x(i,j) = Min( Min( x_x, x_y ), x_z );
               resource.score_y(i,j) = Min( Min( y_x, y_y ), y_z );
               resource.score_z(i,j) = Min( Min( z_x, z_y ), z_z );

               if ( x_x <= x_y )
               {    if ( x_x <= x_z ) resource.x_from(i,j) = 'x';
                    else resource.x_from(i,j) =  'z';    }
               else
               {    if ( x_y <= x_z ) resource.x_from(i,j) = 'y';
                    else resource.x_from(i,j) =  'z';    }

               if ( y_x <= y_y )
               {    if ( y_x <= y_z ) resource.y_from(i,j) = 'x';
                    else resource.y_from(i,j) =  'z';    }
               else
               {    if ( y_y <= y_z ) resource.y_from(i,j) = 'y';
                    else resource.y_from(i,j) =  'z';    }

               if ( z_x <= z_y )
               {    if ( z_x <= z_z ) resource.z_from(i,j) = 'x';
                    else resource.z_from(i,j) =  'z';    }
               else
               {    if ( z_y <= z_z ) resource.z_from(i,j) = 'y';
                    else resource.z_from(i,j) =  'z';    }    }    }

     best_score = Min( resource.score_x(n,N), Min( resource.score_y(n,N), resource.score_z(n,N) ) );

//     vec< vec<unsigned char> > *from;
     vec<unsigned char> *from;
     if ( resource.score_x(n,N) <= resource.score_y(n,N) )
     {    if ( resource.score_x(n,N) <= resource.score_z(n,N) ) from = &resource.x_from();
          else from = &resource.z_from();    }
     else
     {    if ( resource.score_y(n,N) <= resource.score_z(n,N) ) from = &resource.y_from();
          else from = &resource.z_from();    }

     int i = int(n);
     int j = int(N);
     int lcount = 0, g1count = 0, g2count = 0;
     int last_length = 0;
     avector<int> gaps(0), lengths(0);
     while(1)
     {
//          unsigned char dir = (*from)[i][j];
          unsigned char dir = (*from)[i*resource.n2()+j];
          //cout << dir;
          if ( from == &resource.x_from() )
          {    if ( g1count > 0 )
               {    if ( last_length > 0 )
                    {    gaps.Prepend( g1count );
                         lengths.Prepend( last_length );    }
                    g1count = 0;    }
               if ( g2count > 0 )
               {    if ( last_length > 0 )
                    {    gaps.Prepend( -g2count );
                         lengths.Prepend( last_length );    }
                    g2count = 0;    }
               ++lcount;
               --i;
               --j;   }
          else if ( from == &resource.z_from() )  // gap on long sequence
          {    if ( lcount > 0 )
               {    last_length = lcount;
                    lcount = 0;    }
               ForceAssert( g1count == 0 );
               ++g2count;
               --i;    }
          else                           // gap on short sequence
          {    if ( lcount > 0 )
               {    last_length = lcount;
                    lcount = 0;    }
               ForceAssert( g2count == 0 );
               ++g1count;
               --j;    }

          if ( dir == 'x') from = &resource.x_from();
          else if ( dir == 'y') from = &resource.y_from();
          else from = &resource.z_from();

          if( (*from)[i*resource.n2()+j] == 's' ) break;

    }

     //cout << "\n";

     if ( g1count != 0 ) gaps.Prepend( g1count );
     else if ( g2count != 0 ) gaps.Prepend( -g2count );
     else gaps.Prepend(0);

     lengths.Prepend( lcount );

     int pos1 = i;
     int pos2 = j;

     if ( gaps(0) < 0 )
     {   pos2 -= gaps(0);
         gaps(0) = 0;    }

     if ( gaps(0) > 0 )
     {   pos1 += gaps(0);
         gaps(0) = 0;    }

     int errors = best_score;
     a = alignment( pos1, pos2, errors, gaps, lengths );

     return best_score;    }

typedef uint64_t kmer_val_t;
kmer_val_t KmerVal(const basevector&b, size_t K, size_t p = 0) {
    ForceAssert(K <= 32);
    ForceAssert(p + K <= b.size());
    kmer_val_t out = 0;
    for (size_t ii = 0; ii < K; ++ii) {
        out <<= 2;
        out += b[p + ii];
    }
    return out;
}

class kmerizer_t {
public:
    kmerizer_t(const basevector&b, const typename basevector::size_type K)
        : mSeq(b), mK(K), mMask(pow(4, K - 1)), mNext(K) {
        mValid = mK - 1 < mSeq.size();
        if (mValid) mVal = KmerVal(mSeq, mK);
    } ;
    kmerizer_t& operator++() {
        mValid = mNext < mSeq.size();
        if (mValid) {
            mVal %= mMask;
            mVal <<= 2;
            mVal += mSeq[mNext];
//            ForceAssert( mVal == KmerVal(mSeq,mK,mNext+1-mK));
            ++mNext;
        }
        return *this;
    }
    kmer_val_t val() const {
        ForceAssert(mValid);
        return mVal;
    }
    ;
    typename basevector::size_type pos() const { return mNext - mK; } ;
    bool valid() const { return mValid; } ;
private:
    kmerizer_t();
    const basevector& mSeq;
    const typename basevector::size_type mK;
    const kmer_val_t mMask;
    kmer_val_t mVal;
    typename basevector::size_type mNext;
    bool mValid;
};

class input_t {
public:
    input_t(const String& sReads, const String&sAnchors)
        :mReads(sReads+".fastb")
        ,mQuals(sReads+".qualb")
        {
        FetchReads(mAnchors, 0, sAnchors+".fasta");
        ForceAssert(mAnchors.size() == 2);
    }
    const vecbasevector& reads() const { return mReads; } ;
    vecbasevector& reads() { return mReads; } ;
    const basevector& read(size_t ii) const { return mReads[ii]; } ;
    basevector& read(size_t ii) { return mReads[ii]; } ;

    const vecqualvector& quals() const { return mQuals; } ;
    vecqualvector& quals() { return mQuals; } ;
    const qualvector& qual(size_t ii) const { return mQuals[ii]; } ;
    qualvector& qual(size_t ii) { return mQuals[ii]; } ;

    const vecbasevector& anchors() const { return mAnchors; } ;
    vecbasevector& anchors() { return mAnchors; } ;

private:
    vecbasevector mReads;
    vecqualvector mQuals;
    vecbasevector mAnchors;
};

class analysis_t {
public:
    typedef typename basevector::size_type pos_t;
    typedef std::tuple<int,int,pos_t> align_t;
    analysis_t(const String&in_head_, const String&anchors_, const size_t K)
    : mInput(in_head_,anchors_), mK(K), mIN_HEAD(in_head_), mANCHORS(anchors_), mOUT_HEAD(mIN_HEAD.SafeAfterLast("/")+"."+mANCHORS.SafeAfterLast("/")+"."+ToString(mK))
    {
        const auto& anchors = mInput.anchors();
        ForceAssert(anchors.size() == 2);

        for (kmerizer_t kmerizer(anchors[0], mK); kmerizer.valid(); ++kmerizer) {
            ForceAssert(kmer_anchor.find(kmerizer.val()) == kmer_anchor.end());
            kmer_anchor[kmerizer.val()] = std::make_pair(1,-(int)kmerizer.pos());
        }
        for (kmerizer_t kmerizer(anchors[1], mK); kmerizer.valid(); ++kmerizer) {
            ForceAssert(kmer_anchor.find(kmerizer.val()) == kmer_anchor.end());
            kmer_anchor[kmerizer.val()] = std::make_pair(2,(int)kmerizer.pos());
        }
        basevector btmp=anchors[0];
        for (kmerizer_t kmerizer(btmp.ReverseComplement(), mK); kmerizer.valid(); ++kmerizer) {
            ForceAssert(kmer_anchor.find(kmerizer.val()) == kmer_anchor.end());
            kmer_anchor[kmerizer.val()] = std::make_pair(-1,(int)kmerizer.pos());
        }
        btmp=anchors[1];
        for (kmerizer_t kmerizer(btmp.ReverseComplement(), mK); kmerizer.valid(); ++kmerizer) {
            ForceAssert(kmer_anchor.find(kmerizer.val()) == kmer_anchor.end());
            kmer_anchor[kmerizer.val()] = std::make_pair(-2,-(int)kmerizer.pos());
        }
    }
    void run(vec<double>& results,const vecbasevector& candidates){
        anchorReads();
        vote(results,candidates);
    }
    void vote(vec<double>& results, const vecbasevector& candidates){
        const int64_t nCandidates=candidates.size();
        results.clear();
        results.resize(nCandidates,0);

        vec<basevector> paddedCandidates(nCandidates);
        {
            const basevector& left = mInput.anchors()[0];
//            std::cout << "left:" << std::endl;
//            std::cout << left.ToString() << std::endl;
            const basevector& right = mInput.anchors()[1];
//            std::cout << "right:" << std::endl;
//            std::cout << right.ToString() << std::endl;
            for(int64_t cc=0;cc<nCandidates;++cc){
                paddedCandidates[cc].clear();
                paddedCandidates[cc].reserve(candidates[cc].size()+left.size()+right.size());
                paddedCandidates[cc].append(left);
                paddedCandidates[cc].append(candidates[cc]);
                paddedCandidates[cc].append(right);
            }
        }
//    std::cout << "Padded Candidates:" << std::endl;
//    for(const auto& can:paddedCandidates) std::cout << can.ToString() << std::endl;
//    std::cout << std::endl;

        const size_t nN = bridge_read.size();

        unsigned int longest=0;
        unsigned int shortest = std::numeric_limits<unsigned>:: max();
        mBridgeReads.clear();
        mBridgeReads.resize(nN);
        mBridgeQuals.clear();
        mBridgeQuals.resize(nN);
        for(size_t nn=0;nn<nN;++nn){
            paddedRead(mBridgeReads[nn],mBridgeQuals[nn],bridge_read[nn]);
            longest = max(longest,mBridgeReads[nn].size());
            shortest = min(shortest,mBridgeReads[nn].size());
        }

        {
            mBridgeReads.WriteAll(mOUT_HEAD+".fastb");
            mBridgeQuals.WriteAll(mOUT_HEAD+".qualb");
            ofstream ofs(mOUT_HEAD+".fastb.id");
            for(size_t nn=0;nn<nN;++nn){
                ofs << ">" <<bridge_read[nn]<<"\n";
            }
            ofs.close();
        }

        vec< vec<int> > alignments(nN, vec<int>(nCandidates,-1));


        std::cout << "shortest read length: " << shortest << std::endl;
        std::cout << "longest read length: " << longest << std::endl;
        std::cout << Date() << ": alignment for " << nN << " reads"<< std::endl;
        #pragma omp parallel
        {
            vec<double> results_loc(results.size(),0);
            #pragma omp master
            {
                std::cout << "number of threads " << omp_get_num_threads()<<std::endl;
            }

            SmithWatAffineResource_t resource(longest+1,longest+1);
            #pragma omp for schedule(dynamic,1)
            for(size_t nn=0;nn<nN;++nn){
                if(nn%100==0){
                    #pragma omp critical
                    {
                        std::cout << "nn= " << nn << std::endl;
                    }

                }
                int iBestScore=std::numeric_limits<int>::max();
                vec<int> best_c;
                for(int cc=0;cc<nCandidates;++cc){
                    ForceAssert(mBridgeReads[nn].size()>0);
                    ForceAssert(paddedCandidates[cc].size()>0);
                    alignments[nn][cc]= getAlignment(mBridgeReads[nn],paddedCandidates[cc],resource);
                    if(alignments[nn][cc] < iBestScore){
                        best_c.clear();
                        iBestScore=alignments[nn][cc];
                    }
                    if(alignments[nn][cc]==iBestScore){
                        best_c.push_back(cc);
                    }
                }
                ForceAssert(best_c.size()>0);
                for(const auto&c:best_c){
                    results_loc[c] += 1.0/best_c.size();
                }
            }
            #pragma omp critical
            {
                for(int64_t cc=0;cc<nCandidates;++cc){ results[cc]+=results_loc[cc]; }
            }
        }

        for(size_t nn=0;nn<nN;++nn){
            std::cout << nn << " alignment collection: " ;
            std::copy(alignments[nn].begin(),alignments[nn].end(),ostream_iterator<int>(std::cout," "));
            std::cout << std::endl;
        }

        std::cout << Date() << ": resulting votes:" << std::endl;
        for(const auto& vote:results){ std::cout<<vote<< " " ;}
        std::cout<<std::endl;
        std::cout << Date() << ": done" << std::endl;
    }
private:
    void anchorRead(size_t ii) {
        auto& buffer=read_match[ii];
        buffer.clear();
        for (kmerizer_t kmerizer(mInput.read(ii), mK); kmerizer.valid(); ++kmerizer) {
                auto itr = kmer_anchor.find(kmerizer.val());
                if (itr != kmer_anchor.end()) {
                buffer.push_back( std::make_tuple((*itr).second.first,(*itr).second.second,kmerizer.pos()));
            }
        }
        std::sort(buffer.begin(),buffer.end());
        int last_flag=std::numeric_limits<int>::min();
        size_t n_entry=0;
        for(size_t jj=0;jj<buffer.size();++jj){
            const auto& flag=std::get<0>(buffer[jj]);
            if(flag > last_flag){
                if(jj!=n_entry){ buffer[n_entry]=buffer[jj]; }
                last_flag=flag;
                ++n_entry;
            }
        }
        buffer.resize(n_entry);
    }
    void anchorReads() {
        size_t nReads = mInput.reads().size();
        read_match.clear();
        read_match.resize(nReads);
        std::cout << "Anchoring " << nReads << " reads to " << kmer_anchor.size() << " anchoring sequences." << std::endl;
        #pragma omp parallel for schedule(dynamic,1)
        for (size_t rr = 0; rr < nReads; ++rr) {
            anchorRead(rr);
            const auto& buffer=read_match[rr];
            if(    buffer.size()==2 && get<0>(buffer[0])<0 && get<0>(buffer[1])<0){
                mInput.read(rr).ReverseComplement();
                std::reverse(mInput.qual(rr).begin(),mInput.qual(rr).end());
                anchorRead(rr);
            }

        }
        bridge_read.reserve(read_match.size());
        for (size_t rr = 0; rr < read_match.size(); ++rr) {
            const auto& buffer=read_match[rr];
            if( buffer.size()==2 && (get<0>(buffer[0])>0 == get<0>(buffer[1])>0) && get<2>(buffer[0]) < get<2>(buffer[1])){
                bridge_read.push_back(rr);
            }
        }
        std::cout<<"number of bridges: "<<bridge_read.size()<<std::endl;
    }
    int getAlignment(const basevector&left,const basevector&right,SmithWatAffineResource_t& resource){
        alignment a;
        int i= SmithWatAffine_loc(left,right,a,resource, true,true,1,1,1);
//        PrintVisualAlignment(True,std::cout,left,right,a);
        return i;
    }

    int getAlignment(const basevector&left,const basevector&right){
        alignment a;
        int i= SmithWatAffine(left,right,a, true,true,1,1,1);
//        PrintVisualAlignment(True,std::cout,left,right,a);
        return i;
    }
    void paddedRead(basevector&base,qualvector& qual,size_t rr){
        const auto& specs=read_match[rr];
        const basevector& left = mInput.anchors()[0];
        const basevector& right = mInput.anchors()[1];
        ForceAssert( specs.size()==2 && get<0>(specs[0])==1 && get<0>(specs[1])==2 && get<2>(specs[0])<get<2>(specs[1]));
        {
            const basevector& org = mInput.read(rr);
            base.clear();
            base.reserve(org.size()+left.size()+right.size());
            base.append(left.begin(),left.begin()+abs(get<1>(specs[0])));
            base.append(org.begin()+get<2>(specs[0]) ,org.begin()+get<2>(specs[1]));
            base.append(right.begin()+abs(get<1>(specs[1])),right.end());
        }
        {
            const typename qualvector::value_type pad_qual=0;
            const qualvector& org = mInput.qual(rr);
            qual.clear();
            qual.reserve(org.size()+left.size()+right.size());
            qual.append(abs(get<1>(specs[0])),pad_qual);
            qual.append(org.begin()+get<2>(specs[0]) ,org.begin()+get<2>(specs[1]));
            qual.append(right.size()-abs(get<1>(specs[1])),pad_qual);
        }
        ForceAssert(qual.size()==base.size());
//        std::cout << get<2>(specs[1])-get<2>(specs[0]) << " " << out.size() << std::endl;
    }
    input_t mInput;
    const size_t mK;
    std::unordered_map<kmer_val_t, std::pair<int,int>> kmer_anchor;
    vec<vec<align_t>> read_match;

    vec<uint64_t> bridge_read;
    vecbasevector mBridgeReads;
    vecqualvector mBridgeQuals;

    String mIN_HEAD;
    String mANCHORS;
    String mOUT_HEAD;
};

}

#include "MainTools.h"
int main(int argc, char** argv) {
    RunTime();
    BeginCommandArguments;

    CommandArgument_String_Doc(IN_HEAD, "Looks for IN_HEAD.fastb for the reads");
    CommandArgument_String_Doc(ANCHORS, "Looks for ANCHORS.fasta for the anchors");
    CommandArgument_String_Doc(CANDIDATES, "candidates participating with the voting, will be padded with anchors");
    CommandArgument_Int_OrDefault(CUTOFF, 0);
    CommandArgument_UnsignedInt_OrDefault(K, 20u);

    EndCommandArguments;

    ForceAssert(K <= 32);

    analysis_t analysis(IN_HEAD,ANCHORS,K);

    vecbasevector candidates;


    FetchReads(candidates, 0, CANDIDATES+".fasta");
    std::cout << "Candidate lengths:" << std::endl;
    for(const auto& can:candidates) std::cout << can.size() << std::endl;
    std::cout << std::endl;

    vec<double> results;
    analysis.run(results,candidates);

    std::cout << Date() << ": Done." << std::endl;
}
