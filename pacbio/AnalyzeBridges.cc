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
    input_t(const String& sReads, const String&sAnchors,const String&sRefs)
        :mReads(sReads+".fastb")
        ,mQuals(sReads+".qualb")
        {
        FetchReads(mAnchors, 0, sAnchors+".fasta");
        ForceAssert(mAnchors.size() == 2);
        FetchReads(mRefs, 0, sRefs+".fasta");
        ForceAssert(mRefs.size() == 1);
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

    const basevector& ref()const { return mRefs[0]; } ;
    basevector& ref() { return mRefs[0]; } ;
private:
    vecbasevector mReads;
    vecqualvector mQuals;
    vecbasevector mAnchors;
    vecbasevector mRefs;
};

class analysis_t {
public:
    typedef typename basevector::size_type pos_t;
    typedef std::tuple<int,int,pos_t> align_t;
    analysis_t(const String&in_head_, const String&anchors_, const String&refs_,const size_t K)
    : mInput(in_head_,anchors_,refs_), mK(K), mIN_HEAD(in_head_), mANCHORS(anchors_), mOUT_HEAD(mIN_HEAD.SafeAfterLast("/")+"."+mANCHORS.SafeAfterLast("/")+"."+ToString(mK))
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
    void run(size_t nSamples,int cutoff,bool bUseRef) {
        anchorReads();
        buildAlignmentMatrix(nSamples);
        if(bUseRef){
            calcConsensusRef(cutoff);
        }
        else{
            calcConsensus(cutoff);
        }
    }
    void report() {
        size_t nReads = mInput.reads().size();
        std::cout<<std::endl;
        std::cout << "Results of anchoring " << nReads << " reads to " << kmer_anchor.size() << " anchoring sequences:" << std::endl;
        std::map<size_t, size_t> counter;
        size_t nGoodMatch=0;
        size_t nBadMatch=0;
        for (size_t rr = 0; rr < nReads; ++rr) {
            const auto& buffer=read_match[rr];
            if (buffer.size() > 0) {
                bool bGoodMatch = buffer.size()==2 && (get<0>(buffer[0])>0 == get<0>(buffer[1])>0);
                for(const auto& entry: buffer){
                    int a = get<0>(entry)-1;
                    int b = abs(get<1>(entry));
                    int c = get<2>(entry);
                    if(a>=0){
                        for(size_t ss=0;ss<mK;++ss){
                            ForceAssert( mInput.read(rr)[c+ss]==mInput.anchors()[a][b+ss]   );
                        }
                    }
                }
                if(bGoodMatch){
                    if ( get<2>(buffer[0]) < get<2>(buffer[1])){
                        ++nGoodMatch;
                    }
                    else{
                        ++nBadMatch;
                    }
                }
            }
            counter[buffer.size()]+=1;
        }
        std::cout<<std::endl;
        std::cout<<"anchor_match number_of_reads"<<std::endl;
        for (const auto& entry : counter) {
            std::cout << entry.first << " " << entry.second << std::endl;
        }
        std::cout<<std::endl;
        std::cout <<"good match: " << nGoodMatch << std::endl;
        std::cout <<"bad match:  " << nBadMatch << std::endl;

        std::cout<<std::endl;
        ForceAssert( mBridgeReads.size() == size_t(sqrt(mAlignmentMatrix.size())+0.5));
        const size_t nN = mBridgeReads.size();
        std::cout <<"number of all-to-all candidates is " << nN << "/" << bridge_read.size() <<  std::endl;
        for(size_t nn=0;nn<nN*nN;++nn){
            const size_t row=nn/nN;
            const size_t col=nn%nN;
            if(col==0) std::cout << "matrix ";
            std::cout << mAlignmentMatrix[nn];
            if(col == nN-1){ std::cout << std::endl; }
            else{ std::cout << " "; }
        }
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
    void buildAlignmentMatrix(size_t nSample=std::numeric_limits<size_t>::max()){
        const String sCacheFile=mOUT_HEAD+".cache";
        std::cout << "operating with cache file: " << sCacheFile << std::endl;

        struct loc_hash{
            size_t operator()(const std::pair<size_t,size_t>& in)const {
                return std::hash<size_t>()(in.first*in.second);
            };
        };
        std::unordered_map<std::pair<size_t,size_t>,int,loc_hash> cache;
        if(IsRegularFile(sCacheFile)){
            std::cout << Date() << ": loading cache file " << sCacheFile << std::endl;
            ifstream ifs(sCacheFile);
            size_t i,j;
            int k;
            for(ifs>>i>>j>>k ; ifs.good() ; ifs>>i>>j>>k){
                cache[ make_pair(i,j) ] = k;
            }
            std::cout << Date() << ": done" << std::endl;
        }

        const size_t nN = min(bridge_read.size(),nSample);
        if(nN != bridge_read.size()){
            std::random_shuffle(bridge_read.begin(),bridge_read.end());
        }

        mAlignmentMatrix.clear(); mAlignmentMatrix.resize(nN*nN);
        unsigned int longest=0;
        mBridgeReads.clear();
        mBridgeReads.resize(nN);
        mBridgeQuals.clear();
        mBridgeQuals.resize(nN);
        for(size_t nn=0;nn<nN;++nn){
            paddedRead(mBridgeReads[nn],mBridgeQuals[nn],bridge_read[nn]);
            longest = max(longest,mBridgeReads[nn].size());
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


        std::cout << "longest read length: " << longest << std::endl;
        std::cout << Date() << ": all-to-all alignment between " << nN << " reads"<< std::endl;
        #pragma omp parallel
        {
            #pragma omp master
            {
                std::cout << "number of threads " << omp_get_num_threads()<<std::endl;
            }

            SmithWatAffineResource_t resource(longest+1,longest+1);
            #pragma omp for
            for(size_t nn=0;nn<nN*nN;++nn){

                const size_t row=nn/nN;
                const size_t col=nn%nN;
                //just do everything for now for measurement
    //            if( col < row) continue;

                auto itr = cache.find(std::make_pair(bridge_read[row],bridge_read[col]));
                if(itr == cache.end()){
                    const basevector& left=mBridgeReads[row];
                    const basevector& right=mBridgeReads[col];

//                    mAlignmentMatrix[nn]=getAlignment(left,right);
                    mAlignmentMatrix[nn]=getAlignment(left,right,resource);
                }
                else{
                    mAlignmentMatrix[nn]=(*itr).second;
                }
//                mAlignmentMatrix[col*nN+row]=mAlignmentMatrix[nn];
            }
        }

        {
            ofstream ofs(sCacheFile,std::ofstream::out);//|std::ofstream::app);
            for(size_t nn=0;nn<nN*nN;++nn){
                const size_t row=nn/nN;
                const size_t col=nn%nN;
                ofs << bridge_read[row] << " " << bridge_read[col] << " " << mAlignmentMatrix[nn] << std::endl;
            }
            ofs.close();
        }

        std::cout << Date() << ": done" << std::endl;
    }
    void calcConsensusRef(int iCutOff){
        const int64_t nN = mBridgeReads.size();
        unsigned int longest=mInput.ref().size();

        vec<int> vMinScore(nN,std::numeric_limits<int>::max());
        vec<size_t> vIndex(nN,std::numeric_limits<int>::max());
        #pragma omp parallel
        {
            SmithWatAffineResource_t resource(longest+1,longest+1);
            #pragma omp for
            for(int64_t row=0;row<nN;++row){
                vIndex[row]=row;
                vMinScore[row] = getAlignment(mBridgeReads[row],mInput.ref(),resource);
            }
        }
        SortSync(vMinScore,vIndex);
        std::cout << "sorted score:";
        for(int64_t row=0;row<nN;++row){
            std::cout << " " << vMinScore[row] << "," << vIndex[row];
        }
        std::cout << std::endl;


        size_t iEnd = std::distance(vMinScore.begin(),std::upper_bound(vMinScore.begin(),vMinScore.end(),iCutOff));
        std::cout << "cutoff score: " << iCutOff << std::endl;
        std::cout << "cutoff element: " << iEnd << std::endl;

        vecbasevector consensus(1,mInput.ref());

        ofstream ofs(mOUT_HEAD+".consensus.fasta",std::ofstream::out|std::ofstream::app);

        const size_t nCandidates=iEnd;
        vecbasevector candidates;candidates.reserve(nCandidates);

        for( size_t ii=0 ; ii< iEnd ; ++ii){
            candidates.push_back( mBridgeReads[vIndex[ii]] );
        }
        std::cout << "nReads= " << candidates.size() << std::endl;
//        consensus_compute_padded(candidates,&consensus[0],1,1,1,1,"");
        consensus_compute(candidates,&consensus[0],30,1,"");
        ofs << ">consensus_ref" << std::endl;
        ofs << consensus[0].ToString() << std::endl;
        ofs.close();
    }

    void calcConsensus(int iCutOff,const size_t maxCandidates=100){
        const size_t nN = mBridgeReads.size();
        const size_t nCandidates=min(maxCandidates,nN);

        for(size_t row=0; row<nN;++row){
			for(size_t col=row+1; col<nN;++col){
			    const int tmp = (mAlignmentMatrix[row*nN+col] + mAlignmentMatrix[col*nN+row])/2;
			    mAlignmentMatrix[row*nN+col] = tmp;
			    mAlignmentMatrix[col*nN+row] = tmp;
			}
        }

        vec<int> vMinScore(nN,std::numeric_limits<int>::max());
        vec<size_t> vIndex(nN,std::numeric_limits<int>::max());
        for(size_t row=0;row<nN;++row){
            vIndex[row]=row;

            vec<int> tmp(mAlignmentMatrix.begin() + row*nN , mAlignmentMatrix.begin() + row*nN + nN);
            std::sort(tmp.begin(),tmp.end());
            vMinScore[row]=std::accumulate(tmp.begin(),tmp.begin()+nCandidates,0);

/*
            for(int64_t col=0; col<row;++col){
                vMinScore[row] = min(vMinScore[row],mAlignmentMatrix[row*nN+col]);
            }
            for(int64_t col=row+1; col<nN;++col){
                vMinScore[row] = min(vMinScore[row],mAlignmentMatrix[row*nN+col]);
            }
            */

        }
        auto vMinScoreOrg=vMinScore;
        SortSync(vMinScore,vIndex);
        std::cout << "sorted score:";
        for(size_t row=0;row<nN;++row){
            std::cout << " " << vMinScore[row] << "," << vIndex[row];
        }
        std::cout << std::endl;



//        size_t iEnd = std::distance(vMinScore.begin(),std::upper_bound(vMinScore.begin(),vMinScore.end(),iCutOff));
//        std::cout << "cutoff score: " << iCutOff << std::endl;
//        std::cout << "cutoff element: " << iEnd << std::endl;

        size_t nThreads=0;
        #pragma omp parallel
        { if(omp_get_thread_num()==0){ nThreads=omp_get_num_threads(); } }

        //mBridgeReads.WriteAll(mOUT_HEAD+".consensus_bridges.fastb");

        vecbasevector consensus(nThreads,mInput.ref());
        /*
        vecbasevector consensus;
        consensus.reserve(nThreads);
        for(size_t tt=0;tt<nThreads;++tt){
            consensus.push_back(mBridgeReads[vIndex[tt]]);
        }
        */

        ofstream ofs(mOUT_HEAD+".consensus.fasta",std::ofstream::out|std::ofstream::app);
        #pragma omp parallel for schedule(dynamic,1)
        for(size_t cc=0;cc<vIndex.size();++cc)
        {
			size_t row = vIndex[cc];

			vec< std::pair<int,int> > to_sort; to_sort.reserve(nN);
            for(size_t col=0; col<nN;++col){ to_sort.push_back( make_pair(mAlignmentMatrix[row*nN+col],col) ); }
            std::sort(to_sort.begin(),to_sort.end());
            ForceAssert( to_sort[0].first == 0 );
            ForceAssert( to_sort[0].second == int(row) );

            {
                ofstream ofs(mOUT_HEAD+".consensus."+ToString(cc)+".ranked_neighbors",std::ofstream::out|std::ofstream::app);
                for(const auto& entry: to_sort){ ofs << entry.first << " " << entry.second << "\n"; }
                ofs.close();
            }

            if( cc >= nThreads) continue;

			vecbasevector candidates;candidates.reserve(nCandidates);
			int iLast=0;
            #pragma omp critical
            {
                std::cout << "candidates: " ;
                for( size_t ii=0 ; ii< nCandidates && to_sort[ii].first <= iCutOff; ++ii){
                    iLast=to_sort[ii].first;
                    std::cout << " " << to_sort[ii].second;
                    candidates.push_back( mBridgeReads[to_sort[ii].second] );
                }
                std::cout << std::endl;
                std::cout << "consensus " << cc << " getting consensus from a batch of " << candidates.size() << " iLast " << iLast<< std::endl;
                std::cout << "consensus " << cc << " starting consensus length " << consensus[cc].size() << std::endl;
            }


/*

            for(auto& entry: candidates){
                entry.ReverseComplement();
            }
            consensus[id].ReverseComplement();
            consensus_compute(candidates,&consensus[id],30,1,"");
//            consensus_compute_padded(candidates,&consensus[id],1,1,1,1,"");

            for(auto& entry: candidates){
                entry.ReverseComplement();
            }
            consensus[id].ReverseComplement();
*/
            consensus_compute(candidates,&consensus[cc],30,0,"");
//            consensus_compute_padded(candidates,&consensus[id],1,1,1,1,"");
            {
                ofstream ofs(mOUT_HEAD+".consensus."+ToString(cc)+".fasta",std::ofstream::out|std::ofstream::app);
                ofs << ">consensus_" << cc << std::endl;
                ofs << consensus[cc].ToString() << std::endl;
                ofs.close();
            }

            #pragma omp critical
            {
                ofs << ">consensus_" << cc << std::endl;
                ofs << consensus[cc].ToString() << std::endl;
            }
        }
        ofs.close();

    }
    input_t mInput;
    const size_t mK;
    std::unordered_map<kmer_val_t, std::pair<int,int>> kmer_anchor;
    vec<vec<align_t>> read_match;

    vec<uint64_t> bridge_read;
    vecbasevector mBridgeReads;
    vecqualvector mBridgeQuals;
    vec<int> mAlignmentMatrix;

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
    CommandArgument_String_Doc(REF, "Looks for REF.fasta for reference");
    CommandArgument_Int_OrDefault(CUTOFF, 0);
    CommandArgument_UnsignedInt_OrDefault(K, 20u);
    CommandArgument_Bool_OrDefault(USE_REF, False);
    CommandArgument_UnsignedInt_OrDefault(NSAMPLES, std::numeric_limits<unsigned int>::max());

    EndCommandArguments;

    ForceAssert(K <= 32);

    analysis_t analysis(IN_HEAD,ANCHORS,REF, K);

    analysis.run(NSAMPLES,CUTOFF,USE_REF);
//    analysis.report();

    std::cout << Date() << ": Done." << std::endl;
}
