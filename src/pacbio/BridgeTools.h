/*
 * BridgeTools.h
 *
 *  Created on: Dec 20, 2013
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

#include "Histogram.h"

#include "Basevector.h"
#include "FetchReads.h"
#include "PrintAlignment.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "paths/LongReadPatchOptimizer.h"


namespace pacbio_bridge_tools{
typedef uint64_t kmer_val_t;

typedef std::unordered_map<int64_t,int64_t> pos_count_t;
typedef std::unordered_map<kmer_val_t,pos_count_t> kmer_pos_count_t;

typedef std::unordered_map<kmer_val_t,int64_t> kmer_count_t;
typedef std::unordered_map<int64_t,kmer_count_t> pos_kmer_count_t;

kmer_val_t KmerVal(const basevector&b, size_t K, size_t p = 0);

class kmerizer_t {
public:
    kmerizer_t(const basevector&b, const typename basevector::size_type K)
        : mSeq(b), mK(K), mMask(pow(4, K - 1)), mVal(std::numeric_limits<kmer_val_t>::max()),mNext(K) {
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
                               , bool penalize_left_gap, bool penalize_right_gap, const int mismatch_penalty, const int gap_open_penalty, const int gap_extend_penalty );
void pos1_to_pos2(vec<size_t>& map, const alignment&a);

class kp_analysis_t{
public:
    kp_analysis_t(uint64_t a):mK(a){};


    void processRead(const basevector& read,const qualvector& qual,const vec<size_t>&pos_map
                    ,const typename qualvector::value_type threshold,unsigned int flank){
        for (kmerizer_t kmerizer(read, mK); kmerizer.valid(); ++kmerizer) {
            const auto pos = kmerizer.pos();
            ForceAssert( pos < qual.size());
            ForceAssert( mK>= 2*flank);
            if( *std::min_element(qual.begin()+pos+flank,qual.begin()+pos+mK-flank) >= threshold){
                const auto val = kmerizer.val();
                pkc[pos_map[pos]][val] += 1;
                kpc[val][pos_map[pos]] += 1;
            }
        }
    }
    void processReads(const vecbasevector&reads,const vecqualvector&quals, const vec<vec<size_t>>& pos_maps
                     ,const typename qualvector::value_type threshold,unsigned int flank){
        ForceAssert(reads.size()==quals.size());
        ForceAssert(reads.size()==pos_maps.size());
        for(size_t rr=0;rr<reads.size();++rr){ processRead(reads[rr],quals[rr],pos_maps[rr],threshold,flank); }
    }

    void printKmerDistribution( kmer_val_t v, const String& out_head ="hist" )const{
        auto itr = kpc.find(v);
        if(itr == kpc.end()) return;

        vec< std::pair<int64_t,int64_t>> pos_count;
        for( const auto& p_c: (*itr).second){
            pos_count.emplace_back(p_c.first,p_c.second);
        }
        std::sort(pos_count.begin(),pos_count.end());

        histogram<int> histo;
        int MIN=0,BIN_SIZE=50;
        int MAX = pos_count.back().first+BIN_SIZE;
        vec<int> bins;
        int step = 0;
        while( 1 ) {
            int bin = MIN + step * BIN_SIZE;
            if ( bin > MAX ) break;
            bins.push_back( bin );
            step++;
        }
        histo.AddBins( bins );

        for(const auto& entry: pos_count){
            histo.AddIdenticalData(entry.first,entry.second);
            std::cout << "spatial " << entry.first << " " << entry.second << std::endl;
        }

        {
            using namespace ns_psplot;

            vec<ns_psplot::freetext> labels;

            ofstream eps_out( out_head+".eps");
            histo.PrintAsEps( eps_out, labels, 0 );
            eps_out.close( );
        }

    }

    const kmer_pos_count_t& getKPC()const{return kpc;}
    uint64_t K()const{return mK;};

    friend ostream& operator<<(ostream&os,const kp_analysis_t&in){
        vec< std::tuple<int64_t,int64_t,kmer_val_t> > pos_mcount_kmer;
        for(const auto& p_kc: in.pkc){
            const auto& pos=p_kc.first;
            for(const auto& k_c: p_kc.second){
                pos_mcount_kmer.emplace_back( pos, -k_c.second, k_c.first);
            }
        }
        std::sort(pos_mcount_kmer.begin(),pos_mcount_kmer.end());
        for( const auto& entry: pos_mcount_kmer){
            std::cout << get<0>(entry) << " " << -get<1>(entry) << " " << get<2>(entry) << std::endl;
        }
        return os;
    }
private:
    kp_analysis_t();
    pos_kmer_count_t pkc;
    kmer_pos_count_t kpc;
    const uint64_t mK;
};

class alternatives_t
{
public:
    struct elem_t{
        elem_t(size_t a, size_t b, size_t c, size_t d
              , basevector::const_iterator e, basevector::const_iterator f
              ):qfront(a),qback(b),rfront(c),rback(d),qseq(e+a,e+b+1),rseq(f+c,f+d+1),qcount(0),rcount(0),mcount(0){};
        std::map<basevector,size_t> other_count;
        size_t qfront;
        size_t qback;
        size_t rfront;
        size_t rback;
        basevector qseq;
        basevector rseq;
        size_t qcount;
        size_t rcount;
        size_t mcount;
        friend ostream& operator<<(ostream&os,const elem_t&in){
            os << "r(" << /*in.rfront << "," << in.rback << "," <<*/ in.rseq.ToString() << "," << in.rcount << ") q("
                       << /*in.qfront << "," << in.qback << "," <<*/ in.qseq.ToString() << "," << in.qcount << ") ("
                       << in.mcount << ")" ;
            for(const auto& entry: in.other_count){
                os<< " (" << entry.first.ToString()<<","<<entry.second <<")";
            }
            return os;
        }
    };
    alternatives_t(){};
    alternatives_t(const basevector&qseq,const basevector&rseq,const alignment&a,const size_t flank=1){load(qseq,rseq,a,flank);}

    void load(const basevector& qseq, const basevector& rseq, const alignment&a,const size_t flank=1){
        clear();

        int pos1,pos2,errors;
        avector<int> gaps, lengths;
        a.Unpack(pos1,pos2,errors,gaps,lengths);

        /*
        std::cout << pos1 << " " << pos2 << " " << errors << std::endl;
        std::cout << lengths.length << " " << gaps.length << std::endl;
        for(size_t ii=0;ii<lengths.length;++ii){
            std::cout << " " << lengths.x[ii];
        }
        std::cout << std::endl;
        for(size_t ii=0;ii<gaps.length;++ii){
            std::cout << " " << gaps.x[ii];
        }
        std::cout << std::endl;
        */

        ForceAssert(lengths.length==gaps.length);
        const size_t nL=lengths.length;

        if(nL<2) return;

        ForceAssert(gaps.x[0]==0);

        int qnext=pos1;
        int rnext=pos2;
        for(size_t ll=1;ll<nL;++ll){
            for( int ss=0;ss<lengths.x[ll-1];++ss){
                if( qseq[qnext+ss] != rseq[rnext+ss]){
                    mElements.emplace_back(max(qnext+ss-int(flank),0),qnext+ss+int(flank)
                                          ,max(rnext+ss-int(flank),0),rnext+ss+int(flank),qseq.begin(),rseq.begin());
                }
            }
            qnext += lengths.x[ll-1];
            rnext += lengths.x[ll-1];

            size_t gap = abs(gaps.x[ll]);
            if( gaps.x[ll] > 0){
                size_t rfront=rnext;
                size_t rback=rnext+gap-1;
                size_t qfront=qnext;
                size_t qback=qnext-1;
                for(;rfront>0 && rseq[rfront-1]==rseq[rfront];--rfront,--qfront){}
                for(;rback+1<rseq.size() && rseq[rback+1]==rseq[rback];++rback,++qback){}

                qfront-=flank;
                rfront-=flank;

                qback+=flank;
                rback+=flank;

                mElements.emplace_back(qfront,qback,rfront,rback,qseq.begin(),rseq.begin());

                rnext+=gap;
            }
            else{
                ForceAssert(false);
                size_t rfront=rnext;
                size_t rback=rnext-1;
                size_t qfront=qnext;
                size_t qback=qnext+gap-1;

                for(;qfront>0 && qseq[qfront-1]==qseq[qfront];--qfront,--rfront){}
                for(;qback+1<qseq.size() && qseq[qback+1]==qseq[qback];++qback,++rback){}

                mElements.emplace_back(qfront,qback,rfront,rback,qseq.begin(),rseq.begin());

                qnext+=gap;
            }
        }
        for( int ss=0;ss<lengths.x[nL-1];++ss){
            if( qseq[qnext+ss] != rseq[rnext+ss]){
                mElements.emplace_back(max(qnext+ss-int(flank),0),qnext+ss+int(flank)
                                      ,max(rnext+ss-int(flank),0),rnext+ss+int(flank),qseq.begin(),rseq.begin());
            }
        }
    };
    void pileup(const vecbasevector& reads,const vec<alignment>& read_alignments){
        ForceAssert(reads.size()==read_alignments.size());
        for(size_t rr=0;rr<reads.size();++rr){
            pileup(reads[rr],read_alignments[rr]);
        }
    }
    void pileup(const vecbasevector& reads,const vecqualvector& quals, const vec<vec<size_t>>& pos_map
               ,const typename qualvector::value_type threshold,unsigned int flank,int64_t pos_dev);
    void clear(){mElements.clear();}
    friend ostream& operator<<(ostream&os,const alternatives_t&in){
        for(size_t ee=0;ee<in.mElements.size();++ee){
            os << ee << " " << in.mElements[ee] << std::endl;
        }
        return os;
    }
private:
    void getQuerySequence(basevector&out,size_t cfront, size_t cback, const basevector&q,const alignment&a,int verbosity){
        out.clear();
        if(cback<cfront){return;}
        out.reserve(cback-cfront+1);

        int pos1,pos2,errors;
        avector<int> gaps, lengths;
        a.Unpack(pos1,pos2,errors,gaps,lengths);
        ForceAssert(lengths.length==gaps.length);
        const size_t nL=lengths.length;
        if(nL<1) return;
        ForceAssert(gaps.x[0]==0);

        size_t qfront=pos1;
        size_t qback=std::numeric_limits<size_t>::max();
        size_t rfront=pos2;
        size_t rback=std::numeric_limits<size_t>::max();

        size_t lstart=0;
        size_t cpos=cfront;
        for(;lstart<nL;++lstart){
            if( gaps.x[lstart] > 0){
                if( cpos >= rfront && cpos < rfront+gaps.x[lstart] ){
                    out.push_back(q[qfront-1]);
                    cpos = rfront+gaps.x[lstart];
                }
                rfront += gaps.x[lstart];
            }
            else{
                qfront += -gaps.x[lstart];
            }
            qback = qfront + lengths.x[lstart] - 1;
            rback = rfront + lengths.x[lstart] - 1;
//            if( rfront <= cfront && cfront <= rback){
            if( rfront <= cpos && cpos <= rback){
                break;
            }
            if(verbosity){
                std::cout << "pileup seek: " << lstart << "("
                                             << qfront << " " << qback << ")("
                                             << cfront << " " << cback << ")( "
                                             << rfront << " " << rback << ") "
                                             << basevector(q.begin()+qfront,q.begin()+qback+1).ToString()
                                             << std::endl;
            }
            qfront = qback+1;
            rfront = rback+1;
        }
        if(verbosity){
            std::cout << "pileup seek: " << lstart << "("
                                         << qfront << " " << qback << ")("
                                         << cfront << " " << cback << ")( "
                                         << rfront << " " << rback << ") "
                                         << basevector(q.begin()+qfront,q.begin()+qback+1).ToString()
                                         << std::endl;
        }
        if(lstart==nL) return;
        ForceAssert( rfront <= cpos && cpos <=rback);


        for(/*size_t cpos=cfront*/; cpos<=cback; ++cpos){
            if(cpos > rback){
                qfront = qback+1;
                rfront = rback+1;
if(verbosity){
    std::cout << "pileup skip: " << cpos << "("
                                 << qfront << " " << qback << ")("
                                 << cfront << " " << cback << ")( "
                                 << rfront << " " << rback << ")" << std::endl;
}
                ++lstart;
                if( gaps.x[lstart] > 0){
if(verbosity){
    std::cout << "pileup plus: " << cpos << "("
                                 << qfront << " " << qback << ")("
                                 << cfront << " " << cback << ")( "
                                 << rfront << " " << rback << ")" << std::endl;
}
                    rfront += gaps.x[lstart];
                    if(cpos+gaps.x[lstart] > cback){
                        out.push_back( q[qfront]  );
                    }
                    cpos+=gaps.x[lstart];
                }
                else{
if(verbosity){
    std::cout << "pileup minus: " << cpos << "("
                                 << qfront << " " << qback << ")("
                                 << cfront << " " << cback << ")( "
                                 << rfront << " " << rback << ")" << std::endl;
}
                    if(cpos<cback) out.append(q.begin()+qfront,q.begin()+qfront-gaps.x[lstart]);
                    qfront += -gaps.x[lstart];
                }
                qback = qfront + lengths.x[lstart] - 1;
                rback = rfront + lengths.x[lstart] - 1;
if(verbosity){
    std::cout << "pileup skip ends: " << cpos << "("
                                 << qfront << " " << qback << ")("
                                 << cfront << " " << cback << ")( "
                                 << rfront << " " << rback << ")" << std::endl;
}
            }
            if(cpos>cback) break;

            int qpos = cpos+qfront-rfront;
            out.push_back(q[qpos]);
        }
    }
    void pileup(const basevector&read,const alignment& a){
        basevector r_chunk;
        size_t nPass=0;
        for(auto& elem: mElements){
            bool bMicroVerbose=false;
            ++nPass;
            getQuerySequence(r_chunk,elem.rfront,elem.rback,read,a,bMicroVerbose&&false);
            bool bR = r_chunk.size()==elem.rseq.size() && std::equal(r_chunk.begin()+1,r_chunk.end()-1,elem.rseq.begin()+1);
            bool bQ = r_chunk.size()==elem.qseq.size() && std::equal(r_chunk.begin()+1,r_chunk.end()-1,elem.qseq.begin()+1);
            if(bR==bQ){
                if(bMicroVerbose){
                    std::cout << "pileup tie: "
                              << elem.rseq.ToString() << " "
                              << elem.qseq.ToString() << " "
                              << r_chunk.ToString() << std::endl;
                }
                ++elem.mcount;
                elem.other_count[r_chunk]+=1;
            }
            else if(bR){
                if(bMicroVerbose){
                    std::cout << "pileup ref: "
                              << elem.rseq.ToString() << " "
                              << elem.qseq.ToString() << " "
                              << r_chunk.ToString() << std::endl;
                }
                ++elem.rcount;
            }
            else if(bQ){
                if(bMicroVerbose){
                    std::cout << "pileup con: "
                              << elem.rseq.ToString() << " "
                              << elem.qseq.ToString() << " "
                              << r_chunk.ToString() << std::endl;
                }
                ++elem.qcount;
            }
        }
    }
    vec<elem_t> mElements;
};



}
