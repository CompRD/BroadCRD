///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// DisplayNhood: display a neighborhood of a SupportedHyperBasevector, seeded on 
// given edges (to be rendered green) and extended to given depth in the graph.  
// Boundary vertices are labeled red.

// MakeDepend: dependency QueryLookupTable
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "MainTools.h"
#include "ParseSet.h"
#include "lookup/LookAlign.h"
#include "paths/long/CleanEfasta.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/SupportedHyperBasevector.h"
#include <fstream>
#include "ParseSet.h"
#include "paths/reporting/ReftigUtils.h"

namespace
{
    class bubble_t;
    int64_t isHomoSNP(const bubble_t&bubble,const SupportedHyperBasevector&shbv);
    
    class bubble_t{
    public:
        bubble_t():iInEdge(-1),iOutEdge(-1),vEdges(),vWeights(),vReadCounts(){};
        bubble_t(int i,int o,const vec<int>&e,vec<fix64_6>w):iInEdge(i),iOutEdge(o),vEdges(e),vWeights(w),vReadCounts(e.size(),0){};
        
        void logRead(size_t edge){++vReadCounts[edge];};
        
        int In()const{return iInEdge;};
        int Out()const{return iOutEdge;};
        const vec<int>& Edges()const{return vEdges;};
        const vec<fix64_6>& Weights()const{return vWeights;};
        const vec<size_t>& ReadCounts()const{return vReadCounts;};
        friend ostream& operator<< ( ostream &os, const bubble_t &bubble ){
            os<< bubble.iInEdge  << " -> (";
            for(size_t ee=0;ee<bubble.vEdges.size();++ee){
                os << " " << bubble.vEdges[ee] <<":"<<bubble.vReadCounts[ee]
                                               <<":"<<std::scientific<<std::setprecision(5)<<bubble.vWeights[ee].ToDouble();
            }
            os<< " ) -> " << bubble.iOutEdge;
            return os;
        } ;
    private:
        int iInEdge;
        int iOutEdge;
        vec<int> vEdges;
        vec<fix64_6> vWeights;
        vec<size_t> vReadCounts;
    };
    
    void CenteredSequence(String&seq,const int64_t flank, const size_t left, const size_t mid, const size_t right,const SupportedHyperBasevector&shbv){
        const int64_t pad = shbv.K()-1;
        seq.clear();
        const auto& seqL=shbv.EdgeObject(left);
        const int64_t lenL=seqL.size()-pad;
        
        const auto& seqM=shbv.EdgeObject(mid);
        
        const auto& seqR=shbv.EdgeObject(right);
        const int64_t lenR=seqR.size()-pad;
        
        if(flank>0){
            if(lenL <= flank){
                std::transform(seqL.begin()    ,seqL.end()-pad,std::back_inserter(seq),BaseToCharMapper());
            }
            else{
                std::transform(seqL.end()-pad-flank    ,seqL.end()-pad,std::back_inserter(seq),BaseToCharMapper());
            }
        }
        std::transform(seqM.begin()    ,seqM.end()    ,std::back_inserter(seq),BaseToCharMapper());
        if(flank>0){
            if(lenR <= flank){
                std::transform(seqR.begin()+pad,seqR.end(),std::back_inserter(seq),BaseToCharMapper());
            }
            else{
                std::transform(seqR.begin()+pad,seqR.begin()+pad+flank,std::back_inserter(seq),BaseToCharMapper());
            }
        }
        
    }

    void DetectBubbles(vec<bubble_t>& bubbles, const SupportedHyperBasevector& shbv,const String& READS_FASTB, const String& ALIGN_HEAD
                      ,const int64_t max_bubble_length
                      ,const int64_t maximum_target_range
                      ){
        const auto nEdges=shbv.EdgeObjectCount();
        bubbles.clear();
        bubbles.reserve(nEdges/2);
        
        //for read alignments
        const bool bAlignReads=READS_FASTB.size()>0;
        const String FASTA = ALIGN_HEAD+".fasta";
        const String FASTB = ALIGN_HEAD+".fasta.fastb";
        const String LOOKUP = ALIGN_HEAD+".fasta.lookup";
        const String ALIGN = ALIGN_HEAD+".qlt";
        vec< std::tuple<size_t,size_t,int,int,int> > target_bubble_branch;
        target_bubble_branch.reserve(nEdges);
        
        
        vec< vec<std::pair<int,int>> > edge_path(nEdges);
        for(int pp=0;pp<shbv.NPaths();++pp){
            const auto& path=shbv.Path(pp);
            for(int ee=0;ee<path.isize();++ee){
                edge_path[path[ee]].push_back(make_pair(pp,ee));
            }
        }
        
        
        size_t out_count=0;
        std::ofstream fasta(FASTA,ios::out);
        std::ofstream fasta_rev[2];
        fasta_rev[0].open(FASTA+"0",ios::out);
        fasta_rev[1].open(FASTA+"1",ios::out);
        
        for( auto vid = 0 ; vid<shbv.N(); ++vid){
            if( shbv.ToSize(vid) == 1 && shbv.To(vid)[0] != vid){
                int iChoice=-1;
                vec<int> edges;
                for( auto ee=0;ee<shbv.FromSize(vid) ; ++ee){
                    auto dest = shbv.From(vid)[ee];
                    if(dest==vid) continue;
                    if(iChoice == -1) iChoice=dest;
                    if(iChoice == dest){
                        edges.push_back(shbv.EdgeObjectIndexByIndexFrom(vid,ee));
                    }
                    else{
                        iChoice=-2;
                    }
                }
                if(iChoice>=0 && shbv.FromSize(iChoice)==1&& shbv.From(iChoice)[0] != iChoice && shbv.ToSize(iChoice)==edges.isize() && edges.size()>1){
                    auto inEdge=shbv.EdgeObjectIndexByIndexTo(vid,0);
                    auto outEdge=shbv.EdgeObjectIndexByIndexFrom(iChoice,0);
                    
                    vec<fix64_6> weights(edges.size(),0);
                    unsigned int shortest=std::numeric_limits<unsigned>::max();
                    unsigned int longest=std::numeric_limits<unsigned>::min();
                    for(size_t ee=0;ee<edges.size();++ee){
                        shortest=std::min(shortest,shbv.EdgeObject(edges[ee]).size());
                        longest=std::max(longest,shbv.EdgeObject(edges[ee]).size());
                        for(const auto& path_coor:edge_path[edges[ee]]){
                            const auto& path = shbv.Path(path_coor.first);
                            ForceAssert( path[path_coor.second] == edges[ee]);
                            if(   (path_coor.second == 0             || path[path_coor.second-1] == inEdge || path[path_coor.second-1] == outEdge) 
                               && (path_coor.second == path.isize()-1 || path[path_coor.second+1] == inEdge || path[path_coor.second+1] == outEdge) 
                               && ( path[path_coor.second+1] != path[path_coor.second-1] || inEdge==outEdge)
                              ){
                                weights[ee]+=shbv.Weight(path_coor.first);
                            }
                        }
                    }
                    for(size_t ee=0;ee<edges.size();++ee){
                        auto pad = shbv.K()-1;
                        bvec a(shbv.EdgeObject(edges[ee]).begin(),shbv.EdgeObject(edges[ee]).begin()+pad);
                        bvec b(shbv.EdgeObject(edges[ee]).end()-pad,shbv.EdgeObject(edges[ee]).end());
                        
                        bvec c(shbv.EdgeObject(inEdge).end()-pad,shbv.EdgeObject(inEdge).end());
                        bvec d(shbv.EdgeObject(outEdge).begin(),shbv.EdgeObject(outEdge).begin()+pad);
                        ForceAssert(a==c );
                        ForceAssert(b==d );
                        for(size_t ff=ee+1;ff<edges.size();++ff){
                            ForceAssert(edges[ee]!=edges[ff]);
                            ForceAssert(shbv.EdgeObject(edges[ee])!=shbv.EdgeObject(edges[ff]));
                        }
                    }
//                    bool bTake= true || isHomoSNP(bubble_t(inEdge,outEdge,edges,weights),shbv)>=10;
//                    if(bTake){
                        if(longest <= max_bubble_length){
                            const int64_t flank = (maximum_target_range - shortest)/2;
                            if(longest+2*(flank>=0?flank:0)<=maximum_target_range){
                                auto base_count=out_count;
                                for(size_t ee=0;ee<edges.size();++ee){
                                    target_bubble_branch.push_back(std::make_tuple(bubbles.size(),ee,inEdge,edges[ee],outEdge));
                                    String seq;
                                    CenteredSequence(seq,flank,inEdge,edges[ee],outEdge,shbv);
                                    fasta << ">" << out_count++ << "_" << bubbles.size() << "_" << ee 
                                          << "_" << weights[ee]
                                          << "_" << inEdge << "_" << edges[ee] << "_" << outEdge 
                                          << std::endl;
                                    fasta << seq << std::endl;
                                    
                                    if(ee<2){
                                        fasta_rev[ee%2] << ">" << base_count << "_" << bubbles.size() << "_" << ee 
                                              << "_" << weights[ee]
                                              << "_" << inEdge << "_" << edges[ee] << "_" << outEdge 
                                              << std::endl;
                                        fasta_rev[ee%2] << seq << std::endl;
                                    }
                                }
                            }
                        }
                        bubbles.emplace_back(inEdge ,outEdge ,edges,weights);
//                    }
                }
            }
        }
        fasta.close();
        fasta_rev[0].close();
        fasta_rev[1].close();
        
        if(bAlignReads){
            const int K = 26;
            vec<look_align> hits;
            if(! IsRegularFile( ALIGN )) {
                SystemSucceedQuiet( "Fasta2Fastb IN="+FASTA+ " OUT=" +FASTB );
                SystemSucceedQuiet( "MakeLookupTable SOURCE="+FASTA+ " OUT_HEAD=" +FASTA +" LO=True ");
                Mkdir777(ALIGN+".dir");
                GetAlignsFast( K, READS_FASTB, LOOKUP, ALIGN, hits, True, ALIGN+".dir" );
                Rmdir(ALIGN+".dir");
            }
            else{
                GetAlignsFast( K, READS_FASTB, LOOKUP, ALIGN, hits, True, ALIGN+".dir" );
            }
            sort(hits.begin(),hits.end()
                ,[](const look_align&L,const look_align&R) { if(L.query_id<R.query_id) return true;
                                                             if(L.query_id>R.query_id) return false;
                                                             if(L.target_id<R.target_id) return true;
                                                             if(L.target_id>R.target_id) return false;
                                                             if(L.Errors()<R.Errors()) return true;
                                                             if(L.Errors()>R.Errors()) return false;
                                                             if(L.pos2()<R.pos2()) return true;
                                                             if(L.pos2()>R.pos2()) return false;
                                                             return false;
                                                           });
            for(size_t next_align=0;next_align<hits.size();){
                size_t start_align=next_align;
                for(++next_align
                   ;next_align<hits.size() && hits[start_align].query_id==hits[next_align].query_id
                                           && std::get<0>(target_bubble_branch[hits[start_align].target_id])==
                                              std::get<0>(target_bubble_branch[hits[next_align].target_id])
                   ;++next_align){ }
                
                vec<size_t> good_a;
                for(size_t aa=start_align;aa<next_align;++aa){
                    const auto& align=hits[aa];
                    if( align.pos1()==0 && align.Pos1()==int(align.query_length) && align.Errors()<10){
                        good_a.push_back(aa);
                        auto tmp = std::get<0>(target_bubble_branch[align.target_id]);
                        auto tmpe = std::get<1>(target_bubble_branch[align.target_id]);
//                        std::cout << "hit: " << align.query_id << ">" << align.target_id << " "<< align.Errors() << "=" << align.Mutations()<<"+"<<align.Indels()<<" "
//                                  << tmp << "," << tmpe <<" "<< bubbles[tmp] << std::endl;
                    }
                }
//                std::cout<<"section"<<std::endl;
                if(good_a.size()==0) continue;
                if( good_a.size()==1 || hits[good_a[0]].Errors()<hits[good_a[1]].Errors()){
                    const auto& tmp = target_bubble_branch[hits[good_a[0]].target_id];
                    bubbles[std::get<0>(tmp)].logRead( std::get<1>(tmp));
                }
            }
        }
    }
    
    class report_t{
    private:
        vec<size_t> _vNEdges;
        vec<int64_t> _vLengths;
        vec<size_t> _vFlankLengths;
        vec<double> _vWeightsRatio;
        vec<double> _vCountsRatio;
        size_t _nBubbles;
        size_t _nZero;
        size_t _nZeroRead;
        size_t _nHomoSNP;
    public:
        void sort(){
            std::sort(_vNEdges.begin(),_vNEdges.end());
            std::sort(_vLengths.begin(),_vLengths.end());
            std::sort(_vFlankLengths.begin(),_vFlankLengths.end());
            std::sort(_vWeightsRatio.begin(),_vWeightsRatio.end());
            std::sort(_vCountsRatio.begin(),_vCountsRatio.end());
        };
        report_t():_nBubbles(0),_nZero(0),_nZeroRead(0),_nHomoSNP(0){};
        void Consider(const vec<bubble_t>& bubbles, const SupportedHyperBasevector& shbv){
            _nBubbles+=bubbles.size();
            
            for(const auto& bubble: bubbles){
                if(isHomoSNP(bubble,shbv)>=10){
                    ++_nHomoSNP;
                }
            }
            
            _vNEdges.reserve(_vNEdges.size()+bubbles.size());
            for(const auto& entry: bubbles){
                _vNEdges.push_back( entry.Edges().size());
            }
            
            _vLengths.reserve(_vLengths.size()+bubbles.size());
            for(size_t ii=0;ii<bubbles.size();++ii){
                int len=0;
                for(const auto& eid: bubbles[ii].Edges()){
                    len=max( shbv.EdgeObject(eid).isize() , len);
                }
                _vLengths.push_back( len -2*(shbv.K()-1));
            }
            
            _vFlankLengths.reserve(_vFlankLengths.size()+bubbles.size()*2);
            for(const auto& entry: bubbles){
                _vFlankLengths.push_back( shbv.EdgeObject(entry.In()).size());
                _vFlankLengths.push_back( shbv.EdgeObject(entry.Out()).size());
            }
            
            _vWeightsRatio.reserve(_vWeightsRatio.size()+bubbles.size());
            for(const auto& entry: bubbles){
                auto weights = entry.Weights();
                std::sort(weights.begin(),weights.end());
                double max=weights.back().ToDouble();
                double nextmax=weights[weights.size()-2].ToDouble();
                ForceAssert(weights.size()>1);
                ForceAssert(max>=nextmax);
                if(max==0){
                    ++_nZero;
                }
                else{
                    _vWeightsRatio.push_back(nextmax/max);
                }
            }
            _vCountsRatio.reserve(_vCountsRatio.size()+bubbles.size());
            for(const auto& entry: bubbles){
                auto weights = entry.ReadCounts();
                std::sort(weights.begin(),weights.end());
                double max=double(weights.back());
                double nextmax=double(weights[weights.size()-2]);
                ForceAssert(weights.size()>1);
                ForceAssert(max>=nextmax);
                if(max==0){
                    ++_nZeroRead;
                }
                else{
                    _vCountsRatio.push_back(nextmax/max);
                }
            }
        }
        friend ostream& operator<< ( ostream &os, const report_t &report ){
            os << std::endl;
            if(report._nBubbles==0){
                os << "no bubbles in record" << std::endl;
            }
            os << "                       (min 1st-qtile median 3rd-qtile max)" << std::endl;
            os << "# edges in bubble:      " << report._vNEdges.front() << " " 
                                             << report._vNEdges[report._vNEdges.size()/4]<< " "
                                             << report._vNEdges[report._vNEdges.size()/2]<< " "
                                             << report._vNEdges[report._vNEdges.size()/4*3]<< " "
                                             << report._vNEdges.back() << std::endl;
            
            os << "bubble len - 2*(K-1):   " << report._vLengths.front() << " "
                                             << report._vLengths[report._vLengths.size()/4]<< " "
                                             << report._vLengths[report._vLengths.size()/2]<< " "
                                             << report._vLengths[report._vLengths.size()/4*3]<< " "
                                             << report._vLengths.back() << std::endl;
            
            os << "flanking edge length:   " << report._vFlankLengths.front() << " "
                                             << report._vFlankLengths[report._vFlankLengths.size()/4]<< " "
                                             << report._vFlankLengths[report._vFlankLengths.size()/2]<< " "
                                             << report._vFlankLengths[report._vFlankLengths.size()/4*3]<< " "
                                             << report._vFlankLengths.back() << std::endl;
            
            os << "weight:(2nd max)/(max): " << report._vWeightsRatio.front() << " "
                                             << report._vWeightsRatio[report._vWeightsRatio.size()/4]<< " "
                                             << report._vWeightsRatio[report._vWeightsRatio.size()/2]<< " "
                                             << report._vWeightsRatio[report._vWeightsRatio.size()/4*3]<< " "
                                             << report._vWeightsRatio.back() << std::endl;
            
            os << "count:(2nd max)/(max):  " << report._vCountsRatio.front() << " "
                                             << report._vCountsRatio[report._vCountsRatio.size()/4]<< " "
                                             << report._vCountsRatio[report._vCountsRatio.size()/2]<< " "
                                             << report._vCountsRatio[report._vCountsRatio.size()/4*3]<< " "
                                             << report._vCountsRatio.back() << std::endl;
            
            os << "# bubbles:              " << report._nBubbles << std::endl;
            os << "# zero-weight bubbles:  " << report._nZero << std::endl;
            os << "# zero-read bubbles:    " << report._nZeroRead << std::endl;
            os << "# homopolymer SNP:      " << report._nHomoSNP << std::endl;
            
            os << std::endl;
            return os;
        };
    };
    int64_t isHomoSNP(const bubble_t&bubble,const SupportedHyperBasevector&shbv){
        int64_t len=-1;
        const auto& edges=bubble.Edges();
        for( const auto& edge: edges){
            const bvec& seq=shbv.EdgeObject(edge);
            if(len<0) len=seq.size();
            if(seq.size()!=len) return -1;
        }
        int64_t first_diff=-1;
        for(int64_t bb=0;bb<len;++bb){
            auto base=shbv.EdgeObject(edges[0])[bb];
            for(size_t ee=1;ee<edges.size();++ee){
                if(base != shbv.EdgeObject(edges[ee])[bb]){
                    if(first_diff<0){
                        first_diff=bb;
                    }
                    if(first_diff!=bb){
                        return -1;
                    }
                }
            }
        }
        ForceAssert(first_diff>=shbv.K()-1);
        
        int64_t homo_length=0;
        for( const auto& edge: edges){
            const bvec& seq=shbv.EdgeObject(edge);
            auto base=seq[first_diff];
            int64_t loc_length=1;
            for( int64_t bb = first_diff+1 ; bb<seq.size() && seq[bb]==base ; ++bb,++loc_length){}
            for( int64_t bb = first_diff-1 ; bb>=0 && seq[bb]==base ; --bb,++loc_length){}
            homo_length=max(homo_length,loc_length);
        }
        return homo_length;
    };
}

int main(int argc, char *argv[])
{
    RunTime( );
     
    BeginCommandArguments;
    CommandArgument_String_Doc(IN_HEADS, "a list of IN_HEAD, INHEAD.shbv will be parsed")
    CommandArgument_String_OrDefault_Doc(FASTBS, "", "a list of fastb files storing uncorrected reads for the assemblies in IN_HEAD")
    CommandArgument_Int_OrDefault_Doc(MAX_BUBBLE_LEN, std::numeric_limits<int>::max(), "maximum bubble len under consideration")
    CommandArgument_Int_OrDefault_Doc(MAX_TARGET_LEN, std::numeric_limits<int>::max(), "maximum target range for read alignments")
    EndCommandArguments;
    
    
    vec<String> heads;
    ParseStringSet( IN_HEADS, heads );
    
    vec<String> fastbs;
    if(FASTBS.size()>0){
        ParseStringSet( FASTBS, fastbs );
    }
    else{
        fastbs.assign(heads.size(),"");
    }
    
    report_t report;
    for(size_t hh=0;hh<heads.size();++hh){
        const String& head=heads[hh];
        SupportedHyperBasevector shbv;
        const String fSHBV=head+".shbv";
        std::cout << Date( ) << ": loading " << fSHBV << std::endl;
        if ( isReadable(fSHBV) ) BinaryReader::readFile( fSHBV, &shbv );
        else                     FatalErr(fSHBV);
        
        std::cout << Date( ) << ": detecting bubbles" << std::endl;
        vec<bubble_t> bubbles;
        DetectBubbles(bubbles,shbv,fastbs[hh],ToString(hh),MAX_BUBBLE_LEN,MAX_TARGET_LEN);
        
        std::cout << Date( ) << ": filling in report" << std::endl;
        for(const auto&bubble:bubbles){
            std::cout << bubble << std::endl;
        }
        report.Consider(bubbles,shbv);
    }
    report.sort();
    std::cout << report << std::endl;
    std::cout << Date( ) << ": done" << std::endl;
}

