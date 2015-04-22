///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * SummarizeBubblesGapToy.cc
 *
 *  Created on: Mar 21, 2014
 *      Author: blau
 */

#include <array>

#include "MainTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"

class bubble_logger{
public:
    typedef std::tuple<int,int,int> support_t;
    class bubble_data{
    public:
        bubble_data(const vec<int>&branch1_edges,const vec<int>&branch2_edges)
            :branch_edges({branch1_edges,branch2_edges})
            ,branch_supports(2){}
        void addSupport(size_t branch, support_t score){ branch_supports[branch].push_back(score);}
        vec<support_t>const& getSupport(size_t branch)const{ return branch_supports[branch];};
        vec<int>const& getEdges(size_t branch)const{ return branch_edges[branch];};
    private:
        const vec< vec<int> > branch_edges;
        vec< vec<support_t> > branch_supports;
    };
    bubble_logger(const HyperBasevector& hb, const vec<int>& inv)
        :edge_alt_(hb.EdgeObjectCount(),-1)
        ,edge_bubble_branch_(hb.EdgeObjectCount(),std::make_pair(-1,-1))
        ,bubble_data_()
    {
        for(int vv=0;vv<hb.N();++vv){
            if(    hb.ToSize(vv) == 1
                && hb.FromSize(vv) == 2
                && hb.From(vv).front() != vv
                && hb.From(vv).front() == hb.From(vv).back()
                && hb.FromSize( hb.From(vv).front() ) == 1
                && hb.From( hb.From(vv).front() ).front() != hb.From(vv).front()
              ){
                const int edge_0 = hb.EdgeObjectIndexByIndexFrom(vv,0);
                const int edge_0_rc = inv[edge_0];
                const int edge_1 = hb.EdgeObjectIndexByIndexFrom(vv,1);
                const int edge_1_rc = inv[edge_1];

                if( (edge_0_rc <0 && edge_1_rc >=0) || (edge_0_rc >=0 && edge_1_rc <0 ))
                    std::cout << "WARNING: inv-data is inconsistent rc("<<edge_0<<")=" << edge_0_rc << " rc("<<edge_1<<")="<<edge_1_rc << std::endl;

                if(edge_0_rc<0 || edge_1_rc<0) continue;


                if(edge_alt_[edge_0] < 0 || edge_alt_[edge_1]<0){
                    if( edge_alt_[edge_0]>=0 || edge_alt_[edge_1]>=0){
                        FatalErr("inconsistent bubble, probably due to inversion properties");
                    }

                    const int bubble_idx = bubble_data_.size();

                    edge_alt_[edge_0]=edge_1;
                    edge_bubble_branch_[edge_0] = std::make_pair(bubble_idx,0);
                    vec<int> branch0(1,edge_0);

                    edge_alt_[edge_1]=edge_0;
                    edge_bubble_branch_[edge_1] = std::make_pair(bubble_idx,1);
                    vec<int> branch1(1,edge_1);
                    if( edge_0_rc >=0 || edge_1_rc >=0){
                        if( edge_0_rc<0 || edge_1_rc<0){
                            FatalErr("inconsistent bubble, probably due to inversion properties");
                        }

                        if(edge_alt_[edge_0_rc] < 0 || edge_alt_[edge_1_rc]<0){
                            if( edge_alt_[edge_0_rc]>=0 || edge_alt_[edge_1_rc]>=0){
                                FatalErr("inconsistent bubble, probably due to inversion properties");
                            }
                            edge_alt_[edge_0_rc]=edge_1_rc;
                            edge_bubble_branch_[edge_0_rc] = std::make_pair(bubble_idx,0);
                            branch0.push_back(edge_0_rc);
                            edge_alt_[edge_1_rc]=edge_0_rc;
                            edge_bubble_branch_[edge_1_rc] = std::make_pair(bubble_idx,1);
                            branch1.push_back(edge_1_rc);
                        }
                    }
                    bubble_data_.emplace_back(branch0,branch1);
                }
                else{
                    if( edge_alt_[edge_0]!=edge_1 && edge_alt_[edge_1]!=edge_0){
                        FatalErr("inconsistent bubble, probably due to inversion properties");
                    }
                }
            }
        }
//        std::cout << "alt: " << std::endl;
        for(int ee=0;ee<edge_alt_.isize();++ee){
//            std::cout << ee << "->" << edge_alt_[ee] << " : "
//                      << edge_bubble_branch_[ee].first << " " << edge_bubble_branch_[ee].second << " "
//                      << std::endl;
            ForceAssert(edge_bubble_branch_[ee].first<0 || edge_bubble_branch_[ee].first == edge_bubble_branch_[edge_alt_[ee]].first);
            ForceAssert(
                  ( edge_bubble_branch_[ee].second == -1 && edge_bubble_branch_[edge_alt_[ee]].second == -1)
               || ( edge_bubble_branch_[ee].second != edge_bubble_branch_[edge_alt_[ee]].second )
                       );
            ForceAssert(ee!=edge_alt_[ee]);
        }
    }

    int alt(int edge)const{ return edge_alt_[edge];}
    bool inBubble(int edge)const{ return edge_bubble_branch_[edge].first>=0;}
    void addWeight(int edge,support_t score){
        const auto& b_b=edge_bubble_branch_[edge];
        bubble_data_[ b_b.first ].addSupport(b_b.second,score);
    }
    void report(std::ostream&os)const{
        os << "bubble_logger report: "<< std::endl;
        os << "support_t = (qsum_of_read_against_branch, diff_of_qsum_between_two_branch, n_diff_of_read_against_branch" << std::endl;
        os << "(edges0:edges1) (branch0_qsum:branch1_qsum) { per-read-(qsum,qdif)-list-0 }:{ per-read-(qsum,qdif)-list-1 }" << std::endl;
        for( const auto& data: bubble_data_){
            vec<support_t> support0 = data.getSupport(0);
            vec<support_t> support1 = data.getSupport(1);
            vec<int> edges0 = data.getEdges(0);
            vec<int> edges1 = data.getEdges(1);

            auto cmp_fcn = [](support_t const&l,support_t const&r)
                    { if(std::get<0>(l)==std::get<0>(r)) return std::get<1>(l) > std::get<1>(r);
                      else return std::get<0>(l)<std::get<0>(r);};
            std::sort(support0.begin(),support0.end(),cmp_fcn);
            std::sort(support1.begin(),support1.end(),cmp_fcn);

            auto add_fcn = [](int const&i, support_t const&s){return i+std::get<1>(s);};
            int q0 = std::accumulate(support0.begin(),support0.end(),0,add_fcn);
            int q1 = std::accumulate(support1.begin(),support1.end(),0,add_fcn);

            if( q1 < q0){
                std::swap(q0,q1);
                std::swap(support0,support1);
                std::swap(edges0,edges1);
            }

            os << "( ";
            for(const auto&entry: edges0){ os << entry << " "; }
            os << ": ";
            for(const auto&entry: edges1){ os << entry << " "; }
            os << ") (";

            os << q0 << ":" << q1 << ") { ";
            for(const auto&entry: support0){ os << "(" << std::get<0>(entry) << "," << std::get<1>(entry) << "," << std::get<2>(entry) << ") ";}
            os << "}:{ ";
            for(const auto&entry: support1){ os << "(" << std::get<0>(entry) << "," << std::get<1>(entry) << "," << std::get<2>(entry) << ") ";}
            os << "}" << std::endl;
        }
    }
private:
    vec<int> edge_alt_;
    vec< std::pair<int,int> > edge_bubble_branch_;
    vec< bubble_data > bubble_data_;
};

//do an gap-free alignment of the read against the graph according to rp
//returns the sum of read-quality-score at the position the sequence mismatches
std::tuple<int,int,int>
getQ(basevector const&read, qualvector const&qual, ReadPath const&rp, HyperBasevector const&hb){
//std::cout << "getQ: " << rp.getOffset() << " ";
//for(const auto&entry: rp){ std::cout  << entry << "(" << hb.EdgeObject(entry).size() << ") " ; }
//std::cout<<std::endl;
    int out=0;
    int bp=0;
    int shift=rp.getOffset();
    if( shift<0 ){
        bp = -shift;
        shift=0;
    }
//    std::cout << "read " << bp << "("<<read.size() << ")"<<std::endl;

    const int bp_start = bp;
    int nDiff=0;

    for(size_t ee=0;bp<read.isize() && ee<rp.size();++ee){
        basevector const& edge = hb.EdgeObject(rp[ee]);
        for(size_t ep=shift ; bp < read.isize() && ep < edge.size() ; ++bp , ++ep){
            if( read[bp] != edge[ep] && qual[bp]>3){
//std::cout << bp << " vs " << rp[ee] << ":" <<ep << " " << BaseToCharMapper()(read[bp]) << " " << BaseToCharMapper()(edge[ep]) << "+" << int(qual[bp]) << std::endl;;
                ++nDiff;
                out += qual[bp];
            }
            else{
//std::cout << bp << " vs " << rp[ee] << ":" <<ep << " " << BaseToCharMapper()(read[bp]) << " " << BaseToCharMapper()(edge[ep]) << " " << int(qual[bp]) << std::endl;;
            }
        }
        shift=hb.K()-1;
    }
    return std::make_tuple(out,nDiff,bp-bp_start);
}

// if one or more edges in rp is part of a bubble, perform gap-free alignment on the path as the alternate path
// and collect the result to bubble_logger
void log_read(bubble_logger&bl, basevector const&read, qualvector const&qual, ReadPath const&rp, HyperBasevector const&hb){
    if(rp.size()==0) return;
    for(size_t rr=0;rr<rp.size();++rr){
        const int edge = rp[rr];
        const int other_edge = bl.alt(edge);
        if( other_edge >= 0
               && hb.EdgeObject(edge).size() < read.size()*2
               && hb.EdgeObject(other_edge).size() < read.size()*2 ){
            ReadPath other_rp=rp;
            other_rp[rr]=other_edge;
            if( rr == 0){
                other_rp.setOffset(   hb.EdgeObject(other_edge).isize()
                                    - hb.EdgeObject(edge).size()
                                    + rp.getOffset()
                                  );
            }
            const auto cur = getQ(read,qual,rp,hb);
            const auto alt = getQ(read,qual,other_rp,hb);
            const int q_cur = std::get<0>(cur);
            const int q_alt = std::get<0>(alt);
            const int nd_cur = std::get<1>(cur);
            const int nd_alt = std::get<1>(alt);
            if( q_cur > q_alt){
                std::cout << "WARNING: real-path of this read is not the lowest error path, q_alt < q_cur " << q_alt << " " << q_cur << std::endl;
//std::cout<< "adding " <<q_cur << " " << q_alt << std::endl;
                bl.addWeight(other_edge,std::make_tuple(q_alt,q_cur-q_alt,nd_alt));
            }
            else{
//std::cout<< "adding " <<q_alt << " " << q_cur << std::endl;
                bl.addWeight(edge,std::make_tuple(q_cur,q_alt-q_cur,nd_cur));
            }
        }
    }
}

class analyzer_t{
public:
    explicit analyzer_t(const String&b)
        :base_dir_(b)
        ,reads_(b+"/data/frag_reads_orig.fastb")
        ,quals_(b+"/data/frag_reads_orig.qualb")
        { ForceAssert(reads_.size() == quals_.size());}
    void analyze(const String& sub_dir);
private:
    const String base_dir_;
    const vecbasevector reads_;
    const vecqualvector quals_;
};

bool ValidPath(ReadPath const& rp, HyperBasevector const& hb,vec<int> const& to_right){
    for(size_t rr=1;rr<rp.size();++rr){
        int prev_edge = rp[rr-1];
        int curr_edge = rp[rr];
        int vv = to_right[prev_edge];
        bool bFound=false;
        for(int jj=0; jj < hb.FromSize(vv) ; ++jj){
            if( hb.EdgeObjectIndexByIndexFrom(vv,jj) == curr_edge){
                bFound=true;
                break;
            }
        }
        if(!bFound) return false;
    }
    return true;
}



void analyzer_t::analyze(const String& sub_dir){
     const String assembly_dir = base_dir_ + "/" + sub_dir;
     std::cout << Date() << ": analyzing" << assembly_dir << std::endl;
     if( ! IsDirectory(assembly_dir)){ FatalErr(assembly_dir + " does not exist."); }


     if( ! IsRegularFile(assembly_dir+"/a.hbv")){ FatalErr(assembly_dir+"/a.hbv" + " does not exist."); }

     HyperBasevector hb(assembly_dir+"/a.hbv");
     vec<int> inv2;  BinaryReader::readFile(assembly_dir+"/a.inv",&inv2);
     ForceAssert( hb.EdgeObjectCount() == inv2.isize());

     bubble_logger logger( hb , inv2 );

     ReadPathVec paths2(assembly_dir+"/a.paths");
     ForceAssert( paths2.size() == reads_.size() );
     vec<int> to_right;
     hb.ToRight(to_right);

     for(size_t rr=0;rr<reads_.size();++rr){
         if( !ValidPath(paths2[rr],hb,to_right)){
             std::cout << "WARNING: invalid path: offset="<<paths2[rr].getOffset() << " edges=";
             for(const auto&entry:paths2[rr]){ std::cout << entry << " "; }
             std::cout<<std::endl;
         }
         log_read(logger, reads_[rr], quals_[rr], paths2[rr], hb);
     }

     logger.report(std::cout);
}

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault_Doc(INSTANCE, "1",
          "to allow multiple concurrent runs");
     CommandArgument_String_OrDefault_Doc(USER, "",
          "use this instead of Getenv(USER)");
     CommandArgument_String_OrDefault_Doc(DIR, "a.fin",
          "directory name of the assembly, under the conventional work_dir");
     EndCommandArguments;

     const String work_dir =   "/wga/scr4/"
                             + ( USER != "" ? USER : Getenv("USER") )
                             + "/GapToy/" + INSTANCE;
     if( ! IsDirectory(work_dir)){ FatalErr(work_dir + " does not exist."); }

     analyzer_t instance(work_dir);
     std::cout << Date() << ": analyzer set up for " << work_dir << std::endl;

     instance.analyze(DIR);

}
