///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// A 'connection assembly' consists of a collection of DNA sequences, together
// with perfect overlaps between some of them.  These are either between one
// sequence (or its reverse complement, if rc1), and another (or its reverse
// complement, if rc2).

#ifndef CONNECTION_ASSEMBLY
#define CONNECTION_ASSEMBLY

#include "Map.h"
#include "Set.h"
#include "Basevector.h"
#include "CoreTools.h"
#include "graph/Digraph.h"
#include "math/HoInterval.h"
#include "paths/HyperBasevector.h"

class connection {
     public:
     connection( );
     connection( const Bool rc1, const Bool rc2, const int overlap )
          : rc1(rc1), rc2(rc2), overlap(overlap) { }
     connection( const Bool rc1, const Bool rc2, const int overlap ,const int f)
          : rc1(rc1), rc2(rc2), overlap(overlap),flag(f) { }
     Bool rc1, rc2;
     int overlap;
     int flag;
};

inline std::ostream& operator<< (std::ostream& out, const connection& c) {
  return out << (c.rc1 ? "rc" : "fw") << "-(" << c.overlap  << ")-" 
	     << (c.rc2 ? "rc" : "fw"); 
}

/*
void back_trace_fermi_shrink(int src, int iRight, const vec<std::pair<int,int>>& ancestor){
    if(iRight<0){

        return;
    }
    std::pair<int,int> p = ancestor[2*src+int(iRight)];
    if(iRight){
        if(p.first>=0){
          std::cout<<"TRACE\t" << src << "\t<-\t" << p.first << " " << p.second << std::endl;
        }
    }
    else{
        if(p.first>=0){
          std::cout<<"TRACE\t" << p.first << " " << p.second << "\t->\t" << src << std::endl;
        }
    }

}
*/


class connection_assembly {
     public:
     struct node_t{
          node_t():bases(),status(),bConnected(false),bValid(true),bProcessed(false)
                  ,bRC(false),bAsRC(false),bAsFW(false)
                  ,variants(),variant_of(make_tuple(-1,-1,-1))
                  {std::cout<<"WARNING node_t default constructor used.";};
          node_t(const basevector& b):bases(b)
                                     ,status(b.size(),0)
                                     ,mapped_loc(b.size())
//                                     ,mapped_branch(b.size()+1)
                                     ,vertices(b.size()+1,-1)
                                     ,clones(b.size())
                                     ,clones_tgt(b.size())
                                     ,clones_src(b.size())
                                     ,bConnected(false),bValid(true),bProcessed(false)
                                     ,bRC(false),bAsRC(false),bAsFW(false)
                                     ,variants(),variant_of(make_tuple(-1,-1,-1))
                                     {};
          node_t(const node_t& in):bases(in.bases),incoming(in.incoming),outgoing(in.outgoing),zero_length_destination(in.zero_length_destination)
                                  ,mapped_loc(in.mapped_loc)
//                                  ,mapped_branch(in.mapped_branch)
                                  ,vertices(in.vertices)
                                  ,clones(in.clones),clones_tgt(in.clones_tgt),clones_src(in.clones_src)
                                  ,bConnected(in.bConnected),bValid(in.bValid),bProcessed(in.bProcessed),bRC(in.bRC),bAsRC(in.bAsRC),bAsFW(in.bAsFW)
                                  ,variants(in.variants),variant_of(in.variant_of)
                                  {};
          basevector bases;
          vec<short> status;;
          std::map<unsigned int, vec<std::pair<int,unsigned int>>> incoming;
          std::map<unsigned int, vec<std::pair<int,unsigned int>>> outgoing;
          std::map<unsigned int, std::set<std::pair<int,int> >> zero_length_destination;

          vec<std::pair<int,int>> mapped_loc;
//          vec<std::pair<int,int>> mapped_branch;
          vec<int> vertices;

          vec<StdSet<std::pair<int,int>>> clones;

          vec<StdSet<std::pair<int,int>>> clones_tgt; // clones of this location
          vec<StdSet<std::pair<int,int>>> clones_src; // the "original" location of this clone
          bool bConnected;
          bool bValid;
          bool bProcessed;

          bool bRC;
          bool bAsRC;
          bool bAsFW;

          vec<std::tuple<int,int,int> > variants;//vertex number, front, back
          std::tuple<int,int,int> variant_of;
     };
     vecbasevector bases;
//     digraphE<connection> G;
     digraphVE<node_t,connection> G;

  void PrintConnections(std::ostream& out, int v) {
    const vec<int> & from = G.From(v);
    const vec<int> & to = G.To(v);
    for (uint32_t i = 0; i < to.size(); ++i) {
      const connection & c = G.EdgeObjectByIndexTo(v, i);
      out << to[i] << "--" << c << "--" << v << endl;
    }
    for (uint32_t i = 0; i < from.size(); ++i) {
      const connection & c = G.EdgeObjectByIndexFrom(v, i);
      out << v << "--" << c << "--" << from[i] << endl;
    }
    
  }

};

class branch_log{
public:
//    branch_log(size_t n):nBases(n),incoming(n+1),outgoing(n+1){};
//    size_t nBases;
//    vec<std::list<std::tuple<int,int,connection>> > incoming;
//    vec<std::list<std::tuple<int,int,connection>> > outgoing;
//    std::set<int> appearances;

    std::map<int,int> appearances;
    std::set<int> stack_appearances;
//    std::map< std::pair<std::pair<int,int>,int> , size_t > edge_log;// make_pair(minmax(src,des),flag), number
    struct stack_element{
        explicit stack_element(int i):src(i){};
        int src;
//        vec<int> edge_indices_to_childs;
//        vec<std::pair<int,int>> edge_overlap_indices_to_self;
//        vec<std::pair<int,int>> edge_flags_indices_to_ancestors;

        vec<std::tuple<int,int,int>> edge_overlap_index_dir_to_self;
        vec<std::tuple<int,int,int>> edge_flags_index_dir_to_ancestors;
        vec<std::tuple<int,int,int>> edge_flags_index_dir_to_descendants;

    private:
        stack_element();
    };
//    std::stack<stack_element> trail;

    vec<stack_element> trail;

};

void mark_log(branch_log& log, int src, connection_assembly& CA, int offset, int match_front, int match_back);

std::pair<std::pair<int,int>,std::pair<int,int> > getOverlapBaseRange(int iLength1, int iLength2,int overlap,int dir);
//std::pair<std::pair<int,int>,std::pair<int,int> > getOverlapBranchRange(int iLength1, int iLength2,int overlap,int dir);

//std::pair<std::tuple<int,int,int>,std::tuple<int,int,int> > getOverlapBaseRange(bool bRC1,int iLength, bool bRC2, int iLength2,int overlap);
//std::pair<std::tuple<int,int,int>,std::tuple<int,int,int> > getOverlapBranchRange(bool bRC1,int iLength, bool bRC2, int iLength2,int overlap);
//void MergeBranch( int target_idx, int target_pos, int source_idx, int source_pos, connection_assembly& ca, branch_log& log);

void MergeLoc(std::pair<int,int>&tgt_loc,std::pair<int,int>&src_loc,connection_assembly& ca,branch_log& log, bool bStrict=true);

void DetermineTargetVertices(std::set<int>& vertices, int& nVertices, connection_assembly&, branch_log&, int src, int src_pos);

void CollectTargetBranches(vec<vec<std::pair<int,int>>>&branches, const connection_assembly& ca, int src, int src_pos);

std::pair<int,int>& SyncMappedLoc(connection_assembly&, std::pair<int,int>);
//std::pair<int,int>& SyncMappedBranch(connection_assembly&, std::pair<int,int>);

// MatchRef: Find edge sequences that match a reference sequence perfectly.
// Answer is vec of ( ref_id, interval_on_ref_id, seq ), where seq entries are
// -n-1 if rc of edge.

void MatchRef( const connection_assembly& A, const vecbasevector& G,
     vec< triple< int, ho_interval, vec<int> > >& matches );

// ToHbv: for fixed K, attempt to convert a connection assembly into a 
// HyperBasevector assembly, together with an involution.  This can fail.  
// Return True upon success.
//
// Not fully implemented, and probably will never work.

Bool ToHbv( connection_assembly& A, const int K, HyperBasevector& hb, 
     vec<int>& inv, const Bool verbose = False );

void SetOrientation(size_t src, connection_assembly& ca);
void PreprocessGraph(connection_assembly& ca);

void AccommodateRC(connection_assembly& ca);

class BaseVec2KmerPath: public basevector
{
public:
    BaseVec2KmerPath():basevector(),bLeftKmerized(false),bRightKmerized(false){};
    BaseVec2KmerPath(const basevector&in):basevector(in),bLeftKmerized(false),bRightKmerized(false){};
    BaseVec2KmerPath& operator=(const basevector&in){
        basevector::operator=(in);
        bLeftKmerized=false;
        bRightKmerized=false;
        return *this;

    };
    BaseVec2KmerPath& operator=(const BaseVec2KmerPath&in){
        basevector::operator=(in);
        bLeftKmerized=in.bLeftKmerized;
        bRightKmerized=in.bLeftKmerized;
        return *this;
    };
    bool bLeftKmerized;
    bool bRightKmerized;
};

int MakeSequenceGraph(digraphE<BaseVec2KmerPath>&sequence_graph, connection_assembly CA);
bool MakeSequenceGraphConnected(digraphE<BaseVec2KmerPath>& sequence_graph, connection_assembly& CA);
bool MakeSequenceGraphConnectedVertex(digraphE<BaseVec2KmerPath>& sequence_graph, int vertex_index,connection_assembly& CA);

void SimplifySequenceGraph(digraphE<BaseVec2KmerPath>&sequence_graph);

void CollectSimpleTandems(vec<int>& log, int src, const digraphE<BaseVec2KmerPath>& sequence_graph);

int overlapLeft(const basevector&left, const basevector& right );
int overlapRight(const basevector&left, const basevector& right);

std::pair<int,int> DetermineVertexK(int vv, digraphE<BaseVec2KmerPath>&sequence_graph);
void HyperKmerize(digraphE<BaseVec2KmerPath>&sequence_graph,size_t target_k);

#endif

