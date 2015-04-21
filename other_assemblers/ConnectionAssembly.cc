///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <deque>
#include "CoreTools.h"
#include "Equiv.h"
#include "graph/DigraphTemplate.h"
#include "math/HoInterval.h"
#include "other_assemblers/ConnectionAssembly.h"
#include "paths/BigMapTools.h"
#include "paths/long/CreateGenome.h"

template void digraphE<connection>::Initialize( vec<vec<int> > const&, 
     vec<vec<int> > const&, vec<connection> const&, vec<vec<int> > const&, 
     vec<vec<int> > const&, const Bool );

template const connection & digraphE<connection>::EdgeObject(int) const;
template const connection & digraphE<connection>::EdgeObjectByIndexTo( int v, int j ) const;
template const connection & digraphE<connection>::EdgeObjectByIndexFrom( int v, int j ) const;

template void digraphVE<connection_assembly::node_t,connection>::Initialize( vec<vec<int> > const&,
     vec<vec<int> > const&, vec<connection_assembly::node_t>const&,vec<connection> const&, vec<vec<int> > const&,
     vec<vec<int> > const& );

template void digraphE<BaseVec2KmerPath>::RemoveDuplicateEdges() ;


// Breath-first exploration to look for minimal orientation needed to represent overlap graph assembly
void SetOrientation(int root, connection_assembly& ca){
    deque< std::pair<int,bool> > todo; // each entry is node number and bRC
    connection_assembly::node_t& root_node = ca.G.VertMutable(root);

    if( root_node.bAsFW || root_node.bAsRC) return;

    if( ca.G.FromSize(root) > 0){
        const connection& conn = ca.G.EdgeObjectByIndexFrom( root, 0 );
        if ( conn.rc1 ){
            if( !root_node.bAsRC ){
                root_node.bAsRC=true;
                todo.push_back(std::make_pair(root,true));
            }
        }
        else{
            if( !root_node.bAsFW ){
                root_node.bAsFW=true;
                todo.push_back(std::make_pair(root,false));
            }
        }
    }
    else if( ca.G.ToSize(root) > 0){
        const connection& conn = ca.G.EdgeObjectByIndexTo( root, 0 );
        if ( conn.rc2 ){
            if( !root_node.bAsRC ){
                root_node.bAsRC=true;
                todo.push_back(std::make_pair(root,true));
            }
        }
        else{
            if( !root_node.bAsFW ){
                root_node.bAsFW=true;
                todo.push_back(std::make_pair(root,false));
            }
        }
    }
    else{
        if( !root_node.bAsFW ){
            root_node.bAsFW=true;
            todo.push_back(std::make_pair(root,false));
        }
    }

    while(!todo.empty()){
        int src = todo.front().first;
        bool srcRC = todo.front().second;
        todo.pop_front();

        const connection_assembly::node_t& src_node = ca.G.Vert(src);
        if( srcRC ){ ForceAssert( src_node.bAsRC ); }
        else{ ForceAssert( src_node.bAsFW ); }

        for(int conn_index=0;conn_index<ca.G.FromSize(src);++conn_index){
            const connection& conn = ca.G.EdgeObjectByIndexFrom( src, conn_index );
            if( conn.rc1 == srcRC){
                int des = ca.G.From(src)[conn_index];
                connection_assembly::node_t& des_node = ca.G.VertMutable(des);
                if( conn.rc2 ){
                    if( !des_node.bAsRC ){
                        des_node.bAsRC=true;
                        todo.push_back(std::make_pair(des,true));
                    }
                }
                else{
                    if( !des_node.bAsFW ){
                        des_node.bAsFW=true;
                        todo.push_back(std::make_pair(des,false));
                    }
                }
            }

        }
        for(int conn_index=0;conn_index<ca.G.ToSize(src);++conn_index){
            const connection& conn = ca.G.EdgeObjectByIndexTo( src, conn_index );
            if( conn.rc2 == srcRC){
                int des = ca.G.To(src)[conn_index];
                connection_assembly::node_t& des_node = ca.G.VertMutable(des);
                if( conn.rc1 ){
                    if( !des_node.bAsRC ){
                        des_node.bAsRC=true;
                        todo.push_back(std::make_pair(des,true));
                    }
                }
                else{
                    if( !des_node.bAsFW ){
                        des_node.bAsFW=true;
                        todo.push_back(std::make_pair(des,false));
                    }
                }

            }
        }
    }
}

void AccommodateRC(connection_assembly& ca){
    int nExtra=0;
    const int nOrg=ca.G.N();
    int nVertices=ca.G.N();
    //duplicate nodes as needed
    for( int src=0;src<nOrg;++src){

//std::cout << "checking " << src << std::endl;
        if( std::get<0>(ca.G.VertMutable(src).variant_of )>=0){
            continue;
        }
        if( ca.G.VertMutable(src).bAsRC && ca.G.VertMutable(src).bAsFW){
            const int new_idx = nVertices;

//std::cout << "duplicating " << src <<" to " << new_idx<< std::endl;
            ca.G.VertMutable(src).bAsFW=true;
            ca.G.VertMutable(src).bAsRC=false;
            ForceAssert ( ca.G.Vert(src).bAsFW != ca.G.Vert(src).bAsRC);

            const size_t nVariants = ca.G.VertMutable(src).variants.size();


            connection_assembly::node_t new_node = ca.G.Vert(src);

            for(size_t vv=0; vv<nVariants;++vv){
                ForceAssert( std::get<0>(ca.G.VertMutable(src).variants[vv]) < nOrg    );
                ForceAssert( std::get<0>(ca.G.VertMutable(src).variants[vv]) >= 0   );
                std::get<0>(new_node.variants[vv]) = new_idx+1+vv;
            }

            new_node.bAsFW=false;
            new_node.bAsRC=true;
            ca.G.AddVertex(new_node);
            ++nExtra;
            ++nVertices;

            for(size_t vv=0; vv<nVariants;++vv){
                ForceAssert( std::get<0>(ca.G.VertMutable(src).variants[vv]) < nOrg    );
                ForceAssert( std::get<0>(ca.G.VertMutable(src).variants[vv]) >= 0   );
                connection_assembly::node_t new_variant_node = ca.G.Vert(std::get<0>(ca.G.VertMutable(src).variants[vv]));
                std::get<0>(new_variant_node.variant_of) = new_idx;
                ca.G.AddVertex(new_variant_node);
                ++nVertices;
                ++nExtra;
            }



            ForceAssert ( ca.G.Vert(src).bAsFW != ca.G.Vert(src).bAsRC);
            ForceAssert ( ca.G.Vert(new_idx).bAsFW != ca.G.Vert(new_idx).bAsRC);
            ForceAssert ( ca.G.Vert(new_idx).bAsFW != ca.G.Vert(src).bAsFW);
            ForceAssert (new_idx < ca.G.N());
            for(int conn_index=0;conn_index<ca.G.FromSize(src);++conn_index){
                const connection& conn = ca.G.EdgeObjectByIndexFrom( src, conn_index );
                if(conn.rc1){
                    int des = ca.G.From(src)[conn_index];
                    ca.G.AddEdge(new_idx,des,conn);
                }
            }

            for(int conn_index=0;conn_index<ca.G.ToSize(src);++conn_index){
                const connection& conn = ca.G.EdgeObjectByIndexTo( src, conn_index );
                if(conn.rc2){
                    int des = ca.G.To(src)[conn_index];
                    ca.G.AddEdge(des,new_idx,conn);
                }
            }
        }
    }
    //delete edges
    for( int src=0;src<ca.G.N();++src){
        connection_assembly::node_t& src_node = ca.G.VertMutable(src);
        if( std::get<0>(src_node.variant_of )>=0){
            continue;
        }
        ForceAssert( src_node.bAsRC || src_node.bAsFW);
        ForceAssert( src_node.bAsRC != src_node.bAsFW);
        src_node.bRC = src_node.bAsRC;
        for(int conn_index=ca.G.FromSize(src)-1;conn_index>=0;--conn_index){
            connection& conn = ca.G.EdgeObjectByIndexFromMutable( src, conn_index );
            if(conn.rc1!=src_node.bRC){
                ca.G.DeleteEdgeFrom(src,conn_index);
            }
        }
        for(int conn_index=ca.G.ToSize(src)-1;conn_index>=0;--conn_index){
            connection& conn = ca.G.EdgeObjectByIndexToMutable( src, conn_index );
            if(conn.rc2!=src_node.bRC){
                ca.G.DeleteEdgeTo(src,conn_index);
            }
        }
    }
    //flip RC to FW
    for( int src=0;src<ca.G.N();++src){
        connection_assembly::node_t& src_node = ca.G.VertMutable(src);
        if( std::get<0>(src_node.variant_of )>=0){
            continue;
        }
        ForceAssert( src_node.bAsRC || src_node.bAsFW);
        ForceAssert( src_node.bAsRC != src_node.bAsFW);
        ForceAssert( src_node.bAsRC == src_node.bRC);
        for(int conn_index=0;conn_index<ca.G.FromSize(src);++conn_index){
            connection& conn = ca.G.EdgeObjectByIndexFromMutable( src, conn_index );
            ForceAssert(conn.rc1==src_node.bRC);
            conn.rc1=false;
        }
        for(int conn_index=0;conn_index<ca.G.ToSize(src);++conn_index){
            connection& conn = ca.G.EdgeObjectByIndexToMutable( src, conn_index );
            ForceAssert(conn.rc2==src_node.bRC);
            conn.rc2=false;
        }
        if(src_node.bRC){
std::cout << src << " RC" << std::endl;
            src_node.bases.ReverseComplement();
            for(size_t vv=0; vv<src_node.variants.size();++vv){
                connection_assembly::node_t& variant_node = ca.G.VertMutable(std::get<0>(src_node.variants[vv]));

                ForceAssert( std::get<0>(variant_node.variant_of) == src);
                ForceAssert( std::get<1>(variant_node.variant_of) == std::get<1>(src_node.variants[vv]));
                ForceAssert( std::get<2>(variant_node.variant_of) == std::get<2>(src_node.variants[vv]));

                int new_front = src_node.bases.isize()-1-std::get<2>(src_node.variants[vv]);
                int new_back = src_node.bases.isize()-1-std::get<1>(src_node.variants[vv]);

                std::get<1>(src_node.variants[vv]) = new_front;
                std::get<2>(src_node.variants[vv]) = new_back;
                variant_node.bases.ReverseComplement();

                std::get<1>(variant_node.variant_of) = new_front;
                std::get<2>(variant_node.variant_of) = new_back;
            }
        }
        else{
std::cout << src << " FW" << std::endl;
        }
//std::cout << src_node.bases.ToString()<<std::endl;
        src_node.bRC=false;
        for(size_t vv=0; vv<src_node.variants.size();++vv){
            connection_assembly::node_t& variant_node = ca.G.VertMutable(std::get<0>(src_node.variants[vv]));

            ForceAssert( std::get<0>(variant_node.variant_of) == src);
            ForceAssert( std::get<1>(variant_node.variant_of) == std::get<1>(src_node.variants[vv]));
            ForceAssert( std::get<2>(variant_node.variant_of) == std::get<2>(src_node.variants[vv]));
            ForceAssert( variant_node.bases[0] == src_node.bases[std::get<1>(src_node.variants[vv])]);
            ForceAssert( variant_node.bases[variant_node.bases.size()-1] == src_node.bases[std::get<2>(src_node.variants[vv])]);
            variant_node.bRC=false;
            variant_node.bAsFW=true;
        }
    }
    ForceAssert( nOrg+nExtra == ca.G.N());
std::cout << nExtra << "/" << nOrg << " extra nodes to accommodate RC" << std::endl;
}

void PreprocessGraph(connection_assembly& ca){
    for( int src = 0 ; src < ca.G.N(); ++src){
        connection_assembly::node_t& src_node = ca.G.VertMutable(src);
        std::cout << "Preprocess vertex " << src << " " << src_node.bases.size() << std::endl;
        ForceAssert(src_node.bValid);

        const int nBases = src_node.bases.size();

        src_node.mapped_loc.resize(nBases);
//        src_node.mapped_branch.resize(nBases+1);
        src_node.vertices.resize(nBases+1);
        src_node.clones.resize(nBases);
        src_node.clones_tgt.resize(nBases);
        src_node.clones_src.resize(nBases);

        for(int ii=0;ii<nBases;++ii){
            src_node.mapped_loc[ii]=std::make_pair(src,ii);
        }
        if( std::get<0>(src_node.variant_of) >=0){
            std::cout<<"preping variant" << std::endl;
            src_node.mapped_loc[0]=std::make_pair(std::get<0>(src_node.variant_of),std::get<1>(src_node.variant_of));
            src_node.mapped_loc[nBases-1]=std::make_pair(std::get<0>(src_node.variant_of),std::get<2>(src_node.variant_of));

        }
        for(int ii=0;ii<nBases+1;++ii){
//            src_node.mapped_branch[ii]=std::make_pair(src,ii);
            src_node.vertices[ii]=-1;
        }


        if( ca.G.From(src).size()>0){
            src_node.bConnected=true;
            int last_neighbor = ca.G.From(src)[0];
            for(int conn_index=1;conn_index<ca.G.FromSize(src);++conn_index){
                int neighbor = ca.G.From(src)[conn_index];
                ForceAssert( last_neighbor <= neighbor);
                if(last_neighbor == neighbor){
//                    ForceAssert(  ca.G.EdgeObjectByIndexFrom(src,conn_index-1).flag
//                                !=ca.G.EdgeObjectByIndexFrom(src,conn_index).flag);
                    ForceAssert(ca.G.EdgeObjectByIndexFrom(src,conn_index-1).flag<=2);
                    ForceAssert(ca.G.EdgeObjectByIndexFrom(src,conn_index).flag<=2);
                }

                last_neighbor=neighbor;
            }
        }
        if( ca.G.To(src).size()>0){
            src_node.bConnected=true;
            int last_neighbor = ca.G.To(src)[0];
            for(int conn_index=1;conn_index<ca.G.ToSize(src);++conn_index){
                int neighbor = ca.G.To(src)[conn_index];

                ForceAssert( last_neighbor <= neighbor);

                if(last_neighbor == neighbor){
//                    ForceAssert(ca.G.EdgeObjectByIndexTo(src,conn_index-1).flag!=ca.G.EdgeObjectByIndexTo(src,conn_index).flag);
                    ForceAssert(ca.G.EdgeObjectByIndexTo(src,conn_index-1).flag<=2);
                    ForceAssert(ca.G.EdgeObjectByIndexTo(src,conn_index).flag<=2);
                }
                last_neighbor=neighbor;
            }
        }
    }
}

bool MakeSequenceGraphConnected(digraphE<BaseVec2KmerPath>& sequence_graph, connection_assembly& ca){

    for(int vi=0;vi<ca.G.N() ; ++vi){
        const connection_assembly::node_t& node = ca.G.Vert(vi);
        bool bLoopBack=false;


        for(size_t conn_index=0;conn_index<ca.G.To(vi).size();++conn_index){
            if( vi == ca.G.To(vi)[conn_index]){
                bLoopBack=true;
                break;
            }
        }
        if( bLoopBack && node.bValid && !node.bProcessed){
            MakeSequenceGraphConnectedVertex(sequence_graph,vi,ca);
        }
    }

    for(int vi=0;vi<ca.G.N() ; ++vi){
        const connection_assembly::node_t& node = ca.G.Vert(vi);
        if( node.bValid && !node.bProcessed){
            MakeSequenceGraphConnectedVertex(sequence_graph,vi,ca);
        }
    }
    return false;
}
bool MakeSequenceGraphConnectedVertex(digraphE<BaseVec2KmerPath>& sequence_graph, int vertex_index,connection_assembly& ca){
std::cout << std::endl
          <<" starting with vertex " << vertex_index << std::endl;;
    bool bModified=false;
    connection_assembly::node_t& root_node = ca.G.VertMutable(vertex_index);


//    branch_log log(root_node.bases.size());
    branch_log log;

//std::cout << root_node.bases.ToString() << std::endl;

    mark_log(log,vertex_index,ca,0,0,root_node.bases.size()-1);

    std::set<size_t> order;
    for(auto& entry: log.appearances){
        ForceAssert( order.find(entry.second) == order.end());
        order.insert(entry.second);
    }
//    for(auto& entry: order){ std::cout<< entry << " "; } std::cout << std::endl;
    for( size_t oo = 0 ; oo < order.size() ; oo++){
        ForceAssert( order.find(oo) != order.end());
    }

std::cout<< "-----------------------" << std::endl;
    for(size_t ii=0;ii<log.trail.size();++ii){
//        branch_log::stack_element& loc_stack = log.trail[log.trail.size()-1-ii];
        branch_log::stack_element& loc_stack = log.trail[ii];

        int src = loc_stack.src;
        std::cout << "contracting " << src << std::endl;
        connection_assembly::node_t& src_node = ca.G.VertMutable(src);
        vec<std::pair<int,int>> &src_mapped_loc = src_node.mapped_loc;
//        vec<std::pair<int,int>> &src_mapped_branch = src_node.mapped_branch;

        if(loc_stack.edge_overlap_index_dir_to_self.size()>0 ){
            std::cout <<"N edges to self: "<< loc_stack.edge_overlap_index_dir_to_self.size() << std::endl;
            Sort(loc_stack.edge_overlap_index_dir_to_self);

            bool b1st=true;
            int lastoverlap=0;

            for( const auto& entry: loc_stack.edge_overlap_index_dir_to_self){
                int conn_index = std::get<1>(entry);
                int dir = std::get<2>(entry);
                const connection& conn = (dir)?(ca.G.EdgeObjectByIndexFrom(src,conn_index))
                                              :(ca.G.EdgeObjectByIndexTo  (src,conn_index));
                int des = (dir)?(ca.G.From(src)[conn_index])
                               :(ca.G.  To(src)[conn_index]);
                ForceAssert(des==src);

                if( !b1st){
                    if( conn.overlap != lastoverlap){
                        FatalErr("Vertex overlaping with FW of itself at multiple overlaps are not implemented, but it can be done");
                    }
                    ForceAssert( conn.overlap == lastoverlap);
                }
                lastoverlap = conn.overlap;

                ForceAssert(!conn.rc1);
                ForceAssert(!conn.rc2);

                connection_assembly::node_t& des_node = ca.G.VertMutable(des);
                ForceAssert(!des_node.bProcessed);

//                if( log.edge_log.find( std::make_pair(std::minmax(src,des),conn.flag) ) != log.edge_log.end()) continue;
//                size_t tmp = log.edge_log.size();
//                log.edge_log[ std::make_pair(std::minmax(src,des),conn.flag) ] = tmp;


                if( conn.overlap*2 <= src_node.bases.isize()){
                    auto base_range=getOverlapBaseRange(src_node.bases.isize(),des_node.bases.isize(),conn.overlap,dir);

                    for( int i1=base_range.first.first, e1=base_range.first.second, i2=base_range.second.first, e2=base_range.second.second
                       ; i1 != e1 && i2!=e2 ; ++i1,++i2){
                        std::pair<int,int>& base1= SyncMappedLoc(ca,std::make_pair(src,i1));
                        std::pair<int,int>& base2= SyncMappedLoc(ca,std::make_pair(des,i2));

                        if( log.appearances[base2.first] < log.appearances[base1.first]){
//                            ForceAssert( des_node.clones[i2].size()==0);
                            ForceAssert( des_node.clones_src[i2].size()==0);
                            ForceAssert( des_node.clones_tgt[i2].size()==0);
                            MergeLoc( base2, base1,ca,log );
                        }
                        else if( log.appearances[base1.first] < log.appearances[base2.first]){
//                            ForceAssert( des_node.clones[i2].size()==0);
                            ForceAssert( des_node.clones_src[i2].size()==0);
                            ForceAssert( des_node.clones_tgt[i2].size()==0);
                            MergeLoc( base1, base2,ca,log );
                        }
                        else if (base1.first!=src){
                            if( base1.second != base2.second){
    std::cout << "cloning tie-breaker invoked for " << base1.first << " " << base1.second <<" "
                                                    << base2.first << " " << base2.second <<" "
                                                    << i1 << " " << e1 << " "
                                                    << i2 << " " << e2 << " (during loop back)"
                                                    << std::endl;
                                connection_assembly::node_t& clone_node = ca.G.VertMutable(base1.first);
//                                clone_node.clones[base2.second].insert(base1);
//                                clone_node.clones[base1.second].insert(base2);
                                if( base2.second < base1.second){
                                    clone_node.clones_tgt[base2.second].insert(base1);
                                    clone_node.clones_src[base1.second].insert(base2);
                                }
                                else{
                                    clone_node.clones_tgt[base1.second].insert(base2);
                                    clone_node.clones_src[base2.second].insert(base1);
                                }
                            }

//                            ForceAssert( src_node.clones[i1].size()==0);
//                            ForceAssert( des_node.clones[i2].size()==0);
                            ForceAssert( src_node.clones_tgt[i1].size()==0);
                            ForceAssert( src_node.clones_src[i1].size()==0);
                            ForceAssert( des_node.clones_tgt[i2].size()==0);
                            ForceAssert( des_node.clones_src[i2].size()==0);
                            base2=std::make_pair(src,i1);
                        }
                        else{
                            if( base1.second != base2.second){
                                ForceAssert( base1==std::make_pair(src,i1) && base2==std::make_pair(des,i2));
//                                for( auto& entry: des_node.clones[i2] ){
//                                    ForceAssert( entry != base2);
//
//                                    auto& clone_list = ca.G.VertMutable(entry.first).clones[entry.second];
//                                    ForceAssert(clone_list.find(base2)!=clone_list.end());
//                                    clone_list.erase(base2);
//
//                                    src_node.clones[base1.second].insert(entry);
//                                    clone_list.insert(base1);
//                                }
//                                des_node.clones[base2.second].clear();
//                                ForceAssert( des_node.clones[i2].size()==0);
                                for( auto& entry: des_node.clones_src[i2] ){
                                    ForceAssert( entry != base2);

                                    auto& clone_list = ca.G.VertMutable(entry.first).clones_tgt[entry.second];
                                    ForceAssert(clone_list.find(base2)!=clone_list.end());
                                    clone_list.erase(base2);

                                    clone_list.insert(base1);
                                    src_node.clones_src[base1.second].insert(entry);
                                }
                                des_node.clones_src[base2.second].clear();
                                ForceAssert( des_node.clones_src[i2].size()==0);

                                for( auto& entry: des_node.clones_tgt[i2] ){
                                    ForceAssert( entry != base2);

                                    auto& clone_list = ca.G.VertMutable(entry.first).clones_src[entry.second];
                                    ForceAssert(clone_list.find(base2)!=clone_list.end());
                                    clone_list.erase(base2);

                                    clone_list.insert(base1);
                                    src_node.clones_tgt[base1.second].insert(entry);
                                }
                                des_node.clones_tgt[base2.second].clear();
                                ForceAssert( des_node.clones_tgt[i2].size()==0);

                                base2=std::make_pair(src,i1);
                            }
                        }
                    }
                }
                else{
                    int trimmed_length = src_node.bases.isize() - conn.overlap;
                    ForceAssert(trimmed_length>0);
                    for( int i2=trimmed_length ; i2 < src_node.bases.isize() ; ++i2){
                        int i1 = (i2-trimmed_length) % trimmed_length;
                        ForceAssert(i1>=0);

                        std::pair<int,int>& base1= SyncMappedLoc(ca,std::make_pair(src,i1));
                        std::pair<int,int>& base2= SyncMappedLoc(ca,std::make_pair(des,i2));
                        if( log.appearances[base2.first] < log.appearances[base1.first]){
//                            ForceAssert( des_node.clones[i2].size()==0);
                            ForceAssert( des_node.clones_src[i2].size()==0);
                            ForceAssert( des_node.clones_tgt[i2].size()==0);
                            MergeLoc( base2, base1,ca,log );
                        }
                        else if( log.appearances[base1.first] < log.appearances[base2.first]){
//                            ForceAssert( des_node.clones[i2].size()==0);
                            ForceAssert( des_node.clones_src[i2].size()==0);
                            ForceAssert( des_node.clones_tgt[i2].size()==0);
                            MergeLoc( base1, base2,ca,log );
                        }
                        else if (base1.first!=src){
                            if( base1.second != base2.second){
    std::cout << "cloning tie-breaker invoked for " << base1.first << " " << base1.second <<" "
                                                    << base2.first << " " << base2.second <<" "
                                                    << i1 << " " 
                                                    << i2 << "(during loop back)"
                                                    << std::endl;
                                connection_assembly::node_t& clone_node = ca.G.VertMutable(base1.first);
//                                clone_node.clones[base2.second].insert(base1);
//                                clone_node.clones[base1.second].insert(base2);
                                if( base2.second < base1.second){
                                    clone_node.clones_tgt[base2.second].insert(base1);
                                    clone_node.clones_src[base1.second].insert(base2);
                                }
                                else{
                                    clone_node.clones_tgt[base1.second].insert(base2);
                                    clone_node.clones_src[base2.second].insert(base1);
                                }
                            }

//                            ForceAssert( src_node.clones[i1].size()==0);
//                            ForceAssert( des_node.clones[i2].size()==0);
                            ForceAssert( src_node.clones_src[i1].size()==0);
                            ForceAssert( src_node.clones_tgt[i1].size()==0);
                            ForceAssert( des_node.clones_src[i2].size()==0);
                            ForceAssert( des_node.clones_tgt[i2].size()==0);
                            base2=std::make_pair(src,i1);
                        }
                        else{
                            if( base1.second != base2.second){
                                ForceAssert( base1==std::make_pair(src,i1) && base2==std::make_pair(des,i2));
//                                for( auto& entry: des_node.clones[i2] ){
//                                    ForceAssert( entry != base2);
//
//                                    auto& clone_list = ca.G.VertMutable(entry.first).clones[entry.second];
//                                    ForceAssert(clone_list.find(base2)!=clone_list.end());
//                                    clone_list.erase(base2);
//                                    src_node.clones[base1.second].insert(entry);
//
//                                    clone_list.insert(base1);
//                                }
//                                des_node.clones[base2.second].clear();
//                                ForceAssert( des_node.clones[i2].size()==0);
                                for( auto& entry: des_node.clones_src[i2] ){
                                    ForceAssert( entry != base2);

                                    auto& clone_list = ca.G.VertMutable(entry.first).clones_tgt[entry.second];
                                    ForceAssert(clone_list.find(base2)!=clone_list.end());
                                    clone_list.erase(base2);
                                    clone_list.insert(base1);
                                    src_node.clones_tgt[base1.second].insert(entry);

                                }
                                des_node.clones_src[base2.second].clear();
                                ForceAssert( des_node.clones_src[i2].size()==0);

                                for( auto& entry: des_node.clones_tgt[i2] ){
                                    ForceAssert( entry != base2);

                                    auto& clone_list = ca.G.VertMutable(entry.first).clones_src[entry.second];
                                    ForceAssert(clone_list.find(base2)!=clone_list.end());
                                    clone_list.erase(base2);
                                    clone_list.insert(base1);
                                    src_node.clones_src[base1.second].insert(entry);

                                }
                                des_node.clones_tgt[base2.second].clear();
                                ForceAssert( des_node.clones_tgt[i2].size()==0);

                                base2=std::make_pair(src,i1);
                            }
                        }
                    }
                }
/*
                auto branch_range=getOverlapBranchRange(src_node.bases.isize(),des_node.bases.isize(),conn.overlap,dir);
                for( int i1=branch_range.first.first, e1=branch_range.first.second, i2=branch_range.second.first, e2=branch_range.second.second
                   ; i1 != e1 && i2!=e2 ; ++i1,++i2){
                    std::pair<int,int>& branch1= SyncMappedBranch(ca,std::make_pair(src,i1));
                    std::pair<int,int>& branch2= SyncMappedBranch(ca,std::make_pair(des,i2));
                    if (log.appearances[branch1.first]<=log.appearances[branch2.first]){
//                        MergeBranch( src, i1, des, i2, ca,log);
                        branch2=branch1;
                    }
                    else{
//                        MergeBranch( des, i2, src, i1, ca,log);
                        branch1=branch2;
                    }

                }
*/
                b1st=false;
            }
        }//if(loc_stack.edge_overlap_index_dir_to_self.size()>0 )


        if(loc_stack.edge_flags_index_dir_to_descendants.size()>0 ){
            std::cout <<"N descendants: "<< loc_stack.edge_flags_index_dir_to_descendants.size() << std::endl;
            Sort(loc_stack.edge_flags_index_dir_to_descendants);

            int last_priority=-1;
            int last_flag = -1;


            for( const auto& entry: loc_stack.edge_flags_index_dir_to_descendants){
                int des_priority = std::get<0>(entry);
                int conn_index = std::get<1>(entry);
                int dir = std::get<2>(entry);
                const connection& conn = (dir)?(ca.G.EdgeObjectByIndexFrom(src,conn_index))
                                              :(ca.G.EdgeObjectByIndexTo  (src,conn_index));
                ForceAssert(!conn.rc1);
                ForceAssert(!conn.rc2);
                int des = (dir)?(ca.G.From(src)[conn_index])
                               :(ca.G.  To(src)[conn_index]);

                connection_assembly::node_t& des_node = ca.G.VertMutable(des);
std::cout<<"Ancestor: " << des << std::endl;
                ForceAssert(!des_node.bProcessed);
                ForceAssert( des_priority == log.appearances[des]);

//                if( log.edge_log.find( std::make_pair(std::minmax(src,des),conn.flag) ) != log.edge_log.end()) continue;
//                size_t tmp = log.edge_log.size();
//                log.edge_log[ std::make_pair(std::minmax(src,des),conn.flag) ] = tmp;

                std::cout << "des flag priority " << des<< " " << conn.flag << " " << des_priority << std::endl;
                ForceAssert( last_priority <= des_priority );
                ForceAssert( std::get<0>(entry) > log.appearances[src] );
                last_priority=std::get<0>(entry);
                last_flag=conn.flag;
//std::cout << "connection: " << int(conn.rc1) << " " << int(conn.rc2) << conn.overlap << std::endl;
                auto base_range=getOverlapBaseRange(src_node.bases.isize(),des_node.bases.isize(),conn.overlap,dir);
                int nConflict=0;
                ForceAssert( base_range.first.second-base_range.first.first==base_range.second.second-base_range.second.first) ;
std::cout << "Loop: "
          << base_range.first.first << " " << base_range.first.second <<" "
          << base_range.second.first << " " << base_range.second.second <<" "
          << std::endl;
                for( int i1=base_range.first.first, e1=base_range.first.second, i2=base_range.second.first, e2=base_range.second.second
                   ; i1 != e1 && i2!=e2 ; ++i1,++i2){

                    ForceAssert(src_node.bases[i1]==des_node.bases[i2]);

                    std::pair<int,int>& mapped_base1= SyncMappedLoc(ca,std::make_pair(src,i1));
                    std::pair<int,int>& mapped_base2= SyncMappedLoc(ca,std::make_pair(des,i2));
    ForceAssert(  ca.G.Vert(mapped_base1.first).bases[mapped_base1.second]
                ==ca.G.Vert(mapped_base2.first).bases[mapped_base2.second]
               ) ;

                    // if the real bases are on the same sequence
                    if (log.appearances[mapped_base1.first]==log.appearances[mapped_base2.first]){
                        ForceAssert(mapped_base1.first==mapped_base2.first);
                        if( mapped_base1.second != mapped_base2.second){
std::cout << "cloning tie-breaker invoked for " << mapped_base1.first << " " << mapped_base1.second <<" "
                                                << mapped_base2.first << " " << mapped_base2.second <<" "
                                                << i1 << " " << e1 << " "
                                                << i2 << " " << e2 << " "
                                                << std::endl;
                            connection_assembly::node_t& clone_node = ca.G.VertMutable(mapped_base1.first);

//                            clone_node.clones[mapped_base2.second].insert(mapped_base1);
//                            clone_node.clones[mapped_base1.second].insert(mapped_base2);

                            if( mapped_base2.second < mapped_base1.second){
                                clone_node.clones_tgt[mapped_base2.second].insert(mapped_base1);
                                clone_node.clones_src[mapped_base1.second].insert(mapped_base2);
                            }
                            else{
                                clone_node.clones_tgt[mapped_base1.second].insert(mapped_base2);
                                clone_node.clones_src[mapped_base2.second].insert(mapped_base1);
                            }

                                nConflict++;
                        }
                    }
                    else if (log.appearances[mapped_base1.first]<log.appearances[mapped_base2.first]){
                        const connection_assembly::node_t& real_src_node = ca.G.Vert(mapped_base2.first);
//                        if( real_src_node.clones[mapped_base1.second].size() >0 )
                        if( real_src_node.clones_src[mapped_base2.second].size() >0  || real_src_node.clones_tgt[mapped_base2.second].size() >0 ){
std::cout << "part of clone propagated 1<-2 " << mapped_base1.first << " " << log.appearances[mapped_base1.first] << " " << mapped_base1.second
          << " to " << mapped_base1.first << " " << log.appearances[mapped_base1.first] << " " << mapped_base1.second << std::endl;
                            ForceAssert( log.appearances[mapped_base1.first] < log.appearances[mapped_base2.first] );
//                            ForceAssert( ca.G.Vert(mapped_base2.first).clones[mapped_base2.second].size() == 0);
                        }
                        MergeLoc(mapped_base1,mapped_base2,ca,log);
//                        mapped_base2=mapped_base1;
//std::cout << " modifying des with " << mapped_base1.first << " " << mapped_base1.second <<" " << std::endl;
                    }
                    else{
                        const connection_assembly::node_t& real_src_node = ca.G.Vert(mapped_base1.first);
//                        if( real_src_node.clones[mapped_base2.second].size() >0 )
                        if( real_src_node.clones_src[mapped_base1.second].size() >0  || real_src_node.clones_tgt[mapped_base1.second].size() >0 ){
std::cout << "part of clone propagated 2<-1 " << mapped_base1.first<< " " << log.appearances[mapped_base1.first] << " " << mapped_base1.second
          << " to " << mapped_base2.first << " " << log.appearances[mapped_base2.first] << " " << mapped_base2.second << std::endl;
                            ForceAssert( log.appearances[mapped_base2.first] < log.appearances[mapped_base1.first] );
//                            ForceAssert( ca.G.Vert(mapped_base1.first).clones[mapped_base1.second].size() == 0);
                        }
                        MergeLoc(mapped_base2,mapped_base1,ca,log);
//                        mapped_base1=mapped_base2;
//std::cout << " modifying src with  " << mapped_base2.first << " " << mapped_base2.second <<" " << std::endl;
                    }

                }
                std::cout << "conflicts: " << nConflict << "/" << base_range.first.second-base_range.first.first << std::endl;
/*
                auto branch_range=getOverlapBranchRange(src_node.bases.isize(),des_node.bases.isize(),conn.overlap,dir);
                for( int i1=branch_range.first.first, e1=branch_range.first.second, i2=branch_range.second.first, e2=branch_range.second.second
                   ; i1 != e1 && i2!=e2 ; ++i1,++i2){
                    std::pair<int,int>& mapped_branch1= SyncMappedBranch(ca,std::make_pair(src,i1));
                    std::pair<int,int>& mapped_branch2= SyncMappedBranch(ca,std::make_pair(des,i2));
                    if (log.appearances[mapped_branch1.first]==log.appearances[mapped_branch2.first]){
                        if( mapped_branch1.second != mapped_branch1.second){
                            FatalErr("Sketchy branch match not implemented");
                        }
//                        MergeBranch( src, i1, des, i2, ca,log);
                    }
                    else if (log.appearances[mapped_branch1.first]<log.appearances[mapped_branch2.first]){
//                        MergeBranch( src, i1, des, i2, ca,log);
                        mapped_branch2=mapped_branch1;
                    }
                    else{
//                        MergeBranch( des, i2, src, i1, ca,log);
                        mapped_branch1=mapped_branch2;
                    }

                }
*/
            }
std::cout<<"end of stack"<<std::endl;
        }
std::cout<<"done with " << src <<std::endl;
        src_node.bProcessed=true;
    }

    int nVertices=sequence_graph.N();
    int nEdges=0;
std::cout<< "-----------------------" << std::endl;

    for(size_t ii=0;ii<log.trail.size();++ii){
        int src = log.trail[log.trail.size()-1-ii].src;
std::cout << "tracing sequence " << src << std::endl;
        connection_assembly::node_t& src_node = ca.G.VertMutable(src);

        int front = -1, back=-1;
/*
if( src == 0 & false){
        for(int bb=0 ; bb < src_node.bases.isize(); ++bb){
            if(src_node.mapped_loc[bb]==std::make_pair(src,bb)){
              std::cout << 1;
            }
            else{
              std::cout << 0;
            }
        }
        std::cout << std::endl;
        std::cout << std::endl;

        for(int bb=0 ; bb <= src_node.bases.isize(); ++bb){
            if(src_node.mapped_branch[bb].first==std::make_pair(src,bb).first){
              std::cout << 1;
            }
            else{
              std::cout << 0;
            }
        }
        std::cout << std::endl;
        std::cout << std::endl;

        for(int bb=0 ; bb <= src_node.bases.isize(); ++bb){
            if(src_node.mapped_branch[bb].second==std::make_pair(src,bb).second){
              std::cout << 'o';
            }
            else{
              std::cout << src_node.mapped_branch[bb].second << " ";
            }
        }
        std::cout << std::endl;
        std::cout << std::endl;
        for(int bb=0 ; bb <= src_node.bases.isize(); ++bb){
            if(src_node.vertices[bb]==-1){
              std::cout << '-';
            }
            else{
              std::cout << '|';
            }
        }
        std::cout << std::endl;
}
*/
        std::set<int> src_vertices;
        std::set<int> des_vertices;
        for(int bb=0 ; bb < src_node.bases.isize(); ++bb){
            if(src_node.mapped_loc[bb]==std::make_pair(src,bb)){ //if location is unmapped
                if(front<0){ //if there is no unmapped loc before, edge starts here
                    front = bb;

                    src_vertices.clear();
                    DetermineTargetVertices( src_vertices, nVertices,ca,log,src,bb);
                    ForceAssert( src_vertices.find(src_node.vertices[bb]) != src_vertices.end());
                }

                if(  src_node.vertices[bb+1]>=0){ // if there's a vertex after loc bb, edge ends here
                    back=bb;

                    // double check vertices
                    des_vertices.clear();
                    DetermineTargetVertices( des_vertices, nVertices,ca,log,src,bb+1);

                    ForceAssert( des_vertices.find(src_node.vertices[bb+1]) != des_vertices.end());

                    front = -1;
                }
            }
            else{//if location is mapped
                if(front>=0){ // if front is set, the edge ends before the loc bb
                  // edge from front to bb-1
                    back=bb-1;

                    des_vertices.clear();
                    DetermineTargetVertices( des_vertices, nVertices,ca,log,src,bb);
                    ForceAssert( des_vertices.find(src_node.vertices[bb]) != des_vertices.end());

                    front = -1;
                }
            }
        }
        // this is for the last branch point, after the last loc
        if( front >= 0){
            int bb=src_node.bases.isize();
            back = bb-1;

            des_vertices.clear();
            DetermineTargetVertices( des_vertices, nVertices,ca,log,src,bb);
            ForceAssert( des_vertices.find(src_node.vertices[bb]) != des_vertices.end());

            front = -1;
        }
    }
    sequence_graph.AddVertices( nVertices-sequence_graph.N());
    std::cout << "determining zero-length local branches: " << std::endl;
    for(size_t ii=0;ii<log.trail.size();++ii){
        int src = log.trail[log.trail.size()-1-ii].src;
        connection_assembly::node_t& src_node = ca.G.VertMutable(src);
        int front=-1;
        vec< std::pair<int,int> > pattern,last_pattern;

        int nMismatch=0;
        int nCloned=0;
        int nRegion=0;
        set<int> forward_clones;
        int forward_left=0, forward_right=0;
        for(int bb=0 ; bb < src_node.bases.isize(); ++bb){
            if(src_node.mapped_loc[bb]==std::make_pair(src,bb)){ //if location is unmapped
//                if( src_node.clones[bb].size()!=0){
                if( src_node.clones_tgt[bb].size()!=0){
std::cout << bb << " " << src_node.bases.ToString()[bb] << std::endl;
                    ++nCloned;
                    pattern.clear();

                    if( forward_clones.size()==0){
                        bool bOfInterest=false;
                        int tmpL=0;
                        int tmpR=0;
//                        for(auto entry: src_node.clones[bb]){
                        for(auto entry: src_node.clones_tgt[bb]){
                            if( entry.first == src ){
                                if(entry.second>bb){
                                    bOfInterest=true;
                                }
                            }
                            else{
                                const connection_assembly::node_t& destined_node = ca.G.Vert(entry.first);
                                if( destined_node.vertices[entry.second] >= 0){
                                    tmpL++;
                                }
                                if( destined_node.vertices[entry.second+1] >= 0){
                                    tmpR++;
                                }
                            }
                        }
                        if(bOfInterest){
                            forward_left+=tmpL;
                            forward_right+=tmpL;
                            if( src_node.vertices[bb] >= 0&& forward_clones.find(bb)==forward_clones.end()){
                                forward_left++;
                            }
                            if( src_node.vertices[bb+1] >= 0&& forward_clones.find(bb)==forward_clones.end()){
                                forward_right++;
                            }
                            forward_clones.insert(bb);
                        }
                    }
                    if(forward_clones.find(bb) != forward_clones.end()){
//                        for(auto entry: src_node.clones[bb]){
                        for(auto entry: src_node.clones_tgt[bb]){
                            if( entry.first==src){
                                if ( entry.second-bb > 0){
                                    if( src_node.vertices[entry.second] >= 0 && forward_clones.find(entry.second)==forward_clones.end()){
                                        forward_left++;
                                    }
                                    if( src_node.vertices[entry.second+1] >= 0 && forward_clones.find(entry.second)==forward_clones.end()){
                                        forward_right++;
                                    }
                                    forward_clones.insert(entry.second);
                                }
                            }
                            else{
                                const connection_assembly::node_t& destined_node = ca.G.Vert(entry.first);
                                if( destined_node.vertices[entry.second] >= 0){
                                    forward_left++;
                                }
                                if( destined_node.vertices[entry.second+1] >= 0){
                                    forward_right++;
                                }
                            }
                        }
                    }


                    int nValidClones=0;

//                    for(auto entry: src_node.clones[bb]){
                    for(auto entry: src_node.clones_tgt[bb]){
                        if( entry.first != src || entry.second != bb){
                            nValidClones++;
                            pattern.push_back( make_pair(entry.first,entry.second-bb) );

                        }
                    }

                    if(nValidClones==0){
                        front=-1;
                        continue;
                    }

                    Sort(pattern);
                    ForceAssert( pattern.size() > 0);

                    if( front >= 0){
                        bool bMatches=true;
                        if( pattern.size() != last_pattern.size()) bMatches=false;
                        for(size_t ii=0;ii<pattern.size()&&bMatches;++ii){
                            if( pattern[ii] != last_pattern[ii] ){
                                bMatches=false;
                            }
                        }
                        if( !bMatches ){
for(auto entry: pattern){ std::cout << "(" << entry.first << "," << entry.second<<")"; } std::cout<<endl;
                            nMismatch++;
                        }
                        //compare
                    }
                    else{
for(auto entry: pattern){ std::cout << "(" << entry.first << "," << entry.second<<")"; } std::cout<<endl;
                        front=bb;
                        ++nRegion;
                    }

                    last_pattern=pattern;
                }
                else{
                    /*
                    if(forward_clones.size()>0){
                        std::cout<<"clones ";
                        for(auto entry:forward_clones){std::cout<< entry <<",";};
                        std::cout<< "(" << forward_left << " " << forward_right<<")"<<std::endl;
                    }
                    forward_clones.clear();
                    forward_left=0;
                    forward_right=0;
                    */
                    front = -1;
                }
            }
            else{
//                ForceAssert( src_node.clones[bb].size()==0);
                ForceAssert( src_node.clones_tgt[bb].size()==0);
                ForceAssert( src_node.clones_src[bb].size()==0);
                /*
                if(forward_clones.size()>0){
                    std::cout<<"clones ";
                    for(auto entry:forward_clones){std::cout<< entry <<",";};
                    std::cout<< "(" << forward_left << " " << forward_right<<")"<<std::endl;
                }
                forward_clones.clear();
                forward_left=0;
                forward_right=0;
                */


                front =-1;
            }
        }
        if(nCloned>0){
          std::cout<<"sequence " << src << " has " << nMismatch << "/" << nCloned << " mismatches in " << nRegion << " regions" << std::endl;
        }
        if(forward_clones.size()>0){
            std::cout<<"clones ";
            for(auto entry:forward_clones){std::cout<< entry <<",";};
            std::cout<< "(" << forward_left << " " << forward_right << ")"<<std::endl;
        }
        set<int> reverse_clones;
        nMismatch=0;
        nCloned=0;
        nRegion=0;
        int reverse_left=0,reverse_right=0;
        front=-1;
//        for(int bb=0 ; bb < src_node.bases.isize(); ++bb)
        for(int bb=src_node.bases.isize()-1 ; bb >=0 ; --bb){
            if(src_node.mapped_loc[bb]==std::make_pair(src,bb)){ //if location is unmapped
//                if( src_node.clones[bb].size()!=0){
                if( src_node.clones_src[bb].size()!=0){
std::cout << bb << " " << src_node.bases.ToString()[bb] << std::endl;
                    ++nCloned;
                    pattern.clear();

                    if( reverse_clones.size()==0){
                        bool bOfInterest=false;
                        int tmpL=0;
                        int tmpR=0;
//                        for(auto entry: src_node.clones[bb]){
                        for(auto entry: src_node.clones_src[bb]){
                            if( entry.first == src ){
                                if(entry.second <bb){
                                    bOfInterest=true;
                                }
                            }
                            else{
                                const connection_assembly::node_t& destined_node = ca.G.Vert(entry.first);
                                if( destined_node.vertices[entry.second] >= 0){
                                    tmpL++;
                                }
                                if( destined_node.vertices[entry.second+1] >= 0){
                                    tmpR++;
                                }
                            }
                        }
                        if(bOfInterest){
                            reverse_left+=tmpL;
                            reverse_right+=tmpL;
                            if( src_node.vertices[bb] >= 0&& reverse_clones.find(bb)==reverse_clones.end()){
                                reverse_left++;
                            }
                            if( src_node.vertices[bb+1] >= 0&& reverse_clones.find(bb)==reverse_clones.end()){
                                reverse_right++;
                            }
                            reverse_clones.insert(bb);
                        }
                    }
                    if(reverse_clones.find(bb) != reverse_clones.end()){
//                        for(auto entry: src_node.clones[bb]){
                        for(auto entry: src_node.clones_src[bb]){
                            if( entry.first==src){
                                if ( entry.second-bb < 0){
                                    if( src_node.vertices[entry.second] >= 0&& reverse_clones.find(entry.second)==reverse_clones.end()){
                                        reverse_left++;
                                    }
                                    if( src_node.vertices[entry.second+1] >= 0&& reverse_clones.find(entry.second)==reverse_clones.end()){
                                        reverse_right++;
                                    }
                                    reverse_clones.insert(entry.second);
                                }
                            }
                            else{
                                const connection_assembly::node_t& destined_node = ca.G.Vert(entry.first);
                                if( destined_node.vertices[entry.second] >= 0){
                                    reverse_left++;
                                }
                                if( destined_node.vertices[entry.second+1] >= 0){
                                    reverse_right++;
                                }
                            }
                        }
                    }


                    int nValidClones=0;

//                    for(auto entry: src_node.clones[bb]){
                    for(auto entry: src_node.clones_src[bb]){
                        if( entry.first != src || entry.second != bb){
                            nValidClones++;
                            pattern.push_back( make_pair(entry.first,entry.second-bb) );

                        }
                    }

                    if(nValidClones==0){
                        front=-1;
                        continue;
                    }

                    Sort(pattern);
                    ForceAssert( pattern.size() > 0);

                    if( front >= 0){
                        bool bMatches=true;
                        if( pattern.size() != last_pattern.size()) bMatches=false;
                        for(size_t ii=0;ii<pattern.size()&&bMatches;++ii){
                            if( pattern[ii] != last_pattern[ii] ){
                                bMatches=false;
                            }
                        }
                        if( !bMatches ){
for(auto entry: pattern){ std::cout << "(" << entry.first << "," << entry.second<<")"; } std::cout<<endl;
                            nMismatch++;
                        }
                        //compare
                    }
                    else{
for(auto entry: pattern){ std::cout << "(" << entry.first << "," << entry.second<<")"; } std::cout<<endl;
                        front=bb;
                        ++nRegion;
                    }

                    last_pattern=pattern;
                }
                else{
                    /*
                    if(reverse_clones.size()>0){
                        std::cout<<"clones ";
                        for(auto entry:reverse_clones){std::cout<< entry <<",";};
                        std::cout<< "(" << reverse_left << " " << reverse_right<<")"<<std::endl;
                    }
                    reverse_clones.clear();
                    reverse_left=0;
                    reverse_right=0;
                    */
                    front = -1;
                }
            }
            else{
//                ForceAssert( src_node.clones[bb].size()==0);
                ForceAssert( src_node.clones_src[bb].size()==0);
                ForceAssert( src_node.clones_tgt[bb].size()==0);
                /*
                if(reverse_clones.size()>0){
                    std::cout<<"clones ";
                    for(auto entry:reverse_clones){std::cout<< entry <<",";};
                    std::cout<< "(" << reverse_left << " " << reverse_right<<")"<<std::endl;
                }
                reverse_clones.clear();
                reverse_left=0;
                reverse_right=0;
                */


                front =-1;
            }
        }
        if(nCloned>0){
          std::cout<<"sequence " << src << " has " << nMismatch << "/" << nCloned << " mismatches in " << nRegion << " regions" << std::endl;
        }
        if(reverse_clones.size()>0){
            std::cout<<"clones ";
            for(auto entry:reverse_clones){std::cout<< entry <<",";};
            std::cout<< "(" << reverse_left << " " << reverse_right << ")"<<std::endl;
        }

//        ForceAssert( reverse_left==0 || reverse_right==0);
//        ForceAssert( forward_left==0 || forward_right==0);

        bool bFoward = forward_clones.size() > reverse_clones.size();// std::max(forward_left,forward_right) >= std::max(reverse_left,reverse_right);
        if(bFoward){
            int vertex_shift = (forward_right > forward_left)?(1):(0);
            for(int bb=0;bb<src_node.bases.isize();++bb){
                if( forward_clones.find(bb) != forward_clones.end()){
                    if( src_node.vertices[bb+vertex_shift]<0){
                        src_node.vertices[bb+vertex_shift]=nVertices;
                        sequence_graph.AddVertices(1);
                        ++nVertices;
                    }
//                    for(auto entry: src_node.clones[bb]){
                    for(auto entry: src_node.clones_tgt[bb]){
                        if(entry.first==src && entry.second>bb){
                            ForceAssert(forward_clones.find(entry.second)!=forward_clones.end());
                            if( src_node.vertices[entry.second+vertex_shift]<0){
                                src_node.vertices[entry.second+vertex_shift]=nVertices;
                                sequence_graph.AddVertices(1);
                                ++nVertices;
                            }
                            sequence_graph.AddEdge(src_node.vertices[bb+vertex_shift], src_node.vertices[entry.second+vertex_shift],basevector());
std::cout << "Skip Link " << src_node.vertices[bb+vertex_shift] << "->" << src_node.vertices[entry.second+vertex_shift] << std::endl;
                            ++nEdges;
                        }
                    }

                }

            }
        }
        else{
            int vertex_shift = (reverse_right > reverse_left)?(1):(0);
            for(int bb=src_node.bases.isize()-1;bb>=0;--bb){
                if( reverse_clones.find(bb) != reverse_clones.end()){
                    if( src_node.vertices[bb+vertex_shift]<0){
                        src_node.vertices[bb+vertex_shift]=nVertices;
                        sequence_graph.AddVertices(1);
                        ++nVertices;
                    }
//                    for(auto entry: src_node.clones[bb]){
                    for(auto entry: src_node.clones_src[bb]){
                        if(entry.first==src && entry.second<bb){
                            ForceAssert(reverse_clones.find(entry.second)!=reverse_clones.end());
                            if( src_node.vertices[entry.second+vertex_shift]<0){
                                src_node.vertices[entry.second+vertex_shift]=nVertices;
                                sequence_graph.AddVertices(1);
                                ++nVertices;
                            }
                            sequence_graph.AddEdge(src_node.vertices[entry.second+vertex_shift],src_node.vertices[bb+vertex_shift],basevector());
std::cout << "Skip Link " << src_node.vertices[entry.second+vertex_shift] << "->" << src_node.vertices[bb+vertex_shift] << std::endl;
                            ++nEdges;
                        }
                    }

                }

            }
        }
    }
    std::cout << "determining zero-length non-local branches: " << std::endl;
    for(size_t ii=0;ii<log.trail.size();++ii){
        int src = log.trail[log.trail.size()-1-ii].src;
        connection_assembly::node_t& src_node = ca.G.VertMutable(src);
        int front = -1;
        std::pair<int,int> last_pattern(-1,-1), pattern,last_remote_loc;
        for(int bb=0 ; bb < src_node.bases.isize(); ++bb){
            if(src_node.mapped_loc[bb]==std::make_pair(src,bb)){ //if location is unmapped
                std::set<std::pair<int,int>> set_of_interest;
//                for(auto entry: src_node.clones[bb]){
//                    if( entry.first != src && log.appearances[entry.first] > log.appearances[src]){
//                        set_of_interest.insert(entry);
//                    }
//                }
                for(auto entry: src_node.clones_tgt[bb]){
                    if( entry.first != src ){
                        set_of_interest.insert(entry);
                    }
                }
                ForceAssert(set_of_interest.size()<=1);

                if(set_of_interest.size()>0){
                    std::pair<int,int> remote_loc = *set_of_interest.begin();
                    pattern = remote_loc;
                    pattern.second-=bb;

                    if(front == -1){
                        front = bb;
                        if( src_node.vertices[bb]<0){
                            src_node.vertices[bb]=nVertices;
                            sequence_graph.AddVertices(1);
                            ++nVertices;
                        }
                        connection_assembly::node_t& remote_node = ca.G.VertMutable(remote_loc.first);
                        if( remote_node.vertices[remote_loc.second]<0){
                            remote_node.vertices[remote_loc.second]=nVertices;
                            sequence_graph.AddVertices(1);
                            ++nVertices;
                        }
                        sequence_graph.AddEdge(src_node.vertices[front], remote_node.vertices[remote_loc.second],basevector());
std::cout << "Remote Link " << src_node.vertices[front] << "->" << remote_node.vertices[remote_loc.second] << std::endl;
                        ++nEdges;
//                        sequence_graph.AddEdge(remote_node.vertices[remote_loc.second],src_node.vertices[front],basevector());
//std::cout << "Remote Link " << remote_node.vertices[remote_loc.second] << "->" << src_node.vertices[front] << std::endl;
//                        ++nEdges;

                    }
                    else if( pattern != last_pattern){
                        ForceAssert(src_node.vertices[front]>=0);
                        {
                            int back = bb-1;
                            if( src_node.vertices[back+1]<0){
                                src_node.vertices[back+1]=nVertices;
                                sequence_graph.AddVertices(1);
                                ++nVertices;
                            }
                            connection_assembly::node_t& remote_node = ca.G.VertMutable(last_remote_loc.first);
                            if( remote_node.vertices[last_remote_loc.second+1]<0){
                                remote_node.vertices[last_remote_loc.second+1]=nVertices;
                                sequence_graph.AddVertices(1);
                                ++nVertices;
                            }
                            sequence_graph.AddEdge(src_node.vertices[back+1],remote_node.vertices[last_remote_loc.second+1],basevector());
    std::cout << "Remote Link " << src_node.vertices[back+1] << "->" << remote_node.vertices[last_remote_loc.second+1] << std::endl;
                            ++nEdges;
                        }

                        {
                            front=bb;
                            if( src_node.vertices[front]<0){
                                src_node.vertices[front]=nVertices;
                                sequence_graph.AddVertices(1);
                                ++nVertices;
                            }
                            connection_assembly::node_t& remote_node = ca.G.VertMutable(remote_loc.first);
                            if( remote_node.vertices[remote_loc.second]<0){
                                remote_node.vertices[remote_loc.second]=nVertices;
                                sequence_graph.AddVertices(1);
                                ++nVertices;
                            }
                            sequence_graph.AddEdge(src_node.vertices[front], remote_node.vertices[remote_loc.second],basevector());
    std::cout << "Remote Link " << src_node.vertices[front] << "->" << remote_node.vertices[remote_loc.second] << std::endl;
                            ++nEdges;
                        }
                    }

                    last_pattern=pattern;
                    last_remote_loc=remote_loc;
                }
                else if(front >=0){
                    ForceAssert(src_node.vertices[front]>=0);
                    int back = bb-1;
                    if( src_node.vertices[back+1]<0){
                        src_node.vertices[back+1]=nVertices;
                        sequence_graph.AddVertices(1);
                        ++nVertices;
                    }
                    connection_assembly::node_t& remote_node = ca.G.VertMutable(last_remote_loc.first);
                    if( remote_node.vertices[last_remote_loc.second+1]<0){
                        remote_node.vertices[last_remote_loc.second+1]=nVertices;
                        sequence_graph.AddVertices(1);
                        ++nVertices;
                    }
                    sequence_graph.AddEdge(src_node.vertices[back+1],remote_node.vertices[last_remote_loc.second+1],basevector());
std::cout << "Remote Link " << src_node.vertices[back+1] << "->" << remote_node.vertices[last_remote_loc.second+1] << std::endl;
                    ++nEdges;
                    front = -1;
                }
            }
            else{
//                ForceAssert( src_node.clones[bb].size() == 0);
                ForceAssert( src_node.clones_src[bb].size() == 0);
                ForceAssert( src_node.clones_tgt[bb].size() == 0);
                if( front >=0 ){
                    ForceAssert(src_node.vertices[front]>=0);
                    int back = bb;
                    if( src_node.vertices[back+1]<0){
                        src_node.vertices[back+1]=nVertices;
                        sequence_graph.AddVertices(1);
                        ++nVertices;
                    }
                    connection_assembly::node_t& remote_node = ca.G.VertMutable(last_remote_loc.first);
                    if( remote_node.vertices[last_remote_loc.second+1]<0){
                        remote_node.vertices[last_remote_loc.second+1]=nVertices;
                        sequence_graph.AddVertices(1);
                        ++nVertices;
                    }
                    sequence_graph.AddEdge(src_node.vertices[back+1],remote_node.vertices[last_remote_loc.second+1],basevector());
std::cout << "Remote Link " << src_node.vertices[back+1] << "->" << remote_node.vertices[last_remote_loc.second+1] << std::endl;
                    ++nEdges;
//                    sequence_graph.AddEdge(remote_node.vertices[last_remote_loc.second+1],src_node.vertices[back+1],basevector());
//std::cout << "Remote Link " << remote_node.vertices[last_remote_loc.second+1] << "->" << src_node.vertices[back+1] << std::endl;
//                    ++nEdges;
                    front=-1;
                }
            }
        }
        if( front >=0 ){
            ForceAssert(src_node.vertices[front]>=0);
            int back = src_node.bases.isize()-1;
            if( src_node.vertices[back+1]<0){
                src_node.vertices[back+1]=nVertices;
                sequence_graph.AddVertices(1);
                ++nVertices;
            }
            connection_assembly::node_t& remote_node = ca.G.VertMutable(last_remote_loc.first);
            if( remote_node.vertices[last_remote_loc.second+1]<0){
                remote_node.vertices[last_remote_loc.second+1]=nVertices;
                sequence_graph.AddVertices(1);
                ++nVertices;
            }
            sequence_graph.AddEdge(src_node.vertices[back+1],remote_node.vertices[last_remote_loc.second+1],basevector());
std::cout << "Remote Link " << src_node.vertices[back+1] << "->" << remote_node.vertices[last_remote_loc.second+1] << std::endl;
            ++nEdges;
//            sequence_graph.AddEdge(remote_node.vertices[last_remote_loc.second+1],src_node.vertices[back+1],basevector());
//std::cout << "Remote Link " << remote_node.vertices[last_remote_loc.second+1] << "->" << src_node.vertices[back+1] << std::endl;
//            ++nEdges;
            front=-1;
        }
    }
    for(size_t ii=0;ii<log.trail.size();++ii){
        int src = log.trail[log.trail.size()-1-ii].src;
std::cout << "tracing sequence " << src << std::endl;
        connection_assembly::node_t& src_node = ca.G.VertMutable(src);

        int front = -1, back=-1;
/*
if( src == 0 & false){
        for(int bb=0 ; bb < src_node.bases.isize(); ++bb){
            if(src_node.mapped_loc[bb]==std::make_pair(src,bb)){
              std::cout << 1;
            }
            else{
              std::cout << 0;
            }
        }
        std::cout << std::endl;
        std::cout << std::endl;

        for(int bb=0 ; bb <= src_node.bases.isize(); ++bb){
            if(src_node.mapped_branch[bb].first==std::make_pair(src,bb).first){
              std::cout << 1;
            }
            else{
              std::cout << 0;
            }
        }
        std::cout << std::endl;
        std::cout << std::endl;

        for(int bb=0 ; bb <= src_node.bases.isize(); ++bb){
            if(src_node.mapped_branch[bb].second==std::make_pair(src,bb).second){
              std::cout << 'o';
            }
            else{
              std::cout << src_node.mapped_branch[bb].second << " ";
            }
        }
        std::cout << std::endl;
        std::cout << std::endl;
        for(int bb=0 ; bb <= src_node.bases.isize(); ++bb){
            if(src_node.vertices[bb]==-1){
              std::cout << '-';
            }
            else{
              std::cout << '|';
            }
        }
        std::cout << std::endl;
}
*/
        std::set<int> src_vertices;
        std::set<int> des_vertices;
        for(int bb=0 ; bb < src_node.bases.isize(); ++bb){
            if(src_node.mapped_loc[bb]==std::make_pair(src,bb)){ //if location is unmapped
                if(front<0){ //if there is no unmapped loc before, edge starts here
                    front = bb;

                    src_vertices.clear();
                    DetermineTargetVertices( src_vertices, nVertices,ca,log,src,bb);
                    sequence_graph.AddVertices( nVertices-sequence_graph.N());
                    ForceAssert( src_vertices.find(src_node.vertices[bb]) != src_vertices.end());
                }

                if(  src_node.vertices[bb+1]>=0){ // if there's a vertex after loc bb, edge ends here
                    back=bb;

                    // double check vertices
                    des_vertices.clear();
                    DetermineTargetVertices( des_vertices, nVertices,ca,log,src,bb+1);
                    sequence_graph.AddVertices( nVertices-sequence_graph.N());

                    ForceAssert( des_vertices.find(src_node.vertices[bb+1]) != des_vertices.end());

                    for( auto start_vertex: src_vertices){
                        for( auto end_vertex: des_vertices){
                        sequence_graph.AddEdge(start_vertex,end_vertex,basevector(src_node.bases,front,back-front+1));
std::cout << "Edge " << start_vertex << "->" << end_vertex << " " << front << " " << back << std::endl;
                            ++nEdges;
                        }
                    }

                    front = -1;
                }
            }
            else{//if location is mapped
                if(front>=0){ // if front is set, the edge ends before the loc bb
                  // edge from front to bb-1
                    back=bb-1;

                    des_vertices.clear();
                    DetermineTargetVertices( des_vertices, nVertices,ca,log,src,bb);
                    sequence_graph.AddVertices( nVertices-sequence_graph.N());
                    ForceAssert( des_vertices.find(src_node.vertices[bb]) != des_vertices.end());

                    for( auto start_vertex: src_vertices){
                        for( auto end_vertex: des_vertices){
                        sequence_graph.AddEdge(start_vertex,end_vertex,basevector(src_node.bases,front,back-front+1));
std::cout << "Edge " << start_vertex << "->" << end_vertex << " " << front << " " << back << std::endl;
                            ++nEdges;
                        }
                    }

                    front = -1;
                }
            }
        }
        // this is for the last branch point, after the last loc
        if( front >= 0){
            int bb=src_node.bases.isize();
            back = bb-1;

            des_vertices.clear();
            DetermineTargetVertices( des_vertices, nVertices,ca,log,src,bb);
            sequence_graph.AddVertices( nVertices-sequence_graph.N());
            ForceAssert( des_vertices.find(src_node.vertices[bb]) != des_vertices.end());
            for( auto start_vertex: src_vertices){
                for( auto end_vertex: des_vertices){
                    sequence_graph.AddEdge(start_vertex,end_vertex,basevector(src_node.bases,front,back-front+1));
std::cout << "Edge " << start_vertex << "->" << end_vertex << " " << front << " " << back << std::endl;
                    ++nEdges;
                }
            }

            front = -1;
        }
    }



    ForceAssert(nVertices==sequence_graph.N());


std::cout<< "number of vertices: " << nVertices << std::endl;
std::cout<< "number of edges: " << nEdges << std::endl;

    return bModified;
}


std::pair<int,int>& SyncMappedLoc(connection_assembly& ca, std::pair<int,int> loc){
    ForceAssert(loc.second < ca.G.Vert(loc.first).bases.isize());
    auto& mapped=ca.G.VertMutable(loc.first).mapped_loc[loc.second];

    ForceAssert(  ca.G.Vert(loc.first).bases[loc.second]
                ==ca.G.Vert(mapped.first).bases[mapped.second]
               ) ;
    if(mapped!=loc){
//        if( ca.G.Vert(loc.first).clones[loc.second].size() != 0){
        if( ca.G.Vert(loc.first).clones_tgt[loc.second].size() != 0 || ca.G.Vert(loc.first).clones_src[loc.second].size() != 0){
            std::cout << loc.first << " " << loc.second << "->" << mapped.first << " " << mapped.second << std::endl;
        }
        ForceAssert( ca.G.Vert(loc.first).clones_src[loc.second].size() == 0);
        ForceAssert( ca.G.Vert(loc.first).clones_tgt[loc.second].size() == 0);
        mapped = SyncMappedLoc(ca,mapped);
    }
    return mapped;
}
/*
std::pair<int,int>& SyncMappedBranch(connection_assembly& ca, std::pair<int,int> loc){
    ForceAssert(loc.second <= ca.G.Vert(loc.first).bases.isize());
    auto& mapped=ca.G.VertMutable(loc.first).mapped_branch[loc.second];
    if(mapped!=loc){
        mapped = SyncMappedBranch(ca,mapped);
    }
    return mapped;
}
*/



void mark_log( branch_log& log, int root, connection_assembly& ca
             , int offset, int match_front, int match_back){
    connection_assembly::node_t& root_node = ca.G.VertMutable(root);
    ForceAssert(!root_node.bProcessed);

    ForceAssert(log.stack_appearances.find(root)==log.stack_appearances.end());

    ForceAssert(log.appearances.find(root)==log.appearances.end());

    if( log.appearances.find(root) == log.appearances.end()){
        const size_t tmp = log.appearances.size();
        log.appearances[root] = tmp;
    }
    log.stack_appearances.insert(root);

    size_t loc_trail_index=log.trail.size();
    log.trail.push(branch_log::stack_element(root));
    std::cout << " pushing " << root << std::endl;

    //go left
    for(size_t conn_index=0;conn_index<ca.G.To(root).size();++conn_index){
        int des = ca.G.To(root)[conn_index];
        const connection_assembly::node_t& des_node = ca.G.Vert(des);
        ForceAssert(des_node.bConnected);
        ForceAssert(!des_node.bProcessed);
        const connection& conn = ca.G.EdgeObjectByIndexTo( root, conn_index );
        ForceAssert(conn.overlap > 0);

        if( des==root){
            log.trail[loc_trail_index].edge_overlap_index_dir_to_self.push_back(make_tuple(conn.overlap,conn_index,0));
        }
        else if(log.appearances.find(des)==log.appearances.end()){
//            log.trail[loc_trail_index].edge_indices_to_childs.push_back(conn_index);

std::cout
          << des << "(" << ca.G.Vert(des).bases.isize() << ") L\t"
          << root << "(" << root_node.bases.isize() << ")\t"
          << conn.overlap << "\t"
          << std::endl;

            mark_log(log,des,ca,0,0,0);
        }

        ForceAssert(log.appearances.find(des)!=log.appearances.end());
        if(des!=root&&log.appearances[root]>log.appearances[des]){
            log.trail[loc_trail_index].edge_flags_index_dir_to_ancestors.push_back(make_tuple(log.appearances[des],conn_index,0));
        }
        if(des!=root&&log.appearances[root]<log.appearances[des]){
            log.trail[loc_trail_index].edge_flags_index_dir_to_descendants.push_back(make_tuple(log.appearances[des],conn_index,0));
        }
    }

    //go right
    for(size_t conn_index=0;conn_index<ca.G.From(root).size();++conn_index){
        int des = ca.G.From(root)[conn_index];
        const connection_assembly::node_t& des_node = ca.G.Vert(des);
        ForceAssert(des_node.bConnected);
        ForceAssert(!des_node.bProcessed);
        const connection& conn = ca.G.EdgeObjectByIndexFrom( root, conn_index );
        ForceAssert(conn.overlap > 0);

        if( des==root){
            log.trail[loc_trail_index].edge_overlap_index_dir_to_self.push_back(make_tuple(conn.overlap,conn_index,1));
        }
        else if(log.appearances.find(des)==log.appearances.end()){
//            log.trail[loc_trail_index].edge_indices_to_childs.push_back(conn_index);
std::cout << root << "(" << root_node.bases.isize() << ") R\t"
          << des << "(" << ca.G.Vert(des).bases.isize() << ")\t"
          << conn.overlap << "\t"
          << std::endl;

            mark_log(log,des,ca,0,0,0);
        }

        ForceAssert(log.appearances.find(des)!=log.appearances.end());
        if(des!=root&&log.appearances[root]>log.appearances[des]){
            log.trail[loc_trail_index].edge_flags_index_dir_to_ancestors.push_back(make_tuple(log.appearances[des],conn_index,1));
        }
        if(des!=root&&log.appearances[root]<log.appearances[des]){
            log.trail[loc_trail_index].edge_flags_index_dir_to_descendants.push_back(make_tuple(log.appearances[des],conn_index,1));
        }
    }
    for(size_t vv=0; vv<root_node.variants.size();++vv){
        int variant_idx= std::get<0>(root_node.variants[vv]);
        std::cout << " pushing variant " << variant_idx << std::endl;
        log.trail.push(branch_log::stack_element(variant_idx));
    }
//    std::cout << "done" << std::endl;
    log.stack_appearances.erase(root);
}

int MakeSequenceGraph(digraphE<BaseVec2KmerPath>&sequence_graph, connection_assembly CA){
    for( int src=0 ; src < CA.G.N() ; ++src){
        const connection_assembly::node_t& src_node = CA.G.Vert(src);
        if( std::get<0>(src_node.variant_of )<0){
            SetOrientation(src,CA);

        }
    }
    AccommodateRC(CA);
    PreprocessGraph(CA);
    MakeSequenceGraphConnected(sequence_graph,CA);

    return 0;

}

Bool ToHbv( connection_assembly& A, const int K, HyperBasevector& hb, 
     vec<int>& inv, const Bool verbose )
{    


     // First define the edge objects of hb.

     vec<basevector> bases( 2 * A.bases.size( ) );
     for ( size_t j = 0; j < A.bases.size( ); j++ )
     {    bases[j] = A.bases[j];
          bases[ j + A.bases.size( ) ] = A.bases[j];
          bases[ j + A.bases.size( ) ].ReverseComplement( );    }
     int n = A.bases.size( );

     {
     Ofstream( out, "xxx.fasta" );
     for ( size_t j = 0; j < A.bases.size( ); j++ )
          A.bases[j].Print( out, j );
     }

     const int max_length = 10000000;
     vec<int> left_trim(n, -1), right_trim(n, -1);
     vec<pair<int, int> > left(n, make_pair(0,max_length)), right(n, make_pair(0, max_length));
     

     // diagnostic
     vec< std::pair<int,int> > ancestor(n*2,make_pair(-3,-3));

     enum {BUILD, COMPUTE, VALIDATE, COMPLETE, ERROR} status = BUILD;
     bool updated = false;
     size_t nError = 0;
     do {
       updated = false;
       for ( int edge = 0; edge < n; edge++ ) {   // for each edge
         const int edge_length = A.bases[edge].size();
         for ( int conn_index = 0; conn_index < A.G.From(edge).isize( ); conn_index++ ) {  // for each connection
           const int dest_edge = A.G.From(edge)[conn_index];
           const int dest_edge_length = A.bases[dest_edge].size();
           const connection& c = A.G.EdgeObjectByIndexFrom( edge, conn_index );
           const int overlap = c.overlap;
           pair<int, int>& l = (c.rc1 ? left[edge] : right[edge] );
           pair<int, int>& r = (c.rc2 ? right[dest_edge] : left[dest_edge] );
           int& trim_l = (c.rc1 ? left_trim[edge] : right_trim[edge]);
           int& trim_r = (c.rc2 ? right_trim[dest_edge] : left_trim[dest_edge]);

           int l_index = 2*edge + int(!c.rc1)   ;
           int r_index = 2*dest_edge + int(c.rc2);

           if (status == BUILD ) {  // round 1, build overlap table
             const pair<int, int> l_original = l;
             const pair<int, int> r_original = r;
             int l_max = edge_length - (c.rc1 ? right[edge].first : left[edge].first );
             int r_max = dest_edge_length - (c.rc2 ? left[dest_edge].first : right[dest_edge].first );
             l.first = Max(l.first, overlap - r.second);
             r.first = Max(r.first, overlap - l.second);
             l.second = Min(overlap - r.first, l_max);
             r.second = Min(overlap - l.first, r_max);
             if (l_original != l || r_original != r)
               updated = true;
           } else if (status == COMPUTE) {  // round 2, compute edge trims
             if (trim_l == -1 && trim_r == -1) {
std::cout<< "TRIM " << edge << "  \t  \t" << dest_edge << " " << int(c.rc1) << " " << int(c.rc2) << std::endl;
               trim_l = (l.first + l.second) / 2 - (K -1) / 2;
               trim_r = overlap - trim_l - (K - 1);

               ForceAssert( ancestor[l_index] == make_pair(-3,-3));
               ForceAssert( ancestor[r_index] == make_pair(-3,-3));
               ancestor[l_index] = make_pair(dest_edge,int(c.rc2));
               ancestor[r_index] = make_pair(edge,int(!c.rc1));

             }  else if (trim_l == -1 && trim_r != -1) {
std::cout<< "TRIM " << edge << "  \t<-\t" << dest_edge << " " << int(c.rc1) << " " << int(c.rc2) << std::endl;
               trim_l = overlap - trim_r - (K - 1);

               ForceAssert( ancestor[l_index] == make_pair(-3,-3));
               ForceAssert( ancestor[r_index] != make_pair(-3,-3));
               ancestor[l_index] = make_pair(dest_edge,int(c.rc2));
             } else if (trim_l != -1 && trim_r == -1) {
std::cout<< "TRIM " << edge << "  \t->\t" << dest_edge << " " << int(c.rc1) << " " << int(c.rc2) << std::endl;
               trim_r = overlap - trim_l - (K - 1);
               ForceAssert( ancestor[l_index] != make_pair(-3,-3));
               ForceAssert( ancestor[r_index] == make_pair(-3,-3));
               ancestor[r_index] = make_pair(edge,int(!c.rc1));
             }

             if( trim_r + trim_l != overlap - (K-1)){
std::cout<< "TRIM " << edge << "  \t!!\t" << dest_edge << " " << int(c.rc1) << " " << int(c.rc2) << std::endl;
             }
             if (left_trim[edge] + right_trim[edge] > edge_length - K+1 ) {
std::cout<< "TRIML" << edge << "  \t!!\t" << left_trim[edge] << " " << right_trim[edge] << " " << edge_length << " " << K << " " << int(c.rc1) << " " << int(c.rc2) << std::endl;
             }
             if (left_trim[dest_edge] + right_trim[dest_edge] > dest_edge_length - K+1) {
std::cout<< "TRIMR" << dest_edge << "  \t!!\t" << left_trim[dest_edge] << " " << right_trim[dest_edge] << " " << dest_edge_length << " " << K << " " << int(c.rc1) << " " << int(c.rc2) << std::endl;
             }
           } else if (status == VALIDATE) {  // round 3, validate edge trims
             if (trim_l + trim_r != overlap - (K-1)) {
               cout << "Edge " << edge << ( c.rc1 ? " left trim" : " right trim" ) << " + "
                    << "Edge " << dest_edge << ( c.rc2 ? " right trim" : " left trim" )
                    << " != overlap - (K-1)" << endl;
               cout << trim_l << " + " << trim_r << " != " << overlap << " - (" << K << " - 1)" << endl;
               ++nError;
//               status = ERROR;
             }
           }
         }
         if (status == COMPUTE) {  // round 2, compute edge trims part 2
           left_trim[edge] = (left_trim[edge] == -1 ? 0 : left_trim[edge]);
           right_trim[edge] = (right_trim[edge] == -1 ? 0 : right_trim[edge]);
         } else if (status == VALIDATE)  {  // round 3, validate edge trims part 2
           size_t old_error=nError;
           if ( left_trim[edge] < 0) {
             cout << "Left Trim < 0 Edge " << edge << endl;
             ++nError;
           }
           if ( right_trim[edge] < 0) {
             cout << "right Trim < 0 Edge " << edge << endl;
             ++nError;
           }
           if (left_trim[edge] + right_trim[edge] > edge_length - K+1) {
             cout << "Left Trim + Right Trim > Edge Size - K +1 for Edge " << edge << endl;
             cout << left_trim[edge] << " + " << right_trim[edge] << " > "
                  << A.bases[edge].size() << " - " << K-1  << endl;
             if(old_error==nError){
                 ++nError;
             }
//             status = ERROR;
           }
         }
       }
       
       // Update status
       if (status == BUILD && updated == false)
         status = COMPUTE;
       else if (status == COMPUTE)
         status = VALIDATE;
       else if (status == VALIDATE)
         status = COMPLETE;
     } while (status != COMPLETE && status != ERROR);


     if (nError){
       std::cout << "number of errors " << nError << std::endl;
       FatalErr("A problem occured whilst computing overlaps. See above.");
     }
     






     ////////////////////////////////////////////////////////////////////////////////////
     // Various bits of old code - no idea what it all does, so I'll leave it be for now.
     
     
     // Solve overlap equations.  Trim is left_trim[0], right_trim[0], ... .

     // vec<int> trim(2*n, -1);

     // for ( int e = 0; e < n; e++ )
     // for ( int l = 0; l < A.G.From(e).isize( ); l++ )
     // {    int f = A.G.From(e)[l];
     //      const connection& c = A.G.EdgeObjectByIndexFrom( e, l );
     //      cout << ( c.rc1 ? "L" : "R" ) << e << " + " << ( c.rc2 ? "R" : "L" ) << f
     //           << " = " << c.overlap - (K-1) << "\n";    
     //      int v1 = 2*e + ( c.rc1 ? 0 : 1 ), v2 = 2*f + ( c.rc2 ? 1 : 0 );
     //      int len1 = bases[e].size( ), len2 = bases[f].size( );
     //      int rhs = c.overlap - (K-1);
     //      if ( trim[v1] >= 0 && trim[v2] >= 0 )
     //      {    if ( trim[v1] + trim[v2] == rhs ) continue;
     //           PRINT4( v1, v2, trim[v1], trim[v2] ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXX
     //           FatalErr("Unsatisfiable constraint (1)."); }
     //      else if ( trim[v1] < 0 & trim[v2] >= 0 )
     //      {    int t1 = rhs - trim[v2];
     //           if ( t1 < 0 )
     //           {    FatalErr("Unsatisfiable constraint (2).");  }
     //           cout << "setting trim[" << v1 << "] to " << t1 << endl; // XXXXXXXXXX
     //           trim[v1] = t1;    }
     //      else if ( trim[v2] < 0 & trim[v1] >= 0 )
     //      {    int t2 = rhs - trim[v1];
     //           if ( t2 < 0 )
     //           {    PRINT(t2);
     //                FatalErr("Unsatisfiable constraint (3).");   }
     //           cout << "setting trim[" << v2 << "] to " << t2 << endl; // XXXXXXXXXX
     //           trim[v2] = t2;    }
     //      else if ( trim[v1] < 0 && trim[v2] < 0 )
     //      {    int t1 = rhs/2;
     //           int t2 = rhs - t1;
     //           if ( t1 < 0 )
     //           {    FatalErr("Unsatisfiable constraint (4).");    }
     //           cout << "setting trim[" << v1 << "] to " << t1 << endl; // XXXXXXXXXX
     //           cout << "setting trim[" << v2 << "] to " << t2 << endl; // XXXXXXXXXX
     //           trim[v1] = t1, trim[v2] = t2;    }    }
     // for ( int e = 0; e < n; e++ )
     // {    int &left = trim[ 2*e ], &right = trim[ 2*e + 1 ];
     //      if ( left < 0 ) left = 0;
     //      if ( right < 0 ) right = 0;
     //      if ( left + right >= bases[e].isize( ) )
     //      {    FatalErr("Illegal trimming.");    }    }

     // cout << "\ntrims\n";
     // for ( int e = 0; e < n; e++ )
     // {    int left = trim[ 2*e ], right = trim[ 2*e + 1 ];
     //      PRINT3( e, left, right );    }
     // cout << "\n";
     // flush(cout);
     // return True; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     // }

     // Trim basevectors.

     for ( size_t j = 0; j < A.bases.size( ); j++ ) { 
       PRINT3(A.bases[j].size(), left_trim[j], right_trim[j]); 
       bases[j].assign( A.bases[j].begin() + left_trim[j],
			A.bases[j].end() - right_trim[j]);
       PRINT(bases[j].size()); 
       //       bases[ j + A.bases.size( ) ] = A.bases[j];
       //       bases[ j + A.bases.size( ) ].ReverseComplement( ); 
     }

     // // Now define vertices, which are equivalence classes or seq or rc-seq ends.

     // equiv_rel e( 4*n );
     // for ( int v1 = 0; v1 < n; v1++ )
     // for ( int j = 0; j < A.G.From(v1).isize( ); j++ )
     // {    int v2 = A.G.From(v1)[j];
     //      const connection& c = A.G.EdgeObjectByIndexFrom( v1, j );
     //      Bool rc1 = c.rc1, rc2 = c.rc2;
     //      int over = c.overlap;
     //      int e1 = v1 + ( rc1 ? n : 0 ), e2 = v2 + ( rc2 ? n : 0 );
     //      e.Join( 2*e1 + 1, 2*e2 );    }
     // vec<int> reps;
     // e.OrbitRepsAlt(reps);
     // int N = reps.size( );

     // // Now define edges.

     // vec< vec<int> > from(N), to(N), from_edge_obj(N), to_edge_obj(N);
     // for ( int x = 0; x < 2*n; x++ )
     // {    int v1 = BinPosition( reps, e.ClassId( 2*x ) );
     //      int v2 = BinPosition( reps, e.ClassId( 2*x + 1 ) );
     //      from[v1].push_back(v2), to[v2].push_back(v1);
     //      from_edge_obj[v1].push_back(x);
     //      to_edge_obj[v2].push_back(x);    }

     // // For every vertex, define trimming.

     // for ( int v = 0; v < N; v++ )
     // {    vec<int> o;
     //      e.Orbit( v, o );
     //      vec<int> ins, outs;
     //      for ( int j = 0; j < o.isize( ); j++ )
     //      {    int r = o[j];
     //           if ( r % 2 == 1 ) ins.push_back(r/2);
     //           else outs.push_back(r/2);    }
     //      if ( ins.empty( ) || outs.empty( ) ) continue;
     
     //      vec<int> all;
     //      for ( int j = 0; j < ins.isize( ); j++ )
     //      {    if ( ins[j] < n ) all.push_back( ins[j] );
     //           else all.push_back( ins[j]-n );    }
     //      for ( int j = 0; j < outs.isize( ); j++ )
     //      {    if ( outs[j] < n ) all.push_back( outs[j] );
     //           else all.push_back( outs[j]-n );    }
     //      UniqueSort(all);
     //      cout << "\n";
     //      for ( int j = 0; j < all.isize( ); j++ )
     //           cout << all[j] << "[l=" << bases[ all[j] ].size( ) << "]\n";

     //      cout << "ins:";
     //      for ( int j = 0; j < ins.isize( ); j++ )
     //      {    if ( ins[j] < n ) cout << " " << ins[j] << "fw";
     //           else cout << " " << ins[j]-n << "rc";    }
     //      cout << "; outs:";
     //      for ( int j = 0; j < outs.isize( ); j++ )
     //      {    if ( outs[j] < n ) cout << " " << outs[j] << "fw";
     //           else cout << " " << outs[j]-n << "rc";    }
     //      cout << "\n";    
     //      for ( int j = 0; j < ins.isize( ); j++ )
     //      {    int e = ins[j];
     //           int e0 = e;
     //           if ( e >= n ) e -= n;
     //           for ( int l = 0; l < A.G.From(e).isize( ); l++ )
     //           {    int f = A.G.From(e)[l];
     //                const connection& c = A.G.EdgeObjectByIndexFrom( e, l );
     //                if ( ( e0 >= n ) != c.rc1 ) continue;
     //                if ( !Member( outs, f + ( c.rc2 ? n : 0 ) ) ) continue;
     //                cout << e << ( c.rc1 ? "rc" : "fw" ) << " overlaps "
     //                     << f << ( c.rc2 ? "rc" : "fw" ) << " by "
     //                     << c.overlap << " bases, "
     //                     << "so " << ( c.rc1 ? "L" : "R" ) << e << " + "
     //                     << ( c.rc2 ? "R" : "L" ) << f
     //                     << " = " << c.overlap - (K-1) << "\n";    }    }    }
     return True;  
}

void MatchRef( const connection_assembly& A, const vecbasevector& G,
     vec< triple< int, ho_interval, vec<int> > >& matches )
{
     // Find perfect placements of edges or their reverse complements.

     int n = A.bases.size( );
     vec< vec< pair<int,int> > > P(2*n);
     {    const int LG = 12;
          VecIntPairVec Glocs;
          CreateGlocs( G, LG, Glocs );
          for ( size_t e = 0; e < A.bases.size( ); e++ )
          {    vec<placementy> places 
                    = FindGenomicPlacementsY( e, A.bases[e], LG, G, Glocs );
               for ( int j = 0; j < places.isize( ); j++ )
               {    const placementy& p = places[j];
                    P[ 2*e + ( p.fw ? 0 : 1 ) ].push( p.g, p.pos );    }    }    }

     // Form a digraph whose vertices are perfect placements of fw/rc edges and
     // whose edges are induced by connections consistent with the placements.

     int N = 0;
     vec< triple<int,int,int> > V;
     for ( int j = 0; j < P.isize( ); j++ )
     {    N += P[j].size( );
          for ( int l = 0; l < P[j].isize( ); l++ )
               V.push( j, P[j][l].first, P[j][l].second );    }
     vec< vec<int> > from(N), to(N);
     for ( int v1 = 0; v1 < n; v1++ )
     {    for ( int j = 0; j < A.G.From(v1).isize( ); j++ )
          {    int v2 = A.G.From(v1)[j];
               const connection& c = A.G.EdgeObjectByIndexFrom(v1, j);
               int x1, x2;
               for ( int pass = 1; pass <= 2; pass++ )
               {    if ( pass == 1 )
                    {    x1 = 2*v1 + ( c.rc1 ? 1 : 0 ); 
                         x2 = 2*v2 + ( c.rc2 ? 1 : 0 );    }
                    else
                    {    x2 = 2*v1 + ( c.rc1 ? 0 : 1 ); 
                         x1 = 2*v2 + ( c.rc2 ? 0 : 1 );    }
                    for ( int l1 = 0; l1 < P[x1].isize( ); l1++ )
                    for ( int l2 = 0; l2 < P[x2].isize( ); l2++ )
                    {    int g = P[x1][l1].first;
                         if ( g != P[x2][l2].first ) continue;
                         int start1 = P[x1][l1].second, start2 = P[x2][l2].second;
                         int stop1 = start1 + A.bases[x1/2].isize( );
                         int stop2 = start2 + A.bases[x2/2].isize( );
                         if ( stop1 - start2 != c.overlap ) continue;
                         int m1 = l1, m2 = l2;
                         for ( int u = 0; u < x1; u++ )
                              m1 += P[u].size( );
                         for ( int u = 0; u < x2; u++ )
                              m2 += P[u].size( );
                         from[m1].push_back(m2); 
                         to[m2].push_back(m1);    }    }    }    }
     digraph H( from, to );

     // Find all source-sink paths through the graph.

     matches.clear( );
     vec< vec<int> > paths;
     H.AllPaths( -1, -1, paths );
     for ( int j = 0; j < paths.isize( ); j++ )
     {    const vec<int>& p = paths[j];
          int g = V[ p.front( ) ].second, pos = V[ p.front( ) ].third;
          int Pos = V[ p.back( ) ].third
               + A.bases[ V[ p.back( ) ].first / 2 ].isize( );
          vec<int> q;
          for ( int l = 0; l < paths[j].isize( ); l++ )
          {    int e = V[ paths[j][l] ].first;
               if ( e % 2 == 0 ) q.push_back(e/2);
               else q.push_back( -e/2-1 );    }
          matches.push( g, ho_interval( pos, Pos ), q );    }
     Sort(matches);    }
/*
std::pair<std::tuple<int,int,int>,std::tuple<int,int,int> > getOverlapBranchRange(bool bRC1,int iLength1, bool bRC2, int iLength2,int overlap){
    ForceAssert(iLength1>=overlap);
    int sStep,sFront,sBack;
    if(bRC1){
        sStep=-1;
        sFront=overlap;
        sBack= 0;
    }
    else{
        sStep=1;
        sFront=iLength1-overlap;
        sBack = iLength1;
    }

   int dStep,dFront,dBack;
   if(bRC2){
       dStep=-1;
       dFront=iLength2;
       dBack=iLength2-overlap;
   }
   else{
       dStep=1;
       dFront=0;
       dBack = overlap;
   }
   return std::make_pair(  std::make_tuple(sFront,sBack,sStep) , std::make_tuple(dFront,dBack,dStep));
}
*/
/*
std::pair<std::tuple<int,int,int>,std::tuple<int,int,int> > getOverlapBaseRange(bool bRC1,int iLength1, bool bRC2, int iLength2,int overlap){
    ForceAssert(iLength1>=overlap);
    int sStep,sFront,sBack;
    if(bRC1){
        sStep=-1;
        sFront= overlap-1;
        sBack= 0 ;
    }
    else{
        sStep=1;
        sFront= iLength1-overlap;
        sBack = iLength1-1;
    }

   int dStep,dFront,dBack;
   if(bRC2){
       dStep=-1;
       dFront=iLength2-1;
       dBack=iLength2-overlap;
   }
   else{
       dStep=1;
       dFront=0;
       dBack = overlap-1;
   }
   return std::make_pair(  std::make_tuple(sFront,sBack,sStep) , std::make_tuple(dFront,dBack,dStep));
}
*/

std::pair<std::pair<int,int>,std::pair<int,int> > getOverlapBaseRange(int iLength1, int iLength2,int overlap,int dir){
    ForceAssert(iLength1>=overlap);
    ForceAssert(iLength2>=overlap);

    int sBegin,sEnd,dBegin,dEnd;

    if( dir ){ // 1-2
        sBegin=iLength1-overlap;
        sEnd=iLength1;
        dBegin=0;
        dEnd=overlap;
    }
    else{ // 2-1
        sBegin=0;
        sEnd=overlap;
        dBegin=iLength2-overlap;
        dEnd=iLength2;
    }
    return std::make_pair(  std::make_pair(sBegin,sEnd) , std::make_pair(dBegin,dEnd));
}
/*
std::pair<std::pair<int,int>,std::pair<int,int> > getOverlapBranchRange(int iLength1, int iLength2,int overlap,int dir){
    ForceAssert(iLength1>=overlap);
    ForceAssert(iLength2>=overlap);

    int sBegin,sEnd,dBegin,dEnd;
    if( dir ){ // 1-2
        sBegin=iLength1-overlap;
        sEnd=iLength1+1;
        dBegin=0;
        dEnd=overlap+1;
    }
    else{ // 2-1
        sBegin=0;
        sEnd=overlap+1;
        dBegin=iLength2-overlap;
        dEnd=iLength2+1;
    }
    return std::make_pair(  std::make_pair(sBegin,sEnd) , std::make_pair(dBegin,dEnd));
}
*/
/*
void MergeBranch( int tgt_idx, int tgt_pos, int src_idx, int src_pos, connection_assembly& ca, branch_log& log){
    if( tgt_idx==src_idx && tgt_pos == src_pos) return;
//std::cout << "Merging " << tgt_idx << " " << tgt_pos << " " << src_idx <<" " << src_pos << std::endl;

    ForceAssert(tgt_idx!=src_idx || tgt_pos !=src_pos);
    connection_assembly::node_t& tgt_node = ca.G.VertMutable(tgt_idx);
    std::pair<int,int>& tgt_branch=tgt_node.mapped_branch[tgt_pos];

    connection_assembly::node_t& src_node = ca.G.VertMutable(src_idx);
    std::pair<int,int>& src_branch=src_node.mapped_branch[src_pos];

    ForceAssert(log.appearances.find(tgt_branch.first) != log.appearances.end());
    ForceAssert(log.appearances.find(src_branch.first) != log.appearances.end());
    ForceAssert(log.appearances[tgt_branch.first] <= log.appearances[src_branch.first]);

    if( tgt_branch.first != tgt_idx || tgt_branch.second != tgt_pos){ // if tgt's location has been mapped
        ForceAssert(  tgt_node.zero_length_destination.find(tgt_pos)==tgt_node.zero_length_destination.end()
                    ||tgt_node.zero_length_destination[tgt_pos].size()==0);
    }

    if( src_branch.first == src_idx && src_branch.second == src_pos){ // if src's location has not been mapped
        auto sSeek = src_node.zero_length_destination.find(src_pos);
        if( sSeek!= src_node.zero_length_destination.end()){
            std::stack< std::pair<int,int>> to_delete;
            for(auto& entry: (*sSeek).second){ // for each connection from src
                int neighbor_of_src = entry.first;
                int nos_pos = entry.second;
                connection_assembly::node_t& nos_node = ca.G.VertMutable(neighbor_of_src);

                auto& real_nos_loc = src_node.mapped_branch[nos_pos];
                int real_neighbor_of_src = real_nos_loc.first;
                int real_nos_pos = real_nos_loc.second;

//std::cout << "nos  " << neighbor_of_src << " " << nos_pos << std::endl;
//std::cout << "rnos " << real_neighbor_of_src << " " << real_nos_pos << std::endl;
                ForceAssert(neighbor_of_src==real_neighbor_of_src && nos_pos==real_nos_pos);

                connection_assembly::node_t& real_nos_node = ca.G.VertMutable(real_neighbor_of_src);

                ForceAssert(  real_nos_node.zero_length_destination[real_nos_pos].find(std::make_pair(src_idx,src_pos))
                            !=real_nos_node.zero_length_destination[real_nos_pos].end()
                           );

                int real_tgt = tgt_branch.first;
                int real_tgt_pos = tgt_branch.second;
                connection_assembly::node_t& real_tgt_node = ca.G.VertMutable(real_tgt);
                if(real_tgt!=tgt_idx || real_tgt_pos!=tgt_pos){
                    ForceAssert(  tgt_node.zero_length_destination.find(tgt_pos)==tgt_node.zero_length_destination.end()
                                ||tgt_node.zero_length_destination[tgt_pos].size()==0);
                }

                real_nos_node.zero_length_destination[real_nos_pos].erase(std::make_pair(src_idx,src_pos));
                to_delete.push(std::make_pair(real_neighbor_of_src,real_nos_pos));

                real_nos_node.zero_length_destination[real_nos_pos].insert(std::make_pair(real_tgt,real_tgt_pos));

                real_tgt_node.zero_length_destination[real_tgt_pos].insert(std::make_pair(real_neighbor_of_src,real_nos_pos));
            }
            for(;!to_delete.empty();to_delete.pop()){
                src_node.zero_length_destination[src_pos].erase(to_delete.top());
            }
        }
    }
    else{ //there should be no entries if the node has been mapped
        ForceAssert(  src_node.zero_length_destination.find(src_pos)==src_node.zero_length_destination.end()
                    ||src_node.zero_length_destination[src_pos].size()==0);
    }
    src_branch=tgt_branch;
}
*/

void MergeLoc( std::pair<int,int>& tgt_loc, std::pair<int,int>& src_loc, connection_assembly& ca, branch_log& log, bool bStrict){
    connection_assembly::node_t& tgt_node = ca.G.VertMutable(tgt_loc.first);
    connection_assembly::node_t& src_node = ca.G.VertMutable(src_loc.first);
    ForceAssert( tgt_node.mapped_loc[tgt_loc.second] == tgt_loc );
    ForceAssert( src_node.mapped_loc[src_loc.second] == src_loc );

    if( tgt_loc == src_loc ) return;
    ForceAssert(log.appearances.find(tgt_loc.first) != log.appearances.end());
    ForceAssert(log.appearances.find(src_loc.first) != log.appearances.end());
    if(bStrict){
        ForceAssert(log.appearances[tgt_loc.first] < log.appearances[src_loc.first]);
    }

//    if( src_node.clones[src_loc.second].size() > 0){
//        for( auto& entry: src_node.clones[src_loc.second] ){
//            ForceAssert( entry != src_loc);
//
//            auto& clone_list = ca.G.VertMutable(entry.first).clones[entry.second];
//            ForceAssert(clone_list.find(src_loc)!=clone_list.end());
//            clone_list.erase(src_loc);
//
//            clone_list.insert(tgt_loc);
//            tgt_node.clones[tgt_loc.second].insert(entry);
//        }
//        src_node.clones[src_loc.second].clear();
//    }
    if( src_node.clones_src[src_loc.second].size() > 0){
        for( auto& entry: src_node.clones_src[src_loc.second] ){
            ForceAssert( entry != src_loc);

            auto& clone_list = ca.G.VertMutable(entry.first).clones_tgt[entry.second];
            ForceAssert(clone_list.find(src_loc)!=clone_list.end());
            clone_list.erase(src_loc);

            clone_list.insert(tgt_loc);
            tgt_node.clones_src[tgt_loc.second].insert(entry);
        }
        src_node.clones_src[src_loc.second].clear();
    }
    if( src_node.clones_tgt[src_loc.second].size() > 0){
        for( auto& entry: src_node.clones_tgt[src_loc.second] ){
            ForceAssert( entry != src_loc);

            auto& clone_list = ca.G.VertMutable(entry.first).clones_src[entry.second];
            ForceAssert(clone_list.find(src_loc)!=clone_list.end());
            clone_list.erase(src_loc);

            clone_list.insert(tgt_loc);
            tgt_node.clones_tgt[tgt_loc.second].insert(entry);
        }
        src_node.clones_tgt[src_loc.second].clear();
    }
    src_loc = tgt_loc;
}

void CollectTargetBranches(vec<vec<std::pair<int,int>>>& branches, const connection_assembly&ca, int src, int src_pos){
    ForceAssert(src_pos>=0);
    if( src_pos > 0){ // if there's a base location to the left
        std::pair<int,int> left_loc = std::make_pair(src,src_pos-1);
        std::pair<int,int> mapped_left = ca.G.Vert(left_loc.first).mapped_loc[left_loc.second];

        ForceAssert(  ca.G.Vert(left_loc.first).bases[left_loc.second]
                   == ca.G.Vert(mapped_left.first).bases[mapped_left.second]
                   );
        if( mapped_left == left_loc){
//            ForceAssert( ca.G.Vert(left_loc.first).clones[left_loc.second].size() == 0);
//            const std::set<std::pair<int,int>>& clones = ca.G.Vert(left_loc.first).clones[left_loc.second];
            const StdSet<std::pair<int,int>>& clones_src = ca.G.Vert(left_loc.first).clones_src[left_loc.second];
            const StdSet<std::pair<int,int>>& clones_tgt = ca.G.Vert(left_loc.first).clones_tgt[left_loc.second];

            vec<std::pair<int,int>> loc_branches;
            loc_branches.push_back( std::make_pair( left_loc.first, left_loc.second+1 ));
//if( clones.size() >0 ){
//std::cout << left_loc.first <<" " << left_loc.second << " Lclones: ";
//}
//            for(auto entry: clones){
//                if(entry!=left_loc){
//                    loc_branches.push_back( std::make_pair( entry.first, entry.second+1 ));
//                }
//std::cout << "(" << entry.first<<"," << entry.second<<")";
//            }
//if( clones.size() >0 ){
//std::cout << std::endl;
//}

if( clones_src.size() >0 ){
std::cout << left_loc.first <<" " << left_loc.second << " Lclones: ";
}
            for(auto entry: clones_src){
                if(entry!=left_loc){
                    loc_branches.push_back( std::make_pair( entry.first, entry.second+1 ));
                }
std::cout << "(" << entry.first<<"," << entry.second<<")";
            }
if( clones_src.size() >0 ){
std::cout << std::endl;
}
if( clones_tgt.size() >0 ){
std::cout << left_loc.first <<" " << left_loc.second << " Lclones: ";
}
            for(auto entry: clones_tgt){
                if(entry!=left_loc){
                    loc_branches.push_back( std::make_pair( entry.first, entry.second+1 ));
                }
std::cout << "(" << entry.first<<"," << entry.second<<")";
            }
if( clones_tgt.size() >0 ){
std::cout << std::endl;
}

            bool bNew=true;
            for(auto entry: branches){
                if( entry[0] == loc_branches[0]){
                    bNew=false;
                    break;
                }
            }
            if(bNew){
                branches.push_back(loc_branches);
            }

//            branches.insert(std::make_pair(src,src_pos));
        }
        else {
            CollectTargetBranches(branches, ca, mapped_left.first, mapped_left.second+1);
        }
    }
    ForceAssert( src_pos <= ca.G.Vert(src).bases.isize());
    if( src_pos < ca.G.Vert(src).bases.isize()){ //if there's a base location to the right
        std::pair<int,int> right_loc = std::make_pair(src,src_pos);

        ForceAssert( src_pos <= ca.G.Vert(src).bases.isize());

        std::pair<int,int> mapped_right = ca.G.Vert(right_loc.first).mapped_loc[right_loc.second];
        ForceAssert(  ca.G.Vert(right_loc.first).bases[right_loc.second]
                   == ca.G.Vert(mapped_right.first).bases[mapped_right.second]
                   );
        if( mapped_right == right_loc){
//            const std::set<std::pair<int,int>>& clones = ca.G.Vert(right_loc.first).clones[right_loc.second];
            const StdSet<std::pair<int,int>>& clones_src = ca.G.Vert(right_loc.first).clones_src[right_loc.second];
            const StdSet<std::pair<int,int>>& clones_tgt = ca.G.Vert(right_loc.first).clones_tgt[right_loc.second];

            vec<std::pair<int,int>> loc_branches;
            loc_branches.push_back( std::make_pair( right_loc.first, right_loc.second ));
//if( clones.size() >0 ){
//std::cout << right_loc.first <<" " << right_loc.second << " Rclones: ";
//}
//            for(auto entry: clones){
//                if(entry!=right_loc){
//                    loc_branches.push_back( std::make_pair( entry.first, entry.second ));
//                }
//std::cout << "(" << entry.first<<"," << entry.second<<")";
//            }
//if( clones.size() >0 ){
//std::cout << std::endl;
//}
if( clones_src.size() >0 ){
std::cout << right_loc.first <<" " << right_loc.second << " Rclones: ";
}
            for(auto entry: clones_src){
                if(entry!=right_loc){
                    loc_branches.push_back( std::make_pair( entry.first, entry.second ));
                }
std::cout << "(" << entry.first<<"," << entry.second<<")";
            }
if( clones_src.size() >0 ){
std::cout << std::endl;
}
if( clones_tgt.size() >0 ){
std::cout << right_loc.first <<" " << right_loc.second << " Rclones: ";
}
            for(auto entry: clones_tgt){
                if(entry!=right_loc){
                    loc_branches.push_back( std::make_pair( entry.first, entry.second ));
                }
std::cout << "(" << entry.first<<"," << entry.second<<")";
            }
if( clones_tgt.size() >0 ){
std::cout << std::endl;
}
            bool bNew=true;
            for(auto entry: branches){
                if( entry[0] == loc_branches[0]){
                    bNew=false;
                    break;
                }
            }
            if(bNew){
                branches.push_back(loc_branches);
            }


//            branches.insert(std::make_pair(src,src_pos));
        }
        else {
            CollectTargetBranches(branches, ca, mapped_right.first, mapped_right.second);
        }
    }
}
void DetermineTargetVertices(std::set<int>& vertices, int& nVertices, connection_assembly& ca, branch_log& log, int src, int src_pos){
//std::cout << "Determining vertices for " << src << " " << src_pos << std::endl;
    ForceAssert(src_pos>=0);
    ForceAssert(ca.G.Vert(src).bases.isize() >= src_pos);
    vec<vec<std::pair<int,int>>> tgt_branches;
    CollectTargetBranches(tgt_branches, ca, src, src_pos);


    int nClonedVertices=0;

    std::set<std::pair<int,int>> stackable_branch;
    std::set<std::pair<int,int>> cloned_branch;;

    for( auto& entry: tgt_branches){
        for(size_t ii=1;ii<entry.size();++ii){
            cloned_branch.insert(entry[ii]);
        }
    }
    for( auto& entry: tgt_branches){
        if( cloned_branch.find(entry[0]) == cloned_branch.end()){
            stackable_branch.insert(entry[0]);
        }
    }

    vertices.clear();
    int lowest_rank_discovered = std::numeric_limits<int>::max();
    int lowest_rank_vertex_idx = -1;


    for( auto& entry: stackable_branch){
        int vertex_idx = ca.G.Vert(entry.first).vertices[entry.second] ;
        if( vertex_idx >=0 && lowest_rank_discovered > log.appearances[entry.first] ){
            lowest_rank_discovered = log.appearances[entry.first] ;
            lowest_rank_vertex_idx = vertex_idx;
        }
    }
    if( lowest_rank_vertex_idx < 0 ){
        lowest_rank_vertex_idx = nVertices;
        vertices.insert(lowest_rank_vertex_idx);
        ++nVertices;
    }
    for( auto& entry: stackable_branch){
        int& vertex_idx = ca.G.VertMutable(entry.first).vertices[entry.second] ;
        if( vertex_idx <0) vertex_idx = lowest_rank_vertex_idx;

        //skip edge to local sequence
        if( entry.first != src || entry.second==src_pos){
            vertices.insert(vertex_idx);
        }

    }



    for( auto& entry: cloned_branch){
        int& vertex_idx = ca.G.VertMutable(entry.first).vertices[entry.second] ;
        if( vertex_idx <0){
          vertex_idx = nVertices;
          ++nVertices;
        }
        if( entry.first != src || entry.second==src_pos){
            vertices.insert(vertex_idx);
        }
    }


}
void CollectSimpleTandems(vec<int>& log, int src, const digraphE<BaseVec2KmerPath>& sequence_graph){
    if(   log.size()!=0
       && (sequence_graph.ToSize(src) != 2 || sequence_graph.To(src)[0]!=sequence_graph.To(src)[1])
      ){
        return;
    }
    if(   sequence_graph.FromSize(src) == 2
       && sequence_graph.From(src)[0]==sequence_graph.From(src)[1]
       && (sequence_graph.EdgeObjectByIndexFrom(src,0).size()==0 || sequence_graph.EdgeObjectByIndexFrom(src,1).size()==0)
       && (sequence_graph.EdgeObjectByIndexFrom(src,0).size() !=  sequence_graph.EdgeObjectByIndexFrom(src,1).size())
      ){
        log.push_back(src);
        int des = sequence_graph.From(src)[0];
        if(   sequence_graph.FromSize(des) == 2
           && sequence_graph.From(des)[0]==sequence_graph.From(des)[1]
           && (sequence_graph.EdgeObjectByIndexFrom(des,0).size()==0 || sequence_graph.EdgeObjectByIndexFrom(des,1).size()==0)
           && (sequence_graph.EdgeObjectByIndexFrom(des,0).size() !=  sequence_graph.EdgeObjectByIndexFrom(des,1).size())
          ){
            const basevector& org_seq = ( sequence_graph.EdgeObjectByIndexFrom(src,0).size()==0)?(sequence_graph.EdgeObjectByIndexFrom(src,1))
                                                                                                :(sequence_graph.EdgeObjectByIndexFrom(src,0));
            const basevector& des_seq = ( sequence_graph.EdgeObjectByIndexFrom(des,0).size()==0)?(sequence_graph.EdgeObjectByIndexFrom(des,1))
                                                                                                :(sequence_graph.EdgeObjectByIndexFrom(des,0));
            if(des_seq==org_seq){
                CollectSimpleTandems(log,des,sequence_graph);
            }
        }
    }
}

void SimplifySequenceGraph(digraphE<BaseVec2KmerPath>& sequence_graph){
    sequence_graph.RemoveDuplicateEdges();
    sequence_graph.RemoveUnneededVertices();
    sequence_graph.RemoveDeadEdgeObjects();
    sequence_graph.RemoveEdgelessVertices();
    bool bModified=false;
    do{
        bModified=false;
        for(int src = 0 ; src<sequence_graph.N() && !bModified; ++src){
            bool bIsRoot=true;
            for(int conn_idx = 0 ; conn_idx < sequence_graph.ToSize(src) ; ++ conn_idx){
                if(sequence_graph.EdgeObjectByIndexTo(src,conn_idx).size()==0){
                    bIsRoot=false;
                    break;
                }
            }
            if(!bIsRoot) continue;

            vec<int> log;
            CollectSimpleTandems(log,src,sequence_graph);
            size_t nSimpleTandems = log.size();
            if( nSimpleTandems > 0){
                ForceAssert(log[0]==src);
                int des = sequence_graph.From( log.back())[0];
                if( nSimpleTandems > 1){
                    for(size_t vv=1 ; vv< log.size() ; ++vv){
                        sequence_graph.DeleteEdgesAtVertex(log[vv]);
                    }
                    basevector bCurr;
                    const basevector& org_seq = ( sequence_graph.EdgeObjectByIndexFrom(src,0).size()==0)?(sequence_graph.EdgeObjectByIndexFrom(src,1))
                                                                                                        :(sequence_graph.EdgeObjectByIndexFrom(src,0));
                    for(size_t rr=0 ; rr <= nSimpleTandems; ++rr){
                        sequence_graph.AddEdge(src,des,bCurr);
                        bCurr.append(org_seq);

                    }
                    sequence_graph.RemoveDeadEdgeObjects();
                    sequence_graph.RemoveEdgelessVertices();
                    bModified=true;
                }
            }
        }

    }while(bModified);
    /*
    do{
        bModified=false;
        for(int src = 0 ; src<sequence_graph.N() && !bModified; ++src){
            int des=-1;
            int conn_idx=0;
            for(; conn_idx < sequence_graph.FromSize(src) ; ++conn_idx){
                if( sequence_graph.EdgeObjectByIndexFrom(src,conn_idx).size()==0){
                    des=sequence_graph.From(src)[conn_idx];
                    break;
                }
            }
            if(des>=0){
                ForceAssert(sequence_graph.EdgeObjectByIndexFrom(src,conn_idx).size()==0);
                ForceAssert(sequence_graph.FromSize(src)>1);

                sequence_graph.DeleteEdgeFrom(src,conn_idx);
                for(int conn_idx=0; conn_idx < sequence_graph.ToSize(src) ; ++conn_idx){

                    int srcsrc=sequence_graph.To(src)[conn_idx];
                    auto new_edge=sequence_graph.EdgeObjectByIndexTo(src,conn_idx);
                    sequence_graph.AddEdge(srcsrc,des,new_edge);
                }
                bModified=true;
            }
        }

    }while(bModified);
    */
    sequence_graph.RemoveDuplicateEdges();
    sequence_graph.RemoveUnneededVertices();
    sequence_graph.RemoveDeadEdgeObjects();
    sequence_graph.RemoveEdgelessVertices();
}
int overlapLeft(const basevector&left, const basevector& right){
    int n = std::min(left.isize(),right.isize());
    int out = 0;
    for(;out<n;++out){
        if(left[out]!=right[out]){
            break;
        }
    }
    return out;
}
int overlapRight(const basevector&left, const basevector& right){
    int n = std::min(left.isize(),right.isize());
    if(n==0)return 0;

    int out = 0;
    int lshift=left.isize()-1;
    int rshift=right.isize()-1;
    for(;out<n;++out){
        if(left[lshift-out]!=right[rshift-out]){
            break;
        }
    }
    return out;
}

std::pair<int,int> DetermineVertexK(int vv, digraphE<BaseVec2KmerPath>&sequence_graph){
  std::pair<int,int> out = std::make_pair(std::numeric_limits<int>::max(),std::numeric_limits<int>::max());
  const int nTo = sequence_graph.ToSize(vv);
  for(int ii=0;ii<nTo;++ii){
      const basevector& ii_bases = sequence_graph.EdgeObjectByIndexTo(vv,ii);

      out.first = std::min(out.first,ii_bases.isize());
      for(int jj=ii+1;jj<nTo;++jj){
          const basevector& jj_bases = sequence_graph.EdgeObjectByIndexTo(vv,jj);
          out.first = std::min(out.first,overlapRight(ii_bases,jj_bases));
      }
  }

  const int nFrom = sequence_graph.FromSize(vv);
  for(int ii=0;ii<nFrom;++ii){
      const basevector& ii_bases = sequence_graph.EdgeObjectByIndexFrom(vv,ii);

      out.second = std::min(out.second,ii_bases.isize());
      for(int jj=ii+1;jj<nFrom;++jj){
          const basevector& jj_bases = sequence_graph.EdgeObjectByIndexFrom(vv,jj);
          out.second = std::min(out.second,overlapLeft(ii_bases,jj_bases));
      }
  }
  out.first += 1;
  out.second += 1;
  if( nTo==0){
      out.first=-1;
  }
  if( nFrom==0){
      out.second=-1;
  }
  return out;
}

void HyperKmerize(digraphE<BaseVec2KmerPath>& sequence_graph,size_t target_k){
    for(int vv = 0 ; vv < sequence_graph.N() ; ++vv){
        std::pair<int,int> ks = DetermineVertexK(vv,sequence_graph);
        std::cout << vv << " " << ks.first << " " << ks.second << std::endl;
    }
}
