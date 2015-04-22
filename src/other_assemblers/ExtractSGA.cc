///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Extract a fermi assembly.

#include <map>

#include "FastIfstream.h"
#include "Fastavector.h"
#include "CoreTools.h"
#include "other_assemblers/ConnectionAssembly.h"
#include "other_assemblers/ExtractSGA.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"

template connection_assembly::node_t& digraphVE<connection_assembly::node_t,connection>::VertMutable(int);

void ExtractSGA( const String& IN_DIR, connection_assembly& A, bool bVariants )
{
     String sDelim=" \t";
     // Check for file.
std::cout << ">>>ExtractSGA" << std::endl;
     if ( !IsDirectory(IN_DIR) )
     {    FatalErr("Can't find input directory.");    }
     String fn = IN_DIR + "/primary-graph.asqg.gz";
     if ( !IsRegularFile(fn) )
     {    FatalErr("Failed to find " << fn << ".");    }

     //should be using unordered_map, but String has no stock hash function
     std::map<String,std::pair<size_t,basevector>> vertices;
     size_t nv=0;

     // This is for dealing with SGA graph's ED line
     enum {C1,C2,S1,E1,L1,S2,E2,L2,RC,DIFF};
     typedef std::tuple<String,String,int64_t,int64_t,int64_t,int64_t,int64_t,int64_t,bool,int64_t> edge_log_t;
     vec< edge_log_t > edge_log;
     {
         size_t nHeaderLines=0;
         String sBuffer,sLocBuffer;
         fast_pipe_ifstream in( "zcat " + fn );
         for(getline(in,sBuffer);!in.fail();getline(in,sBuffer)){
             // 1st field - entry type
             auto start=sBuffer.find_first_not_of(sDelim,0);
             if (start==String::npos) continue;
             auto end=sBuffer.find_first_of(sDelim,start);
             if (end==String::npos) continue;

             String sMarker=sBuffer.substr(start,end-start);

             if(sMarker=="VT") { // if it is a vertex  line
                 //2nd field - vertex name
                 start=sBuffer.find_first_not_of(sDelim,end); ForceAssert(start!=String::npos);
                 end=sBuffer.find_first_of(sDelim,start);     ForceAssert(end!=String::npos && end>start);
                 String sVertexName = sBuffer.substr(start,end-start);
                 ForceAssert(vertices.find(sVertexName)==vertices.end());

                 //3rd field - sequence
                 start=sBuffer.find_first_not_of(sDelim,end); ForceAssert(start!=String::npos);
                 end=sBuffer.find_first_of(sDelim,start);     if(end==String::npos) end=sBuffer.size(); ForceAssert(end>=start);

                 vertices[sVertexName]=std::pair<size_t,basevector>(nv,basevector(end-start));
std::cout << nv << " " << sVertexName << std::endl;
                 ++nv;
                 auto& bases = vertices[sVertexName].second;
                 for( size_t ii=start ; ii<end ; ++ii){
                     bases.Set(ii-start,as_char(sBuffer[ii]));
                 }
             }
             else if(sMarker=="ED") { //edge line
                 //2nd field - contig 1 name
                 start=sBuffer.find_first_not_of(sDelim,end); ForceAssert(start!=String::npos);
                 end=sBuffer.find_first_of(sDelim,start);     ForceAssert(end!=String::npos && end>start);
                 String sContig1 = sBuffer.substr(start,end-start);

                 //3rd field - contig 2 name
                 start=sBuffer.find_first_not_of(sDelim,end); ForceAssert(start!=String::npos);
                 end=sBuffer.find_first_of(sDelim,start);     ForceAssert(end!=String::npos && end>start);
                 String sContig2 = sBuffer.substr(start,end-start);

                 //4th field - contig 1 overlap starts
                 start=sBuffer.find_first_not_of(sDelim,end); ForceAssert(start!=String::npos);
                 end=sBuffer.find_first_of(sDelim,start);     ForceAssert(end!=String::npos && end>start);
                 sLocBuffer = sBuffer.substr(start,end-start);ForceAssert(sLocBuffer.IsInt());
                 int64_t start1=sLocBuffer.Int();

                 //5th field - contig 1 overlap ends
                 start=sBuffer.find_first_not_of(sDelim,end); ForceAssert(start!=String::npos);
                 end=sBuffer.find_first_of(sDelim,start);     ForceAssert(end!=String::npos && end>start);
                 sLocBuffer = sBuffer.substr(start,end-start);ForceAssert(sLocBuffer.IsInt());
                 int64_t end1=sLocBuffer.Int()+1;

                 //6th field - contig 1 length
                 start=sBuffer.find_first_not_of(sDelim,end); ForceAssert(start!=String::npos);
                 end=sBuffer.find_first_of(sDelim,start);     ForceAssert(end!=String::npos && end>start);
                 sLocBuffer = sBuffer.substr(start,end-start);ForceAssert(sLocBuffer.IsInt());
                 int64_t length1=sLocBuffer.Int();

                 //7th field - contig 2 overlap starts
                 start=sBuffer.find_first_not_of(sDelim,end); ForceAssert(start!=String::npos);
                 end=sBuffer.find_first_of(sDelim,start);     ForceAssert(end!=String::npos && end>start);
                 sLocBuffer = sBuffer.substr(start,end-start);ForceAssert(sLocBuffer.IsInt());
                 int64_t start2=sLocBuffer.Int();

                 //8th field - contig 2 overlap ends
                 start=sBuffer.find_first_not_of(sDelim,end); ForceAssert(start!=String::npos);
                 end=sBuffer.find_first_of(sDelim,start);     ForceAssert(end!=String::npos && end>start);
                 sLocBuffer = sBuffer.substr(start,end-start);ForceAssert(sLocBuffer.IsInt());
                 int64_t end2=sLocBuffer.Int()+1;

                 //9th field - contig 2 length
                 start=sBuffer.find_first_not_of(sDelim,end); ForceAssert(start!=String::npos);
                 end=sBuffer.find_first_of(sDelim,start);     ForceAssert(end!=String::npos && end>start);
                 sLocBuffer = sBuffer.substr(start,end-start);ForceAssert(sLocBuffer.IsInt());
                 int64_t length2=sLocBuffer.Int();

                 //10th field - contig 2 is RC (1) or not (0)
                 start=sBuffer.find_first_not_of(sDelim,end); ForceAssert(start!=String::npos);
                 end=sBuffer.find_first_of(sDelim,start);     ForceAssert(end!=String::npos && end>start);
                 sLocBuffer = sBuffer.substr(start,end-start);ForceAssert(sLocBuffer.IsInt());
                 int64_t iRC=sLocBuffer.Int();
                 ForceAssert(iRC==0 || iRC==1);
                 bool bRC = iRC==1;

                 //11th field - number of differences in overlap (seems to be depricated and can take on negative values)
                 start=sBuffer.find_first_not_of(sDelim,end); ForceAssert(start!=String::npos);
                 end=sBuffer.find_first_of(sDelim,start);     if(end==String::npos) end=sBuffer.size(); ForceAssert(end>start);
                 sLocBuffer = sBuffer.substr(start,end-start);ForceAssert(sLocBuffer.IsInt());
                 int64_t iDiff=sLocBuffer.Int();



                 ForceAssert(end2-start2==end1-start1);
                 edge_log.push_back(edge_log_t(sContig1,sContig2,start1,end1,length1,start2,end2,length2,bRC,iDiff));
             }
             else if (sMarker=="HT"){ //header line
                 std::cout << "Header line: " << sBuffer << std::endl;
                 ++nHeaderLines;
             }
             else{
                 std::cerr << "Warning: unrecognized line marker: " << sBuffer << std::endl;
             }

         }
         if(nHeaderLines!=1){
           std::cerr << "Warning: "<< nHeaderLines << " header lines detected." << std::endl;
         }
     }
     ForceAssert(nv==vertices.size());
     vecbasevector bases(nv);
     for( auto& entry: vertices){
         ForceAssert(bases[entry.second.first].size()==0);
         bases[entry.second.first] = entry.second.second;
     }


     //parse variants

     std::map<String,fastavector> variants_dict;

     std::map<String,std::tuple<basevector,int,int>> base_sequence;

     std::map<int,vec<std::tuple<int,int,int>>> vertex_variants_list;
     if(bVariants)
     {
         vec<fastavector> variants_list;
         vec<String> variants_tags;
         String variant_fasta="primary-variants.fa";
         String variant_query="variants_contigs.query";

         LoadFromFastaFile(variant_fasta,variants_list,variants_tags);
         for(size_t ii=0;ii<variants_list.size();++ii){
             int end=0;
             int start=variants_tags[ii].find_first_not_of(sDelim,end);
             end = variants_tags[ii].find_first_of(sDelim);
             variants_tags[ii]=variants_tags[ii].substr(start,end-start);
//             std::cout << variants_tags[ii] << std::endl;
             ForceAssert( variants_dict.find(variants_tags[ii])==variants_dict.end());
             variants_dict[variants_tags[ii]]=variants_list[ii];
         }

         string sBuffer,sLocBuffer;

         vec<String> row;

         std::map<std::string,std::tuple<int,int,int>> var2vertex;

         enum{VAR,ALT,ALT_START,ALT_END,ALT_L,ALT_RC,CONTIG,CONTIG_START,CONTIG_END,CONTIG_L};
         {
             ifstream in(variant_query);
             for(getline(in,sBuffer);!in.fail();getline(in,sBuffer)){
                 if(sBuffer[0]!='#'){
                     string sLocBuffer;
                     istringstream stream(sBuffer);
                     row.clear();
                     while( getline(stream, sLocBuffer, '\t') ){
                         row.push_back(sLocBuffer);
                     }

                     int iAlt = atoi(row[ALT].c_str());
                     if( iAlt == 0){
                         ForceAssert(vertices.find(row[CONTIG])!=vertices.end());
                         int orgVertex = vertices[row[CONTIG]].first;

                         ForceAssert(var2vertex.find(row[VAR])==var2vertex.end());
                         var2vertex[row[VAR]]=make_tuple(orgVertex
                                                        ,atoi(row[CONTIG_START].c_str())
                                                        ,atoi(row[CONTIG_END].c_str())-1
                                                        );
                         basevector tmp = variants_dict[row[VAR]+"/"+row[ALT]].ToBasevector();
                         if( atoi(row[ALT_RC].c_str())!=0){
                             tmp.ReverseComplement();
                         }
                         ForceAssert(base_sequence.find(row[VAR])==base_sequence.end());
                         base_sequence[row[VAR]] = std::make_tuple(tmp,std::numeric_limits<int>::max(),std::numeric_limits<int>::max());
                     }
                     else{
                         ForceAssert(var2vertex.find(row[VAR])!=var2vertex.end());

                         ForceAssert(variants_dict.find(row[VAR]+"/"+row[ALT])!=variants_dict.end());

                         ForceAssert(atoi(row[ALT_START].c_str())==0);
                         ForceAssert(atoi(row[ALT_END].c_str())==atoi(row[ALT_L].c_str()));

                         basevector tmp = variants_dict[row[VAR]+"/"+row[ALT]].ToBasevector();
                         if( atoi(row[ALT_RC].c_str())!=0){
                             tmp.ReverseComplement();
                         }


                         ForceAssert(base_sequence.find(row[VAR])!=base_sequence.end());
                         const auto& org_sequence=std::get<0>(base_sequence[row[VAR]]);

                         int firstDiff=-1,lastDiff=-1;

                         int nmin = std::min( org_sequence.size(),tmp.size());

                         for( int ii=0;ii<nmin;++ii){
                             if( org_sequence[ii] != tmp[ii]){
                                 firstDiff=ii;
                                 break;
                             }
                         }
                         ForceAssert( firstDiff>0 && firstDiff < nmin-1);
                         for( int ii=0;ii<nmin;++ii){
                             if( org_sequence[org_sequence.isize()-1-ii] != tmp[tmp.isize()-1-ii]){
                                 lastDiff=ii;
                                 break;
                             }
                         }
                         ForceAssert( lastDiff>0 && lastDiff < nmin-1);

                         int front = std::min(firstDiff-1,tmp.isize()-1-lastDiff)-1;
                         int back = std::max(firstDiff-1,tmp.isize()-1-lastDiff)+1;

                         int org_back_shift = tmp.isize()-1-back;

                         std::get<1>(base_sequence[row[VAR]]) = std::min( std::get<1>(base_sequence[row[VAR]]) , front    );
                         std::get<2>(base_sequence[row[VAR]]) = std::min( std::get<2>(base_sequence[row[VAR]]) , org_back_shift    );

                         ForceAssert(front<back);

#if 0
                         for( int ii = 0 ; ii <= front ;++ii){
                             ForceAssert( tmp[ii] == org_sequence[ii]);
                         }

                         for( int ii = 0 ; ii <= org_back_shift ;++ii){
                             ForceAssert( tmp[tmp.isize()-1-ii] == org_sequence[org_sequence.isize()-1-ii]);
                         }

                         ForceAssert( tmp[back] == org_sequence[org_sequence.size()-1 - org_back_shift ]);

                         tmp=basevector(tmp,front,back-front+1);

                         ForceAssert( tmp[0] == org_sequence[front]);
                         ForceAssert( tmp[tmp.isize()-1] == org_sequence[org_sequence.size()-1 - org_back_shift ]);


                         int new_idx = bases.size();
    std::cout<< new_idx << " " << row[VAR]<<"/"<<row[ALT] << std::endl;
                         bases.push_back(tmp);
                         ++nv;

                         int org_idx = std::get<0>(var2vertex[row[VAR]]);
                         int org_front = std::get<1>(var2vertex[row[VAR]])+front;
                         int org_back = std::get<2>(var2vertex[row[VAR]]) - org_back_shift;
                         vertex_variants_list[org_idx].push_back(make_tuple(new_idx,org_front,org_back));
#endif
                     }
                 }
             }
        }
         {
             ifstream in(variant_query);
             for(getline(in,sBuffer);!in.fail();getline(in,sBuffer)){
                 if(sBuffer[0]!='#'){
                     string sLocBuffer;
                     istringstream stream(sBuffer);
                     row.clear();
                     while( getline(stream, sLocBuffer, '\t') ){
                         row.push_back(sLocBuffer);
                     }

                     int iAlt = atoi(row[ALT].c_str());
                     if( iAlt == 0){
                         ForceAssert(vertices.find(row[CONTIG])!=vertices.end());

                         ForceAssert(var2vertex.find(row[VAR])!=var2vertex.end());
                         ForceAssert(base_sequence.find(row[VAR])!=base_sequence.end());
                     }
                     else{
                         ForceAssert(var2vertex.find(row[VAR])!=var2vertex.end());

                         ForceAssert(variants_dict.find(row[VAR]+"/"+row[ALT])!=variants_dict.end());

                         ForceAssert(atoi(row[ALT_START].c_str())==0);
                         ForceAssert(atoi(row[ALT_END].c_str())==atoi(row[ALT_L].c_str()));

                         basevector tmp = variants_dict[row[VAR]+"/"+row[ALT]].ToBasevector();
                         if( atoi(row[ALT_RC].c_str())!=0){
                             tmp.ReverseComplement();
                         }


                         ForceAssert(base_sequence.find(row[VAR])!=base_sequence.end());
                         const auto& org_sequence=std::get<0>(base_sequence[row[VAR]]);



                         int org_back_shift = std::get<2>(base_sequence[row[VAR]]);

                         int front = std::get<1>(base_sequence[row[VAR]]);
                         int back  = tmp.isize()-1-std::get<2>(base_sequence[row[VAR]]);

//                         std::cout << std::get<1>(base_sequence[row[VAR]]) << " " << std::get<2>(base_sequence[row[VAR]]) << std::endl;
//                         std::cout << front << " " << back << std::endl;

                         ForceAssert(front<back);


                         for( int ii = 0 ; ii <= front ;++ii){
                             ForceAssert( tmp[ii] == org_sequence[ii]);
                         }

                         for( int ii = 0 ; ii <= org_back_shift ;++ii){
                             ForceAssert( tmp[tmp.isize()-1-ii] == org_sequence[org_sequence.isize()-1-ii]);
                         }

                         ForceAssert( tmp[back] == org_sequence[org_sequence.size()-1 - org_back_shift ]);

                         tmp=basevector(tmp,front,back-front+1);

                         ForceAssert( tmp[0] == org_sequence[front]);
                         ForceAssert( tmp[tmp.isize()-1] == org_sequence[org_sequence.size()-1 - org_back_shift ]);


                         int new_idx = bases.size();
    std::cout<< new_idx << " " << row[VAR]<<"/"<<row[ALT] << std::endl;
                         bases.push_back(tmp);
                         ++nv;

                         int org_idx = std::get<0>(var2vertex[row[VAR]]);
                         int org_front = std::get<1>(var2vertex[row[VAR]])+front;
                         int org_back = std::get<2>(var2vertex[row[VAR]]) - org_back_shift;
                         vertex_variants_list[org_idx].push_back(make_tuple(new_idx,org_front,org_back));
                     }
                 }
             }
        }
    }


     // Build graph, this is a modified version of the last part of ExtractFermi
     std::map< std::pair<int,int> , int> edge_count; // in case there are multiple edges between two nodes


     vec< pair< pair<int,int>, std::tuple<int,int,int,int> > > X;

     size_t nNonEndToEnd=0;

     for( auto& entry: edge_log){
         auto itr = vertices.find(std::get<C1>(entry)); ForceAssert(itr!=vertices.end());
         auto id1 = (*itr).second.first;

         itr = vertices.find(std::get<C2>(entry)); ForceAssert(itr!=vertices.end());
         auto id2 = (*itr).second.first;

         auto overlap = std::get<E1>(entry)-std::get<S1>(entry);

         if (  std::get<E1>(entry) != std::get<L1>(entry)
             &&std::get<S1>(entry) != 0
            ){
             std::cout<< "WARNING: In " << std::get<C1>(entry) << " to " << std::get<C2>(entry) << " the 1st vertex is not an end-to-end overlap." << std::endl;
             ++nNonEndToEnd;
         }
         ForceAssert (  std::get<E1>(entry) == std::get<L1>(entry) || std::get<S1>(entry) == 0);

         if (  std::get<E2>(entry) != std::get<L2>(entry)
             &&std::get<S2>(entry) != 0
            ){
             std::cout<< "WARNING: In " << std::get<C1>(entry) << " to " << std::get<C2>(entry) << " the 2nd vertex is not an end-to-end overlap." << std::endl;
             ++nNonEndToEnd;
         }
         ForceAssert (  std::get<E2>(entry) == std::get<L2>(entry) || std::get<S2>(entry) == 0);

         ForceAssert (  std::get<DIFF>(entry) <= 0);

         basevector b1 = bases[id1], b2 = bases[id2];

         ForceAssert(b1.size()==std::get<L1>(entry));
         ForceAssert(b2.size()==std::get<L2>(entry));
         int& loc_edge_count = edge_count[std::minmax(id1,id2)];

         if( std::get<RC>(entry) ){
             b2.ReverseComplement();
             for(int64_t ii=0;ii<overlap;++ii){ ForceAssert(b1[std::get<S1>(entry)+ii] == b2[b2.isize()-std::get<E2>(entry)+ii]); }

             bool b1Left = std::get<E1>(entry) == std::get<L1>(entry);
             bool b2Left = std::get<S2>(entry) == 0;
             ForceAssert(b1Left!=b2Left);


             if(b1Left){
                 X.push( make_pair(id1,id2), make_tuple(0,1,overlap,loc_edge_count) );
                 X.push( make_pair(id2,id1), make_tuple(0,1,overlap,loc_edge_count) );
             }
             else{
                 X.push( make_pair(id2,id1), make_tuple(1,0,overlap,loc_edge_count) );
                 X.push( make_pair(id1,id2), make_tuple(1,0,overlap,loc_edge_count) );
             }
         }
         else{
             basevector b1 = bases[id1], b2 = bases[id2];
             for(int64_t ii=0;ii<overlap;++ii){ ForceAssert(b1[std::get<S1>(entry)+ii] == b2[std::get<S2>(entry)+ii]); }

             bool b1Left = std::get<E1>(entry) == std::get<L1>(entry);
             bool b2Left = std::get<E2>(entry) == std::get<L2>(entry);
             ForceAssert(b1Left!=b2Left);
             if(b1Left){
                 X.push( make_pair(id1,id2), make_tuple(0,0,overlap,loc_edge_count) );
                 X.push( make_pair(id2,id1), make_tuple(1,1,overlap,loc_edge_count) );
             }
             else{
                 X.push( make_pair(id2,id1), make_tuple(0,0,overlap,loc_edge_count) );
                 X.push( make_pair(id1,id2), make_tuple(1,1,overlap,loc_edge_count) );
             }

         }
         ++loc_edge_count;

     }


     vec< vec<int> > from(nv), to(nv), from_edge_obj(nv), to_edge_obj(nv);
     vec<connection> edges;
     // this is a similar version of the last portion of ExtractFermi, note that consistency checks are done in the previous section
     int min_overlap = 1000000000;
     UniqueSort(X);
     for ( int i = 0; i < X.isize( ); i++ ){
          int id1 = X[i].first.first, id2 = X[i].first.second;
          int pos1 = std::get<0>(X[i].second), pos2 = std::get<1>(X[i].second);
          int over = std::get<2>(X[i].second);
          int loc_edge_count = std::get<3>(X[i].second);
          min_overlap = Min( min_overlap, over );

          basevector b1 = bases[id1], b2 = bases[id2];

          Bool rc1 = pos1 ,rc2 =  pos2 ;

          from[id1].push_back(id2), to[id2].push_back(id1);
          from_edge_obj[id1].push_back( edges.size( ) );
          to_edge_obj[id2].push_back( edges.size( ) );
          edges.push( rc1, rc2, over, loc_edge_count );
     }
     A.bases = bases;
     vec<connection_assembly::node_t> graph_vertices(bases.size());
     for(size_t ii=0;ii<bases.size();++ii){
         graph_vertices[ii]=connection_assembly::node_t(bases[ii]);
     }
     A.G.Initialize( from, to, graph_vertices, edges, to_edge_obj, from_edge_obj );

     //put in variants

     for( const auto& entry: vertex_variants_list){
         int root = entry.first;
         connection_assembly::node_t& root_node = A.G.VertMutable(root);

         for(auto tgt_entry: entry.second){
             int tgt_idx = std::get<0>(tgt_entry);
             int front = std::get<1>(tgt_entry);
             int back = std::get<2>(tgt_entry);
             root_node.variants.push_back( make_tuple(tgt_idx,front,back) );
             connection_assembly::node_t& tgt_node = A.G.VertMutable(tgt_idx);
             tgt_node.variant_of = std::make_tuple(root,front,back);
         }
     }
std::cout << "Number of non end-to-end / edges " << nNonEndToEnd << "/" << edge_log.size()*2 << std::endl;
std::cout << "<<<ExtractSGA" << std::endl;
}
