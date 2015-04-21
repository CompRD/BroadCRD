///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "STLExtensions.h"
#include "TokenizeString.h"
#include "PrintAlignment.h"
#include "FetchReads.h"
#include "paths/long/fosmid/Fosmids.h"
#include <string>

namespace{


class fosmid_alignment_t{
public:
    fosmid_alignment_t(const vecbasevector& refs, const String& sDir, int id):m_references(refs){ load(sDir,id); }
    void load(const String& sDir, int id){
        String ID = ToString(id),gid,loc;
        GetRegionInfo(ID,m_g,m_start,m_stop,gid,loc);

        BinaryReader::readFile( sDir + "/" + ID + "/region.align" , &m_alignment);
        const String fos_fasta = "/wga/dev/references/Homo_sapiens/NA12878_Fosmid_Pool.regions.fin/fos."+ID+".fasta";
        vecbasevector buffer;
        FetchReads(buffer,0,fos_fasta);
        m_query=buffer[0];
    };
    void PrintVisualAlignment(ostream&os)const{ ::PrintVisualAlignment( True, os, m_query, m_references[m_g], m_alignment ); }
    void PrintVisualAlignment(ostream&os,int min_dist_from_indel,int max_dist_from_indel, int flank)const;
private:
    fosmid_alignment_t();
    align m_alignment;
    const vecbasevector& m_references;
    basevector m_query;
    int m_g, m_start, m_stop;
};

class fosmid_alignments_collection{
public:
    typedef std::pair<int,fosmid_alignment_t> fid_alignment_t;
    fosmid_alignments_collection(const String& sDir)
        :m_references("/wga/scr4/bigrefs/human19/genome.fastb"){
        auto fids = AllFosmids();
        m_fid_alignment.reserve(fids.size());
        for(const auto&fid: fids){ m_fid_alignment.emplace_back(fid,fosmid_alignment_t(m_references,sDir,fid)); }
    }
    vec<fid_alignment_t>::const_iterator begin()const{return m_fid_alignment.cbegin();};
    vec<fid_alignment_t>::const_iterator end()const{return m_fid_alignment.cend();};
private:
    fosmid_alignments_collection();
    const vecbasevector m_references;

    vec<fid_alignment_t> m_fid_alignment;

};

void fosmid_alignment_t::PrintVisualAlignment(ostream& os,const int min_dist_from_indel, const int max_dist_from_indel,const int flank)const{
    struct print_coor_t{
        int L1b,L1e,L2b,L2e,M1b,M1e,M2b,M2e,R1b,R1e,R2b,R2e;
        bool bHasSub;
        print_coor_t():L1b(-1),L1e(-1),L2b(-1),L2e(-1)
                      ,M1b(-1),M1e(-1),M2b(-1),M2e(-1)
                      ,R1b(-1),R1e(-1),R2b(-1),R2e(-1),bHasSub(false){};
    };
    const basevector& ref=m_references[m_g];
    int p1=m_alignment.pos1();
    int p2=m_alignment.pos2();
    vec<print_coor_t> print_coors(m_alignment.Nblocks());
    for ( int j = 0; j < m_alignment.Nblocks( ); j++ )
    {   if ( m_alignment.Gaps(j) > 0 ) {
            print_coors[j].M2b=p2;
            p2 += m_alignment.Gaps(j);
            print_coors[j].M2e=p2;
        }
        if ( m_alignment.Gaps(j) < 0 ) {
            print_coors[j].M1b=p1;
            p1 -= m_alignment.Gaps(j);
            print_coors[j].M1e=p1;
        }
        int limit_R1e=p1+m_alignment.Lengths(j);
        int limit_R2e=p2+m_alignment.Lengths(j);
        int limit_L1b=p1;
        int limit_L2b=p2;

        print_coors[j].R1b = p1;
        print_coors[j].R1e = min(p1+flank,limit_R1e);
        print_coors[j].R2b = p2;
        print_coors[j].R2e = min(p2+flank,limit_R2e);


        if( j+1 < m_alignment.Nblocks()){
            print_coors[j+1].L1b = max(p1+m_alignment.Lengths(j)-flank,limit_L1b);
            print_coors[j+1].L1e = p1+m_alignment.Lengths(j);
            print_coors[j+1].L2b = max(p2+m_alignment.Lengths(j)-flank,limit_L2b);
            print_coors[j+1].L2e = p2+m_alignment.Lengths(j);
        }

        for ( int x = 0; x < m_alignment.Lengths(j); x++ ){
            if ( m_query[p1] != ref[p2] ) {
                if( x+1 > min_dist_from_indel && x < max_dist_from_indel ){
                    print_coors[j].R1e = max(min(p1+flank+1,limit_R1e),print_coors[j].R1e);
                    print_coors[j].R2e = max(min(p2+flank+1,limit_R2e),print_coors[j].R2e);
                    print_coors[j].bHasSub=true;
                }
                if( j+1<m_alignment.Nblocks() && x+max_dist_from_indel >= m_alignment.Lengths(j)
                   && min_dist_from_indel+x <= m_alignment.Lengths(j)
                  ){
                    print_coors[j+1].L1b = min(max(p1-flank,limit_L1b),print_coors[j+1].L1b);
                    print_coors[j+1].L2b = min(max(p2-flank,limit_L2b),print_coors[j+1].L2b);
                    print_coors[j+1].bHasSub=true;
                }
            }
            ++p1; ++p2;
        }

    }
    String top,middle,bottom;
    for(int bb=1 ; bb < m_alignment.Nblocks() ; ++bb){
        const auto& entry= print_coors[bb];
        if(entry.bHasSub){
            top.clear();
            middle.clear();
            bottom.clear();
            std::transform(m_query.begin()+entry.L1b,m_query.begin()+entry.L1e
                          ,ref.begin()+entry.L2b,std::back_inserter(top)
                          ,[](unsigned char a,unsigned char b){  return (a==b)?' ':'*';}
                          );
            std::transform(m_query.begin()+entry.L1b,m_query.begin()+entry.L1e,std::back_inserter(middle),BaseToCharMapper());
            std::transform(ref.begin()+entry.L2b,ref.begin()+entry.L2e,std::back_inserter(bottom),BaseToCharMapper());
            ForceAssert(top.size()==middle.size() && top.size()==bottom.size());
            ForceAssert( (entry.M1b==entry.M1e) !=  (entry.M2b==entry.M2e) );
            if(entry.M1b!=entry.M1e){
                top.append(entry.M1e-entry.M1b,'|');
                std::transform(m_query.begin()+entry.M1b,m_query.begin()+entry.M1e,std::back_inserter(middle),BaseToCharMapper());
                bottom.append(entry.M1e-entry.M1b,' ');
                ForceAssert(entry.R2b==entry.L2e);
                ForceAssert(entry.M1b==entry.L1e);
                ForceAssert(entry.M1e==entry.R1b);
            }
            else{
                top.append(entry.M2e-entry.M2b,'|');
                middle.append(entry.M2e-entry.M2b,' ');
                std::transform(ref.begin()+entry.M2b,ref.begin()+entry.M2e,std::back_inserter(bottom),BaseToCharMapper());
                ForceAssert(entry.R1b==entry.L1e);
                ForceAssert(entry.M2b==entry.L2e);
                ForceAssert(entry.M2e==entry.R2b);
            }
            ForceAssert(top.size()==middle.size() && top.size()==bottom.size());

            std::transform(m_query.begin()+entry.R1b,m_query.begin()+entry.R1e
                          ,ref.begin()+entry.R2b,std::back_inserter(top)
                          ,[](unsigned char a,unsigned char b){  return (a==b)?' ':'*';}
                          );
            std::transform(m_query.begin()+entry.R1b,m_query.begin()+entry.R1e,std::back_inserter(middle),BaseToCharMapper());
            std::transform(ref.begin()+entry.R2b,ref.begin()+entry.R2e,std::back_inserter(bottom),BaseToCharMapper());
            os << "query:  " << entry.L1b << "-" << entry.L1e << " gap "
                            << entry.R1b << "-" << entry.R1e << '\n';
            os << "target: " << entry.L2b << "-" << entry.L2e << " gap "
                             << entry.R2b << "-" << entry.R2e << '\n';
            ForceAssert(top.size()==middle.size() && top.size()==bottom.size());
            const size_t width=80;
            for(size_t cc=0 ; cc<top.size() ;cc+=width ){
                for(size_t xx=0;xx<width&&cc+xx<top.size();++xx){ os << top[cc+xx]; }
                os<<'\n';
                for(size_t xx=0;xx<width&&cc+xx<top.size();++xx){ os << middle[cc+xx]; }
                os<<'\n';
                for(size_t xx=0;xx<width&&cc+xx<top.size();++xx){ os << bottom[cc+xx]; }
                os<<"\n\n";
            }
        }
    }
};

}

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(BASE_DIR,"base directory of CompareVars data");
     CommandArgument_Int_OrDefault_Doc(MIN_DIST, 0, "minimum distance from indel");
     CommandArgument_Int_OrDefault_Doc(MAX_DIST, -1, "maximum distance from indel");
     CommandArgument_Int_OrDefault_Doc(FLANKS, 80, "print flanking region");
     EndCommandArguments;

     fosmid_alignments_collection collection(BASE_DIR);

     for(const auto& entry: collection){
         std::cout << "fid=" << entry.first<<std::endl;
         if(MAX_DIST==-1){
             entry.second.PrintVisualAlignment(std::cout);
         }
         else{
             entry.second.PrintVisualAlignment(std::cout,MIN_DIST,MAX_DIST,FLANKS);
         }
         std::cout << std::endl;
     }
}
