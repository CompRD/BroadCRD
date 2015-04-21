///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * CompareSNPSets.cc
 *
 *  Created on: Feb 6, 2014
 *      Author: blau
 */

#include <map>
#include <unordered_set>

#include "MainTools.h"
#include "FastIfstream.h"
#include "Vec.h"
#include "ParallelVecUtilities.h"

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "STLExtensions.h"
#include "TokenizeString.h"
#include "util/TextTable.h"

typedef vec< pair< pair<String,int>, pair<String,String> > > dbsnp_t;


void LoadDBSNP(dbsnp_t& dbsnp){
    String line,chr,junk,ref,alt;
    int64_t pos;
    const String cachefile = "dbSNP.cache";
    if ( IsRegularFile( cachefile )){
         BinaryReader::readFile( cachefile , &dbsnp );
    }
    else{
         fast_ifstream din( "/humgen/gsa-hpprojects/GATK/bundle/current/b37/"
              "dbsnp_138.b37.excluding_sites_after_129.vcf" );
         vec<String> dlines;
         dlines.reserve(14000000);
         while(1)
         {    getline( din, line );
              if ( din.fail( ) ) break;
              if ( line.Contains( "#", 0 ) ) continue;
              if ( line[0] != 'X' && line[0] != 'Y'
                   && ( line[0] < '0' || line[0] > '9' ) )
              {     continue;     }
              dlines.push_back(line);    }
         dbsnp.reserve( dlines.size( ) );
         for ( int z = 0; z < dlines.isize( ); z++ )
         {    const String& line = dlines[z];
              istrstream iline( line.c_str( ) );
              iline >> chr >> pos >> junk >> ref >> alt;
              if ( ref.size( ) > 1 || alt.size( ) > 1 ) continue;
              dbsnp.push( make_pair(chr,pos), make_pair(ref,alt) );    }
         ParallelSort(dbsnp);
         BinaryWriter::writeFile( cachefile , dbsnp);}
};

typedef std::pair<int64_t,int64_t> interval_t;
typedef std::tuple<int64_t,String,String> sub_t;
namespace std
{
    template<>
    struct hash<sub_t>
    {
        size_t operator()(sub_t const& in)const{
            return std::hash<size_t>()(std::get<0>(in));
        }
    };

}
typedef std::unordered_set<sub_t> sub_set_t;

class vcf_coors_t{
public:
    vcf_coors_t(const String& file, int64_t indel_pad,const dbsnp_t& dbsnp, const double dSubMinQual, const double dIdsMinQual){LoadVCF(file,indel_pad,dbsnp,dSubMinQual,dIdsMinQual);};
    vcf_coors_t(){};
    void LoadVCF(const String& file, int64_t indel_pad, const dbsnp_t& dbsnp,const double dSubMinQual,const double dIdsMinQual);
    void CompressCoors();

    vec<String> ChromList()const{
        vec<String> out;
        for(const auto& entry: chrom_indel_intervals){
            out.push_back(entry.first);
        }
        for(const auto& entry: chrom_ts_sub){
            out.push_back(entry.first);
        }
        for(const auto& entry: chrom_tv_sub){
            out.push_back(entry.first);
        }
        UniqueSort(out);
        return out;
    }
    vec<interval_t>& getChromIndelIntervals(const String& c){  return chrom_indel_intervals[c];};
    sub_set_t& getChromTsPos(const String& c){return chrom_ts_sub[c];};
    sub_set_t& getChromTvPos(const String& c){return chrom_tv_sub[c];};

private:
    std::map<String,vec<interval_t>> chrom_indel_intervals;
    std::map<String,sub_set_t> chrom_ts_sub;
    std::map<String,sub_set_t> chrom_tv_sub;
};

bool bIntersect(int64_t pos, const vec<interval_t>& intervals){
    if(intervals.size() == 0) return false;
    bool bIntersect=false;
    auto lb = std::lower_bound(intervals.begin(),intervals.end(), make_pair(pos,pos));
    if( lb == intervals.end()){
        if(   pos >= intervals.back().first && pos <= intervals.back().second ){ return true;}
    }
    else{
        if(   pos >= (*lb).first && pos <= (*lb).second ){ return true;}

        if( std::distance(intervals.begin(),lb) >0){
            --lb;
            if(   pos >= (*lb).first && pos <= (*lb).second ){ return true;}
        }
    }
    auto ub = std::upper_bound(intervals.begin(),intervals.end(), make_pair(pos,pos));
    if( ub == intervals.end()){
        if(   pos >= intervals.back().first && pos <= intervals.back().second ){ return true;}
    }
    else{
        if(   pos >= (*ub).first && pos <= (*ub).second ){ return true;}
        if( std::distance(intervals.begin(),ub) >0){
            --ub;
            if(   pos >= (*ub).first && pos <= (*ub).second ){ return true;}
        }
    }
    return bIntersect;

}


std::pair<int64_t,int64_t> calcTsTv(const sub_set_t& ts, const sub_set_t& tv, const vec<interval_t>& banned_intervals){
    int64_t ns=0,nv=0;

    for(const auto& entry: ts){
        if(!bIntersect(std::get<0>(entry),banned_intervals)){
            ++ns;
        }
    }
    for(const auto& entry: tv){
        if(!bIntersect(std::get<0>(entry),banned_intervals)){
            ++nv;
        }
    }
    return make_pair(ns,nv);
}




void vcf_coors_t::LoadVCF(const String& file, int64_t indel_pad, const dbsnp_t& dbsnp
                         ,const double dSubMinQual,const double dIdsMinQual){
     fast_ifstream in(file);
     String line,junk,chr,ref,alts,gen,filter;
     vec<String> alt_list;
     int pos;
     double qual;
     vec<String> lines;
     std::cout << Date() << ": loading " << file << std::endl;
     std::cout << Date() << ": SubMinQual=" << dSubMinQual << " IdsMinQual=" << dIdsMinQual << std::endl;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          if ( line.Contains( "#", 0 ) ) continue;
          if ( line[0] != 'X' && line[0] != 'Y'
               && ( line[0] < '0' || line[0] > '9' ) )
          {     continue;     }
          istrstream iline( line.c_str( ) );
          iline >> chr >> pos >> junk >> ref >> alts >> qual >> filter
               >> junk >> junk >> gen;

          if(!filter.Contains("PASS")) continue;

          ForceAssert(alts.size()>0);
          ForceAssert(ref.size()>0);
          if(ref.size()==1 && ref[0]=='.') continue;
          if(alts.size()==1 && alts[0]=='.') continue;

          Tokenize(alts,',',alt_list);
          if(alts[0]==',' || alts.back()==','){
              alt_list.push_back("");
          }


          vec<int> alleles;
          if ( gen.Contains( ":" ) ) gen = gen.Before( ":" );

          if( gen[0] == '.' ) continue;
          gen.GlobalReplaceBy( "|", "/" );
          //copy-and-pasted from TiTv.cc
          ForceAssert(gen.Contains("/"));
          {
              int hi=gen.Before("/").Int();
              int hi2=gen.After("/").Int();
              alleles.push_back(min(hi,hi2));
              if(hi!=hi2) alleles.push_back(max(hi,hi2));
          }
          if ( alleles.empty( ) || alleles.size( ) > 2 ) continue;
          if ( alleles.solo( ) && alleles[0] == 0 ) continue;
          if ( alleles.size( ) > 2 && alleles[0] != 0 ) continue;
          //end of copy-and-pasted from TiTv.cc

          String alt;
          if(alleles.size()==1){
              ForceAssert(alleles.front()!=0);
              alt=alt_list[alleles.front()-1];
          }
          else if(alleles.size()==2){
              if( alleles.front() != 0 ) continue; //TiTv doesn't have this
              alt=alt_list[alleles[1]-1];
          }
          ForceAssert(alt.size()>0);

//          for(const auto& alt: alt_list){
          {
              ref.ToUpper();
              alt.ToUpper();
              if( ref.size() == alt.size() ){
                  if(    ref.size()==1
                      && qual >= dSubMinQual
                      && !BinMember(dbsnp, make_pair(make_pair(chr,pos),make_pair(ref,alt))  )){
                      if(  (ref=="A" && alt == "G")
                         ||(ref=="G" && alt == "A")
                         ||(ref=="C" && alt == "T")
                         ||(ref=="T" && alt == "C")
                        )
                      {
                          auto info = chrom_ts_sub[chr].insert(make_tuple(int64_t(pos),ref,alt));
                          if(!info.second){ std::cout << "WARNING: dupped ts " << line << std::endl; }
                      }
                      else{
                          auto info = chrom_tv_sub[chr].insert(make_tuple(int64_t(pos),ref,alt));
                          if(!info.second){ std::cout << "WARNING: dupped tv " << line << std::endl; }
                      }
                  }
              }
              else if(qual >= dIdsMinQual){
//the second line is used by TiTv
                  chrom_indel_intervals[chr].emplace_back( max(int64_t(0),pos - indel_pad) , pos+ref.size()-1 + indel_pad );
//                  chrom_indel_intervals[chr].emplace_back( max(int64_t(0),pos - indel_pad) , pos + indel_pad );
              }
          }
     }
     CompressCoors();
}

void CompressCoors(vec<interval_t>& in){
    std::sort(in.begin(),in.end());
    int64_t last_back = in.front().first - 1;
    int64_t last_front = last_back+1;

    vec<interval_t> new_intervals;
    new_intervals.reserve(in.size());

    for(const auto& interval: in){
        if( interval.first > last_back){
            if(last_front <= last_back){
                new_intervals.emplace_back(last_front,last_back);
            }
            last_front = interval.first;
            last_back  = interval.second;
        }
        last_back = max(last_back,interval.second);
    }
    if(last_front <= last_back){
        new_intervals.emplace_back(last_front,last_back);
    }
    std::sort(new_intervals.begin(),new_intervals.end());
    for(const auto& interval: new_intervals){
        ForceAssert( interval.first <= interval.second);
    }
    for(size_t vv=1; vv< new_intervals.size() ;++vv){
        ForceAssert( new_intervals[vv-1].second < new_intervals[vv].first);

    }
    using std::swap;
    swap(in,new_intervals);
};

void vcf_coors_t::CompressCoors(){
    for(auto& entry: chrom_indel_intervals){
        ::CompressCoors(entry.second);
    }
}




int main(int argc, char *argv[])
{
    RunTime( );
    BeginCommandArguments;

    CommandArgument_String_OrDefault_Doc(VCF_L,
         "/wga/scr4/human_assemblies/1/v9/v9_combined.filtered.vcf",
         "VCF file to use as L)");
    CommandArgument_Double_OrDefault_Doc(L_SUB_MIN_QUAL, std::numeric_limits<double>::lowest(), "min qual for sub for L");
    CommandArgument_Double_OrDefault_Doc(L_IDS_MIN_QUAL, std::numeric_limits<double>::lowest(), "min qual for ids for L");

    CommandArgument_String_OrDefault_Doc(VCF_R,
         "/wga/scr4/NA12878_calls/mem4/haplotype-caller-pcr-none/reverted.12.aligned.wholegenome.sorted.indel_cleaned_local.recal.unfiltered.recal_snp_recal_indel.vcf",
         "VCF file to use as R");
    CommandArgument_Double_OrDefault_Doc(R_SUB_MIN_QUAL, 180. , "min qual for sub for R");
    CommandArgument_Double_OrDefault_Doc(R_IDS_MIN_QUAL, 100. , "min qual for ids for R");

    CommandArgument_Int_OrDefault_Doc(INDEL_PADDING, 10, "pad indel region by this many bp");

    EndCommandArguments;

    dbsnp_t dbsnp;
    LoadDBSNP(dbsnp);

    std::cout << Date() << ": loading L" << std::endl;
    vcf_coors_t left(VCF_L,INDEL_PADDING ,dbsnp,L_SUB_MIN_QUAL,L_IDS_MIN_QUAL);

    std::cout << Date() << ": loading R" << std::endl;
    vcf_coors_t right(VCF_R,INDEL_PADDING,dbsnp,R_SUB_MIN_QUAL,R_IDS_MIN_QUAL);

    vec<String> chrom_list=left.ChromList();
    chrom_list.append(right.ChromList());
    UniqueSort(chrom_list);
    //TiTv doesn't have this filter
    chrom_list.erase( std::remove(chrom_list.begin(),chrom_list.end(), "Y"), chrom_list.end());
    std::cout << "chromosomes: " ;
    std::copy(chrom_list.begin(),chrom_list.end(),ostream_iterator<String>(std::cout," "));
    std::cout << std::endl;
////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        int64_t ns=0,nv=0;
        for(const auto& chr: chrom_list){
            auto tmp = calcTsTv(left.getChromTsPos(chr),left.getChromTvPos(chr)
                               ,vec<interval_t>());
            ns+=tmp.first;
            nv+=tmp.second;
        }
        std::cout << "novel tstv L all: " << ns << "/" << nv << "=" << double(ns)/double(nv) << std::endl;;
    }

    {
        int64_t ns=0,nv=0;
        for(const auto& chr: chrom_list){
            auto tmp = calcTsTv(right.getChromTsPos(chr),right.getChromTvPos(chr)
                               ,vec<interval_t>());
            ns+=tmp.first;
            nv+=tmp.second;
        }
        std::cout << "novel tstv R all: " << ns << "/" << nv << "=" << double(ns)/double(nv) << std::endl;;
    }

    {
        int64_t ns=0,nv=0;
        for(const auto& chr: chrom_list){
            auto tmp = calcTsTv(left.getChromTsPos(chr),left.getChromTvPos(chr)
                               ,left.getChromIndelIntervals(chr));
            ns+=tmp.first;
            nv+=tmp.second;
        }
        std::cout << "novel tstv L " << INDEL_PADDING << "bp from L indel: " << ns << "/" << nv << "=" << double(ns)/double(nv) << std::endl;;
    }

    {
        int64_t ns=0,nv=0;
        for(const auto& chr: chrom_list){
            auto tmp = calcTsTv(right.getChromTsPos(chr),right.getChromTvPos(chr)
                               ,right.getChromIndelIntervals(chr));
            ns+=tmp.first;
            nv+=tmp.second;
        }
        std::cout << "novel tstv R " << INDEL_PADDING << "bp from R indel: " << ns << "/" << nv << "=" << double(ns)/double(nv) << std::endl;;
    }
////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        int64_t ns=0,nv=0;
        for(const auto& chr: chrom_list){

            vec<interval_t> bad_intervals;
            bad_intervals.append(left.getChromIndelIntervals(chr));
            bad_intervals.append(right.getChromIndelIntervals(chr));
            CompressCoors(bad_intervals);

            auto tmp = calcTsTv(left.getChromTsPos(chr),left.getChromTvPos(chr)
                               ,bad_intervals);

            ns+=tmp.first;
            nv+=tmp.second;
        }
        std::cout << "novel tstv L "<< INDEL_PADDING<<"bp from L+R indel: " << ns << "/" << nv << "=" << double(ns)/double(nv) << std::endl;;
    }

    {
        int64_t ns=0,nv=0;
        for(const auto& chr: chrom_list){
            vec<interval_t> bad_intervals;
            bad_intervals.append(left.getChromIndelIntervals(chr));
            bad_intervals.append(right.getChromIndelIntervals(chr));
            CompressCoors(bad_intervals);
            auto tmp = calcTsTv(right.getChromTsPos(chr),right.getChromTvPos(chr)
                               ,bad_intervals);
            ns+=tmp.first;
            nv+=tmp.second;
        }
        std::cout << "novel tstv R " << INDEL_PADDING <<"bp from L+R indel: " << ns << "/" << nv << "=" << double(ns)/double(nv) << std::endl;;
    }
////////////////////////////////////////////////////////////////////////////////////////////////////
    {
        int64_t ns=0,nv=0;
        for(const auto& chr: chrom_list){
            const auto& ts_L = left.getChromTsPos(chr);
            const auto& tv_L = left.getChromTvPos(chr);
            const auto& ts_R = right.getChromTsPos(chr);
            const auto& tv_R = right.getChromTvPos(chr);

            sub_set_t ts;

            for( const auto& lpos: ts_L){
                if( ts_R.find(lpos) == ts_R.end() ){
                    ts.insert(lpos);
                }
            }

            sub_set_t tv;
            for( const auto& lpos: tv_L){
                if( tv_R.find(lpos) == tv_R.end() ){
                    tv.insert(lpos);
                }
            }
            vec<interval_t> bad_intervals;
            bad_intervals.append(left.getChromIndelIntervals(chr));
            bad_intervals.append(right.getChromIndelIntervals(chr));
            CompressCoors(bad_intervals);

            auto tmp = calcTsTv(ts,tv,bad_intervals);
            ns+=tmp.first;
            nv+=tmp.second;
        }
        std::cout << "novel tstv L-R " <<INDEL_PADDING<< "bp from L+R indel: " <<ns << "/" << nv << "=" <<  double(ns)/double(nv) << std::endl;;
    }

    {
        int64_t ns=0,nv=0;
        for(const auto& chr: chrom_list){
            const auto& ts_L = left.getChromTsPos(chr);
            const auto& tv_L = left.getChromTvPos(chr);
            const auto& ts_R = right.getChromTsPos(chr);
            const auto& tv_R = right.getChromTvPos(chr);

            sub_set_t ts;
            for( const auto& rpos: ts_R){
                if( ts_L.find(rpos) == ts_L.end() ){
                    ts.insert(rpos);
                }
            }

            sub_set_t tv;
            for( const auto& rpos: tv_R){
                if( tv_L.find(rpos) == tv_L.end() ){
                    tv.insert(rpos);
                }
            }
            vec<interval_t> bad_intervals;
            bad_intervals.append(left.getChromIndelIntervals(chr));
            bad_intervals.append(right.getChromIndelIntervals(chr));
            CompressCoors(bad_intervals);

            auto tmp = calcTsTv(ts,tv,bad_intervals);
            ns+=tmp.first;
            nv+=tmp.second;
        }
        std::cout << "novel tstv R-L "<<INDEL_PADDING<<"bp from L+R indel: " << ns << "/" << nv << "=" << double(ns)/double(nv) << std::endl;;
    }
////////////////////////////////////////////////////////////////////////////////////////////////////
    std::cout << Date() << ": done." << std::endl;
}
