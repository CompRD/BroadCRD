///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include "system/Assert.h"
#include "paths/simulation/Regions.h"

namespace Regions{

bool operator<(const Region&left,const Region&right){
  return   (left.llStart<right.llStart)
        || (left.llStart==right.llStart&&left.llEnd<right.llEnd);
};


bool region_records::isInRecord (const std::string& sChromosome)const {
  return (dictionary.find(sChromosome)!=dictionary.end());
}

bool region_records::isInRange (const std::string& sChromosome
                               ,const ssize_t llRangeStart
                               ,const ssize_t llRangeEnd)const{

  const_iterator itr_chr = dictionary.find(sChromosome);
  if(itr_chr==dictionary.end()){
      return false;
  }

  bool bOut = false;

  for( mapped_type::const_iterator itr = itr_chr->second.begin()
                                  ,end = itr_chr->second.end()
     ; !bOut && itr != end && llRangeEnd >= itr->llStart
     ; ++itr){
      if      ( llRangeStart == itr->llStart ) bOut = true;
      else if ( llRangeStart >  itr->llStart ) bOut = llRangeStart < itr->llEnd;
      else                                     bOut = itr->llStart < llRangeEnd;
  }
  return bOut;
}

Region region_records::find(const std::string& sChromosome
                           ,const ssize_t llRangeStart
                           ,const ssize_t llRangeEnd)const{
    const_iterator itr_chr = dictionary.find(sChromosome);
    if(itr_chr==dictionary.end()){ return Region(-1,-2,""); }

    for( mapped_type::const_iterator itr = itr_chr->second.begin()
                                    ,end = itr_chr->second.end()
       ; itr != end && llRangeEnd >= itr->llStart
       ; ++itr){
        if ( llRangeStart == itr->llStart ) {
            return *itr;
        }
        else if ( llRangeStart >  itr->llStart ){
            if(llRangeStart < itr->llEnd){
                return *itr;
            }
        }
        else {
            if(itr->llStart < llRangeEnd){
                return *itr;
            }
        }
    }
    return Region(-1,-2,"");
}

void region_records::readRegions(const std::string&sFileName
                                ,const ssize_t pad_length
                                ,const size_t nSkippedColumns){
  std::ifstream infile(sFileName);
  ForceAssert( infile.is_open() );

  std::string sBuffer;
  while(std::getline(infile,sBuffer)){
      size_t info=addRecord(sBuffer,pad_length,nSkippedColumns);
      if(info){std::cerr<<"error adding record from line (error code "<<info<<")"<<std::endl;}
  }
  infile.close();

  // this is important
  compact();
};

int region_records::addRecord(const std::string& sLine
                             ,const ssize_t pad_length
                             ,const size_t nSkippedColumns
                             ,const std::string sDelim){
  if(sLine.size()==0) return 0;
  //ignore all space/tabs before any entry
  size_t uStart=sLine.find_first_not_of(sDelim);
  size_t uEnd=uStart;

  if(uStart>=sLine.size()) return 1;

  std::string sTag;

  //skip the first nSkippedColumns column
  for(size_t ii=0;ii<nSkippedColumns;++ii){
    uEnd=sLine.find_first_of(sDelim,uStart);
    if(uStart>=sLine.size()) return 2;

    if(ii>0) sTag.append(" ");
    sTag.append(sLine,uStart,uEnd-uStart);

    uStart=sLine.find_first_not_of(sDelim,uEnd);
    if(uStart>=sLine.size()) return 3;
  }

  //read the header
  uEnd = sLine.find_first_of(":",uStart);
  std::vector<Region>& ref=dictionary[sLine.substr(uStart,uEnd-uStart)];

  //read the begining range
  uStart = uEnd+1;
  if(uStart>=sLine.size()) return 4;
  uEnd = sLine.find_first_of("-",uStart);
  ssize_t range_start = atoll(sLine.substr(uStart,uEnd-uStart).c_str())-pad_length;
  if(range_start<0) range_start=0;


  //read the end range
  uStart = uEnd+1;
  if(uStart>=sLine.size()) return 5;
  uEnd = sLine.find_first_of(sDelim,uStart);
  ssize_t range_end = atoll(sLine.substr(uStart,uEnd-uStart).c_str())+pad_length;
  if( range_start <= range_end)
  {
    ref.push_back(Region(range_start,range_end,sTag));
  }
  else{
    std::cerr << "WARNING: ignoring empty/reverse range: [" <<range_start<<","<<range_end<<")"<<std::endl;
  }
  return 0;
}

void region_records::compact(){
    for( auto itr=dictionary.begin();itr!=dictionary.end(); ++itr){ // for each chromosome
        mapped_type& ref = itr->second;
        std::sort(ref.begin(),ref.end(),Regions::operator<); //sort

        for( auto jj = ref.begin(); jj!=ref.end();++jj){//go through the ranges
            //merge the overlaping ranges
            auto range_end = jj+1;
            while( range_end != ref.end() && range_end->llStart <= jj->llEnd){
                jj->llEnd = std::max(jj->llEnd,range_end->llEnd);
                jj->sTag.append(",");
                jj->sTag.append(range_end->sTag);
                ++range_end;
            }
            if( range_end!=jj+1) ref.erase(jj+1,range_end);
        }
    }
};

}//namespace Regions
