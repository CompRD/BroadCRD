///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
#include "MainTools.h"

#include <algorithm>
#include <vector>
#include <fstream>
#include <string>
#include <unordered_map>
*/

#ifndef REGIONS_H
#define REGIONS_H

#include<string>
#include<unordered_map>
#include<unordered_map>
#include<vector>

#include <sys/types.h>

namespace Regions{
//specifier of a range
struct Region{
    //range specified by [llStart,llEnd), semi-open
    ssize_t llStart;
    ssize_t llEnd;
    std::string sTag;

    Region(ssize_t a, ssize_t b,const std::string&c):llStart(a),llEnd(b),sTag(c){};
    Region():llStart(0l),llEnd(0l){};

    bool isValid()const{return llStart<=llEnd;};
};
//comparison operator
bool operator<(const Region&left,const Region&right);


// a dictionary with string key and a vector of ranges 
// unordered_map doesn't have virtual destructor -unsafe to inherit
class region_records{
  public:
    typedef std::vector<Region> mapped_type;
    typedef std::string key_type;
    typedef typename std::unordered_map<key_type,mapped_type> dictionary_t;
    typedef typename dictionary_t::iterator iterator;
    typedef typename dictionary_t::const_iterator const_iterator;

    //constructor
    region_records(const std::string& sFileName
                  ,const ssize_t pad_length=0     // positive number to enlarge region, negative to shrink
                  ,const size_t nSkippedColumns=1){
      readRegions(sFileName,pad_length,nSkippedColumns);
    };

    // return a copy of the region if there's overlap, Region::llStart=-1 > Region::llEnd=-2 if no overlap
    Region find(const std::string& sChromosome
               ,const ssize_t uRangeStart
               ,const ssize_t uRangeEnd)const;

    bool isInRange(const std::string& sChromosome
                  ,const ssize_t uRangeStart
                  ,const ssize_t uRangeEnd)const;

    // is there a record of this chromosome?
    bool isInRecord(const std::string& sChromosome)const;

    // load regions of interest from a file
    void readRegions(const std::string&sFileName
                    ,const ssize_t pad_length=0     // positive number to enlarge region, negative to shrink
                    ,const size_t nSkippedColumns=1);


    // sort and compact the data structure
    void compact();

// unordered_map doesn't have virtual destructor -unsafe to inherit
    const_iterator begin()const {return dictionary.begin();};
    const_iterator end()const {return dictionary.end();};
    iterator begin(){return dictionary.begin();};
    iterator end(){return dictionary.end();};
    size_t size()const {return dictionary.size();};

//  there is no such thing as const operater[], because hash-map creates null value by default
//    const mapped_type& operator[](const key_type&k)const{return dictionary[k];};
    mapped_type& operator[](const key_type&k){return dictionary[k];};
    const_iterator find(const key_type&k)const{return dictionary.find(k);};
    iterator find(const key_type&k){return dictionary.find(k);};
    void Print()const{
      std::cout<<"dictionary content:"<<std::endl;
      for(auto ii=dictionary.begin();ii!=dictionary.end();++ii){
          std::cout << ii->first << std::endl;
          for(auto jj=ii->second.begin();jj!=ii->second.end();++jj){
              std::cout <<"\t"<<jj->sTag<<"\t"<<jj->llStart<<"\t"<<jj->llEnd<<std::endl;
          }
      }
    }

  private:
    dictionary_t dictionary; // a dictionary with key (string) and value which is a vector of Region

    // add a region record from a string, note that it doesn't do "compact"
    int addRecord(const std::string& sLine,const ssize_t pad_length=0
                 , const size_t nSkippedColumns=1
                 , const std::string sDelim= " \t");
};

};

#endif
