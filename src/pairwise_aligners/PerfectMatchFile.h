/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef PERFECT_MATCH_FILE_H
#define PERFECT_MATCH_FILE_H

#include <fstream>

/** This class provides a file output iterator for PerfectMatches.

\class PerfectMatchFile

Putting the PerfectMatches into a file is helpful if they are going to
overwhelm  memory, or if memory is already full because we have a huge
vector of kmers, for example.

*/
class PerfectMatchFile: private std::ofstream {
private:
  longlong size_;

public:
  PerfectMatchFile(const char * name): 
    ofstream(name, std::ios::binary), 
    size_(0) {}
  
  longlong size() const { return size_; }
   
  void WriteToFile(const PerfectMatch & m) {
    write((char *)(&m), sizeof(PerfectMatch));
    ++size_;
  }

  class Iterator {
  private:
    PerfectMatchFile * f;

  public:
    Iterator(PerfectMatchFile * f): f(f) {}
    Iterator & operator*() { return *this; }
    void operator=(const PerfectMatch & m) { f->WriteToFile(m); }  
    Iterator & operator++() { return *this; }
    Iterator & operator++(int) { return *this; }
  };

  Iterator iterator() { return Iterator(this); }

};

template<class ForwardIter>
void copy(ForwardIter begin, ForwardIter end, PerfectMatchFile::Iterator i) {
  while (begin != end) {
    *i = *begin;
    ++i;
    ++begin;
  }
}

///Meant to be used as follows:
///Ifstream (is, matchfilename);
///PerfectMatch m;
///while (NextPerfectMatch(is, m) { //process this match; }
bool NextPerfectMatch(istream & is, PerfectMatch & m) {
  is.read((char *)&m, sizeof(PerfectMatch));
  if (!is) return false;  
  else return true;
}
 
#endif // PERFECT_MATCH_FILE_H

