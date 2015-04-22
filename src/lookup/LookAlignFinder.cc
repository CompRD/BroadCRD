///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/**
   \file
   \copydoc LookAlignFinder
*/

#include "lookup/LookAlignFinder.h"
#include "CoreTools.h"

LookAlignFinder::LookAlignFinder(const String & fname): 
  aligns_(), last_(), found_() {
  ForceAssert(IsRegularFile(fname));
  isp_ = new ifstream(fname.c_str());

  SetEmpty();
  
  String line;
  while (true) {
    getline((*isp_), line);
    if (!isp_->good()) break;
    if (!line.Contains("QUERY",0)) continue;
    last_.ReadParseable(line);
    break;
  }

  if (!empty()) operator++();
}

LookAlignFinder & LookAlignFinder::operator++() {
  ForceAssert(!empty());
  //PRINT2("start", QueryId());
  if(found_.find(QueryId()) != found_.end()) {
    FatalErr("Alignments for read with id " << QueryId()
             << " are not all consecutive.");
  }
  found_.insert(QueryId());

  if (!isp_->good()) {
    SetEmpty();
    return *this; //we have reached the end.
  }

  aligns_.clear();
  aligns_.push_back(last_);
  String line;

  while (true) {
    if (isp_->peek() != 'Q') {
      if (!isp_->good()) break;
      isp_->ignore(numeric_limits<std::streamsize>::max(),'\n');
      continue;
    }

    getline((*isp_), line);
    if (!isp_->good()) break;
    if (!line.Contains("QUERY",0)) continue;
    last_.ReadParseable(line);
    if (last_.query_id == QueryId()) {
      aligns_.push_back(last_);
    }
    else break;
  }
  //PRINT2("end", QueryId());
  //for (int i=0; i != aligns.isize(); ++i) aligns[i].PrintReadableBrief(cout);
  
  return *this;
}
