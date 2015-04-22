/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


#include "polymorphism/BaseCountHolder.h"

int LoadCoverageContig(istream & is, BaseCountHolder & raw,
BaseCountHolder & accepted, basevector * ref) {
  String line;
  vec<String> data;
  vec<int> idata(16);
  int pos=0;

  //read in the Contig line
  getline(is, line);
  if (is.fail()) return -1;
  Tokenize(line, data);
  if (data[0] != "Contig") return -1;
  int contig = data[1].Int();
  int length = data[3].SafeBefore(":").Int();
  if (ref) ref->resize(length);
  raw.SetNumLines(length);
  accepted.SetNumLines(length);
  getline(is, line);//get rid of the pos ref line
  for (int pos=0; pos != length; ++pos) {
    getline(is, line);
    Tokenize(line, data);
    AssertEq(pos,data[0].Int()); 
    if (ref) ref->Set(pos, as_base(data[1][0]));
    transform(data.begin()+2, data.end(), idata.begin(), 
	      mem_fun_ref(&String::Int));
    raw.SetLine(pos, idata.begin(), idata.begin() + 8);
    accepted.SetLine(pos, idata.begin()+8, idata.end());
  }
  return contig;
}

int LoadCoverageFile(const String & fname, vec<BaseCountHolder> & raw,
		      vec<BaseCountHolder> & accepted,
		     vecbasevector * ref) {
  Ifstream(is, fname);
  BaseCountHolder r, a;
  basevector * b = 0;
  raw.clear();
  accepted.clear();
  if (ref) {
    ref->destroy();
    b = new basevector();
  }
  int ret;
  while ((ret = LoadCoverageContig(is, r, a, b)) != -1) {
    raw.push_back(r);
    accepted.push_back(a);
    if (ref) ref->push_back_reserve(*b);
  }
  return raw.size();
}
