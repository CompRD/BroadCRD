/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   \file
   
   FastbDiff.cc: reads in two fastb files, A and B, and outputs a *sorted*
   fastb file which is all reads in A that is not in B.

*/

#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <map>

#include "MainTools.h"
#include "Basevector.h"
#include "Qualvector.h"






bool bv_sz_lt(const BaseVec & bva_fw, const BaseVec & bvb_fw)
{
  const size_t nba = bva_fw.size();
  const size_t nbb = bvb_fw.size();
  if (nba < nbb) return true;
  if (nba > nbb) return false;

  BaseVec bva_can = bva_fw;
  BaseVec bvb_can = bvb_fw;
  
  bva_can.Canonicalize();
  bvb_can.Canonicalize();

  return bva_can < bvb_can;
}

bool bv_sz_eq(const BaseVec & bva_fw, const BaseVec & bvb_fw)
{
  const size_t nba = bva_fw.size();
  const size_t nbb = bvb_fw.size();
  if (nba < nbb) return false;
  if (nba > nbb) return false;

  BaseVec bva_can = bva_fw;
  BaseVec bvb_can = bvb_fw;
  
  bva_can.Canonicalize();
  bvb_can.Canonicalize();

  return bva_can == bvb_can;
}



String hieroglyph(unsigned base)
{
  if (base == 0) return "^";
  if (base == 1) return "(";
  if (base == 2) return "-";
  if (base == 3) return ".";
  ForceAssertLt(base, 4u);
  return "?";
}

String hieroglyphs(const BaseVec & bv)
{
  String h = "";
  unsigned nb = bv.size();
  for (unsigned ib = 0; ib < nb; ib++)
    h += hieroglyph(bv[ib]);
  return h;
}




int main( int argc, char *argv[] ) 
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String_OrDefault(HEAD_A, "");
  CommandArgument_String_OrDefault(HEAD_B, "");
  CommandArgument_String_OrDefault(READS_A, "");
  CommandArgument_String_OrDefault(READS_B, "");
  CommandArgument_String_OrDefault(READS_A_OUT, "");
  CommandArgument_String_OrDefault(READS_B_OUT, "");
  CommandArgument_String_OrDefault(OUT_INDICES, "");
  CommandArgument_Bool_OrDefault(PRINT, False);
  CommandArgument_Bool_OrDefault(DO_QUALS, False);
  CommandArgument_Bool_OrDefault(DO_SORT, True);
  EndCommandArguments;

  if (READS_A == "" && HEAD_A == "") 
    FatalErr("You must specify either READS_A or HEAD_A.");

  if (READS_B == "" && HEAD_B == "") 
    FatalErr("You must specify either READS_B or HEAD_B.");

  if (DO_QUALS && (HEAD_A == "" || HEAD_B == "")) 
    FatalErr("You must specify HEAD_A and HEAD_B if DO_QUALS=True.");

  
  if (READS_A == "") READS_A = HEAD_A + ".fastb";
  if (READS_B == "") READS_B = HEAD_B + ".fastb";



  const bool do_sort = (OUT_INDICES == "" && !DO_QUALS && DO_SORT);

  bool equal = true;

  if (!do_sort) {
    BaseVecVec reads_a(READS_A);
    BaseVecVec reads_b(READS_B);
    
    size_t n_reads_a = reads_a.size();
    size_t n_reads_b = reads_b.size();
  
    vec<size_t> iv;

    const size_t n_reads_min = (n_reads_a < n_reads_b) ? n_reads_a : n_reads_b;
    const size_t n_reads_max = (n_reads_a > n_reads_b) ? n_reads_a : n_reads_b;

    for (size_t i = 0; i != n_reads_min; i++)
      if (reads_a[i] != reads_b[i])
	iv.push_back(i);

    if (iv.size() > 0)
      equal = false;
    

    if (OUT_INDICES != "") {
      ofstream os;
      os.open((OUT_INDICES + ".ind").c_str());
      for (size_t i = 0; i != iv.size(); i++)
        os << iv[i] << endl;
      os.close();
    }
    
    cout << "             No sorting done." << endl;
    cout << setw(12) << (n_reads_min - iv.size()) << " reads are the same." << endl;
    cout << setw(12) << iv.size() + (n_reads_max - n_reads_min) << " reads are different." << endl;

    if (DO_QUALS) {
      QualVecVec quals_a(HEAD_A + ".qualb");
      QualVecVec quals_b(HEAD_B + ".qualb");
      
      size_t nq_diff = 0;
      for (size_t i = 0; i != n_reads_min; i++)
        if (quals_a[i] != quals_b[i])
          nq_diff++;
      
      if (nq_diff > 0) 
        equal = false;

      cout << setw(12) << (n_reads_min - nq_diff) << " quals are the same." << endl;
      cout << setw(12) << nq_diff + (n_reads_max - n_reads_min) << " quals are different." << endl;
    }


    

  }
  else { // don't do indices  =>  sorting

    BaseVecVec reads_a(READS_A);
    BaseVecVec reads_b(READS_B);

    sort(reads_a.rbegin(), reads_a.rend(), bv_sz_lt);
    sort(reads_b.rbegin(), reads_b.rend(), bv_sz_lt);
    
    size_t n_reads_a = reads_a.size();
    size_t n_reads_b = reads_b.size();
    
    size_t n_reads_a_and_b = 0;
    size_t n_reads_a_not_b = 0;
    size_t n_reads_b_not_a = 0;

    BaseVecVec reads_a_out;
    BaseVecVec reads_b_out;

    size_t i_a = 0;
    size_t i_b = 0;
    while (i_a < n_reads_a && i_b < n_reads_b) {
      
      BaseVec & bva = reads_a[i_a];
      BaseVec & bvb = reads_b[i_b];

      if (bv_sz_eq(bva, bvb)) {
	++i_a;  
	++i_b;
	++n_reads_a_and_b;
      } 
      else {
        if (bv_sz_lt(bva, bvb)) { // read_a < read_b
          bva.Canonicalize();
          if (READS_A_OUT != "") reads_a_out.push_back_reserve(bva, 0, 2.0);
          ++i_a;
          ++n_reads_a_not_b;
        } 
        else { // reads_a > reads_b  
          bvb.Canonicalize();
          if (READS_B_OUT != "") reads_b_out.push_back_reserve(bvb, 0, 2.0);
          ++i_b;
          ++n_reads_b_not_a;
        }
      }
    }

    while(i_a < n_reads_a) {
      BaseVec & bva = reads_a[i_a];
      bva.Canonicalize();
      if (READS_A_OUT != "") reads_a_out.push_back_reserve(bva, 0, 2.0);
      ++i_a;
      ++n_reads_a_not_b;
    }

    while(i_b < n_reads_b) {
      BaseVec & bvb = reads_b[i_b];
      bvb.Canonicalize();
      if (READS_B_OUT != "") reads_b_out.push_back_reserve(bvb, 0, 2.0);
      ++i_b;
      ++n_reads_b_not_a;
    }



    if (READS_A_OUT != "") reads_a_out.WriteAll(READS_A_OUT);
    if (READS_B_OUT != "") reads_b_out.WriteAll(READS_B_OUT);

    if (n_reads_a_not_b > 0 || n_reads_b_not_a > 0)
      equal = false;

    cout << setw(12) << n_reads_a << " reads in A" << endl;
    cout << setw(12) << n_reads_b << " reads in B" << endl;
    cout << setw(12) << n_reads_a_and_b << " reads in A and in B" << endl;
    cout << setw(12) << n_reads_a_not_b << " reads in A and not in B" << endl;
    cout << setw(12) << n_reads_b_not_a << " reads in B and not in A" << endl;
    cout << endl;
  
  }

  return (equal) ? 0 : 1;
}
