///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "random/Shuffle.h"
#include "feudal/VirtualMasterVec.h"
#include "feudal/IncrementalWriter.h"

typedef VirtualMasterVec<BaseVec> VBaseVecVec;
typedef VirtualMasterVec<QualVec> VQualVecVec;


/** Pick a subset of reads

\file FastbQualbSelect.cc

*/

static inline 
String Tag(String S = "FQS") { return Date() + " (" + S + "): "; } 



int main( int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String_Doc(HEAD_IN, "Take base vectors from '<HEAD_IN>.fastb'."); 
  CommandArgument_String_Doc(HEAD_OUT, "Put selected base vectors in '<HEAD_OUT>.fastb'."); 
  CommandArgument_Bool_OrDefault_Doc(DO_QUALS, False, "Also process quals in '<HEAD_IN>.qualb'."); 
  CommandArgument_IntSet_OrDefault_Doc(READ_IDS, "", "The indexes of the selected base vectors."); 
  CommandArgument_UnsignedInt_OrDefault_Doc(N_RANDOM, 0, "Choose a random subset of N_RANDOM base vectors.");
  CommandArgument_Double_OrDefault_Doc(FRACTION, 0, "Choose a random fraction of the base vectors.");
  CommandArgument_UnsignedInt_OrDefault_Doc(SEED, 666, "Random seed.");
 
  EndCommandArguments;

  const String fn_fastb_in = HEAD_IN + ".fastb";
  const String fn_qualb_in = HEAD_IN + ".qualb";

  const String fn_fastb_out = HEAD_OUT + ".fastb";
  const String fn_qualb_out = HEAD_OUT + ".qualb";


  vec<size_t> i_bv_select;
  size_t n_bv_out;

  if (N_RANDOM > 0 || FRACTION > 0) {
    const size_t n_bv_in = MastervecFileObjectCount(fn_fastb_in);
    n_bv_out = (N_RANDOM ? N_RANDOM : FRACTION * n_bv_in);
    cout << Tag() << "Shuffling " << n_bv_in << " base vector indexes." << endl;
    Shuffle64(n_bv_in, i_bv_select, SEED);
    i_bv_select.resize(n_bv_out);
  }
  else {  // use provided read ids
    n_bv_out = READ_IDS.size();
    for (size_t i = 0; i != n_bv_out; i++) 
      i_bv_select.push_back(READ_IDS[i]);
  }
  
  sort(i_bv_select.begin(), i_bv_select.end());

  if (true) {
    cout << Tag() << "Picking " << n_bv_out << " base vectors from '" << fn_fastb_in << "' and writing them to '" << fn_fastb_out << "'." << endl;
    VBaseVecVec bvv_in(fn_fastb_in.c_str());
    IncrementalWriter<BaseVec> bvv_out(fn_fastb_out.c_str());
    
    for (size_t i = 0; i < n_bv_out; i++) {
      bvv_out.add(bvv_in[i_bv_select[i]]);
      dots_pct(i, n_bv_out);
    }
  }

  if (DO_QUALS) {
    cout << Tag() << "Picking " << n_bv_out << " qual vectors from '" << fn_qualb_in << "' and writing them to '" << fn_qualb_out << "'." << endl;
    VBaseVecVec qvv_in(fn_qualb_in.c_str());
    IncrementalWriter<BaseVec> qvv_out(fn_qualb_out.c_str());
    
    for (size_t i = 0; i < n_bv_out; i++) {
      qvv_out.add(qvv_in[i_bv_select[i]]);
      dots_pct(i, n_bv_out);
    }
  }

}

