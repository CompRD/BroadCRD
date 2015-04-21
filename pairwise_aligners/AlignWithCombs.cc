/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "TaskTimer.h"

#include "pairwise_aligners/CombAligner.h"
#include <fstream>





int main( int argc, char** argv )
{
 
  RunTime();

  BeginCommandArguments;
  CommandArgument_String(REFERENCE_FASTB);
  CommandArgument_String(SEQUENCE_FASTB);
  CommandArgument_String_OrDefault(OUT_ALIGNS, "comb.aligns");
  CommandArgument_Bool_OrDefault(RC, True);
  CommandArgument_Bool_OrDefault(DOUBLE_SPACED, True);
  CommandArgument_Bool_OrDefault(PRINT, True);
  CommandArgument_UnsignedInt_OrDefault(MAX_ERRORS, 4);
  CommandArgument_Bool_OrDefault(TRUNCATE, False);
  CommandArgument_Bool_OrDefault(TRUNCATE_N, 36);
  CommandArgument_Bool_OrDefault(RC_SEQ, False);

  CommandArgument_UnsignedInt_OrDefault(K, 10);
  EndCommandArguments;




  CombAligner ccAligner;

  vecbasevector seq;
 
  bool rc = true;
  if (RC == False)
    rc = false;

  bool bPrint = true;


  if (PRINT == False)
    bPrint = false;

  ccAligner.SetPrintAligns(bPrint);


  ccAligner.SetK((int)K);
  int i;

  if (DOUBLE_SPACED == True) {
    ccAligner.UseDoubleSpacing(true);
  } else {
    ccAligner.UseDoubleSpacing(false);
  }

  vec<String> names;
 
  cout << "Reading sequences " << SEQUENCE_FASTB << endl;
  seq.ReadAll(SEQUENCE_FASTB);
  for (i=0; i<(int)seq.size(); i++) {
    char tmp[128];
    sprintf(tmp, "read_%d", i);
    names.push_back(tmp);    
  }

  if (TRUNCATE == True) {
    cout << "Truncating bases..." << endl;
    for (i=0; i<(int)seq.size(); i++) {
      basevector b;
      b.SetToSubOf(seq[i], 0, (int)TRUNCATE_N);
      seq[i] = b;
    }
    cout << "done!" << endl;
  }
  if (RC_SEQ == True) {
    cout << "rc'ing bases..." << endl;
    for (i=0; i<(int)seq.size(); i++) {
      seq[i].ReverseComplement();
    }
    cout << "done!" << endl;
  }

  ccAligner.ReadFastb(REFERENCE_FASTB, rc);

  

  vec<PrelimMatch> matches;

  for (i=0; i<(int)seq.size(); i++) {
    ccAligner.Align(matches, seq[i], names[i], i);
  }


  cout << "Writing summary to file " << OUT_ALIGNS << endl;

  std::ofstream out(OUT_ALIGNS.c_str());
  out << "read\trefcontig\tstart\torient\terrors\tsnps\n";
  for (i=0; i<matches.isize(); i++) {
    out << matches[i];
  }
  out.close();
  ForceAssert(out);
  
  return 0;
}
