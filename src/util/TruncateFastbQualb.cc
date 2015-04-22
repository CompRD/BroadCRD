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

#include "Basevector.h"
#include "Qualvector.h"

int main( int argc, char** argv )
{
 
  RunTime();

  BeginCommandArguments;
  CommandArgument_String(IN_FASTB);
  CommandArgument_String(IN_QUALB);
  CommandArgument_String(OUT_FASTB);
  CommandArgument_String(OUT_QUALB);
  CommandArgument_UnsignedInt_OrDefault(TRUNCATE_AT, 36);
  EndCommandArguments;



  vecbasevector b;
  vecqualvector q;
  cout << "Reading..." << endl;
  b.ReadAll(IN_FASTB);
  q.ReadAll(IN_QUALB);

 
  int i;
  cout << "Truncating..." << endl;
  for (i=0; i<(int)b.size(); i++) {
    qualvector qv = q[i];
    basevector bv = b[i];
    q[i].SetToSubOf(qv, 0, (int)TRUNCATE_AT);
    b[i].SetToSubOf(bv, 0, (int)TRUNCATE_AT);
  }

  cout << "Writing..." << endl;
  b.WriteAll(OUT_FASTB);
  q.WriteAll(OUT_QUALB);
  
  return 0;
}
