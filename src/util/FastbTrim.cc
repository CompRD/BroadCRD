/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// FastbExtractSubset: read a table of regions "contig start end", 
//                     and extract them from a Fastb.

#include "Basevector.h"
#include "MainTools.h"

int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String(IN);
  CommandArgument_String(OUT);
  CommandArgument_Int_OrDefault(TRIM_START, 0);
  CommandArgument_Int_OrDefault(TRIM_END,   0);
  CommandArgument_Int_OrDefault(TRIM_TO,    0);
  EndCommandArguments;

    vecbasevector in(IN);
    vecbasevector out(in.size());

    for (size_t i = 0; i < in.size(); i++)
    {
        if (TRIM_TO > 0) {
            out[i].SetToSubOf(in[i], TRIM_START, TRIM_TO);
        } else {
            out[i].SetToSubOf(in[i], TRIM_START, (in[i].size()-TRIM_END)-TRIM_START);
        }

        /*
        cout << in[i].ToString() << endl;
        cout << out[i].ToString() << endl;
        cout << endl; cout.flush();
        */
    }

    out.WriteAll(OUT);    

  return 0;
}
