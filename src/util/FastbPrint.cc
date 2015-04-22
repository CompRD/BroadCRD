///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Basevector.h"
#include "Qualvector.h"



inline String ANSI256RGBColorID(int red, int green, int blue)
{  return "8;5;" + ToString(16 + blue + 6 * (green + 6 * red)); }

inline String ANSISetBackground(String color) { return "\x1b[4" + color + "m"; }
inline String ANSISetForeground(String color) { return "\x1b[3" + color + "m"; }
inline String ANSISetDefault() { return "\x1b[0m"; }


void print_base_vec(const BaseVecVec & bvv,
		    const QualVecVec & qvv,
		    const size_t ibv,
		    const size_t ib0,
		    const size_t ib1,
		    const int COLOR,
		    const bool RC,
		    const bool QUALS)
{
  const BaseVec & bv = bvv[ibv];
  for (size_t ib = ib0; ib < ib1; ib++) {
    char base = as_base((RC) ? 3 - bv[ib1 - (ib - ib0) - 1] : bv[ib]);
    
    String color_str;
    if (QUALS) {
      const QualVec & qv = qvv[ibv];
      
      if (COLOR == 1) {
	color_str = ((qv[ib] >= 20) ? "2" :  // green
		     ((qv[ib] >= 10) ? "3" :  // yellow
		      "1")); // red
      }
      else {
	// Color scheme for 256-color case: Red progressing to green
	int color = Min(qv[ib] / 8, 5);
	color_str = ANSI256RGBColorID(5 - color, color, 0);
      }
    }
    else {
      if (base == 'A') color_str = ANSI256RGBColorID(5, 0, 0);
      if (base == 'C') color_str = ANSI256RGBColorID(2, 5, 0);
      if (base == 'G') color_str = ANSI256RGBColorID(3, 0, 5);
      if (base == 'T') color_str = ANSI256RGBColorID(0, 5, 5);
    }
    
    if (COLOR == 3) {
      if (base == 'A') base = '^';
      if (base == 'C') base = '(';
      if (base == 'G') base = '-';
      if (base == 'T') base = '.';
    }
    
    cout << ANSISetForeground(color_str) + base;
  }
  
  cout << ANSISetDefault();
}


/*
 * PrintBases
 *
 * Print the given bases from the specifed entry in a fastb file.
 */
int main(int argc, char *argv[])
{
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_String(HEAD);
  CommandArgument_IntSet(READ_IDS);
  CommandArgument_UnsignedInt_OrDefault(BEGIN, 0);
  CommandArgument_UnsignedInt_OrDefault(END, (INT_MAX));
  CommandArgument_Bool_OrDefault(QUALS, False);
  CommandArgument_Bool_OrDefault(RC, False);
  CommandArgument_Bool_OrDefault(FULL, False);
  CommandArgument_UnsignedInt_OrDefault(COLOR, 3);
  EndCommandArguments;

  ForceAssertGt(END, BEGIN);

  const size_t nb_max = 100;

  BaseVecVec bvv;
  bvv.SparseRead(HEAD + ".fastb", READ_IDS, 0);
  size_t nibv = READ_IDS.size();
    

  QualVecVec qvv;
  if (QUALS) 
    qvv.SparseRead(HEAD + ".qualb", READ_IDS, 0);


  for (size_t iibv = 0; iibv != nibv; iibv++) {
    const size_t ibv = READ_IDS[iibv];

    const BaseVec & bv = bvv[ibv];
  
    const size_t ib0 = BEGIN;
    const size_t ib1 = (END > bv.size() ? bv.size() : END);
    ForceAssertGt(ib1, ib0);

    cout << "read ID = " << ibv << endl;
    cout << "read size = " << bv.size() << endl;
    
    if (!FULL && ib1 - ib0 > 2 * nb_max) { 
      if (RC) {
	print_base_vec(bvv, qvv, ibv, ib0, ib0 + nb_max, COLOR, RC, QUALS);
	cout << "   ...(" << (ib1 - ib0 - 2 * nb_max) << " bases)...   ";
	print_base_vec(bvv, qvv, ibv, ib1 - nb_max, ib1, COLOR, RC, QUALS);
      }
      else {
	print_base_vec(bvv, qvv, ibv, ib1 - nb_max, ib1, COLOR, RC, QUALS);
	cout << "   ...(" << (ib1 - ib0 - 2 * nb_max) << " bases)...   ";
	print_base_vec(bvv, qvv, ibv, ib0, ib0 + nb_max, COLOR, RC, QUALS);
      }
    }
    else {  
      print_base_vec(bvv, qvv, ibv, ib0, ib1, COLOR, RC, QUALS);
    }
    cout << endl << endl;
  
  }
  cout << endl;

}
