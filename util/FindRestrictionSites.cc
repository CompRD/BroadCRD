/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/** Make a FASTA file surrounding all restriction sites
 \file FindRestrictionSites.cc

   - SITE: enzyme restriction site as a string (e.g. GCGGGCCGC)
   - REF: fasta or fastb file for the reference (default human)
   - OUT_PREFIX: prefix for output fasta file
   - MARGIN: number of bases to print before and after the site.
*/

#include "MainTools.h"
#include "Basevector.h"
#include "FastaFileset.h"

void PrintContext(ostream & fasta, const basevector & bases, int pos,
		  int length, int MARGIN, const String & name);

int main(int argc, char ** argv) {
  RunTime();
  BeginCommandArguments;
  CommandArgument_String(SITE);
  CommandArgument_String_OrDefault(REF,"/wga/scr3/LOOKUP/human35/human.fastb");
  CommandArgument_String_Abbr_OrDefault(OUT_PREFIX,O,"");
  CommandArgument_Int_OrDefault(MARGIN,200);
  EndCommandArguments;

  vecbasevector ref;
  if (REF.Contains(".fastb"))    ref.ReadAll(REF);
  else FastFetchReads(ref,0,REF);

  if (OUT_PREFIX.empty()) OUT_PREFIX=SITE;

  basevector site;
  site.SetFromString(SITE);

  Ofstream(fasta,OUT_PREFIX + ".fasta");
  String name;

  for (size_t i=0; i != ref.size(); ++i) {
    int pos = 0; 
    while (true) {
      pos = ref[i].Find(site,pos);
      if (ref[i].isize() == pos) break;
      PrintContext(fasta,ref[i],pos, site.size(), MARGIN,
		   "chr" + ToString(i) + "_" + ToString(pos));
      ++pos;
      PRINT2(i,pos);
    }
  }
  return 0;
}

void PrintContext(ostream & fasta, const basevector & bases, int pos,
		  int length, int MARGIN, const String & name) {
    int start = max(0,int(pos - MARGIN));
    int end = min(pos + MARGIN + length, bases.isize());
    basevector b;
    b.SetToSubOf(bases, start, end-start);
    b.Print(fasta, name);
}
   
