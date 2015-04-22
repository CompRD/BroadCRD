/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

///Select flows from fasta/qual files, random or by name, and put in new file.
/// \file MakeFastaSubset.cc
///
/// Parameters:
/// - FASTA: Input fasta file
/// - QUAL: Input Qual file (optional)
/// - SELECTED_NAMES: should point to a vec file of Strings,
/// which are the ones to select from the big file. If this parameter is
/// empty, then we assume we want a random selection of NREADS reads
/// from the big file, and we will produce that.
/// - PARTIAL_MATCH: if true, then accept any names that contain one of
/// the names in SELECTED_NAMES, even if they are not exactly identical.
/// Works correctly with COMPLEMENT. 
/// This can be  expensive if SELECTED_NAMES is long!
/// - COMPLEMENT: if True, then get all names that are _not_ in SELECTED_NAMES
/// - SHORTEN_NAMES: if true, assume the names in selected names are
/// the short (integer) version of the long 111111_2222_3333 names
/// that flow reads usually have, and do the comparison accordingly.
/// - SHORTEN_OUTPUT_NAMES: if true, save the names as short names in 
/// the output fasta file.
/// - NFLOWS: number of reads to get from FASTA file in the case of a 
/// random selection. Ignored if SELECTED_NAMES is not empty.
/// .

#include "MainTools.h"
#include "FastaFileset.h"
#include "Basevector.h"
#include "random/Shuffle.h"


#include <sstream>
#include <set>

template<class V> 
int PrintIfNameMatches(ostream & os, const V & v, 
		       vecString & names,  const set<String> & good,
		       bool SHORTEN_NAMES, bool SHORTEN_OUTPUT_NAMES, 
		       bool COMPLEMENT, bool PARTIAL_MATCH, bool COUT_INDEX );

vec<String> GetRandomNamesVec(const vecString & names, const int N,
                      const int seed = 0) {
  const int MAX_N=names.size();
  ForceAssertLe(N, MAX_N);
  vec<int> shuffled;
  vec<String> selected;
  Shuffle(MAX_N, shuffled, seed);
  shuffled.resize(N);
  for (int i = 0; i !=shuffled.isize(); ++i) {
    selected.push_back(names[shuffled[i]]);
  }
  return selected;
}

int main(int argc, char ** argv) {

  BeginCommandArguments;

  CommandArgument_String(FASTA);
  CommandArgument_String_OrDefault(QUAL,"");
  CommandArgument_String_OrDefault(OUT_PREFIX,"")
  CommandArgument_String_OrDefault(SELECTED_NAMES,"");
  CommandArgument_Bool_OrDefault(SHORTEN_NAMES,True);
  CommandArgument_Bool_OrDefault(SHORTEN_OUTPUT_NAMES,False);
  CommandArgument_Bool_OrDefault(COMPLEMENT, False);
  CommandArgument_Bool_OrDefault(PARTIAL_MATCH, False);
  // whether to print the selected indices to cout
  CommandArgument_Bool_OrDefault(COUT_INDEX, False); 
  CommandArgument_UnsignedInt_OrDefault(NREADS,2000);
  CommandArgument_UnsignedInt_OrDefault(SEED,0);

  EndCommandArguments;

  if (OUT_PREFIX.empty()) {
    OUT_PREFIX = FASTA.SafeBefore("fa") + "subset.";
    if (SELECTED_NAMES.empty()) {
      OUT_PREFIX += "random.";
    } else {
      OUT_PREFIX += Basename(SELECTED_NAMES);
    }
  }
  if (OUT_PREFIX[OUT_PREFIX.size() -1] != '.') {
    OUT_PREFIX += ".";
  }

  vecbasevector reads;
  vecString names, namesq;
  vecqualvector readsq;
  bool doquals = false;

  FastFetchReads(reads, &names, FASTA);
  if (!QUAL.empty()) {
    FastFetchQuals(readsq, &namesq, QUAL);
    doquals = true;
  }

  //If we didn't get a list of flow names, generate and write out a random one.
  vec<String> selectedNames;
  if (SELECTED_NAMES.empty()) {
    if (COMPLEMENT) {
      InputErr( "COMPLEMENT can only be set to True if "
		<< "SELECTED_NAMES has been set to a filename" );
    }
    selectedNames = GetRandomNamesVec(names, NREADS, SEED);
    SELECTED_NAMES = OUT_PREFIX + "names";
    Ofstream(namestream, SELECTED_NAMES);
    namestream << selectedNames;
    cout << "Wrote file of " << NREADS << " randomly selected names to " 
         << SELECTED_NAMES << endl;
  }

  //Get the names
  Ifstream(namestream, SELECTED_NAMES);
  namestream >> selectedNames;
  NREADS = selectedNames.size();
  cout << "Reading " << NREADS << " selected names from " 
       << SELECTED_NAMES  << endl;
   
  //open the output files.
  String fastaname = OUT_PREFIX + "fasta";
  Ofstream(fasta, fastaname);

  //Put the desired names into a set for easy searching.
  set<String> good;
  String name;
  for (int i = 0; i != selectedNames.isize(); ++i) { 
    name = selectedNames[i];
    DeleteLeadingWhiteSpace(name);
    DeleteTrailingWhiteSpace(name);
    if (SHORTEN_NAMES) {
      name = name.SafeBefore("_");
      int n = name.Int();
      name = ToString(n);
    }
    good.insert(name);
  }
  
  int lookedFor = COMPLEMENT ? reads.size() - good.size() : good.size();

  //Select and print the reads and quals
  int ffound = PrintIfNameMatches(fasta, reads, names, good,
				  SHORTEN_NAMES,  SHORTEN_OUTPUT_NAMES, 
				  COMPLEMENT, PARTIAL_MATCH, COUT_INDEX);
  int qfound = 0;
  if (doquals) {
    String qualname = OUT_PREFIX + "qual";
    Ofstream(qual, qualname);    
    qfound = PrintIfNameMatches(qual, readsq, namesq, good,
				SHORTEN_NAMES,  SHORTEN_OUTPUT_NAMES, 
				COMPLEMENT, PARTIAL_MATCH, False);
    cout << "Put " << qfound << " quals out of " 
	 << lookedFor  << " looked for into " << qualname << endl;
  }

  //Output result.
  cout << "Put " << ffound << " reads ";
  if (!PARTIAL_MATCH) cout << "out of " << lookedFor  << " looked for ";
  cout << "into " << fastaname << endl;
  return 0;
}

template<class V> 
int PrintIfNameMatches(ostream & os, const V & v, 
		       vecString & names,  const set<String> & good,
		       bool SHORTEN_NAMES, bool SHORTEN_OUTPUT_NAMES, 
		       bool COMPLEMENT, bool PARTIAL_MATCH, bool COUT_INDEX ) {
  ForceAssertEq(v.size(), names.size());
  String name;
  int found=0;
  set<String>::iterator iter;
    
  if (COUT_INDEX)
    cout << "Selected indices:\n";
  for (size_t i =0 ; i != names.size(); ++i) {
    bool ok = false;
    name=names[i];
    DeleteLeadingWhiteSpace(name);
    DeleteTrailingWhiteSpace(name);
    if (SHORTEN_NAMES) {
      name = name.SafeBefore("_");
      int n = name.Int();
      name = ToString(n);
      if (SHORTEN_OUTPUT_NAMES) {
	names[i] = name;
      }
    }
    
    iter = good.find(name);
    if (!PARTIAL_MATCH) {
      if ( (iter == good.end() && COMPLEMENT) 
	   || (iter != good.end() && !COMPLEMENT)) { 
	ok = true;
      }
    }
    else { //look at every good name for a partial match
      ok = COMPLEMENT;
      for (iter = good.begin(); iter != good.end(); ++iter) {
	if (name.Contains(*iter)) {
	  ok = !COMPLEMENT;
	  break;
	}
      }
    }
    if (ok) {
      if (COUT_INDEX)
	cout << i << "\n";
      Print(os, v[i], names[i]);
      ++found;
    }
  }
  return found;
}
