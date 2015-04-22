///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

const char *DOC =
  " Merges two fastb files into a single file in which the reads are interleaved."
  " It will also merge the associated qualb files and create an appropriate" 
  " .pairs file it takes as input two fastb files representing each of the two"
  " reads in a set of pairs (e.g., the i'th read in <HEAD1_IN>.fastb and the i'th"
  " read in <HEAD2_IN>.fastb form a pair).";
  

#include "Basevector.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "Qualvector.h"
#include "String.h"
#include "feudal/IncrementalWriter.h"
#include "feudal/VirtualMasterVec.h"


static inline 
String Tag(String S = "MPF") { return Date() + " (" + S + "): "; } 


int main(int argc, char *argv[])
{
  RunTime();
  
  BeginCommandArguments;
  CommandDoc(DOC);
  
  // Files to input/output
  CommandArgument_String_Doc(HEAD1_IN,
    "Input reads 1 from '<HEAD1_IN>.{fastb,qualb}'.");
  CommandArgument_String_Doc(HEAD2_IN,
    "Input reads 2 from '<HEAD2_IN>.{fastb,qualb}'.");
  CommandArgument_String_Doc(HEAD_OUT,
    "Output interleaved reads to '<HEAD_OUT>.{fastb,qualb}'.");
  
  // Read trimming and reversing
  CommandArgument_Int_OrDefault_Doc(TRIM_START, 0,
    "Position of first base to include. If 0, include from the start.");
  CommandArgument_Int_OrDefault_Doc(TRIM_END, 0,
    "Position of last base to include. If 0, include to the end.");
  CommandArgument_Bool_OrDefault_Doc(REVERSE_READS, False,
    "Whether to reverse reads before trimming.");
  CommandArgument_Bool_OrDefault_Doc(WRITE_PAIRS, False,
    "Write .pairs pairing information file.");
  CommandArgument_String_OrDefault_Doc(LIB_NAME, "",
    "Library name used in the .pairs file.");
  CommandArgument_Int_OrDefault_Doc(LIB_SEP, 0,
    "Library mean separation used in the .pairs file.");
  CommandArgument_Int_OrDefault_Doc(LIB_STDEV, 0,
    "Library standard deviation used in the .pairs file.");
  EndCommandArguments;

  
  typedef VirtualMasterVec<qualvector> VmvQv_t;
  typedef VirtualMasterVec<basevector> VmvBv_t;
  VmvBv_t basesA_in((HEAD1_IN + ".fastb").c_str());
  VmvBv_t basesB_in((HEAD2_IN + ".fastb").c_str());
  VmvQv_t qualsA_in = VmvQv_t((HEAD1_IN + ".qualb").c_str());
  VmvQv_t qualsB_in = VmvQv_t((HEAD2_IN + ".qualb").c_str());
  ForceAssertEq(basesA_in.size(), basesB_in.size());
  ForceAssertEq(basesA_in.size(), qualsA_in.size());
  ForceAssertEq(basesA_in.size(), qualsB_in.size());
  
  VmvBv_t::const_iterator bvA_itr = basesA_in.begin();
  VmvBv_t::const_iterator bvB_itr = basesB_in.begin();
  VmvQv_t::const_iterator qvA_itr = qualsA_in.begin();
  VmvQv_t::const_iterator qvB_itr = qualsB_in.begin();
  
  IncrementalWriter<basevector> bases_out((HEAD_OUT + ".fastb").c_str());
  IncrementalWriter<qualvector> quals_out((HEAD_OUT + ".qualb").c_str());
 
  
  
  // Trim reads?
  bool trim_reads = false;
  int start_base = 0, end_base = 0;
  if (TRIM_START != 0 || TRIM_END != 0) {
    trim_reads = true;
    start_base = TRIM_START;
    end_base = (TRIM_END ? TRIM_END + 1 : basesA_in[0].size());
    cout << Tag() << "Trimming, keeping bases " << start_base << " to " << end_base - 1
	 << " (" << end_base - start_base << " bases)" << endl;
  }
  
  
  if (REVERSE_READS)
    cout << Tag() << "Reversing read pairs." << endl;  // Actual reversing done later
  
  
  cout << Tag() << "Merging." << endl;
  
  // Read alternately from each file in order to create the new read order.
  basevector bv1,bv2;
  qualvector qv1,qv2;
  
  size_t pair_count = 0;
  while (bvA_itr != basesA_in.end()) {
    
    if (!trim_reads && !REVERSE_READS) {
      bases_out.add(*bvA_itr);
      bases_out.add(*bvB_itr);
      quals_out.add(*qvA_itr);
      quals_out.add(*qvB_itr);
    } else {
      if (trim_reads) {
	bv1.assign(bvA_itr->begin() + start_base, bvA_itr->begin() + end_base);
	bv2.assign(bvB_itr->begin() + start_base, bvB_itr->begin() + end_base);
	qv1.assign(qvA_itr->begin() + start_base, qvA_itr->begin() + end_base);
	qv2.assign(qvB_itr->begin() + start_base, qvB_itr->begin() + end_base);
	
      } else {
	bv1 = *bvA_itr;
	bv2 = *bvB_itr;
	qv1 = *qvA_itr;
	qv2 = *qvB_itr;
      }
      if (REVERSE_READS) {
	bv1.ReverseComplement();
	bv2.ReverseComplement();
	qv1.ReverseMe();
	qv2.ReverseMe();
      }
      bases_out.add(bv1);
      bases_out.add(bv2);
      quals_out.add(qv1);
      quals_out.add(qv2);
    }
    
    bvA_itr++;
    bvB_itr++;
    qvA_itr++;
    qvB_itr++;
    pair_count++;
  }
  
  bases_out.close();
  quals_out.close();


  cout << Tag() << "Merged " << pair_count << " read pairs." << endl;

  
  if (WRITE_PAIRS) {
    cout << Tag() << "Creating PairsManager." << endl;

    PairsManager pairs(pair_count * 2);
    int lib_id = pairs.addLibrary(LIB_SEP, LIB_STDEV, LIB_NAME);
    for (size_t i = 0; i < pair_count; i++) 
      pairs.addPairToLib( i * 2, i * 2 + 1, lib_id);
    pairs.Write(HEAD_OUT + ".pairs");
  }
  
  // Done!
  cout << Tag() << "Done." << endl;
  return 0;
}
