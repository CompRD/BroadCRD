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
#include "ReadPairing.h"
#include "String.h"
#include "Vec.h"

/**
 * TrimFastb
 *
 * Take in input a set of reads (bases, quals, and pairing info), and
 * trim them to TRIM_SIZE bases, starting from the beginning of the
 * read. 
 *
 * HEAD_IN: it needs <HEAD_IN>.{fastb,qualb,pairto}
 * HEAD_OUT: it saves <HEAD_OUT>.{fastb,qualb,pairto}
 * TRIM_SIZE: how many bases to keep
 */

void ComputeTrimIndices(size_t TRIMMED_SIZE,
                        size_t NUM_TO_TRIM,
                        String TRIM_SIDE,
                        size_t orig_size,
                        size_t * ib0,
                        size_t * ib1)
{
  if (TRIMMED_SIZE > 0) { // NUM_TO_TRIM == 0
    if (TRIMMED_SIZE >= orig_size) {  // no trimming needed
      *ib0 = 0;
      *ib1 = orig_size;
    }
    else {
      if (TRIM_SIDE == "begin") {  // trim at the beginning
        *ib0 = orig_size - TRIMMED_SIZE;
        *ib1 = orig_size;
      }
      else {                       // trim at the end
        *ib0 = 0;
        *ib1 = TRIMMED_SIZE;
      }
    }
  }
  else { // NUM_TO_TRIM > 0
    if (NUM_TO_TRIM >= orig_size) { // trim everything
      *ib0 = *ib1 = 0;
    }
    else {
      if (TRIM_SIDE == "begin") {  // trim at the beginning
        *ib0 = NUM_TO_TRIM;
        *ib1 = orig_size;
      }
      else {                       // trim at the end
        *ib0 = 0;
        *ib1 = orig_size - NUM_TO_TRIM;
      }
    }
  }
}




int main(int argc, char *argv[])
{
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_String_Doc(HEAD_IN, "input from HEAD_IN.fastb");
  CommandArgument_String_OrDefault_Doc(HEAD_OUT, "", "output to HEAD_OUT.fastb");

  CommandArgument_UnsignedInt_OrDefault_Doc(TRIMMED_SIZE, 0, 
    "Final size after optional trimming");
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_TO_TRIM, 0, 
    "number of bases to trim");
  CommandArgument_String_OrDefault_Valid_Doc(TRIM_SIDE, "begin",
    "{begin,end}",
    "for paired reads, which side (AFTER REVERSAL) to trim the first read of the pair");

  EndCommandArguments;
  
  if ((TRIMMED_SIZE == 0 && NUM_TO_TRIM == 0) || 
      (TRIMMED_SIZE != 0 && NUM_TO_TRIM != 0)) {
    cout << "Must specify one (and only one) of TRIMMED_SIZE and NUM_TO_TRIM" << endl;
    return 0;
  }


  // Dir and file names.
  String bases_in_fn = HEAD_IN + ".fastb";
  if (! IsRegularFile(bases_in_fn)) {
    cout << bases_in_fn << " not found. Abort.\n" << endl;
    return 0;
  }

  String quals_in = HEAD_IN + ".qualb";
  bool quals_exist = IsRegularFile(quals_in);

  if (HEAD_OUT == "") 
    HEAD_OUT = HEAD_IN + ".trimmed";

  String pairs_in = HEAD_IN + ".pairto";
  bool pairs_exist = IsRegularFile(pairs_in);

  String out_dir = (HEAD_OUT.Contains("/")) ? HEAD_OUT.RevBefore("/") : ".";
  Mkpath(out_dir);
  
  String info_out_fn = HEAD_OUT + ".log";
  String bases_out_fn = HEAD_OUT + ".fastb";
  String quals_out_fn = HEAD_OUT + ".qualb";
  String base_out_fn = (out_dir == ".") ? HEAD_OUT : HEAD_OUT.After(out_dir);
  

  ofstream out(info_out_fn.c_str());
  PrintCommandPretty(out);
  out.close();

  // Original read lengths will be created (and used) only if TRIM_SIDE="end".
  vec<size_t> bv_size;



  
  // Bases.
  {
    cout << Date() << ": bases from " + bases_in_fn << endl;
    BaseVecVec bvv(bases_in_fn);
    size_t num_bvs = bvv.size();
    cout << Date() << ": num_bvs = " << num_bvs << endl;
    cout << Date() << ": trimming bases" << endl;


    if (TRIM_SIDE == "end") {
      bv_size.resize(num_bvs);
      for (size_t v_ID = 0; v_ID < num_bvs; v_ID ++)
        bv_size[v_ID] = bvv[v_ID].size();
    }
      
    IncrementalWriter<BaseVec> bwriter(bases_out_fn.c_str());

    for (size_t v_ID = 0; v_ID < num_bvs; v_ID ++) {
      const BaseVec &bv = bvv[v_ID];
      
      size_t ib0, ib1;
      ComputeTrimIndices(TRIMMED_SIZE, NUM_TO_TRIM, TRIM_SIDE, 
                         bv.size(), &ib0, &ib1);

      bwriter.add(BaseVec(bv, ib0, ib1 - ib0));
    }
    bwriter.close();
  }  








  // Quals.
  if (quals_exist) {
    cout << Date() << ": trimming quals" << endl;
    QualVecVec qvv(quals_in);
    size_t num_qvs = qvv.size();

    IncrementalWriter<QualVec> qwriter(quals_out_fn.c_str());

    for (size_t v_ID = 0; v_ID < num_qvs; v_ID++) {
      const QualVec &qv = qvv[v_ID];
      
      size_t ib0, ib1;
      ComputeTrimIndices(TRIMMED_SIZE, NUM_TO_TRIM, TRIM_SIDE, 
                         qv.size(), &ib0, &ib1);
    
      QualVec qv_new(ib1 - ib0);
      
      CopyQuals(qv, ib0, qv_new, 0, ib1 - ib0);
      qwriter.add(qv_new);
    }
    qwriter.close();
  }  











  
  // Pairing info.
  
  if (pairs_exist) {  
    cout << Date() << ": adjusting pairing info: NOT IMPLEMENTED!" << endl;
    /*
    vec<read_pairing> pairs;
    ReadPairsFile(pairs_in, pairs);
    if (! TAIL) {
      for(size_t v_ID = 0; v_ID < pairs.size(); v_ID++) {
	int id1 = pairs[v_ID].id1;
	int id2 = pairs[v_ID].id2;
	int amt1 = Max(0, bv_size[id1] - TRIM_SIZE);
	int amt2 = Max(0, bv_size[id2] - TRIM_SIZE);
	pairs[v_ID].sep += amt1 + amt2;
      }
    }
    WritePairs(out_dir, pairs, nbases, False, base_out_fn);
    */
  }
  
  // Done.
  cout << Date() << ": done" << endl;
  
}

