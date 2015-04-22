/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "paths/AlignSeqsToHyper.h"
#include "paths/SeqOnHyper.h"

int main()
{
  RunTime();

  int K = 20;

  HyperBasevector hb(K);

  basevector e0,e1,e2,e3,e4,e5;
  e0.SetFromString(       "GGGGGGGGGGAAAAAAAAAACCCCC");
  e1.SetFromString("GGGG" "AAAAAAAAAACCCCCCCCCCTTTTTTTTTT");
  e2.SetFromString("CCCC" "CCCCCTTTTTTTTTTAAAAAAAAAA");
  e3.SetFromString("GGGG" "AAAAAAAAAACCCCCCCCCC");
  e4.SetFromString("AAAA" "AAAAACCCCCCCCCCTTTTTTTTTT");
  e5.SetFromString("AAAA" "AAAAACCCCCCCCCCTTTTTGGGGGGGGGG");

  hb.AddVertices(6);
  hb.AddEdge(0,1,e0);
  hb.AddEdge(1,2,e1);
  hb.AddEdge(2,3,e2);
  hb.AddEdge(1,4,e3);
  hb.AddEdge(4,2,e4);
  hb.AddEdge(4,5,e5);

  basevector s0,s1,s2;
  s0.SetFromString("GGGGGAAAAAAAAAACCCCCCCCCCTTTTTGGGGG");
  s1.SetFromString("GGGGGGGGGGAAAAAAAAAACCCCCCCCCCTTTTTTTTTTAAAAA");
  s2.SetFromString("GGGGGGGGGGAAAAAAAAAACCCCCCCCCCTTTTT");
  
  vecbasevector seqs;
  seqs.push_back(s0);
  seqs.push_back(s1);
  seqs.push_back(s2);
  // test rc too!
  s0.ReverseComplement();
  s1.ReverseComplement();
  s2.ReverseComplement();
  seqs.push_back(s0);
  seqs.push_back(s1);
  seqs.push_back(s2);

  int errors = VerifyAlignSeqsToHyper( hb, seqs );

  PRINT( CompressedSeqOnHyper::AmbigCount() );
  PRINT( errors );

  ForceAssertEq( errors, 0 );
}
