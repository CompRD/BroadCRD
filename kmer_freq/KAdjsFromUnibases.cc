/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


#include "Basevector.h"
#include "MainTools.h"
#include "graph/Digraph.h"
#include "kmers/KmerRecord.h"
#include "kmer_freq/KAdjGraph.h"
#include "feudal/BinaryStream.h"

template < int K >
void ConvertUnibasesToKmerAdjs(String run_dir, String READS, String UNIBASES_IN,
			       String UNIGRAPH_IN, String KMER_ADJS_OUT,
			       String KADJ_GRAPH_OUT) {

  String KS = ToString(K);
  String K1S = ToString(K + 1);

  // Load unibases and unibase/unipath adj graph
  cout << "Load unibases..." << endl;
  vecbasevector unibases(run_dir + "/" + READS + "." + UNIBASES_IN + ".k" + KS + ".fastb") ;
  cout << "Load unibases adj graph..." << endl;
  digraph unipathAdjGraph;
  BinaryReader::readFile((run_dir+"/"+READS+"."+UNIGRAPH_IN+".k"+ToString(K)).c_str(),&unipathAdjGraph);

  ForceAssertEq(unibases.size(), static_cast<size_t>(unipathAdjGraph.N()));

  cout << "Found << " << unibases.size() << " unibases" << endl;
  
  // Compute number of kmers adjs in unibases
  size_t kmerAdjCount = (unibases.sumSizes() - unibases.size() * K) / 2;
  PRINT(kmerAdjCount);
  
  // Build kmer adj as a vec of kmer_with_counts

  // Kmers from unibases...
  cout << "Computing kmers adjs." << endl;
  vec< kmer_with_count< K+1 > > kadjs_with_count(kmerAdjCount);
  basevector bv;
  int index = 0;
  for (size_t unibaseNo = 0; unibaseNo < unibases.size(); unibaseNo++) {
    basevector& currentUnibase = unibases[unibaseNo];
    for( int i = 0; i < currentUnibase.isize() - K ; i++ ) {
      bv.SetToSubOf(currentUnibase, i, K + 1);
      if (bv.Canonicalize() != CanonicalForm::REV)
	kadjs_with_count[index++].Set( bv, 1 );
    }
  }
  PRINT(index);

  // Kmers from unibase graph...
  for (int vertex = 0; vertex < unipathAdjGraph.N(); vertex++) {
    vec<int> from = unipathAdjGraph.From(vertex);
    basevector currentUnibase = unibases[vertex];
    for (int i = 0; i < from.isize(); i++) {
      bv.SetToSubOf(currentUnibase, currentUnibase.size() - K, K);
      bv.AppendBase(unibases[from[i]][K-1]);
      if (bv.Canonicalize() != CanonicalForm::REV) {
	kmer_with_count< K+1 > kawc(bv,1);
	kadjs_with_count.push_back(kawc );
      }
    }
  }

  UniqueSort(kadjs_with_count);
  PRINT(kadjs_with_count.isize());

  // Write shaved kmer_with_counts for error correction
  cout << "Writing kmer adjs..." << endl;
  BinaryWriter::writeFile((run_dir+"/"+READS+"."+KMER_ADJS_OUT+".k"+K1S).c_str(),
                      kadjs_with_count);

  // Compute annd save KAdjGraph
  vec< kmer< K+1 > > kadjs( kadjs_with_count.size() );
  for( int i = 0; i < kadjs_with_count.isize(); i++ ) {
    kadjs_with_count[ i ].GetBasevector( bv );
    kadjs[ i ].Set( bv );
  }
  
  DPRINT2( "Sorting kadjs", kadjs.size() );
  UniqueSort( kadjs );
  DPRINT( kadjs.size() );

  cout << "Constructing kadj graph from " << kadjs.size() << " kadjs..." << endl;
  KAdjGraph<K> kadjGraph( kadjs );
  String fileName(run_dir+"/"+READS+"."+KADJ_GRAPH_OUT+".k"+ToString( K ));
  BinaryWriter::writeFile(fileName.c_str(),kadjGraph);

  cout << "Conversion Complete" << endl;
}


int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN); 
  CommandArgument_Int(K);
  CommandArgument_String_OrDefault(READS, "reads_orig");
  CommandArgument_String_OrDefault(UNIBASES_IN, "trusted_unibases");
  CommandArgument_String_OrDefault(UNIGRAPH_IN, "unipathAdjGraph");
  CommandArgument_String_OrDefault(KADJ_GRAPH_OUT, "adjGraph");
  CommandArgument_String_OrDefault(KMER_ADJS_OUT, "fastb.trusted_as_strong");

  EndCommandArguments;

  // Define directories.
  
  String run_dir = PRE + "/" + DATA + "/" + RUN;

  ConvertUnibasesToKmerAdjs<20>(run_dir, READS, UNIBASES_IN, UNIGRAPH_IN,
				KMER_ADJS_OUT, KADJ_GRAPH_OUT);
}

