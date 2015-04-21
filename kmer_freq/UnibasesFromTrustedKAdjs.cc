/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   Program: UnibasesFromTrustedKAdjs

   Construct <unibases> directly from trusted <kadjs>.

   By default the output files are named using the same root as the
   input files, i. e. READS_OUT = READS_IN. This behaviour may be
   overridden by explicitly setting the value of READS_OUT.

   See also: MarkTrustedA, MarkTrustedB.

   INPUT/OUTPUT:
   
   - All Input/Output files are in: PRE/DATA/RUN
   - Input files:
        ../READS_IN.KADJ_IN.k(K+1) (vec of kmer_with_count<K+1>)
   - Output files:
        ../READS_OUT.KADJ_GRAPH_OUT.kK (KAdjGraph<K>)
	../READS_OUT.UNIBASES_OUT.kK.fastb (vecbasevector of unibases)
	../READS_OUT.UNIGRAPH_OUT.kK (digraph)
*/


#include "MainTools.h"

#include "math/Functions.h"
#include "kmers/KmerShape.h"
#include "kmers/KmerRecord.h"
#include "Bitvector.h"
#include "Intvector.h"
#include "Vec.h"
#include "paths/KmerPath.h"
#include "paths/KmerBaseBroker.h"
#include "paths/UnibaseUtils.h"
#include "feudal/BinaryStream.h"

#include "kmer_freq/KAdjGraph.h"

template < int K >
void CreateUnibases( String run_dir,
		     String READS_IN,
		     String KADJS_IN,
		     String READS_OUT,
		     String UNIBASES_OUT,
		     String UNIPATH_GRAPH_OUT,
		     String KADJ_GRAPH_OUT,
		     String ALT_EXT)
{
  vec< kmer_with_count< K+1 > > kadjs_with_count;

  BinaryReader::readFile( run_dir + "/" + READS_IN + "." + KADJS_IN + ".k" + ToString( K+1 ),
                          &kadjs_with_count );

  PRINT( kadjs_with_count.isize() );
  
  vec< kmer< K+1 > > kadjs( kadjs_with_count.size() );

  basevector bv;
  for( int i = 0; i < kadjs_with_count.isize(); i++ ) {
    kadjs_with_count[ i ].GetBasevector( bv );
    kadjs[ i ].Set( bv );
  }

  PRINT( kadjs.size() );
  
  DPRINT2( "Sorting kadjs", kadjs.size() );
  UniqueSort( kadjs );
  DPRINT( kadjs.size() );

  cout << "Constructing kadj graph from " << kadjs.size() << " kadjs..." << endl;
  KAdjGraph<K> kadjGraph( kadjs );

  String fileName(run_dir+"/"+READS_OUT+"."+KADJ_GRAPH_OUT+".k"+ToString(K));
  BinaryWriter::writeFile(fileName.c_str(),kadjGraph);

  cout << " Computing unibases from the kadj graph..." << endl;
  vecbasevector unibases;
  kadjGraph.ComputeUnibases( unibases );
  PRINT( unibases.size() );

  ForceAssertGt( unibases.size(), 0u );
  
  vec< size_t > sizes( unibases.size() );
  for ( size_t i = 0; i < unibases.size(); i++ ) {
    sizes[ i ] = unibases[ i ].size() - K + 1;
    ForceAssertGt( sizes[ i ], 0u );
  }
  Sort( sizes );
  PRINT( N50( sizes ) );

  // Write to unibases file (possibly with alternate file extension)
  String unibases_ext = ( ALT_EXT != "" ?
			  ALT_EXT :
			  "k" + ToString( K ) + ".fastb" );
  
  unibases.WriteAll( run_dir + "/" + READS_OUT + "." + UNIBASES_OUT + "." + unibases_ext );

  digraph unipathAdjGraph;
  kadjGraph.ComputeUnipathAdjGraph( unipathAdjGraph );
  BinaryWriter::writeFile( run_dir + "/" + READS_OUT + "." + UNIPATH_GRAPH_OUT + ".k" + ToString( K ), unipathAdjGraph );

}  // CreateUnibases()


int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA); 
  CommandArgument_String(RUN);
  CommandArgument_Int(K);
  CommandArgument_String_OrDefault(READS_IN, "reads_orig");
  CommandArgument_String_OrDefault(READS_OUT, "");
  CommandArgument_String_OrDefault(KADJ_IN, "fastb.trusted_as_strong");
  CommandArgument_String_OrDefault(UNIBASES_OUT, "trusted_unibases");
  CommandArgument_String_OrDefault(UNIGRAPH_OUT, "unipathGraph");
  CommandArgument_String_OrDefault(KADJ_GRAPH_OUT, "adjGraph");
  CommandArgument_String_OrDefault(ALT_EXT, "");
  EndCommandArguments;

  String run_dir = PRE + "/" + DATA + "/" + RUN;

  if (READS_OUT == "")
    READS_OUT = READS_IN;

#define CASE(K) CreateUnibases<K>( run_dir, READS_IN, KADJ_IN, READS_OUT, \
				   UNIBASES_OUT, UNIGRAPH_OUT, KADJ_GRAPH_OUT,\
                                   ALT_EXT )

  DISPATCH_ON_K_WITH_K_PLUS_1( K, CASE );
  
}
