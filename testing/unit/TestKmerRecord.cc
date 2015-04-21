/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "system/RunTime.h"

#include "kmers/KmerRecord.h"

int main( int argc, char** argv ) {

  RunTime();
  
  basevector in;
  in.SetFromString( "ATGCATAGAGAGAAATTACACGAG" );
  
  // kmer_record
  kmer_record<24> kr;
  kr.Set( in, 0, 0 );
  
  basevector kr_out;
  kr.GetBasevector( kr_out );
  
  ForceAssert( kr_out == in );

  // kmer
  kmer<24> k( in );
  
  basevector k_out;
  k.GetBasevector( k_out );

  ForceAssert( k_out == in );

  // kmer_with_count
  kmer_with_count<24> kwc( in, 1 );
  
  basevector kwc_out;
  kwc.GetBasevector( kwc_out );

  ForceAssert( kwc_out == in );

  cout << "All tests passed." << endl;
}
