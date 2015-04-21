/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Vec.h"
#include "feudal/BinaryStream.h"
#include "kmers/KmerRecord.h"

// We have to force these instantiations to make the compiler happy.
Bool lt_kmer( const kmer_with_count<12>&, const kmer_with_count<12>& );
Bool lt_kmer( const kmer_with_count<16>&, const kmer_with_count<16>& );
Bool lt_kmer( const kmer_with_count<20>&, const kmer_with_count<20>& );
Bool lt_kmer( const kmer_with_count<24>&, const kmer_with_count<24>& );


template <int K> 
void Compare( const String& file1, const String& file2 ) {
  vec<kmer_with_count<K> > table1;
  BinaryReader::readFile( file1, &table1 );
  vec<kmer_with_count<K> > table2;
  BinaryReader::readFile( file2, &table2 );

  typename vec< kmer_with_count<K> >::const_iterator iKmer1, iKmer2;
  iKmer1 = table1.begin();
  iKmer2 = table2.begin();

  const int cap = 1000;

  vec< vec<longlong> > transitions( cap+1, vec<longlong>( cap+1, 0 ) );
  while ( iKmer1 != table1.end() && iKmer2 != table2.end() ) {
    int count1, count2;
    if ( lt_kmer( *iKmer1, *iKmer2 ) ) {
      count1 = iKmer1->Count();
      count2 = 0;
      ++iKmer1;
    }
    else if ( lt_kmer( *iKmer2, *iKmer1 ) ) {
      count1 = 0;
      count2 = iKmer2->Count();
      ++iKmer2;
    }
    else { 
      ForceAssert( eq_kmer( *iKmer1, *iKmer2 ) ); 
      count1 = iKmer1->Count();
      count2 = iKmer2->Count();
      ++iKmer1;
      ++iKmer2;
    }
    ++transitions[ min(cap,count1) ][ min(cap,count2) ];
  }

  for ( unsigned int freq1 = 0; freq1 < transitions.size(); ++freq1 ) {
    for ( unsigned int freq2 = 0; freq2 < transitions[freq1].size(); ++freq2 ) {
      if ( freq2 > 0 ) cout << "\t";
      cout << transitions[freq1][freq2];
    }
    cout << "\n";
  }
}


int main( int argc, char** argv ) {
  
  RunTime();

  BeginCommandArguments;
  CommandArgument_String( TABLE1 );
  CommandArgument_String( TABLE2 );
  CommandArgument_Int( K );
  EndCommandArguments;

  // mode can be "diff" "missing" or "extra"

  switch ( K ) {
    case 12:
      Compare<12>( TABLE1, TABLE2 );
      break;
    case 16:
      Compare<16>( TABLE1, TABLE2 );
      break;
    case 20:
      Compare<20>( TABLE1, TABLE2 );
      break;
    case 24:
      Compare<24>( TABLE1, TABLE2 );
      break;
    default:
      cout << "K must be one of 12, 16, 20, or 24." << endl;
      break;
  }
}
      
  
 


