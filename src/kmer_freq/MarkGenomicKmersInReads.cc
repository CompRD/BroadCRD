/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   Program: MarkGenomicKmersInReads

   For each position of each read, determine whether the kmer starting
   at that position is genomic.  More generally, determine that kmer's
   copy number.  Record this info in a vecbitvector and VecIntVec,
   respectively.

   This is useful e.g. for helping develop MarkTrusted algorithms:
   when looking at the instances of a kmer, we can instantly tell whether
   the kmer is genomic or not (and so test when a particular heuristic
   gives the right answer and when it does not).

   Since this info depends only on the particular set of reads,
   it is recorded in the DATA directory.
 */


#include <limits>
#include "Basevector.h"
#include "Bitvector.h"
#include "CoreTools.h"
#include "CommonSemanticTypes.h"
#include "FeudalMimic.h"
#include "kmers/KmerRecord.h"
#include "kmers/KmerShape.h"
#include "MainTools.h"
#include "Qualvector.h"
#include "Charvector.h"
#include "math/Functions.h"
#include "kmers/SortKmers.h"
#include "CommonSemanticTypes.h"
#include "feudal/BinaryStream.h"


template< class KSHAPE >
void DoMarkGenomicKmersInReads( dirname_t data_dir, filenamepart_t READS ) {
  const nbases_t K = KSHAPE::KSIZE;

  vecbasevector genome( data_dir + "/genome.fastb" );
  vecbasevector reads( data_dir + "/" + READS + ".fastb" );

  VecUCharVec kmerCopyNums( reads.size() );
  for ( size_t i = 0; i < reads.size(); i++ )
    kmerCopyNums[ i ].resize( reads[ i ].size(), char(0) );

  vecbitvector kmerGenomic;
  Mimic( reads, kmerGenomic );

  vec< kmer_with_count< K > > genomicKmersNotInReads;

  typedef kmer_record< K, 2 > krec_t;

  typedef typename vec< krec_t >::const_iterator kocc_it_t;

  vec< int > counts( int( numeric_limits< uchar >::max() ) + 1, 0 );

  vec<bvec const*> allReads;
  allReads.reserve(reads.size()+genome.size());
  for ( vecbvec::const_iterator end(reads.cend()), itr(reads.cbegin()); itr != end; ++itr )
      allReads.push_back(&*itr);
  for ( vecbvec::const_iterator end(genome.cend()), itr(genome.cbegin()); itr != end; ++itr )
      allReads.push_back(&*itr);

  FOR_KMER_OCCS_BEG( KSHAPE, krec_t, reads, USE_STABLE_SORT, PALIND_FW_ONLY, from, to ) {

    ForceAssert( to > from );

    kocc_it_t startGenomeOccs = to - 1;
    while ( startGenomeOccs >= from && static_cast<size_t>(startGenomeOccs->GetId()) >= reads.size() )
      startGenomeOccs--;

    startGenomeOccs++;
    
    ForceAssert( startGenomeOccs == to  ||  static_cast<size_t>(startGenomeOccs->GetId()) >= reads.size() );
    if ( !( startGenomeOccs-1 < from  ||  static_cast<size_t>((startGenomeOccs-1)->GetId()) < reads.size() ) ) {
      for ( kocc_it_t z = from; z < to; z++ ) {
	PRINT3( reads.size(), genome.size(), to-from );
	PRINT5( z - from, z->GetId(), z->GetReadPos(), int( z->IsFw() ), z->GetPos() );
      }
    }
    ForceAssert( startGenomeOccs-1 < from  ||  static_cast<size_t>((startGenomeOccs-1)->GetId()) < reads.size() );
    
    copy_num_t trueCopyNum = to - startGenomeOccs;
    uchar trueCopyNumBounded = uchar( min( trueCopyNum, int( numeric_limits< uchar >::max() ) ) );
    counts[ trueCopyNumBounded ]++;
    
    if ( trueCopyNum > 0 )
      for ( kocc_it_t occInReads = from; occInReads != startGenomeOccs; occInReads++ ) {
	kmerCopyNums[ occInReads->GetId() ][ occInReads->GetReadPos() ] = trueCopyNumBounded;
	kmerGenomic[ occInReads->GetId() ].Set( occInReads->GetReadPos(), True );
      }

    if ( startGenomeOccs == from ) {
      basevector b;
      startGenomeOccs->GetBasevector( b );
      genomicKmersNotInReads.push( b, 1 );
    }
    
  }
  FOR_KMER_OCCS_END();

  kmerCopyNums.WriteAll( data_dir + "/" + READS + ".kmerCopyNums.k" + ToString( KSHAPE::getId() ) );
  kmerGenomic.WriteAll( data_dir + "/" + READS + ".kmerGenomic.k" + ToString( KSHAPE::getId() ) );
  UniqueSort( genomicKmersNotInReads );
  PRINT( genomicKmersNotInReads.size() );
  BinaryWriter::writeFile( data_dir + "/" + READS + ".missingGenomic.k" + ToString( KSHAPE::getId() ), genomicKmersNotInReads );
  PRINT( counts );
  
}


int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(K);
  CommandArgument_String_OrDefault(READS, "reads_orig");
  EndCommandArguments;

  cout << Date( ) << ": loading data" << endl;

  dirname_t data_dir = PRE + "/" + DATA;

#define CASE(KS) DoMarkGenomicKmersInReads<KS>( data_dir, READS )
  DISPATCH_ON_KSHAPE( K, CASE );
  //ForceAssertEq( K, ToString( "20" ) );
  //CASE( KmerShapeDefaultClass(20) );
  
}
