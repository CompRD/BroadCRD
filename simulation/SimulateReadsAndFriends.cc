///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Oct 18, 2012
//

#include "MainTools.h"
#include "feudal/BaseVec.h"
#include "simulation/ReadSimulatorSimpleCore.h"
#include "simulation/ReferenceIterator.h"
#include "Qualvector.h"


int main(int argc, char* argv[])
{
  // Commandline args
  int N=250;
  float ERR_SUB=0.10;
  int K=12;
  int F1=50;
  int F2=50;


  // Simulate Random Read of size N
  RandomGen base_random;
  BaseVec read(N);
  QualVec qual(N,40);
  for ( int i = 0; i < N; ++i )
    read.set(i,base_random.unsignedN(4));

  // Add errors to read with mismatch rate M
  ReferenceIterator err_gen(read, 0, false, false, 0., 0., ERR_SUB);
  BaseVec read_err(N);
  QualVec qual_err(N);
  for ( int i = 0; i < N && !err_gen.done(); ++i, ++err_gen ) {
    read_err.set(i, *err_gen);
    qual_err.set(i, ( read[i] == read_err[i] ) ? 40 : 25 );
  }


  // Simulate F1 true friends with minimum overlap K and size N
  BaseVecVec friends(F1+F2);
  QualVecVec friends_quals(F1+F2);

  int start_start = K - N;
  int start_end   = N - K;
  int start = start_start;
  for ( int i = 0; i < F1; ++i ) {
    friends[i].resize(N);
    friends_quals[i].resize(N);

    int perfect = 0;
    for ( int j = 0; j < N; ++j ) {
      if ( j < start ) {
	friends[i].set(j,4);
	friends_quals[i].set(j,1);
      } else if ( j >= start + N ) {
	friends[i].set(j,4);
	friends_quals[i].set(j,1);
      } else {
	if ( perfect < K || base_random.float01() >= ERR_SUB ) {
	  friends[i].set(j,read[j]);
	  friends_quals[i].set(j,40);
	  if ( perfect < K ) perfect++;
	} else {
	  unsigned char err_base = (base_random.unsignedN(3)+1)^read[j];	// create one of the other 3 bases
	  friends[i].set(j, err_base);
	  friends_quals[i].set(j, 25);
	}
      }
    }

    start += K-1;
    if ( start > start_end )		// we add K-1 so we tile a bit
      start = start_start;

  }


  // Simluate F2 false friends with perfect match K and random otherwise
  start_start = K - N;
  start_end   = N - K;
  start = start_start;
  for ( int i = F1; i < F1+F2; ++i ) {
    friends[i].resize(N);
    friends_quals[i].resize(N);

    int perfect = 0;
    for ( int j = 0; j < N; ++j ) {
      if ( j < start ) {
	friends[i].set(j,4);
	friends_quals[i].set(j,1);
      } else if ( j >= start + N ) {
	friends[i].set(j,4);
	friends_quals[i].set(j,1);
      } else {
	if ( perfect < K ) {
	  friends[i].set(j,read[j]);
	  friends_quals[i].set(j,40);
	  if ( perfect < K ) perfect++;
	} else {
	  unsigned char err_base = base_random.unsignedN(4);
	  friends[i].set(j, err_base);
	  friends_quals[i].set(j, base_random.float01() < ERR_SUB ? 25 : 40 );
	}
      }
    }

    start += K-1;
    if ( start > start_end )		// we add K-1 so we tile a bit
      start = start_start;

  }




  // Write read, read_err, friends, and friend_indicator file
  BaseVecVec tmp_read(2);
  tmp_read[0] = read;
  tmp_read[1] = read_err;
  tmp_read.WriteAll("reads.fastb");

  QualVecVec tmp_qual(2);
  tmp_qual[0] = qual;
  tmp_qual[1] = qual_err;
  tmp_qual.WriteAll("reads.qualb");

  friends.WriteAll("friends.fastb");

  ofstream ind("indicator.txt");
  for ( int i = 0; i < F1+F2 && ind.good() ; ++i )
    if ( i < F1 )
      ind << 0 << endl;
    else
      ind << 1 << endl;
  ind.close();

  ForceAssert(ind.good());

  friends_quals.WriteAll("friends.qualb");

  return 0;
}




