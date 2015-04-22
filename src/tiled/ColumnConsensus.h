/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef COLUMN_CONSENSUS_H
#define COLUMN_CONSENSUS_H

#include "Vec.h"
#include "tiled/CSnp.h"
#include "tiled/Haploqual.h"
#include "tiled/PaddedSeq.h"

/**
 * class column_consensus
 *
 * Organizes column data to find the consensus on a set of bases and
 * quality scores. Pads are valid bases.
 */
class column_consensus {
  
public:
  
  column_consensus( );

  ~column_consensus( );
  
  void SetDefaults( );

  void Add( char base, int qual );

  int Coverage( ) const;

  int Coverage( char base ) const;
  
  // This is being superseded by SnpConsensus.
  bool Consensus( char &base, int &qual, bool tag_snp = true ) const;

  bool SnpConsensus( CSnp &snp, int nhap, bool trusted = true ) const;

  void PrintInfo( ostream &out, bool newline = false ) const;

  const haploqual *A( ) const { return A_; }

  const haploqual *C( ) const { return C_; }

  const haploqual *G( ) const { return G_; }

  const haploqual *T( ) const { return T_; }

  const haploqual *Pad( ) const { return pad_; }
  
  
private:

  // Used by Consensus( ) only.
  bool IsSNP( const haploqual &haplo1, const haploqual &haplo2 ) const;
  
  // Used by Consensus( ) only.
  char BaseSNP( char base1, char base2 ) const;

  void PackageSort( vec< haploqual* > &haplo ) const;
  

private:
  
  int good_qual_;       // defines good quality score
  int max_to_exclude_;  // discard third haplotype if qual is low
  int min_to_include_;  // allow second haplotype if qual is high
  
  haploqual *A_;
  haploqual *C_;
  haploqual *G_;
  haploqual *T_;
  haploqual *pad_;

};

#endif
