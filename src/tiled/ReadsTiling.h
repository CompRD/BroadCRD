// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef READS_TILING_H
#define READS_TILING_H

#include "PackAlign.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "Vec.h"
#include "tiled/PaddedSeq.h"



/*
 * class reads_tiling
 *
 * A reads_tiling is a way to store informations about a specific type
 * of contig. Input data are a master sequence and tiling reads, i.e.
 * reads which overlap master without leaving any gap on it.
 */
class reads_tiling {
  
public:
  
  reads_tiling( );
  
  void SetPointers( const basevector *master_bases,
		    const vecbasevector *primate_bases,
		    const vecqualvector *primate_quals );
  
  void SetFromAligns( const vec<int> *read_ids,
		      const vec<align> *aligns,
		      const vec<Bool> *isRC,
		      ostream *log );
  
  void Clear( );

  void Compactify( );

  void GenerateConsensus( bool tag_snp=false );

  int ConsensusBegin( ) const;

  int ConsensusLength( ) const;
  
  int ReadsCount( ) const;
  
  int WinBegin( ) const;
  
  int WinEnd( ) const;

  int BeginOnConsensus( int read_index ) const;

  void ReadsIds( vec<int> &ids ) const;
  
  void ReadsFullNames( vec<String> &names ) const;

  void MasterToFasta( ostream &out, int width=60 ) const;

  void ConsensusBasesToFasta( ostream &out, int width=60 ) const;
  
  void ConsensusQualsToQual( ostream &out, int width=25 ) const;

  const padded_seq &MasterSeq( ) const;

  const padded_seq &PaddedSeq( int ii ) const;
  
  const vec< vec<char> > &RBases( ) const { return rbases_; }
  
  const vec< vec<int> > &RQuals( ) const { return rquals_; }

  void SlimPrint( ostream &out,
		  int width = 150 ) const;
  
  void FoldPrint( ostream &out,
		  bool quals = true,
		  int width = 150,
		  int from = -1,
		  int to = -1 ) const;
  
  friend istream &operator>> ( istream &in, reads_tiling &tiles );
  
  friend ostream &operator<< ( ostream &out, const reads_tiling &tiles );
  
  
private:
  
  void Adjust( const vec<align> *aligns, int read_index, ostream *log );
  
  void PadClone( int &pos, int read_index );
  
  void PadRead( int &pos, int read_index );
  
  void ColumnConsensus( int column, bool tag_snp=false );

  void GenerateRectangle( bool tag_snp=false );
  
  void GenerateNames( vec<String> &read_name ) const;

  char PrettyQual( int n_qual ) const;
  
  
private:
  
  const basevector *master_bases_;      // may be null
  const vecbasevector *primate_bases_;  // may be null
  const vecqualvector *primate_quals_;  // may be null
  
  int win_begin_;                       // begin of window on master sequence
  padded_seq master_pad_;               // master sequence
  vec<String> tile_name_;           // name of each tile read
  vec<padded_seq> tile_pad_;            // tiles
  
  vec< vec<char> > rbases_;             // rectangle bases
  vec< vec<int> > rquals_;              // rectangle quals
  
  mutable bool rectangle_begin_found_;  // if rectangle_begin has been found
  mutable int rectangle_begin_;         // actual rectangle begin (may be <0)

};


// Load tiled contigs (reads_tilings). Not efficient if file_name contains
//  a very large number of read_tilings.
void LoadReadsTilings( const String &file_name,
		       vec<reads_tiling> &tilings );


#endif
