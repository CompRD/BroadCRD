// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef PRIMATE_CONTIG_H
#define PRIMATE_CONTIG_H

#include "Basevector.h"
#include "PackAlign.h"
#include "Qualvector.h"
#include "lookup/LookAlign.h"
#include "tiled/ReadsTiling.h"



/*
 * class primate_contig
 *
 * It generates tiled-type contigs. A tiled contig is built by looking
 * at how reads align to another genome (for instance, chimp reads on
 * human), and by inserting pads where needed.
 */
class primate_contig {
  
public:
  
  primate_contig( );
  
  primate_contig( const basevector *clone_bases,
		  const vecbasevector *primate_bases,
		  const vecqualvector *primate_quals,
		  ostream *log = 0 );
  
  void Set( const basevector *clone_bases,
	    const vecbasevector *primate_bases,
	    const vecqualvector *primate_quals,
	    ostream *log = 0 );
  
  void AddAlign( const look_align &new_align );
  
  void GenerateTiles( );

  const reads_tiling &Tiles( ) const;

  int BulkBeginOnClone( ) const;
  
  int BulkEndOnClone( ) const;
  
  int BulkLength( ) const;
  
  
private:
  
  const basevector *clone_bases_;
  const vecbasevector *primate_bases_;
  const vecqualvector *primate_quals_;
  ostream *log_;
  
  vec< look_align > aligns_;             // original alignments to human
  reads_tiling tiles_;                   // tiling for this contig

};



#endif
