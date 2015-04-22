/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef C_PADDER_H
#define C_PADDER_H

#include "Basevector.h"
#include "Qualvector.h"
#include "tiled/PaddedSeq.h"
#include "tiled/Tiling.h"

/**
 * class CPadder
 *
 * Class to manage gap movements in tilings.
 */
class CPadder {

public:

  CPadder( const vecbasevector &cbases,
	   const vecbasevector &rbases,
	   ostream *log = 0 );
  
  // Set tiling (tiling will be modified!)
  void SetTiling( tiling *theTiling ) { theTiling_ = theTiling; }
  
  // Set streams to print tilings before and after shifts.
  void SetOutStreams( ostream *before_out, ostream *after_out );

  // Justify all reads in the tiling
  void BlockJustifyLeft( tiling *theTiling );

  // Justify pads of given read (must setTiling first! See .cc)
  int BlockJustifyLeft( int pos );

  
private:

  // Print alignment of read (given by position in tiling) on consensus
  void PrintAlign( int pos, ostream &out ) const;
  
  // Print tiling.
  void PrintTiling( const tiling &tiles, ostream *out ) const;
  
  
private:
  
  const vecbasevector &cbases_;  // contig bases
  const vecbasevector &rbases_;  // read bases
  ostream *log_;                 // optional log file
  ostream *before_out_;          // optional, to print tilings before shifts
  ostream *after_out_;           // optional, to print tilings after shifts

  tiling *theTiling_;            // the core tiling info
  
};

#endif
