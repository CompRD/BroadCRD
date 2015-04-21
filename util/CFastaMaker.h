/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef C_FASTA_MAKER_H
#define C_FASTA_MAKER_H

#include "Basevector.h"
#include "Qualvector.h"
#include "String.h"
#include "SupersHandler.h"

/**
 * CFastaMaker
 *
 * Load fastb/qualb and save them into fasta/qual chunks. For example:
 * save the bases in a window [a,b) of a given super as a single fasta/
 * qual entry, with Ns in the gaps. Notice that gaps will be reset at
 * >= min_gap_.
 */
class CFastaMaker {

public:
  
  CFastaMaker( const String &supers_file,
	       const vecbasevector &bases, 
	       const vecqualvector &quals );

  ~CFastaMaker( );

  void SetGapFloor( int min_gap );

  void SetStreams( ostream *fout, ostream *qout = 0, ostream *log = 0 );
  
  void SuperChunk( int super_id, int *begin = 0, int *end = 0 );
  
  
private:

  void SetDefaults( const String &supers_file );

  void SetToPrint( int contig_id, int *begin = 0, int *end = 0 );

  void PrintFasta( );
  
  void PrintQual( );

  void EndFasta( );

  void EndQual( );
  
  
private:

  shandler *supers_;             // supers (with gaps reset >= min_gap_)
  const vecbasevector &bases_;   // contig bases
  const vecqualvector &quals_;   // contig quals

  int fwidth_;      // page width (fasta)
  int qwidth_;      // page width (qual)
  int min_gap_;     // min gap size
  ostream *fout_;   // output fasta stream
  ostream *qout_;   // output qual stream (may be null)
  ostream *log_;    // log stream (may be null)
  
  int fpos_;        // current position on the line (fasta)
  int qpos_;        // current position on the line (qual)
  int contig_id_;   // id of contig to be printed
  int contig_beg_;  // begin base of contig to be printed
  int contig_end_;  // end base of contig to be printed
  int contig_Ns_;   // how many Ns to print (i.e. size of gap after contig)
  
};

#endif
