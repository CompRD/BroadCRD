/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef T_ALIGN_H
#define T_ALIGN_H

#include "PackAlign.h"

/**
 * class t_align
 *
 * A structure for the data needed to describe the alignment onto
 * a master sequence: an align, the id of the aligned sequence, and
 * the orientation of the aligned sequence (the master is always fw).
 */
class t_align {
  
public:

  t_align( ) :
    id_ ( -1 ), isRC_ ( False ) { }
  
  t_align( int id, Bool isRC, const align al ) :
    id_ ( id ), isRC_ ( isRC ), al_ ( al ) { }
  
  void Set( int id, Bool isRC, const align al ) {
    id_ = id; isRC_ = isRC; al_ = al; }
  
  void Write( ostream &out ) const;

  void Read( istream &in );
  
  void WriteBinary( ostream &out ) const;

  void ReadBinary( istream &in );
  
  
 public:

  int id_;
  Bool isRC_;
  align al_;
  
};

/**
 * class t_align_plus
 *
 * It contains also the id of the contig the read aligns to.
 */
class t_align_plus : public t_align {

public:

  t_align_plus( ) : t_align( ), contig_id_ ( -1 ) { }

  t_align_plus( const t_align &t_al, int contig_id ) :
    t_align( t_al ), contig_id_ ( contig_id ) { }
  
  void Set( const t_align &t_al, int contig_id ) {
    t_align::Set( t_al.id_, t_al.isRC_, t_al.al_ ), contig_id_ = contig_id; }
  
  
public:

  int contig_id_;

};

#endif
