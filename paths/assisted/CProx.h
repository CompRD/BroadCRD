/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__ASSISTED__C_PROX__H
#define PATHS__ASSISTED__C_PROX__H

#include "paths/reporting/CLinkBundle.h"
#include "paths/assisted/CLibLinks.h"
#include "paths/assisted/CProxInfo.h"

/**
 * class CProx
 *
 * Container for a proximity based adjacency between two oriented
 * contigs, where the gap between the contigs can be computed from
 * different sources (jump reads, long jump reads, or as evinced from
 * the alignments of the two contigs onto a reference).
 *
 * See also CProxInfo (metadata for paired libraries, and heuristics).
 */
class CProx {
  
public:
  
  CProx( const CProxInfo *info = 0 );
  
  void SetInfo( const CProxInfo *info ) { info_ = info; }

  // Add ref-based gap (if valid).
  void AddRefGap( const int gap );

  // Add link-based gaps from the same lib (specified by tag_id).
  void AddLinkGaps( const int tag_id, const vec<int> &gaps );

  // Add link-based gaps from the same lib (specified by tag_id).
  void AddBundle( CLinkBundle &bundle) { clb_ = bundle; }

  // Sort links, and compute expected gap/dev.
  void ComputeLinkGap( );
  
  // Various const accessors to link-based gaps.
  int NLinksTotal( ) const { return gaps_.isize( ); }

  int NLinksWinner( ) const { return lwin_.second - lwin_.first; }

  int Gap( int ii ) const { return gaps_[ii].first; }

  int LibDev( int ii ) const { return info_->LibDev( gaps_[ii].second ); }

  bool LibType( int ii ) const { return info_->LibType( gaps_[ii].second ); }

  // Return INT_MAX if there is no valid ref-based estimate for the gap.
  int RefGap( ) const { return rgap_; }

  // Return associated CLinkBundle
  CLinkBundle Bundle( ) const { return clb_; }

  // Return (INT_MAX, 0) if there is no valid link-based estimate for the gap.
  pair<int,int> LinkGap( ) const { return make_pair( lgap_, ldev_ ); }
  
  // Estimate a gap size (and stdev), based on all the info available.
  pair<int,int> EstimatedGap( ) const;

  // Is this a Link-based gap?
  bool IsLinkGap() const { return (lgap_ != INT_MAX); }

  // Is this a Ref-based gap?
  bool IsRefGap() const { return (rgap_ != INT_MAX); }

  // Print info for one pair of oriented contigs.
  void Print( const bool brief, const int v, const int w, ostream &out ) const;
  
  // SELF_SERIALIZABLE methods (note that info_ will not be written/read).
  void writeBinary( BinaryWriter& writer ) const;

  void readBinary( BinaryReader& reader );
  
  static size_t externalSizeof( ) { return 0; }
  
  
private:

  const CProxInfo *info_;      // metadata for paired libraries, and heuristics
  vec< pair<int,int> > gaps_;  // link-based gaps, with tags
  pair<int,int> lwin_;         // window in gaps_ used to compute lgap_/ldev_
  int rgap_;                   // ref-based gap (or INT_MAX)
  int lgap_;                   // link-based gap (or INT_MAX)
  int ldev_;                   // dev of link-based gap (or 0)
  CLinkBundle clb_;	       // CLinkBundle for this link
};

SELF_SERIALIZABLE( CProx );    // See feudal/BinaryStreamTraits.h

#endif
