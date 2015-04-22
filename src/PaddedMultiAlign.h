// Copyright (c) 2000-2003 Whitehead Institute of Biomedical Research
// 

#ifndef PADDEDMULTIALIGN_H
#define PADDEDMULTIALIGN_H

#include "Alignment.h"
#include "CompressedSequence.h"
#include "PaddedBasevector.h"
#include "PaddedReadLocation.h"
#include "Vec.h"

#include <map>

#include <iosfwd>

class PaddedMultiAlign {

public:
  PaddedMultiAlign() {};
  PaddedMultiAlign( int contig_id,
		    vec<unsigned int>& contig_pads,
		    vec<padded_read_location>& reads );
  
  friend
  ostream& operator<< ( ostream& out, const PaddedMultiAlign& pma );

  friend
  istream& operator>> ( istream& in, PaddedMultiAlign& pma );

  void setContigId( int contig_id );
  int getContigId() const;

  void setContigPads( const vec<unsigned int>& contig_pads );
  vec<unsigned int> getContigPads() const;
  
  // return a sorted list of read ids
  vec<int> getReadIds() const;

  void setReads( const vec<padded_read_location>& reads );
  const vec<padded_read_location>& getReads() const;

  void print( ostream &out, 
              const basevector &contig_bases, 
              const veccompseq &read_bases,
              const vec<int>& left_trims,
              const vec<int>& right_trims ) const;

private:
  int contig_id_;
  vec<unsigned int> contig_pads_;
  vec<padded_read_location> reads_;
};



class PaddedMultiAlignBuilder {

public:
  PaddedMultiAlignBuilder() { }

  void getPaddedMultiAlign( PaddedMultiAlign& pma,
			    const int contig_id,
			    const PaddedBasevector& padded_contig );
  
  void addPaddedReadLocation( padded_read_location &rc );

  vec<padded_read_location>* getPaddedReadLocations();

  // ---- methods to update the sequences after alignment
  bool CheckIfNeedRealignment( const padded_read_location &realign, 
			       const vec<unsigned int> &pads_in_contig );

  void addNewPadsToRead( const PaddedBasevector &padded_contig, 
			 padded_read_location &realign );

  void UpdatePads( PaddedBasevector &padded_read, 
		   PaddedBasevector &padded_consensus, 
		   vec<char> &padded_consensus_sequence,
		   alignment &the_alignment );

private:
  vec<padded_read_location> realign_candidates_;

};

#endif
