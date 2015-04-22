// Copyright (c) 2000-2003 Whitehead Institute of Biomedical Research
// 

#ifndef PADDEDREADLOCATION_H
#define PADDEDREADLOCATION_H

#include "Vec.h"
#include "String.h"
#include "system/System.h"
#include "Basevector.h"
#include "PaddedBasevector.h"


class padded_read_location {

public:
  padded_read_location();
  
  padded_read_location( const padded_read_location& RC );

  padded_read_location& operator=( const padded_read_location& rhs );

  ~padded_read_location() {}
  
  int getStartOnContig() const;
  int getEndOnContig() const;
  unsigned int getReadNumber() const;
  unsigned int getNumContigPads() const;
  String getReadName() const;
  bool getComplementation() const;
  vec<unsigned int> getReadPads() const;
 
  void setStartOnContig( const int start_on_contig );
  void setEndOnContig( const int end_on_contig );
  void setReadNumber( const unsigned int read_number );
  void setNumContigPads( const unsigned int contig_pads_when_aligned );
  void setReadName( const String read_name );
  void setComplementation( const bool complementation );
  void setReadPads( const vec<unsigned int> &read_pads );

  void addPad( unsigned int pad );
  char returnRC() const;

  PaddedBasevector* ConstructPaddedBasevector( const basevector& the_basevector ) const;

  String ConstructPaddedSequence( const String& the_sequence ) const;

  friend bool operator<(const padded_read_location& lhs, const padded_read_location& rhs);

  friend ostream& operator<< (ostream& out, const padded_read_location& rc );
  friend istream& operator>> (istream& in, padded_read_location& rc );

private:
  int start_on_contig_, end_on_contig_;
  unsigned int contig_pads_when_aligned_, read_number_;
  String name_;
  bool complementation_;
  vec<unsigned int> pads_;

};

inline int padded_read_location::getStartOnContig() const
{  return start_on_contig_; }

#endif

