// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef PROTO_CONTIG_MAP_H
#define PROTO_CONTIG_MAP_H

#include "tiled/ProtoContigId.h"
#include "String.h"
#include "Vec.h"



/*
 * class proto_contig_map
 *
 * Identifies where a contig_id comes from. It is basically a map from
 * (chromosome_id, contig_pos_in_chromosome) to contig_id, to be used
 * with the OnMaster code (chimp on human, etc.).
 */
class proto_contig_map {

public:

  proto_contig_map( );
  
  void Clear( );
  
  void Reserve( int size );
  
  void Add( const proto_contig_id &proto_id );
  
  void Add( int cg_id, int chr_id, int pos );

  int Find( int chr_id, int pos );

  void Save( const String &out_file );

  void Load( const String &in_file );

  int Size( ) const;
  
  const proto_contig_id &operator[] ( int jj ) const;
  
  
private:
  
  void CheckData( );


private:
  
  vec<proto_contig_id> proto_id_;  // contig id
  bool checked_;                   // switch to check data intergrity  

};



#endif
