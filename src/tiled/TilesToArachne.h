// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef TILES_TO_ARACHNE_H
#define TILES_TO_ARACHNE_H

#include "Basevector.h"
#include "Qualvector.h"
#include "ReadLocation.h"
#include "String.h"
#include "Vec.h"
#include "tiled/ChrMapper.h"
#include "tiled/ReadsTiling.h"



/*
 * class tiles_to_arachne
 *
 * Convert a PrimateContig-type assembly into Arachne.
 */
class tiles_to_arachne {
  
public:
  
  tiles_to_arachne( const String &tiles_dir,
		    const String &consensus_dir );

  void SaveChrMap( const String &map_file );
  
  void SaveFastbQualb( const String &bases_file,
		       const String &quals_file,
		       const String &pads_file );

  void SaveLocs( const String &contigs_file,
		 const String &reads_file,
		 const String &pads_file,
		 const String &locs_file );
  
  void SaveSuperInfo( const String &asupers_file );
  
  
private:
  
  void GenerateChrMap( );

  void CountReadsAndBases( );
  
  
private:
  
  String tiles_dir_;        // where contigs (as reads_tilings) are
  String consensus_dir_;    // where consensus fasta and qual are
  
  int n_reads_;                 // total number of reads
  longlong n_bases_;            // total number of valid bases
  vec<chr_mapper> chr_map_;     // contig to chromosome locations map
  
};



#endif
