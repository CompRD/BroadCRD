// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef TILING_H
#define TILING_H

#include "PackAlign.h"
#include "Vec.h"
#include "tiled/PaddedSeq.h"
#include "tiled/TAlign.h"
#include "tiled/Tile.h"



/*
 * class tiling
 *
 * A tiling is a way to store informations about an Arachne contig.
 * Input is basically a set of alignments, either to a known contig,
 * or to the consensus bases.
 *
 * There are three different entities in play: known sequences, consensus
 * contigs, and plain reads. Tiling consists in aligning the reads onto
 * a master sequence (see later), taking pads into account (a simplified
 * ace representation).
 *
 * The master sequence can be either the known sequence, or the consensus
 * contig. Possible scenarios:
 *
 *  1. ( known_id_ >= 0 ): in this case master is the known sequence,
 *     and "aligns" means the aligns onto the known sequence. If an
 *     align for the contig is given, then it is interpreted as the
 *     align onto the known sequence as well.
 *
 *  2. ( known_id_ >= 0, and contig_id >=0 ): pretty much like 1., the
 *     difference being that in this case reads are aligned to the
 *     contig, and then the contig is aligned to known sequence (use
 *     either SetFromAligns or SetFromAlignsAlt to reproduce, respectively,
 *     case 1. or case 2.).
 *
 *  3. ( known_id_ < 0 ): master is the consensus contig, "aligns" is 
 *     interpreted as alignments to the contig, and of course an align
 *     for the contig is meaningless and should not be passed.
 *
 * Both tile known_ and tile contig_ may be empty, but not at the same time.
 *
 * Remark about the alignments. When Setting from Aligns keep in mind that
 * master is always oriented fw, but aligns must be passed so that pos1
 * and Pos1 denote begin and end of the alignment on the sequence aligning
 * to master; and pos2, Pos2 the begin and end of the alignment on the master.
 * This reflects the behavior of how nobbits are saved in EvaluateConsensus,
 * where id1 means Arachne and id2 means known_contig.
 */
class tiling {
  
public:
  
  tiling( );
  
  tiling( int known_id, int contig_id );

  void SetIds( int known_id, int contig_id );
  
  void ResetContigId( int contig_id );

  int SetFromAligns( const vec<t_align> &read_al, const t_align *cg_al = 0 );

  int SetFromAlignsAlt( const vec<t_align> &read_al, const t_align &cg_al );
  
  int KnownId( ) const { return known_id_; }

  int ContigId( ) const { return contig_id_; }

  int ReadsCount( ) const { return (int)reads_.size( ); }
  
  const tile &Known( ) const { return known_; }
  
  const tile &Contig( ) const { return contig_; }

  const tile &Read( int ii ) const { return reads_[ii]; }

  tile &Known( ) { return known_; }
  
  tile &Contig( ) { return contig_; }

  tile &Read( int ii ) { return reads_[ii]; }

  void Shift( int amount );

  friend istream &operator>> ( istream &in, tiling &the_tiles );
  
  friend ostream &operator<< ( ostream &out, const tiling &the_tiles );
    
  
private:
  
  void Adjust( const t_align &tile_align, int tile_id );
  
  void PadClone( int &pos, int tile_id );
  
  void PadRead( int &pos, int tile_id );
  
  void Compactify( );
  
  
private:
  
  int known_id_;      // known sequence id, or -1
  int contig_id_;     // contig id or -1
  int n_tiles_;       // total number of tiles
  tile known_;        // tile for known (may be empty)
  tile contig_;       // tile for contig (may be empty)
  vec<tile> reads_;   // tiles for reads
  
};



#endif
