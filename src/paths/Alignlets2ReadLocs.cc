///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "PairsManager.h"
#include "Vec.h"
#include "VecUtilities.h"
#include "feudal/BinaryStream.h"
#include "paths/Alignlet.h"
#include "paths/Alignlets2ReadLocs.h"
#include "paths/ReadLoc.h"

/**
 * AppendReadLocs
 */
void AppendReadLocs( const SAlignletsData &aligns_data,
		     const int readClass,
		     vec<read_loc> &locs,
		     ostream *log )
{
  // Log stream.
  ofstream devnull ( "/dev/null" );
  ostream &out = log ? *log : devnull;

  // Descriptors.
  String rt = "";
  if ( readClass == 0 ) rt = "frags";
  else if ( readClass == 1 ) rt += "jumps";
  else if ( readClass == 2 ) rt = "Jumps";
  else ForceAssert( 1 == 0 );   // unsupported type.
  String descr = "AppendReadLocs[" + rt + "]";
  
  // Nothing to do.
  if ( aligns_data.IsEmpty( ) ) {
    out << Date( ) << ": " << descr << " - no " << rt << " reads found" << endl;
    return;
  }

  // Core data.
  const PairsManager &pairs = *( aligns_data.pairs_ );
  const vec<alignlet> &aligns = *( aligns_data.aligns_ );
  const vec<int> &index = *( aligns_data.index_ );
    
  // Map reads onto target.
  out << Date( ) << ": " << descr << " - mapping reads to target" << endl;
  vec< triple<int,int,int64_t> > helper;   // cg_id, start_on_cg, read_id
  helper.reserve( aligns.size( ) );
  for (int64_t read_id=0; read_id<(int64_t)index.size( ); read_id++) {
    if ( index[read_id] < 0 ) continue;
    int cid = aligns[ index[read_id] ].TargetId( );
    int start = aligns[ index[read_id] ].pos2( );
    helper.push_back( triple<int,int,int64_t>( cid, start, read_id ) );
  }
  sort( helper.begin( ), helper.end( ) );
  
  // Reserve memory.
  locs.reserve( helper.size( ) );

  // Generate vector of read_locs.
  out << Date( ) << ": " << descr << " - generating read_locs" << endl;
  for (size_t ii=0; ii<helper.size( ); ii++) {
    const alignlet &al = aligns[ index[ helper[ii].third ] ];
    int64_t rid = helper[ii].third;
    int64_t pairid = pairs.getPairID( rid );
    int64_t p_id =
      pairs.ID1( pairid ) == rid
      ? pairs.ID2( pairid )
      : pairs.ID1( pairid );
    const alignlet *p_al = ( index[p_id] < 0 ) ? 0 : &( aligns[ index[p_id] ] );
    int32_t cid = al.TargetId( );
    int32_t p_cid = p_al ? p_al->TargetId( ) : 0;
    int cpos = al.pos2( );
    int p_cpos = p_al ? p_al->pos2( ) : 0;
    Bool fw = al.Fw1( );
    Bool p_fw = p_al ? p_al->Fw1( ) : False;
    int8_t rclass = readClass;
    uint8_t lib_id = pairs.libraryID( pairid );
    uint16_t rlen = al.Pos2( ) - al.pos2( );
    uint16_t p_rlen = p_al ? p_al->Pos2( ) - p_al->pos2( ) : 0;
    uint8_t band = 0;    // bandwidth - undefined!
    uint8_t p_band = 0;    // bandwidth - undefined!
    bool p_placed = p_al;

    locs.push_back( read_loc( rid, p_id, cid, p_cid, cpos, p_cpos,
			      fw, p_fw, rclass, lib_id, rlen, p_rlen,
			      band, p_band, p_placed ) );
  }
  
  // Done.
  out << Date( ) << ": " << descr << " - done" << endl;
  
}

/**
 * Alignlets2ReadLocs
 */
void Alignlets2ReadLocs( const SAlignletsData &frags,
			 const SAlignletsData &jumps,
			 const SAlignletsData &Jumps,
			 vec<read_loc> &locs,
			 ostream *log )
{
  locs.clear( );
  
  AppendReadLocs( frags, 0, locs, log );
  AppendReadLocs( jumps, 1, locs, log );
  AppendReadLocs( Jumps, 2, locs, log );

  if ( log ) *log << Date( ) << ": Alignlets2ReadLocs - sorting locs" << endl;  
  sort( locs.begin( ), locs.end( ) );

  if ( log ) *log << Date( ) << ": Alignlets2ReadLocs done" << endl;  
}
