///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "PairsManager.h"
#include "Superb.h"
#include "Vec.h"
#include "paths/Alignlet.h"
#include "paths/reporting/AllContigLinksCore.h"
#include "paths/reporting/CBundle.h"
#include "paths/reporting/CSuperLinks.h"

/**
 * AllContigLinksCore
 */
void AllContigLinksCore( const int MIN_LINKS,
			 const vec<fastavector> &contigs,
			 const vec<alignlet> &aligns,
			 const vec<int> &index,
			 const PairsManager &pairs,
			 ostream &out,
			 ostream &log )
{
  // Build subset of supers (new super_id == old contig_id).
  log << Date( ) << ": initializing scaffolds" << endl;
  vec<superb> supers;
  supers.resize( contigs.size( ) );
  for (int tig=0; tig<contigs.isize( ); tig++) {
    supers[tig].SetNtigs( 1 );
    supers[tig].SetTig( 0, tig );
    supers[tig].SetLen( 0, contigs[tig].size( ) );
  }
  
  // All links.
  CSuperLinks linker( &pairs, &supers, &aligns, &index );
  
  // Warning for the links file.
  out << "Each bundle (line) represents the offset between the two contigs'\n"
      << "starting points (not than the gap between them).\n" << endl;
  
  // Loop over all supers.
  for (int super_id=0; super_id<supers.isize( ); super_id++) {
    vec<COffset> all_offsets;
    linker.AllLinks( super_id, all_offsets );

    // Loop over all linked pairs (oriented supers linked to super_id).
    for (int set_id=0; set_id<all_offsets.isize( ); set_id++) {
      const COffset &offset = all_offsets[set_id];
      int id1 = offset.Super1( );
      int id2 = offset.Super2( );
      bool rc2 = offset.Rc2( );
      
      // Loop over all clusters of the linked pair.
      for (size_t cid=0; cid<offset.NClusters( ); cid++) {
	int weight = offset.NLinks( cid );
	if ( weight < MIN_LINKS ) continue;
	float score = offset.Score( cid );
	NormalDistribution nd = offset.Offset( cid );
	int offm = int( nd.mu_ );
	int offs = Max( 1, int( nd.sigma_ ) );
	pair<int,int> off_pair = make_pair( offm, offs );
	pair<int,int> h1 = offset.SpreadWinBases1( cid );
	pair<int,int> h2 = offset.SpreadWinBases2( cid );
	
	CBundle bundle( id1, id2, rc2, weight, score, off_pair, h1, h2 );
	bundle.Print( contigs[id1].size( ), contigs[id2].size( ), out );
      }
    }
  } // loop over all supers

  out << "\nEOF" << endl;
}
