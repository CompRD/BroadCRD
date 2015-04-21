///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "lookup/LookAlign.h"
#include "pairwise_aligners/KmerAligner.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "paths/FindVector.h"

/**
 * FindVectorCore
 */
void FindVectorCore( const int v_id,
		     const int r_id,
		     const bool v_rc,
		     const bvec &vbases,
		     const bvec &read,
		     const int MAX_MISMATCHES,
		     const int MAX_INDELS,
		     const double MAX_ER,
		     KmerAligner<8> &aligner8,
		     vec<look_align> &hits )
{
  hits.clear( );

  // Constants.
  const int band = Max( MAX_INDELS, 1 );

  // Find candidate placements.
  vec<int> pos;
  vec<int> ids;
  aligner8.FindPossibleAlignments( vbases, pos, ids );
      
  // Loop over all distinct placements.
  sort( pos.begin( ), pos.end( ) );
  pos.erase( unique( pos.begin( ), pos.end( ) ), pos.end( ) );
  for (size_t pos_id=0; pos_id<pos.size( ); pos_id++) {
    const int off = - pos[pos_id];

    align al;
    int err = 0;
    SmithWatBandedA( vbases, read, off, band, al, err );
    if ( al.Pos1( ) - al.pos1( ) < 8 ) continue;

    unsigned int v_len = vbases.size( );
    unsigned int r_len = read.size( );
    vector<int> muts_gap1_gap2 = al.MutationsGap1Gap2( vbases, read );
    int n_muts = muts_gap1_gap2[0];
    int n_ind =  muts_gap1_gap2[1] + muts_gap1_gap2[2];
    look_align theLA( v_id, r_id, v_len, r_len, v_rc, al, 0, n_muts, n_ind );
    if ( ! theLA.IsProper( ) ) continue;

    bool good1 = ( n_muts <= MAX_MISMATCHES && n_ind <= MAX_INDELS );
    bool good2 = ( theLA.ErrorRate( ) <= MAX_ER );
    if ( ! ( good1 || good2 ) ) continue;
    
    // Add hit to set.
    hits.push_back( theLA );
  }
  
}

/**
 * FindVector
 */
void FindVector( const String &hits_file,
		 const vecbvec &vectors,
		 const vecbvec &reads,
		 const int MAX_MISMATCHES,
		 const int MAX_INDELS,
		 const double MAX_ER,
		 const bool SKIP_RC,
		 ostream *log,
		 vec<look_align> *hits )
{
  if ( hits ) hits->clear( );

  // Output stream.
  String tmp_file = hits_file + ".tmp";
  ofstream out( tmp_file.c_str( ) );
  
  // Local version of vector (both fw and rc copies, if SKIP_RC is false).
  vecbvec all_vectors;
  vec<int> id_vectors;
  vec<bool> rc_vectors;
  all_vectors.reserve( vectors.size( ) );
  id_vectors.reserve( vectors.size( ) );
  rc_vectors.reserve( vectors.size( ) );
  for (size_t v_id=0; v_id<vectors.size( ); v_id++) {
    const bvec &vbases = vectors[v_id];
    bvec vbases_rc = vbases;
    vbases_rc.ReverseComplement( );
    bool palyndromic = ( vbases_rc == vbases );

    all_vectors.push_back( vbases );
    id_vectors.push_back( v_id );
    rc_vectors.push_back( false );
    if ( palyndromic || SKIP_RC ) continue;

    all_vectors.push_back( vbases_rc );
    id_vectors.push_back( v_id );
    rc_vectors.push_back( true );
  }

  // Loop over all reads.
  const size_t n_dotter = 1000000;
  const size_t n_reads = reads.size( );
  size_t n_done = 0;
  if ( log )
    *log << Date( ) << ": digesting "
	 << ToStringAddCommas( n_reads ) << " reads"
	 << endl;
  
  #pragma omp parallel for
  for (size_t r_id=0; r_id<reads.size( ); r_id++) {

    // Log progress.
    #pragma omp critical
    if ( log && n_done > 0 && n_done % n_dotter == 0 ) {
      double ratio = 100. * SafeQuotient( n_done, n_reads );
      *log << " " << ToStringAddCommas( n_done )
	   << " reads done (" << ToString( ratio, 1 )
	   << "%)" << endl;
    }
    
    // The kmer aligner for this read.
    const bvec &read = reads[r_id];
    KmerAligner<8> aligner8;
    aligner8.SetKmerStep( 1 );
    aligner8.SetBases( read );
    
    // Loop over all vectors (both fw and rc).
    for (size_t allv_id=0; allv_id<all_vectors.size( ); allv_id++) {
      const bvec &vbases = all_vectors[allv_id];
      const int v_id = id_vectors[allv_id];
      const int v_rc = rc_vectors[allv_id];
      vec<look_align> lochits;
      
      FindVectorCore( v_id, r_id, v_rc, vbases, read,
		      MAX_MISMATCHES, MAX_INDELS, MAX_ER, aligner8, lochits );
      
      #pragma omp critical
      for (size_t ii=0; ii<lochits.size( ); ii++)
	lochits[ii].PrintParseable( out );
    }

    // Update counter.
    #pragma omp critical
    n_done++;
  }
  
  // Close output stream.
  out.close( );

  // Reload and sort aligns (remove duplicates).
  if ( log ) *log << Date( ) << ": reloading and sorting aligns" << endl;
  vec<look_align> loc_qlts;
  vec<look_align> &qlts = hits ? *hits : loc_qlts;
  LoadLookAligns( tmp_file, qlts );

  order_lookalign_TargetQueryStartEnd sorter;
  sort( qlts.begin( ), qlts.end( ), sorter );
  qlts.erase( unique( qlts.begin( ), qlts.end( ) ), qlts.end( ) );

  if ( log ) *log << Date( ) << ": saving sorted aligns" << endl;
  WriteLookAligns( hits_file, qlts );
  if ( IsRegularFile( tmp_file ) ) Remove( tmp_file );

  // Done.
  if ( log ) *log << Date( ) << ": FindVector done" <<  endl;

}

