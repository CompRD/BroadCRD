///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoverageAnalyzer.h"
#include "Fastavector.h"
#include "math/NStatsTools.h"
#include "reporting/ScaffoldedFastaToNStats.h"

/**
 * ScaffoldedFastaToNStats
 */
void ScaffoldedFastaToNStats( ostream &out,
			      const int GAPN,
			      const vec<fastavector> &fasta )
{
  vec<int> fasta_lens( fasta.size( ), 0 );
  for (int ii=0; ii<fasta.isize( ); ii++)
    fasta_lens[ii] = (int)fasta[ii].size( );

  // Find all gaps (loop over all sequences in fasta).
  vec<seq_interval> gap_si;
  for (int seq_id=0; seq_id<fasta.isize( ); seq_id++) {
    const fastavector &seq = fasta[seq_id];
    
    // Loop over all bases in sequence.
    int pos = 0;
    while( pos < (int)seq.size( ) ) {

      // Look for first N.
      while( pos < (int)seq.size( ) ) {
	if ( seq[pos] == 'n' || seq[pos] == 'N' ) break;
	pos++;
      }
      if ( pos == (int)seq.size( ) ) break;

      // Count Ns.
      int begin = pos;
      while( pos < (int)seq.size( ) ) {
	if ( seq[pos] != 'n' && seq[pos] != 'N' ) break;
	pos++;
      }
      int len = pos - begin;
      
      // Not a gap.
      if ( len < GAPN ) continue;

      // Add stretch to list.
      gap_si.push_back( seq_interval( gap_si.isize( ), seq_id, begin, pos ) );
    }
  }
  
  // Now gaps have coverage 1, contigs have coverage 0.
  CoverageAnalyzer cov( gap_si, &fasta_lens );

  vec<seq_interval> contig_si;
  cov.GetCoveragesExactly( 0, contig_si );

  // Collect contig and gap length.
  vec<int> cg_lens;
  vec<int> gap_lens;
  cg_lens.reserve( contig_si.size( ) );
  for (int ii=0; ii<contig_si.isize( ); ii++)
    cg_lens.push_back( contig_si[ii].Length( ) );
  gap_lens.reserve( gap_si.size( ) );
  for (int ii=0; ii<gap_si.isize( ); ii++)
    gap_lens.push_back( gap_si[ii].Length( ) );

  // Print stats.
  PrintBasicNStats( "scaffolds", fasta_lens, out );
  out << endl;

  PrintBasicNStats( "contigs", cg_lens, out );
  out << endl;

  PrintBasicNStats( "gaps", gap_lens, out );
  out << endl;
  
}
			      

