///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"

#include "Basevector.h"
#include "Histogram.h"
#include "Intvector.h"
#include "lookup/LookAlign.h"
#include "math/NStatsTools.h"
#include "util/RunCommand.h"
  
/**
 * EvalReadsOnRef
 *
 * Assess a set of reads against a reference: count mismatches,
 * insertions, deletions; and estimate the clean stretch distribution,
 * ie the distribution of the lengths of perfect matches.
 *
 * INPUT:
 *   <ALIGNS>      aligns of reads on reference
 *   <READS>       full path name of reads fastb
 *   <REF>         full path name of ref fastb
 *
 * OUTPUT:
 *   <OUT>.log     log file
 *   <OUT>.eps     eps of histogram of clean stretch
 *   
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( ALIGNS );
  CommandArgument_String( READS );
  CommandArgument_String( REF );
  CommandArgument_String( OUT );

  // Worst case scenario distance between adjacent errors (to reserve memory)
  CommandArgument_Double_OrDefault( EDIST, 5.0 );
  
  EndCommandArguments;
  
  
  // Dir and file names.
  String log_file = OUT + ".log";
  String eps_file = OUT + ".eps";
  
  vec<String> needed;
  needed.push_back( ALIGNS );
  needed.push_back( READS );
  needed.push_back( REF );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;

  ofstream log( log_file.c_str( ) );
  PrintCommandPretty( log );
  cout << Date( ) << "Sending log to " << log_file << "\n" << endl;

  // Load.
  log << Date( ) << ": loading reads and reference" << endl;
  vecbvec reads( READS );
  vecbvec ref( REF );
  
  log << Date( ) << ": loading aligns of reads on reference" << endl;
  vec<look_align> hits;
  LoadLookAligns( ALIGNS, hits);
  
  // Counters.
  size_t n_mis = 0;      // number of mismatches
  size_t n_ins = 0;      // number of insertions
  size_t n_del = 0;      // number of deletions
  vec<int> clean;        // vector of clean lengths
  
  // Back of the envelope computation to reserve mem for clean.
  size_t n_clean_estim = 0;
  for (size_t ii=0; ii<reads.size( ); ii++) 
    n_clean_estim += 1 + (int)( (double)reads[ii].size( ) / EDIST );
  String str_estim = ToStringAddCommas( n_clean_estim );
  log << Date( ) << ": estimating " << str_estim << " clean stretches" << endl;
  clean.reserve( n_clean_estim );

  // Keep track of reads that have already been seen.
  vec<bool> seen( reads.size( ), false );
  size_t n_multiple = 0;

  // Loop over all aligns.
  size_t dotter = 1000;
  log << Date( ) << ": parsing "
      << hits.size( ) << " aligns (. = "
      << dotter << " aligns)\n";
  for (size_t hit_id=0; hit_id<hits.size( ); hit_id++) {
    if ( hit_id % dotter == 0 ) Dot( log, hit_id / dotter );

    const look_align &hit = hits[hit_id];
    const align &al = hit.a;
    const int rid = hit.QueryId( );
    const int tid = hit.TargetId( );
    const bool rc = hit.Rc1( );

    // Tag read as seen (or as already analyzed).
    if ( seen[rid] ) {
      n_multiple++;
      continue;
    }
    seen[rid] = true;

    // Walk on align.
    bvec rc_read;
    if ( rc ) {
      rc_read = reads[rid];
      rc_read.ReverseComplement( );
    }
    
    const bvec &b1 = rc ? rc_read : reads[rid];
    const bvec &b2 = ref[tid];
    
    int p1 = al.pos1( );
    int p2 = al.pos2( );
    for (int block=0; block<al.Nblocks( ); block++) {
      int gap = al.Gaps( block );
      int len = al.Lengths( block );
      if ( gap < 0 ) {
	p1 += -gap;
	n_ins += -gap;
      }
      if ( gap > 0 ) {
	p2 += gap;
	n_del += gap;
      }
      int stretch = 0;
      for (int ii=0; ii<len; ii++) {
	if ( b1[p1] == b2[p2] ) stretch++;
	else {
	  if ( stretch > 0 ) clean.push_back( stretch );
	  stretch = 0;
	  n_mis++;
	}
	p1++;
	p2++;
      }
      if ( stretch > 0 ) clean.push_back( stretch );
    }
     
  } // loop over all aligns.
  log << endl;

  String str_nclean = ToStringAddCommas( clean.size( ) );
  log << Date( ) << ": found " << str_nclean << " clean stretches" << endl;

  // Report discarded aligns.
  if ( n_multiple > 0 ) {
    log << "\nWarning: " << n_multiple
	<< " aligns were discarded, as they belong to reads with more than\n"
	<< " one align (in these cases, one random align was selected, and\n"
	<< " the others were discarded)." << endl;
  }
  
  // Report.
  log << "\nCLEAN STRETCH LENGTHS STATISTICS\n\n";
  PrintBasicNStats( "length", clean, log );
  size_t n_matches = BigSum( clean );
  size_t n_tot = n_matches + n_mis + n_ins + n_del;
  
  vec< vec<String> > table;
  
  table.push_back( MkVec( String( "total" ),
			  ToStringAddCommas( n_tot ),
			  ToString( 100.0, 2 ) + "%" ) );

  double ratio = 100. * SafeQuotient( n_matches, n_tot );
  table.push_back( MkVec( String( "matches" ),
			  ToStringAddCommas( n_matches ),
			  ToString( ratio, 2 ) + "%" ) );
		   
  ratio = 100. * SafeQuotient( n_mis, n_tot );
  table.push_back( MkVec( String( "mismatches" ),
			  ToStringAddCommas( n_mis ),
			  ToString( ratio, 2 ) + "%" ) );
  
  ratio = 100. * SafeQuotient( n_ins, n_tot );
  table.push_back( MkVec( String( "insertions" ),
			  ToStringAddCommas( n_ins ),
			  ToString( ratio, 2 ) + "%" ) );

  ratio = 100. * SafeQuotient( n_del, n_tot );
  table.push_back( MkVec( String( "deletions" ),
			  ToStringAddCommas( n_del ),
			  ToString( ratio, 2 ) + "%" ) );

  log << "\nSUMMARY STATISTICS\n\n";
  PrintTabular( log, table, 3, "rrrr" );
  log << "\n"
      << "See " << eps_file << " for an histogram of the clean stretch\n"
      << endl;
  
  // Generate histogram of separations.
  histogram<int> histo;
  vec<int> bins;
  for (int ii=1; ii<51; ii++)
    bins.push_back( ii );
  histo.AddBins( bins );
  histo.AddData( clean.begin( ), clean.end( ) );
  {
    using namespace ns_psplot;
    
    vec<ns_psplot::freetext> labels( 1 );
    labels[0] = freetext( "Histogram of clean stretches", black, 16 );
    ofstream eps_out( eps_file.c_str( ) );
    histo.PrintAsEps( eps_out, labels, 0 );
    eps_out.close( );
  }
  
  // Done.
  String date = Date( );
  cout << Date( ) << ": EvalReadsOnRef done" << endl;
  log << Date( ) << ": EvalReadsOnRef done" << endl;
  log.close( );

}
