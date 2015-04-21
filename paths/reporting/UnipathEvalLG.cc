///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"
#include "CoverageAnalyzer.h"
#include "PrettyPrintTable.h"
#include "SeqInterval.h"
#include "lookup/LookAlign.h"
#include "math/NStatsTools.h"

/**
 * UnipathEvalLG
 *
 * Evaluation took for unipaths, designed for use with the RunAllPathsLG
 * pipeline.
 * 
 * If you set USE_TRUTH=True, UnipathEvalLG will also look for alignments to
 * reference, which you should have already generated.
 *
 * REQUIRED INPUT:
 *   <HEAD>.unibases.k<K>          (fastb of unibases)
 *
 * OPTIONAL INPUT:
 *   <HEAD>.unibases.k<K>.qltout   (alignments of unibases on reference)
 *   <REF_FASTB>                   (fastb of reference genome)
 *
 * OUTPUT:
 *   sent to cout
 *
 ******************************************************************************/


// TODO: add EvalAdjGraph, EvalPairSeparations

void EvalNStats( const int &K, const vecbvec &unibases, ostream &out );

void EvalRefAlignment( const int &K,
		       const vecbvec &unibases,
		       const vecbvec &genome,
		       const vec<look_align_plus> &hits,
		       ostream &out );






int main( int argc, char *argv[] )
{
  RunTime( );
  BeginCommandArguments;
  CommandArgument_Int( K );
  CommandArgument_String_Doc( HEAD, "Look for unibases at HEAD.unibases.k<K>" );
  CommandArgument_Bool_Doc( USE_TRUTH, "If True, alignments must exist at file HEAD.unibases.K<k>.qltout" );
  CommandArgument_String_OrDefault_Doc( REF_FASTB, "", "Genome fastb file.  Must specify if USE_TRUTH=True." );
  EndCommandArguments;

  // Check that inputs make sense.
  if ( USE_TRUTH && REF_FASTB == "" )
    FatalErr( "If USE_TRUTH=True, you must set REF_FASTB." );
  
  // File names.
  String strK = ToString( K );
  String unibases_file = HEAD + ".unibases.k" + strK;
  String hits_file = unibases_file + ".qltout";

  // Load files.
  // Some data structues are only loaded from file if USE_TRUTH=True;
  // otherwise, they remain empty, which is okay because they are not used.
  cout << Date( ) << ": loading unibases" << endl;
  vecbasevector unibases( unibases_file );
  vecbasevector genome;
  vec<look_align_plus> hits;
  
  if ( USE_TRUTH ) {
    cout << Date( ) << ": loading hits" << endl;
    LoadLookAlignPlus( hits_file, hits );
    
    cout << Date( ) << ": loading reference fastb" << endl;
    genome.ReadAll( REF_FASTB );
  }
    
  // Report ranges.
  EvalNStats( K, unibases, cout );

  // Report on reference.
  if ( USE_TRUTH )
    EvalRefAlignment( K, unibases, genome, hits, cout );

  // Done.
  cout << Date( ) << ": Done!" << endl;
}



/************************
 *                      *
 * FUNCTIONS START HERE *
 *                      *
 ***********************/



void PrintSeqInts( const vec<seq_interval> &sints, ostream &out )
{
  vec<int> lens;
  lens.reserve( sints.size( ) );
  for (uint ii=0; ii<sints.size( ); ii++)
    lens.push_back( sints[ii].End( ) - sints[ii].Begin( ) );

  PrintBasicNStats( "length", lens, out );
  out << endl;
}

void EvalNStats( const int &K, const vecbvec &unibases, ostream &out )
{
  int n_unibases = unibases.size( );
  out << "\n" << Date( ) << ": REPORT ON N-STATS" << endl;
  // Translate to kmer space.
  vec<int> lens;
  lens.reserve( n_unibases );
  for (int ii=0; ii<n_unibases; ii++)
    lens.push_back( (int)unibases[ii].size( ) - (int)( K - 1 ) );
  longlong tot_len = BigSum( lens );

  // Nstats of unipath lengths.
  out << "There are " << n_unibases << " unipaths in the data set\n"
      << "Largest unipath: " << Max( lens ) << " kmers\n"
      << "Total kmer length: " << tot_len << "\n\n"
      << "N-stats for unipath lengths (in kmers, where K=" << K
      << "):\n";
  PrintBasicNStats( "length", lens, out );
  out << endl;

  // Bins and data (number of unipaths, number of kmers).
  int zeros = int( ceil( log10( Max( lens ) ) ) );
  vec<int> bins( 1, 10 );
  while ( bins.isize( ) < zeros )
    bins.push_back( 10 * bins.back( ) );

  vec<longlong> n_unipaths( bins.size( ), 0 );
  vec<longlong> n_kmers( bins.size( ), 0 );

  for (uint ii=0; ii<lens.size( ); ii++) {
    int bin_id = bins.isize( ) - 1;
    while ( bin_id > 0 && lens[ii] <= bins[bin_id-1] ) bin_id--;
    n_unipaths[bin_id] += 1;
    n_kmers[bin_id] += lens[ii];
  }
  
  // Table with number of unipaths, and total kmers.
  vec<String> str_ranges;
  {
    vec< vec<String> > ranges;
    vec<String> aline;
    for (int ii=0; ii<bins.isize( ); ii++) {
      aline.clear( );
      aline.push_back( ii == 0 ? "1" : ToString( bins[ii-1]+1 ) );
      aline.push_back( "-" );
      aline.push_back( ToString( bins[ii] ) );
      ranges.push_back( aline );
    }
    vec<Bool> justify( 3, True );
    BeautifyTable( ranges, &justify );

    for (int ii=0; ii<ranges.isize( ); ii++) {
      String range_line = ToString( ranges[ii][0] );
      for (int jj=1; jj<ranges[ii].isize( ); jj++)
	range_line += "  " + ToString( ranges[ii][jj] );
      str_ranges.push_back( range_line );
    }
  }

  vec< vec<String> > table;
  vec<String> aline;
  aline.push_back( "range of sizes" );
  aline.push_back( "# of unipaths" );
  aline.push_back( "total kmers" );
  table.push_back( aline );

  for (int ii=0; ii<bins.isize( ); ii++) {
    aline.clear( );
    aline.push_back( str_ranges[ii] );
    aline.push_back( ToString( n_unipaths[ii] ) );
    aline.push_back( ToString( n_kmers[ii] ) );
    table.push_back( aline );
  }

  vec<Bool> justify( 3, True );
  BeautifyTable( table, &justify );

  for (int ii=0; ii<table.isize( ); ii++) {
    for (int jj=0; jj<table[ii].isize( ); jj++)
      out << table[ii][jj] << "     ";
    out << "\n";
  }
  out << endl;
  
}

void EvalRefAlignment( const int &K,
		       const vecbvec &unibases,
		       const vecbvec &genome,
		       const vec<look_align_plus> &hits,
		       ostream &out )
{
  out << "\n" << Date( ) << ": REPORT ON REFERENCE ALIGNMENT" << endl;
  // Info on reference genome.
  vec<int> gklens;
  longlong genome_totklen = 0;
  longlong genome_totblen = 0;
  for (size_t ii=0; ii<genome.size( ); ii++) {
    gklens.push_back( genome[ii].size( ) - ( K - 1 ) );
    genome_totklen += (int)genome[ii].size( ) - (int)( K - 1 );
    genome_totblen += genome[ii].size( );
  }

  out << "There are " << genome.size( )
      << " contigs in the reference, for a total of " << genome_totklen
      << " kmers,\ncorresponding to " << genome_totblen
      << " bases\n" << endl;

  // Find covered/uncovered regions (in kmer space!)
  vec<seq_interval> si_covered_100;
  vec<seq_interval> si_covered_10;
  vec<seq_interval> si_covered;
  vec<seq_interval> si_uncovered;  
  {
    vec<seq_interval> si_hits;
    vec<seq_interval> si_hits_10;
    vec<seq_interval> si_hits_100;
    si_hits.reserve( hits.size( ) );
    si_hits_10.reserve( hits.size( ) );
    si_hits_100.reserve( hits.size( ) );
    for (uint ii=0; ii<hits.size( ); ii++) {
      const  look_align_plus &hit = hits[ii];
      int seq_id = hit.target_id;
      int int_id = hit.query_id;
      int begin = hit.a.pos2( );
      int end = hit.a.Pos2( ) - ( K - 1 );
      if ( end - begin < 1 ) continue;
      
      seq_interval newsi( int_id, seq_id, begin, end );

      si_hits.push_back( newsi ); 
      if ( end - begin >= 100 ) si_hits_100.push_back( newsi );
      if ( end - begin >= 10 ) si_hits_10.push_back( newsi );
    }
    
    CoverageAnalyzer cov( si_hits, &gklens );
    cov.GetCoveragesAtLeast( 1, si_covered );
    cov.GetCoveragesExactly( 0, si_uncovered );

    CoverageAnalyzer cov10( si_hits_10, &gklens );
    cov10.GetCoveragesAtLeast( 1, si_covered_10 );
    
    CoverageAnalyzer cov100( si_hits_100, &gklens );
    cov100.GetCoveragesAtLeast( 1, si_covered_100 );
  }    

  out << "N-stats of reftigs (in kmers, where K=" << K << "):\n";
  PrintSeqInts( si_covered, cout );

  out << "N-stats of gaps between reftigs (in kmers, where K=" << K << "):\n";
  PrintSeqInts( si_uncovered, cout );

  // Base coverage statistics.
  longlong basecov_100 = 0;
  longlong basecov_10 = 0;
  longlong basecov_all = 0;

  {
    for (int pass=0; pass<3; pass++) {
      const vec<seq_interval> *sints = 0;
      longlong *bcov = 0;
      
      if ( pass == 0 ) {
	sints = &si_covered_100;
	bcov = &basecov_100;
      }
      else if ( pass == 1 ) {
	sints = &si_covered_10;
	bcov = &basecov_10;
      }
      else {
	sints = &si_covered;
	bcov = &basecov_all;
      }

      for (uint ii=0; ii<sints->size( ); ii++)
	*bcov += ((*sints)[ii]).End( ) - ((*sints)[ii]).Begin( );
    }
  }

  // Report base coverage.
  float ratio_100 = 100.0 * SafeQuotient( basecov_100, genome_totklen );
  float ratio_10 = 100.0 * SafeQuotient( basecov_10, genome_totklen );
  float ratio_all = 100.0 * SafeQuotient( basecov_all, genome_totklen );
  String str_100 = ToString( ratio_100, 2 );
  String str_10 = ToString( ratio_10, 2 );
  String str_all = ToString( ratio_all, 2 );

  out << "Overall base coverage statistics (sizes in kmers):\n"
      << " coverage by unipaths of size >= 100: " << str_100 << "%\n"
      << " coverage by unipaths of size >= 10: " << str_10 << "%\n"
      << " coverage by all unipaths: " << str_all << "%\n"
      << endl;
  
  // core rows
  //   0: no errors
  //   1: < 0.1% error rate
  //   2: < 1% error rate
  //   3: < 10% error rate 
  //   4: all
  // core columns
  //   0: >= 100 kmers
  //   1: >= 10 kmers
  //   2: all
  vec< vec<int> > core;
  core.resize( 5 );
  for (int ii=0; ii<5; ii++)
    core[ii].resize( 3 );
  
  // For each unibase, only consider best placement (lowest error rate).
  uint n_unibases = unibases.size( );
  vec<float> pcerr_rates( n_unibases, -1.0 );
  vec<int> klens( n_unibases, 0 );
  
  for (uint ii=0; ii<hits.size( ); ii++) {
    int id = hits[ii].query_id;
    int klen = hits[ii].a.Pos2( ) - hits[ii].a.pos2( ) - ( K - 1 );
    if ( klen < 1 ) continue;
    klens[id] = klen;
    
    float rate = 100.0 * hits[ii].ErrorRate( );
    if ( pcerr_rates[id] < 0 ) pcerr_rates[id] = rate;
    else pcerr_rates[id] = Min( pcerr_rates[id], rate );
  }
  
  // Build core (non cumulative yet).
  for (uint ii=0; ii<pcerr_rates.size( ); ii++) {
    float pcerr = pcerr_rates[ii];
    int klen = klens[ii];
    if ( klen < 1 ) continue;
    
    int row_id = 0;
    if ( pcerr == 0.0 ) row_id = 0;
    else if ( pcerr < 0.1 ) row_id = 1;
    else if ( pcerr < 1.0 ) row_id = 2;
    else if ( pcerr < 10.0 ) row_id = 3;
    else row_id = 4;

    int col_id = 0;
    if ( klen >= 100 ) col_id = 0;
    else if ( klen >= 10 ) col_id = 1;
    else col_id = 2;

    core[row_id][col_id] += 1;
  }
  
  // Adjust columns and rows (cumulative).
  for (int row_id=0; row_id<5; row_id++)
    for (int col_id=1; col_id<3; col_id++)
      core[row_id][col_id] += core[row_id][col_id-1];
  
  for (int col_id=0; col_id<3; col_id++)
    for (int row_id=1; row_id<5; row_id++) 
      core[row_id][col_id] += core[row_id-1][col_id];

  // Make core into a table (as percent of total unibases).
  vec< vec<String> > table;
  vec<String> aline;

  aline.push_back( "" );
  aline.push_back( "length>=100" );
  aline.push_back( "length>=10" );
  aline.push_back( "all lengths" );
  table.push_back( aline );

  for (int row_id=0; row_id<5; row_id++) {
    aline.clear( );

    if ( row_id == 0 ) aline.push_back( "no errors" );
    else if ( row_id == 1 ) aline.push_back( "e.r. <0.1%" );
    else if ( row_id == 2 ) aline.push_back( "e.r. <1%" );
    else if ( row_id == 3 ) aline.push_back( "e.r. <10%" );
    else aline.push_back( "all" );

    for (int col_id=0; col_id<3; col_id++) {
      double ratio = SafeQuotient( (uint)core[row_id][col_id], n_unibases );
      aline.push_back( ToString( 100.0 * ratio, 2 ) );
    }

    table.push_back( aline );
  }
  
  vec<Bool> just( 4, True );

  out << "Quality of alignments onto reference:\n";
  BeautifyAndPrintTable( table, out, "     ", &just );
  out << "Legend:\n"
      << "  [1] e.r. = error rate\n"
      << "  [2] all lengths are in kmers\n"
      << "  [3] results in percent (against a total of " << n_unibases
      << " unibases)\n" << endl;
  
}

