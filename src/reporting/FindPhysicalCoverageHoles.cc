///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/* FindPhysicalCoverageHoles
 *
 * Take a set of read pairs, typically jumping reads, possibly consisting of
 * many libraries.  Find holes in physical coverage of the reference.
 * ("Physical coverage" is defined as coverage by an entire insert, not just
 * the two reads but also the space between them.)
 * The reads' alignments to reference must be pre-computed.
 *
 *
 * INPUTS:
 *   <PRE>/<DATA>/<GENOME>.fastb
 *   <PRE>/<DATA>/<READS>.pairs
 *   <PRE>/<DATA>/<READS>.qltout
 *
 * REPORTS:
 * -- Finds the number of bases in the reference that are uncovered.
 * -- Reports on the size of physical coverage holes.
 * -- Breaks down physical coverage by library and reports how much coverage
 *    is missing from each library.
 * 
 *
 *
 * NOTE: The following dataset sizes will cause integer overflow.
 * -- More than 2G reads (LookAlign)
 * -- A reference genome with a single contig longer than 2GB (CoverageAnalyzer)
 *
 *
 *
 *
 * Josh Burton
 * April 2010
 *
 ******************************************************************************/



#include "MainTools.h"
#include "CoverageAnalyzer.h"
#include "Intvector.h"
#include "PairsManager.h"
#include "lookup/LookAlign.h"


/* Print a pretty chart with information about each library in the PairsManager:
 * its name, its size, and the amount of reference coverage it misses.
 *
 ******************************************************************************/
void
ReportMissingCoverageByLibrary( ostream & out,
				const PairsManager & pairs,
				const vec<size_t> & missing_by_lib,
				const size_t genome_size,
				const size_t genome_uncovered_size )
{
  
  // Make a chart header line.
  vec< vec<String> > rows;
  vec<String> row1, row2;
  row1.push_back( "Library", "N pairs", "Missing", "(%)", "MissingHere", "(%)" );
  row2.push_back( "-------", "-------", "-------", "---", "-----------", "---" );
  rows.push_back( row1, row2 );
  
  
  vec<String> lib_names = pairs.getLibraryNames();
  vec<size_t> lib_sizes = pairs.getLibrarySizes();
  
  // Make a line of info for each library.
  for ( size_t i = 0; i < pairs.nLibraries(); i++ ) {
    vec<String> row;
    size_t missing_here = missing_by_lib[i] - genome_uncovered_size;
    double pctA = double( 100.0 * missing_by_lib[i] ) / genome_size;
    double pctB = double( 100.0 * missing_here ) / ( genome_size - genome_uncovered_size );
    
    row.push_back( lib_names[i],
		   ToString( lib_sizes[i] ),
		   ToString( missing_by_lib[i] ),
		   ToString( pctA, 2 ),
		   ToString( missing_here ),
		   ToString( pctB, 2 ) );
    rows.push_back( row );
  }
  
  
  // Print the chart!
  out << endl << "MISSING COVERAGE BY LIBRARY" << endl << endl;
  PrintTabular( out, rows, 4, "lrrrrr" );
  out << endl;
  
  // Print some explanation.
  out << "'Missing': bases with less than MIN_COVERAGE by pairs in this library." << endl;
  out << "'MissingHere': bases with less than MIN_COVERAGE by pairs in this library but at least MIN_COVERAGE over all libraries." << endl;
  
  out << endl;
}







int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String_OrDefault( GENOME, "genome" );
  CommandArgument_String_OrDefault( READS, "jump_reads_extract" );
  CommandArgument_Int_OrDefault_Doc( MIN_COVERAGE, 3, "A 'coverage hole' is defined as a region of the genome with physical coverage less than MIN_COVERAGE" );
  CommandArgument_Int_OrDefault_Doc( MAX_JUMP, 1e5, "Jumping reads further apart than this are assumed to be mis-aligned" );
  EndCommandArguments;
  
  
  
  
  String data_dir = PRE + "/" + DATA;
  String reads_head = data_dir + "/" + READS;
  
  
  // Load files.
  vec<int> genome_lens;
  cout << Date() << ": Loading genome (to get its contig sizes only)" << endl;
  {
    vecbasevector genome( data_dir + "/" + GENOME + ".fastb" );
    genome.ElementSizes( genome_lens );
  }
  size_t genome_size = BigSum( genome_lens );
  
  
  cout << Date() << ": Loading PairsManager" << endl;
  PairsManager pairs( reads_head + ".pairs" );
  size_t n_libs = pairs.nLibraries();
  cout << endl;
  pairs.printLibraryStats( cout );
  cout << endl;
  
  cout << Date() << ": Loading LookAligns" << endl;
  vec<look_align> aligns;
  LoadLookAligns( reads_head + ".qltout", aligns );
  
  // Make aligns index.
  cout << Date() << ": Making alignment index" << endl;
  VecLongVec aligns_index( pairs.nReads() );
  for ( size_t i = 0; i < aligns.size(); i++ )
    // If the following line causes a seg fault, the input pairs and aligns
    // file are mismatched.
    aligns_index[ aligns[i].QueryId() ].push_back(i);
  
  
  
  size_t n_pairs_total = pairs.nPairs();
  size_t n_pairs_CN1 = 0;
  size_t n_pairs_same_contig = 0;
  size_t n_pairs_well_aligned = 0;
  
  
  
  cout << Date() << ": Finding each read pair's true physical coverage" << endl;
  vec<seq_interval> intervals;
  vec< vec<seq_interval> > intervals_by_lib( n_libs );
  
  // For each pair, create a seq_interval representing that pair's
  // physical coverage of the reference genome.
  for ( size_t i = 0; i < n_pairs_total; i++ ) {
    int64_t read1 = pairs.ID1(i);
    int64_t read2 = pairs.ID2(i);
    
    // Require both reads to align uniquely.  If the aligns file was generated
    // from a BAM file via SAM2CRDDump, the reads will never align multiply.
    if ( aligns_index[read1].size() != 1 ||
	 aligns_index[read2].size() != 1 ) continue;
    n_pairs_CN1++;
    
    // Require both reads to fall on the same contig in the reference genome.
    const look_align & a1 = aligns[ aligns_index[read1][0] ];
    const look_align & a2 = aligns[ aligns_index[read2][0] ];
    if ( a1.TargetId() != a2.TargetId() ) continue;
    n_pairs_same_contig++;
    
    // Find the portion of genome covered by this pair.
    // Note that "innie" jumping reads will cover a smaller region.
    int start = Min( a1.pos2(), a2.pos2() );
    int stop  = Max( a1.Pos2(), a2.Pos2() );
    
    // If the physical coverage region appears too long, one of the reads is
    // probably mis-aligned.  Skip the pair.
    if ( stop - start > MAX_JUMP ) continue;
    n_pairs_well_aligned++;
    
    // We have a good pair!  Create the seq_interval and fill data structures.
    seq_interval s( i, a1.TargetId(), start, stop );
    intervals.push_back( s );
    intervals_by_lib[ pairs.libraryID(i) ].push_back( s );
  }
  
  
  
  // Report on the pair alignments.
  cout << endl;
  cout << n_pairs_total << "\tTotal pairs" << endl;
  cout << n_pairs_CN1 << "\tPairs where both reads align uniquely" << endl;
  cout << n_pairs_same_contig << "\tPairs where both reads align uniquely on the same contig" << endl;
  cout << n_pairs_well_aligned << "\tPairs appearing to be well-aligned, with apparent insert size <= MAX_JUMP (only these pairs are used as coverage)" << endl;
  cout << endl;
  
  
  // Create a CoverageAnalyzer to calculate genome coverage by these intervals.
  cout << Date() << ": Creating CoverageAnalyzer" << endl;
  CoverageAnalyzer analyzer( intervals, &genome_lens );
  
  // Find all holes in coverage.
  cout << Date() << ": Finding coverage holes" << endl;
  vec<seq_interval> holes;
  analyzer.GetCoveragesAtMost( MIN_COVERAGE-1, holes );
  sort( holes.begin(), holes.end(), seq_interval::OrderByLength() );
  holes = Reverse( holes );
  
  vec<size_t> hole_sizes( holes.size() );
  for ( size_t i = 0; i < holes.size(); i++ )
    hole_sizes[i] = holes[i].Length();
  size_t total_hole_size = BigSum( hole_sizes );
  
  // Calculate the amount of missing coverage by library.
  cout << Date() << ": Finding coverage holes by library" << endl;
  vec<size_t> missing_by_lib( n_libs, 0 );
  for ( size_t i = 0; i < n_libs; i++ ) {
    
    CoverageAnalyzer analyzer( intervals_by_lib[i], &genome_lens );
    vec<seq_interval> holes;
    analyzer.GetCoveragesAtMost( MIN_COVERAGE-1, holes );
    for ( size_t j = 0; j < holes.size(); j++ )
      missing_by_lib[i] += holes[j].Length();
  }
  
  
  
  // Basic report.
  cout << Date() << ": COVERAGE REPORT" << endl << endl;
  cout << "Out of " << genome_size << " bases in the reference, " << total_hole_size << " bases (" << ( ( 100.0 * total_hole_size ) / genome_size ) << "%) have less than MIN_COVERAGE=" << MIN_COVERAGE << " physical coverage by read pairs." << endl;
  cout << "These bases occur in " << holes.size() << " separate regions ('coverage holes'), which have N50 size " << N50( hole_sizes ) << "." << endl;
  cout << endl;
  
  // Report coverage holes.
  // This relies on the fact that the holes have been sorted by length.
  const size_t MAX_N_HOLES_TO_PRINT = 10;
  cout << "Largest " << MAX_N_HOLES_TO_PRINT << " coverage holes:" << endl;
  for ( size_t i = 0; i < holes.size() && i < MAX_N_HOLES_TO_PRINT; i++ )
    cout << "size = " << hole_sizes[i] << "\tat " << holes[i].SeqId() << ":" << holes[i].Begin() << "-" << holes[i].End() << endl;
  cout << endl;
  
  // Report missing coverage by library.
  ReportMissingCoverageByLibrary( cout, pairs, missing_by_lib, genome_size, total_hole_size );
  
  
  // Done.
  cout << Date() << ": Done with FindPhysicalCoverageHoles!" << endl;
  return 0;
}
