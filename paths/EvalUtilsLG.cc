///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// For documentation, see EvalUtilsLG.h

#include "Basevector.h"
#include "Bitvector.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "lookup/LookAlign.h"
#include "paths/EvalUtils.h"
#include "paths/EvalUtilsLG.h"
#include "paths/HyperKmerPath.h"




void
EvalSummary( const String & output_file,
	     const HyperKmerPath & hkp,
	     const vecbasevector& genome, const vecbasevector& genome_diploid,
	     const vecbitvector& genome_amb,
	     const vec<look_align> & aligns,
	     const vec< vec<int> > & aligns_index )
{
  cout << Date( ) << ": EVAL - summary" << endl;
  
  // Open output stream 'out'.  If output_file == "stdout", then out = cout;
  // otherwise, out goes to output_file.
  ofstream o;
  if ( output_file != "stdout" )
    o.open( output_file.c_str( ) );
  ostream &out = output_file == "stdout" ? * (ostream *) &cout : o;
  
  
       String horiz_bar(80, '=');
  horiz_bar = "\n" + horiz_bar + "\n";
  out << horiz_bar << "\t\tEVAL SUMMARY: HyperKmerPath aligned to reference"
      << horiz_bar;
  
  
  // Write "SUMMARY:" line.
  // Also, if not writing to cout, report on gaps in coverage.
  longlong total_bases, total_covered;
  SummarizeReferenceCoverage( total_bases, total_covered, out,
			      genome, genome_diploid, genome_amb, aligns,
			      output_file == "stdout" );
  
  // Report mismatches and indels.
  int mismatches = 0, indels = 0;
  for ( int i = 0; i < aligns.isize( ); i++ ) {
    mismatches += aligns[i].mutations;
    indels += aligns[i].indels;
  }
  out << "\nBASE ERRORS: " << mismatches << " mismatches and "
      << indels << " indels.\n";
  
  // OPTIONALLY:
  // Create a second summary line that assumes the longer unaligned edges
  // actually represent correct sequence that is not in the reference.
  bool second_summary = false;
  if ( second_summary ) {
    out << "\nEdges of size >= 10 kb that do not align "
	 << "to reference at all:\n\n";
    vec<Bool> aligned( hkp.EdgeObjectCount( ), False );
    for ( int i = 0; i < aligns.isize( ); i++ )
      aligned[ aligns[i].query_id ] = True;
    longlong total_covered_plus = total_covered;
    for ( int i = 0; i < hkp.EdgeObjectCount( ); i++ ) {
      int n = hkp.EdgeObject(i).KmerCount( );
      if ( aligned[i] ) continue;
      total_covered_plus += n;
      if ( n < 10000 ) continue;
      out << "[" << BaseAlpha(i) << "] - "
	   << n << " kmers\n";
    }
    out << "\nSECOND SUMMARY: "
	<< PERCENT_RATIO(4, total_covered_plus, total_bases)
	<< " (" << total_covered_plus << "/" << total_bases << ") "
	<< "of reference covered (including unaligned edges)\n";
    out << "\nEdges of size >= 10 kb that align imperfectly "
	 << "to reference:\n\n";
    vec<Bool> perf( hkp.EdgeObjectCount( ), False );
    for ( int i = 0; i < aligns.isize( ); i++ ) {
      if ( aligns[i].FullLength( ) && aligns[i].Errors( ) == 0 )
	perf[ aligns[i].query_id ] = True;
    }
    for ( int i = 0; i < perf.isize( ); i++ ) {
      int n = hkp.EdgeObject(i).KmerCount( );
      if ( !aligned[i] || perf[i] || n < 10000 ) continue;
      out << "[" << BaseAlpha(i) << "] - "
	  << n << " kmers\n";
    }
    vec<Bool> short_or_bad( aligns.size(), False );
    for ( int i = 0; i < aligns.isize( ); i++ )
      if ( aligns[i].query_length < 10000 ||
	   aligns[i].Errors() != 0 ||
	   ! aligns[i].FullLength() )
	short_or_bad[i] = True;
    vec<look_align> good_and_long_aligns = aligns;
    EraseIf( good_and_long_aligns, short_or_bad );
    out << "\nCoverage of reference by edges of size >= 10 kb "
	 << "that align perfectly:";
    SummarizeReferenceCoverage( total_bases, total_covered, out,
				genome, genome_diploid, genome_amb, good_and_long_aligns,
				True );
    out << endl;
  }
  
  // Report on graph statistics (n components, component N50, n edges, etc.)
  PrintGraphStatistics( out, hkp, genome, genome_diploid, aligns, aligns_index );
  
  // Report on putative misassemblies.
  // Skip this step if writing to cout.
  if ( output_file != "stdout" )
    ReportMisassemblies( out, hkp, aligns, aligns_index );
  
  out << endl;
}





void
EvalLibraryStats( const String & output_file, const PairsManager & pairs )
{
  cout << Date( ) << ": EVAL - library stats" << endl;
  
  ofstream out( output_file.c_str( ) );
  pairs.printLibraryStats( out );
  out.close( );
}


// MakeDepend: dependency ScaffoldAccuracy
void
EvalScaffoldAccuracy( const String & output_file,
		      const String & genome_file,
		      const String & scaffold_file )
{
  cout << Date( ) << ": EVAL - scaffold accuracy" << endl;
  
  if ( !IsRegularFile( genome_file ) ) {
    cout << "EvalScaffoldAccuracy: Can't find genome file " + genome_file << endl;
    cout << "You need to provide the genome in fasta format, not fastb, to leave room for ambiguities." << endl;
    return;
  }
  
  if ( !IsRegularFile( scaffold_file ) ) {
    cout << "EvalScaffoldAccuracy: Can't find scaffold file "
	 << scaffold_file << endl;
    cout << "Have you run MakeScaffolds yet?" << endl;
    return;
  }
  
  SystemSucceed( "ScaffoldAccuracy"
		 + ARG( REFHEAD, genome_file.Before( ".fasta" ) )
		 + ARG(ASSEMBLY, scaffold_file)
                 + ARG(MIN_CONTIG, 0)
		 + " > " + output_file );
}
