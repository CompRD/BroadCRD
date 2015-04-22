/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"
#include "CoverageAnalyzer.h"
#include "PairsManager.h"
#include "paths/reporting/CHits.h"
#include "paths/reporting/ReftigUtils.h"
#include "util/RunCommand.h"
// MakeDepend: dependency MakeLookupTable

void PrintTable( const vec<seq_interval> &windows,
		 const vec<int> &clens,
		 const String &title,
		 const String &tag,
		 ostream &out );

/**
 * CleanViralContigs
 *
 * Align reads to contigs, and perform various cleaning operations:
 *
 * 1. Remove windows on contigs with no read coverage.
 *
 * 2. Remove windows with no physical coverage (it only applies to
 *    contigs >= MIN_PHYSCOV_LEN).
 *
 * 3. Remove contigs < MIN_CLEN.
 *
 * 4. Sort contigs by length (longer first).
 *
 * 5. Overwrite the pairs file using the new separations computed on
 *    the assembly contigs.  NB: again, the pairs file will be
 *    overwritten!  The original file will be saved as a backup.
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( RUN_DIR );
  CommandArgument_String( SUB_DIR );
  CommandArgument_String( OUT_DIR );
  CommandArgument_String_OrDefault( READS, "reads" );
  CommandArgument_String_OrDefault( CONTIGS, "contigs" );
  CommandArgument_Int_OrDefault( MIN_PHYSCOV_LEN, 1000 );
  CommandArgument_Int_OrDefault( MIN_CLEN, 250 );
  CommandArgument_Bool_OrDefault( FORCE, False );
  EndCommandArguments;
  
  // No overwriting!
  ForceAssertNe( SUB_DIR, OUT_DIR );

  // Dir and file names.
  String sub_dir = RUN_DIR + "/" + SUB_DIR;
  String out_dir = RUN_DIR + "/" + OUT_DIR;

  String reads_file = RUN_DIR + "/" + READS + ".fastb";
  String pairs_file = RUN_DIR + "/" + READS + ".pairs";

  String str_reads_on_contigs = READS + "_on_" + CONTIGS;
  String contigs_head = sub_dir + "/" + CONTIGS;
  String contigs_file = contigs_head + ".fastb";
  String lookup_file = contigs_head + ".lookup";
  String rawqlt_file = sub_dir + "/" + str_reads_on_contigs + ".raw.qlt";
  String qlt_file = sub_dir + "/" + str_reads_on_contigs + ".qlt";
  
  String log_file = out_dir + "/CleanViralContigs.log";
  String contigs_out_file = out_dir + "/" + CONTIGS + ".fastb";
  
  // These are needed.
  vec<String> needed;
  needed.push_back( reads_file );
  needed.push_back( pairs_file );
  needed.push_back( contigs_file );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;
  
  Mkpath( out_dir );
  ofstream log( log_file.c_str( ) );
  PrintCommandPretty( log );

  // Load pairs.
  cout << Date( ) << ": loading pairs" << endl;
  PairsManager pairs( pairs_file );

  // Assume one library only.
  ForceAssertEq( (int)pairs.nLibraries( ), 1 );

  // Generate lookup table.
  if ( FORCE || ! IsRegularFile( lookup_file ) ) {
    cout << Date( ) << ": generating lookup table" << endl;
    String theCommand
      = "MakeLookupTable LOOKUP_ONLY=True QUIET=True SOURCE= " + contigs_file
      + " OUT_HEAD=" + contigs_head;
    RunCommand( theCommand );
  }

  // Generate or load aligns (logging within).
  vec<look_align> hits;
  if ( ! IsRegularFile( qlt_file ) ) {
    int K = 20;
    vec<look_align> raw_hits;
    GetAlignsFast( K, reads_file, lookup_file, rawqlt_file,
		   raw_hits, !FORCE, sub_dir );

    cout << Date( ) << ": generating aligns... " << flush;
    hits.reserve( raw_hits.size( ) );
    for (size_t ii=0; ii<raw_hits.size( ); ii++)
      if ( raw_hits[ii].IsProper( ) ) hits.push_back( raw_hits[ii] );
    sort( hits.begin( ), hits.end( ) );
    WriteLookAligns( qlt_file, hits );
  }
  else { 
    cout << Date( ) << ": loading aligns... " << flush;
    LoadLookAligns( qlt_file, hits );
  }
  cout << hits.size( ) << " aligns found" << endl;

  // Load bases.
  cout << Date( ) << ": loading contigs... " << flush;
  vecbvec contigs( contigs_file );
  vec<int> clens( contigs.size( ), 0 );
  for (size_t ii=0; ii<contigs.size( ); ii++)
    clens[ii] = contigs[ii].size( );
  cout << contigs.size( ) << " contigs found" << endl;

  // Log info.
  cout << "\nSending log to " << log_file << "\n" << endl;

  // Find regions of no read coverage.
  vec<seq_interval> no_readcov;
  {
    vec<seq_interval> read_si;
    read_si.reserve( hits.size( ) );
    for( size_t ii=0; ii<hits.size( ); ii++) {
      const look_align &al = hits[ii];
      int seq_id = al.target_id;
      int begin = al.pos2( );
      int end = al.Pos2( );
      read_si.push_back( seq_interval( (int)ii, seq_id, begin, end ) );
    }
    CoverageAnalyzer analyzer( read_si, &clens );
    analyzer.GetCoveragesAtMost( 0, no_readcov );

    String title = "WINDOWS REMOVED - NO READ COVERAGE";
    String tag = "no_read_cov";
    PrintTable( no_readcov, clens, title, tag, log );
  }
  
  // Find regions of no physical coverage.
  vec<seq_interval> no_physcov;
  {
    int lib_id = 0;   // assuming one library only.
    CHits hitter( hits, pairs, clens );
    hitter.CoveragesAtMost( lib_id, 0, no_physcov );

    cout << Date( ) << ": updating library stats" << endl;
    hitter.UpdatePairsFile( pairs_file, &cout );

    vec<seq_interval> temp;
    temp.reserve( no_physcov.size( ) );
    for (size_t ii=0; ii<no_physcov.size( ); ii++) {
      int clen = clens[ no_physcov[ii].SeqId( ) ];
      if ( clen < MIN_PHYSCOV_LEN ) continue;
      temp.push_back( no_physcov[ii] );
    }
    swap( temp, no_physcov );
    
    String str_min = "(MIN_PHYSCOV_LEN = " + ToString( MIN_PHYSCOV_LEN ) + ")";
    String title = "WINDOWS REMOVED - NO PHYSICAL COVERAGE " + str_min;
    String tag = "no_phys_cov";
    PrintTable( no_physcov, clens, title, tag, log );
  }
  
  // Combine empty regions, and find regions to keep.
  vec<seq_interval> keepers;
  {
    vec<seq_interval> to_del = no_readcov;
    copy( no_physcov.begin( ), no_physcov.end( ), back_inserter( to_del ) );
    
    CoverageAnalyzer analyzer( to_del, &clens );
    analyzer.GetCoveragesAtMost( 0, keepers );
  }

  vecbvec new_contigs;
  new_contigs.reserve( keepers.size( ) );
  for (size_t ii=0; ii<keepers.size( ); ii++) {
    const seq_interval &region = keepers[ii];
    int id = region.SeqId( );
    int begin = region.Begin( );
    int length = region.Length( );
    bvec new_contig( contigs[id], begin, length );
    new_contigs.push_back( new_contig );
  }

  // Sort contigs by length, and save.
  cout << Date( ) << ": renumbering contigs and saving" << endl;
  vec< pair<int,int> > len2id;
  len2id.reserve( new_contigs.size( ) );
  for (int ii=0; ii<(int)new_contigs.size( ); ii++)
    len2id.push_back( make_pair( (int)new_contigs[ii].size( ), ii ) );
  sort( len2id.rbegin( ), len2id.rend( ) );
  vecbvec temp;
  temp.reserve( new_contigs.size( ) );
  for (int ii=0; ii<(int)len2id.size( ); ii++) {
    if ( len2id[ii].first < MIN_CLEN ) break;
    temp.push_back( new_contigs[ len2id[ii].second ] );
  }
  swap( temp, new_contigs );
  new_contigs.WriteAll( contigs_out_file );
  
  // Done.
  cout << Date( ) << ": done" << endl;

}



////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                            FUNCTIONS START HERE                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////



/**
 * PrintTable
 */
void PrintTable( const vec<seq_interval> &windows,
		 const vec<int> &clens,
		 const String &title,
		 const String &tag,
		 ostream &out )
{
  vec< vec<String> > table;
  
  vec<String> line = MkVec( String( "win_range" ),
			    String( "win_size" ),
			    String( "cg_id" ),
			    String( "cg_length" ),
			    String( "tag" ) );
  table.push_back( line );
  
  for (size_t ii=0; ii<windows.size( ); ii++) {
    if ( windows[ii].Length( ) < 1 ) continue;

    const seq_interval &si = windows[ii];
    String str_beg = ToString( si.Begin( ) );
    String str_end = ToString( si.End( ) );
    line[0] = "[" + str_beg + ", " + str_end + ")";
    line[1] = ToString( si.Length( ) );
    line[2] = ToString( si.SeqId( ) );
    line[3] = ToString( clens[ si.SeqId( ) ] );
    line[4] = tag;
    table.push_back( line );
  }
  
  out << title << "\n\n";
  PrintTabular( out, table, 3, "rrrrl" );
  out << endl;
}
