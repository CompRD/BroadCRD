///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <omp.h>

#include "MainTools.h"
#include "Basevector.h"
#include "CoverageAnalyzer.h"
#include "PairsManager.h"
#include "VecUtilities.h"
#include "math/NStatsTools.h"
#include "paths/ReadLoc.h"
#include "util/RunCommand.h"
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

void ReportCoverage( const PairsManager &pairs,
		     const vec<seq_interval> &raw,
		     const vec<int> &clens,
		     const String &type,
		     const int &threshold,
		     ostream &out );

/**
 * ReportLowPhysCov
 *
 * It detects windows of low physical coverage on a set of contigs
 * (these can be a reference or contigs from an assembly, although in
 * this case the super structure will not be taken into account),
 * breaking the output by library type (frags, jumps, long jumps). It
 * needs the output from ReadAligns.
 *
 * WARNING! Make sure to point FRAGS, JUMPS, LONG_JUMPS to the same
 *          sets of reads used by ReadAligns!
 *

 * Remark: the lengths of the reads are not used in the computation,
 *         ie only the separations between the reads are used to
 *         compute the coverage. This means that overlapping fragment
 *         reads will appear as having almost no physical coverage.
 *
 * Remark: if a set of reads (LONG_JUMPS, for example), is not found,
 *         then print a warning and skip it.
 *
 * ASSEMBLY: It needs ../<ASSEMBLY>.{contigs,readlocs}
 * TAG_FRAG, TAG_JUMP, TAG_LONGJUMP: tag windows of frag (etc) coverage <= this
 * OUTDIR: relative to SUBDIR
 * FRAGS, JUMPS, LONG_JUMPS: read sets to be used (see Remark above)
 * MAX_STRETCH: used to define valid inserts
 * NUM_THREADS: use all available if 0
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String( SUBDIR );
  CommandArgument_String( ASSEMBLY );
  CommandArgument_String_OrDefault( OUTDIR, ASSEMBLY + ".ReportLowPhyCov" );
  CommandArgument_String_OrDefault( FRAGS, "" );
  CommandArgument_String_OrDefault( JUMPS, "" );
  CommandArgument_String_OrDefault( LONG_JUMPS, "" );
  CommandArgument_Int_OrDefault( TAG_FRAG, 9 )
  CommandArgument_Int_OrDefault( TAG_JUMP, 0 )
  CommandArgument_Int_OrDefault( TAG_LONGJUMP, 0 )
  CommandArgument_Double_OrDefault( MAX_STRETCH, 3.5 )
  CommandArgument_UnsignedInt_OrDefault( NUM_THREADS, 0 );
  EndCommandArguments;
  
  // Dir and file names.
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;

  String fpairs_file = run_dir + "/" + FRAGS + ".pairs";
  String jpairs_file = run_dir + "/" + JUMPS + ".pairs";
  String Jpairs_file = run_dir + "/" + LONG_JUMPS + ".pairs";
  bool use_frags = ( FRAGS != "" );
  bool use_jumps = ( JUMPS != "" );
  bool use_Jumps = ( LONG_JUMPS != "" );

  String ass_head = sub_dir + "/" + ASSEMBLY;
  String contigs_file = ass_head + ".contigs.fastb";
  String readlocs_file = ass_head + ".readlocs";

  vec<String> needed;
  needed.push_back( readlocs_file );
  needed.push_back( contigs_file );
  if ( use_frags ) needed.push_back( fpairs_file );
  if ( use_jumps ) needed.push_back( jpairs_file );
  if ( use_Jumps ) needed.push_back( Jpairs_file );

  if ( ! CheckFilesExist( needed, &cout ) ) return 1;
  if ( ! ( use_frags || use_jumps || use_Jumps ) ) {
    cout << "Fatal error: no frags, jumps, or long jump reads found\n" << endl;
    return 1;
  }
  
  // Thread control.
  NUM_THREADS = configNumThreads( NUM_THREADS );
  omp_set_num_threads( NUM_THREADS );
  
  // Load.
  PairsManager fpairs;
  if ( use_frags ) {
    cout << Date( ) << ": loading frag pairs" << endl;
    fpairs.Read( fpairs_file );
  }
  else
    cout << Date( ) << ": WARNING - skipping frags" << endl;

  PairsManager jpairs;
  if ( use_jumps ) {
    cout << Date( ) << ": loading jumps pairs" << endl;
    jpairs.Read( jpairs_file );
  }
  else
    cout << Date( ) << ": WARNING - skipping jumps" << endl;
  
  PairsManager Jpairs;
  if ( use_Jumps ) {
    cout << Date( ) << ": loading long jumps pairs" << endl;
    Jpairs.Read( Jpairs_file );
  }
  else
    cout << Date( ) << ": WARNING - skipping long jumps" << endl;
  
  cout << Date( ) << ": loading contigs" << endl;
  vec<int> clens;
  vecbvec contigs( contigs_file );
  clens.reserve( contigs.size( ) );
  for (int ii=0; ii<(int)contigs.size( ); ii++)
    clens.push_back( (int)contigs[ii].size( ) );
  
  read_locs_on_disk locs_parser( ass_head, run_dir );

  // Store coverages by inserts (internal intervals).
  vec<seq_interval> finserts;
  vec<seq_interval> jinserts;
  vec<seq_interval> Jinserts;
  if ( use_frags ) finserts.reserve( fpairs.nPairs( ) );
  if ( use_jumps ) jinserts.reserve( jpairs.nPairs( ) );
  if ( use_Jumps ) Jinserts.reserve( Jpairs.nPairs( ) );
  
  // Loop over all contigs.
  cout << Date( ) << ": loop over " << contigs.size( ) << " contigs" << endl;

  #pragma omp parallel for
  for (int cid=0; cid<(int)contigs.size( ); cid++) {
    vec<read_loc> locs;
    #pragma omp critical
    locs_parser.LoadContig( cid, locs );
    
    for (int loc_id=0; loc_id<locs.isize( ); loc_id++) {
      const read_loc &loc = locs[loc_id];
      if ( loc.Rc( ) ) continue;
      if ( loc.PartnerFw( ) ) continue;
      if ( loc.ContigId( ) != loc.PartnerContigId( ) ) continue;
      int obs_sep = loc.PartnerStart( ) - loc.Stop( );
      double stretch = double( obs_sep - loc.Sep( ) ) / double( loc.Dev( ) );
      if ( Abs( stretch > MAX_STRETCH ) ) continue;
      
      int seq_id = loc.ContigId( );
      int beg = loc.Stop( );
      int end = loc.PartnerStart( );
      seq_interval si( 1, seq_id, beg, end );

      #pragma omp critical
      if ( loc.Frag( ) && use_frags ) finserts.push_back( si );
      else if ( loc.Jump( ) && use_jumps ) jinserts.push_back( si );
      else if ( loc.LongJump( ) && use_Jumps ) Jinserts.push_back( si );
    }
  }
  cout << endl;
  
  // Find low coverage windows.
  ReportCoverage( fpairs, finserts, clens, "f", TAG_FRAG, cout );
  ReportCoverage( jpairs, jinserts, clens, "j", TAG_JUMP, cout );
  ReportCoverage( Jpairs, Jinserts, clens, "J", TAG_LONGJUMP, cout );
  
  // Done.
  cout << Date( ) << ": EvalJumps done" << endl;

}

/**
 * ReportCoverage
 */
void ReportCoverage( const PairsManager &pairs,
		     const vec<seq_interval> &raw,
		     const vec<int> &clens,
		     const String &type,
		     const int &threshold,
		     ostream &out )
{
  if ( raw.size( ) < 1 ) return;
  
  CoverageAnalyzer covs( raw, &clens );
  vec<seq_interval> wins;
  covs.GetCoveragesAtMost( threshold, wins );

  vec<int> wlens;
  wlens.reserve( wins.size( ) );
  for (size_t ii=0; ii<wins.size( ); ii++)
    wlens.push_back( wins[ii].Length( ) );
  
  vec< vec<String> > table;
  vec<String> line;
  for (size_t ii=0; ii<wins.size( ); ii++) {
    const seq_interval &win = wins[ii];
    int cid = win.SeqId( );
    int clen = clens[cid];
    int beg = win.Begin( );
    int end = win.End( );
    
    line.clear( );
    line = MkVec( type,
		  "c" + ToString( cid ),
		  "[" + ToString( beg ) + ", " + ToString( end ) + ")",
		  ToString( clen ),
		  ToString( end - beg ) );
    table.push_back( line );
  }
  
  String str_type;
  if ( type == "f" ) str_type = "FRAG";
  else if ( type == "j" ) str_type = "JUMP";
  else str_type = "LONG_JUMP";

  out << str_type << " READS\n"
      << "  listing windows of coverage <= " << threshold << "\n"
      << "  pairs in input: " << pairs.nPairs( ) << "\n"
      << "  valid pairs:    " << raw.size( ) << "\n"
      << endl;
  
  PrintTabular( out, table, 3, "rrrrr" );
  out << endl;

  PrintBasicNStats( ToLower( str_type ), wlens, out );
  out << endl;
  
}

