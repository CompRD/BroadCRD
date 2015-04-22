///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "TokenizeString.h"
#include "CoverageAnalyzer.h"
#include "FetchReads.h"
#include "lookup/LookAlign.h"
#include "lookup/QueryLookupTableCore.h"
#include "paths/HyperBasevector.h"
#include "paths/reporting/CNhoodEval.h"
#include "system/System.h"
#include "util/RunCommand.h"
// MakeDepend: dependency QueryLookupTable
// MakeDepend: dependency MakeLookupTable

/**
 * CNhoodEval
 * ParseLog
 */
void CNhoodEval::ParseLog( const String &log_file )
{
  seed_id_ = -1;
  edge_id_.clear( );

  begin_.clear( );
  end_.clear( );
  cn_.clear( );

  true_target_id_.clear( );
  true_begin_.clear( );
  true_end_.clear( );
  true_rc_.clear( );
  true_cn_.clear( );
  
  // Set head_. Notice the last char in head_ must be either '.' or '/'.
  head_ = log_file.Before( "log" );

  // Parse log file.
  ifstream in( log_file.c_str( ) );
  while ( 1 ) {
    String line;
    getline( in, line );
    if ( ! in ) break;

    // Move down the log file until the line with "Local ID".
    if ( ! line.Contains( "Local", 0 ) ) continue;
    getline( in, line );
    if ( ! in ) break;
    
    // Core info is saved after a line with dashes.
    if ( line.Contains( "--------" ) ) {
      while ( 1 ) {
	getline( in, line );
	if ( ! in ) break;
	vec<String> tokens;
	Tokenize( line, tokens );
	if ( tokens.size( ) < 7 ) continue;

	// The SEED.
	if ( tokens[0] == "SEED" ) {
	  seed_id_ = edge_id_.size( );
	  tokens.erase( tokens.begin( ) );
	}

	// End of lines with core info.
	if ( tokens[5] != "+/-" ) break;
	ForceAssert( tokens[0].Int( ) == edge_id_.isize( ) );

	String predicted_loc = tokens[4].After( "." );
	vec<String> range;
	Tokenize( predicted_loc, '-', range );
	int ibegin = range[0].Int( );
	int iend = range[1].Int( );
	if ( predicted_loc.Before( range[0] ) == "-" )
	  ibegin *= -1;
	if ( predicted_loc.After( range[0] ).Before( range[1] ) == "--" )
	  iend *= -1;

	edge_id_.push_back( tokens[1].Int( ) );
	begin_.push_back( ibegin );
	end_.push_back( iend );
	cn_.push_back( tokens[3].Int( ) );

	// Keep parsing line, if eval code was run.
	if ( tokens.size( ) < 10 ) continue;
	String str_true_id = tokens[7].Before( "." );
	String str_true_loc = tokens[7].After( "." );
	Tokenize( str_true_loc, '-', range );
	ibegin = range[0].Int( );
	iend = range[1].Int( );
	if ( str_true_loc.Before( range[0] ) == "-" )
	  ibegin *= -1;
	if ( str_true_loc.After( range[0] ).Before( range[1] ) == "--" )
	  iend *= -1;
	
	true_target_id_.push_back( str_true_id.Int( ) );
	true_begin_.push_back( ibegin );
	true_end_.push_back( iend );
	true_rc_.push_back( tokens[8] == "rc" );
	true_cn_.push_back( tokens[9].Int( ) );
      }
    }

    // We can leave the log file now.
    break;

  }
  in.close( );
    
}
 
/**
 * CNhoodEval
 * CloudWin
 */
pair<int,int> CNhoodEval::CloudWin( ) const
{
  pair<int,int> win = make_pair( -1, -1 );

  if ( seed_id_ > -1 ) {
    win.first = begin_[0];
    win.second = end_[0];
  }
  for (int ii=1; ii<begin_.isize( ); ii++) {
    win.first = Min( win.first, begin_[ii] );
    win.second = Max( win.second, end_[ii] );
  }

  return win;
}

/**
 * CNhoodEval
 * TrueWin
 */
pair<int,int> CNhoodEval::TrueWin( ) const
{
  pair<int,int> win = make_pair( -1, -1 );
  if ( true_begin_.size( ) < 1 ) return win;

  if ( seed_id_ > -1 ) {
    win.first = true_begin_[0];
    win.second = true_end_[0];
  }
  for (int ii=1; ii<true_begin_.isize( ); ii++) {
    win.first = Min( win.first, true_begin_[ii] );
    win.second = Max( win.second, true_end_[ii] );
  }

  return win;
}

/**
 * CNhoodEval
 * ImpliedOnRef
 */
pair<int,int> CNhoodEval::ImpliedOnRef( int signed_pos ) const
{
  pair<int,int> win = this->CloudWin( );
  if ( win.first == win.second ) return win;   // An error.

  bool rc = signed_pos < 0;
  int pos = rc ? - signed_pos - 1 : signed_pos;
  if ( rc ) {
    int old_first = win.first;
    int seed_len = this->CloudLength( seed_id_ );
    win.first = pos + seed_len - win.second;
    win.second = pos + seed_len - old_first;
  }
  else {
    win.first += pos;
    win.second += pos;
  }
  
  return win;
}

/**
 * CNhoodEval
 * PrintInfo
 */
void CNhoodEval::PrintInfo( ostream &out ) const
{
  vec< vec<String> > table;
  vec<String> row;
  for (int ii=0; ii<edge_id_.isize( ); ii++) {
    row.clear( );
    row.push_back( ToString( ii ) );
    row.push_back( ToString( edge_id_[ii] ) );
    if ( true_begin_.size( ) > 0 ) {
      String str_rc = true_rc_[ii] ? " -" : " +";
      row.push_back( ToString( true_begin_[ii] ) );
      row.push_back( ToString( true_end_[ii] ) + str_rc );
      row.push_back( ToString( true_cn_[ii] ) );
    }
    else {
      row.push_back( ToString( begin_[ii] ) );
      row.push_back( ToString( end_[ii] ) );
      row.push_back( ToString( cn_[ii] ) );
    }
    row.push_back( ii == seed_id_ ? "seed" : "" );
    
    table.push_back( row );
  }

  PrintTabular( out, table, 2, "rrrrl" );
  
}

/**
 * CNhoodEval
 * Eval
 *
 * The align of seed onto reference is passed in the same format used
 * in SearchFastb2Core (edge_id, target_id, signed_start).
 */
String CNhoodEval::Eval( const vecbvec &bases,
			 const String &unibases_fastb,
			 const triple<int64_t,int64_t,int> &align ) const
{
  // These are used to enlarge windows on reference
  const double PERCENT_EXTRA = 0.05;
  const int MAX_EXTRA = 2000;

  // Overwrite existing files.
  const bool OVERWRITE = False;
  
  // Local file names.
  String local_dir = head_;
  if ( local_dir.Contains( "/" ) ) local_dir = local_dir.RevBefore( "/" );

  String log_file = head_ + "CNhoodEval.log";
  String target_head = head_ + "target";
  String unibases_ids_file = head_ + "unibases.ids";
  String hbv_fasta_file = head_ + "hbv.fasta";
  String qlt_assembly_file = head_ + "hbv.qlt";
  String qlt_unibases_file = head_ + "unibases.qlt";

  String target_fastb_file = target_head + ".fastb";
  String target_lookup_file = target_head + ".lookup";
  String target_range_file = target_head + ".range";
  
  ofstream log( log_file.c_str( ) );
  log << Date( ) << ": starting evaluation of local assembly" << endl;
  
  // Select window on target.
  log << Date( ) << ": generating target window" << endl;
  int target_id = align.second;
  pair<int,int> target_win = this->ImpliedOnRef( align.third );
  target_win.first = Max( 0, target_win.first );
  target_win.second = Min( bases[target_id].isize( ), target_win.second );
  int tbegin = target_win.first;
  int tend = target_win.second;
  double dtlen = double( tend - tbegin );
  int extra_amt = Min( MAX_EXTRA, int( dtlen * PERCENT_EXTRA ) );
  tbegin += - extra_amt;
  tend += extra_amt;
  tbegin = Max( 0, tbegin );
  tend = Min( bases[target_id].isize( ), tend );
  int tlen = tend - tbegin;
  
  log << "\n"
      << "Nhood expected align: " << ( align.third < 0 ? "rc" : "fw" )
      << " on t" << target_id
      << " [" << target_win.first << ", " << target_win.second << ")\n"
      << "Selected window on reference: t" << target_id
      << " [" << tbegin << ", " << tend << ")\n\n";
  
  // Generate lookup table.
  if ( OVERWRITE || ! IsRegularFile( target_lookup_file ) ) {
    vecbvec local_target;
    local_target.push_back( bvec( bases[target_id], tbegin, tlen ) );
    local_target.WriteAll( target_fastb_file );
    {
      ofstream rout( target_range_file.c_str( ) );
      rout << target_id << "\t" << tbegin << "\t" << tend << "\n";
      rout.close( );
    }
    
    log << Date( ) << ": generating lookup table" << endl;
    String theCommand
      = "MakeLookupTable LOOKUP_ONLY=True SOURCE=" + target_fastb_file
      + " OUT_HEAD=" + target_head;
    RunCommandWithLog( theCommand, "/dev/null" );
    
  }      

  // Align assembly.
  String qlt_args_core
    = "AI=True PARSEABLE=True LIST_UNPLACED_BY_PASS=True TMP_DIR=" + local_dir
    + " K=12 L=" + target_lookup_file;
  
  if ( OVERWRITE || ! IsRegularFile( qlt_assembly_file ) ) {
    if ( ! IsRegularFile( hbv_fasta_file ) )
      log << Date( ) << ": ERROR! no assembly found" << endl;
    else {
      log << Date( ) << ": aligning assembly to target" << endl;
      String theCommand
	= "QueryLookupTable " + qlt_args_core
	+ " SEQS=" + hbv_fasta_file;
      RunCommandWithLog( theCommand, qlt_assembly_file );
    }
  }

  // Align unibases.
  if ( OVERWRITE || ! IsRegularFile( qlt_unibases_file ) ) {
    log << Date( ) << ": saving ids of unibases in nhood (sorted!)" << endl;
    vec<int> ids = edge_id_;
    sort( ids.begin( ), ids.end( ) );
    {
      ofstream iout( unibases_ids_file.c_str( ) );
      for (int ii=0; ii<ids.isize( ); ii++) iout << ids[ii] << "\n";
      iout.close( );
    }
    
    log << Date( ) << ": aligning unibases to target" << endl;
    String theCommand
      = "QueryLookupTable " + qlt_args_core
      + " SEQS_IS_FASTB=True SEQS=" + unibases_fastb
      + " SEQS_TO_PROCESS=@" + unibases_ids_file;
    RunCommandWithLog( theCommand, qlt_unibases_file );

  }

  // Generate overall statistics.
  log << Date( ) << ": generating overal statistics" << endl;
  int orig_tlen = target_win.second - target_win.first;
  int offset = target_win.first - tbegin;

  double acov = this->Coverage( offset, orig_tlen, qlt_assembly_file );
  double ucov = this->Coverage( offset, orig_tlen, qlt_unibases_file );

  String result
    = ToString( target_id ) + " "
    + ToString( target_win.first ) + " "
    + ToString( target_win.second ) + " "
    + ToString( orig_tlen ) + " "
    + ToString( acov, 2 ) + " "
    + ToString( ucov, 2 );

  log << "\n"
      << "SUMMARY: " << result
      << "\n\n";
  
  // Done.
  log << Date( ) << ": done" << endl;
  log.close( );
  return result;
  
}

/**
 * CNhoodEval
 * Coverage
 */
double CNhoodEval::Coverage( const int offset,
			     const int tlen,
			     const String &qltout ) const
{
  double ratio = - 1.0;

  if ( ! IsRegularFile( qltout ) ) return ratio;

  vec<int> seq_lens( 1, tlen );
  vec<look_align_plus> hits;
  LoadLookAlignPlus( qltout, hits );

  vec<seq_interval> si_hits;
  si_hits.reserve( hits.size( ) );
  for (int ii=0; ii<hits.isize( ); ii++) {
    int begin = hits[ii].pos2( ) - offset;
    int end = hits[ii].Pos2( ) - offset;
    si_hits.push_back( seq_interval( ii, 0, begin, end ) );
  }

  CoverageAnalyzer cov( si_hits, &seq_lens );
  vec<seq_interval> regions;
  cov.GetCoveragesAtLeast( 1, regions );
  
  int rlen = 0;
  for (int ii=0; ii<regions.isize( ); ii++)
    rlen += regions[ii].Length( );

  if ( tlen > 0 ) ratio = 100.0 * double( rlen ) / double( tlen );
  
  return ratio;
}

