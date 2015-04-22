///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"

#include "Basevector.h"
#include "Qualvector.h"
#include "PrintAlignment.h" 
#include "btl/WalkAlign.h"
#include "lookup/LookAlign.h"
#include "util/RunCommand.h"
  
/**
 * PrintChimney
 *
 * Quick print method, to display a slice of consensus, as confirmed
 * by supporting reads.
 */
void PrintChimney( const int &rpos,
		   const String &str_ref,
		   const String &str_tag,
		   const vec< pair<int,String> > &n2chimney,
		   ostream &out )
{
  out << "Column starting at " << rpos
      << " (" << str_ref.size( )
      << " bp, failing " << str_tag
      << ")\n";
  
  vec< vec<String> > table;  

  table.push_back( MkVec( String( "  ref" ), String( "" ), str_ref ) );
  
  int end_show = Min( 5, n2chimney.isize( ) );
  for (int ii=0; ii<end_show; ii++)
    table.push_back( MkVec( "c" + ToString( ii ),
			    ToString( n2chimney[ii].first ),
			    n2chimney[ii].second ) );

  PrintTabular( out, table, 4, "rrr" );

  int n_remain_c = n2chimney.isize( ) - end_show;
  int n_remain_r = 0;
  for (int ii=end_show; ii<n2chimney.isize( ); ii++)
    n_remain_r += n2chimney[ii].first;
  
  out << "  " << n_remain_c
      << " more lines (" << n_remain_r
      << " reads)\n\n";

}

/**
 * EvalConsensus
 *
 * Esperimental code to evaluate consensus by looking at evidence from
 * read coverage.  The code does not take any action, at this point,
 * it just reports windows if any of the following tests fails:
 *
 * 1) empty_test: there are no reads at all, on this segment of
 *    consensus.
 *
 * 2) delta_test: the ratio between best and second best leading read
 *    is < DELTA.
 *
 * 3) matching_test: the leading (best) read does not match the
 *    consensus.
 *
 * Note that tests are performed in the order above, and only the
 * first failure mode is reported (ie, a window failing both
 * delta_test and matching_test will be tagged as failing the
 * delta_test only).
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( READS_HEAD );
  CommandArgument_String( READS_ALIGNS );
  CommandArgument_String( REF_HEAD );
  CommandArgument_String( OUT_HEAD );
  CommandArgument_Int( REF_ID );

  // Size of sliding window.
  CommandArgument_Int_OrDefault( WIN_SIZE, 13 );
  
  // For the delta test.
  CommandArgument_Double_OrDefault( DELTA, 1.5 );

  EndCommandArguments;
  
  // Dir and file names.
  String reads_fastb_file = READS_HEAD + ".fastb";
  
  String ref_fastb_file = REF_HEAD + ".fastb";

  String log_file = OUT_HEAD + ".log";
  String info_file = OUT_HEAD + ".info";

  vec<String> needed;
  needed.push_back( reads_fastb_file );
  needed.push_back( READS_ALIGNS );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;

  ofstream log( log_file.c_str( ) );
  ofstream info_out( info_file.c_str( ) );
  PrintCommandPretty( log );
  cout << "Sending log to " << log_file << endl;

  // Load.
  log << Date( ) << ": loading reads and reference" << endl;
  vecbvec rbases( reads_fastb_file );
  vecbvec ref( ref_fastb_file );
  
  log << Date( ) << ": loading aligns" << endl;
  vec<look_align> hits;
  LoadLookAligns( READS_ALIGNS, hits);
  const int n_hits = hits.size( );

  // Reference.
  const bvec &b2 = ref[REF_ID];
  const int b2len = b2.size( );

  // Loop over all windows on ref.
  const int dotter = 100;
  log << Date( ) << ": analyzing "
       << b2len << " windows (. = "
       << dotter << " windows)\n";

  int n_valid = 0;
  int n_empty_fails = 0;
  int n_delta_fails = 0;
  int n_match_fails = 0;
  for (int rpos=0; rpos<b2len - WIN_SIZE; rpos++) {
    if ( rpos % dotter == 0 ) Dot( log, rpos / dotter );
    
    // Chunk of ref.
    const bvec chunk2( b2, rpos, WIN_SIZE );
    
    // Build chimney;
    vec<String> chimney;
    for (int ii=0; ii<n_hits; ii++) {
      const look_align &lookal = hits[ii];
      const align &al = lookal.a;
      const int rid = lookal.QueryId( );
      const int tid = lookal.TargetId( );
      if ( tid != REF_ID ) continue;
      
      pair<int,int> walker = WalkOn2( al, rpos );
      int on1 = walker.first;
      if ( on1 == -1 ) continue;
      if ( on1 == numeric_limits<int>::max( ) ) continue;
      if ( on1 + WIN_SIZE > al.Pos1( ) ) continue;
      
      bvec rc_read;
      if ( lookal.Rc1( ) ) {
	rc_read = rbases[rid];
	rc_read.ReverseComplement( );
      }
      const bvec &b1 = lookal.Rc1( ) ? rc_read : rbases[rid];
      
      bvec chunk1( b1, on1, WIN_SIZE );
      chimney.push_back( chunk1.ToString( ) );
    }
    
    // Sort and select winner.
    sort( chimney.begin( ), chimney.end( ) );
    vec< pair<int,String> > n2chimney( 1, make_pair( 1, chimney[0] ) );
    for (int ii=1; ii<chimney.isize( ); ii++) {
      if ( chimney[ii] == n2chimney.back( ).second )
	n2chimney[ n2chimney.size( )-1 ].first += 1;
      else
	n2chimney.push_back( make_pair( 1, chimney[ii] ) );
    }
    sort( n2chimney.rbegin( ), n2chimney.rend( ) );
    
    // It fails empty test.
    const String str_ref = chunk2.ToString( );
    if ( n2chimney.size( ) < 1 ) {
      PrintChimney( rpos, str_ref, "empty_test", n2chimney, info_out );
      n_empty_fails++;
      continue;
    }
    
    // It fails delta test.
    if ( n2chimney.size( ) > 1 ) {
      double n0 = n2chimney[0].first;
      double n1 = n2chimney[1].first;
      if ( n0 < DELTA * n1 ) {
	PrintChimney( rpos, str_ref, "delta_test", n2chimney, info_out );
	n_delta_fails++;
	continue;
      }
    }    
    
    // It fails matching_test.
    const String &str0 = n2chimney[0].second;
    if ( str0 != str_ref ) {
      PrintChimney( rpos, str_ref, "matching_test", n2chimney, info_out );
      n_match_fails++;
      continue;
    }

    // A valid segment.
    n_valid++;

  } // loop over all bases on consensus.
  info_out.close( );
  log << endl;
  
  // Report stats.
  log << "\n"
      << "Tested " << b2len - WIN_SIZE << " windows. Of these, "
      << n_valid << " appear to be well confirmed by reads.\n"
      << "  " << n_empty_fails << " windows failed empty_test,\n"
      << "  " << n_delta_fails << " windows failed delta_test,\n"
      << "  " << n_match_fails << " windows failed matching_test.\n"
      << "\n"
      << "Note that delta = " << ToString( DELTA, 1 )
      << ", and that that tests are performed in the order above:\n"
      << "only the first failure mode is reported (ie, a window failing\n"
      << "both delta_test and matching_test will be tagged as failing the\n"
      << "delta_test only).\n"
      << endl;
  
  // Done.
  String date = Date( );
  log << "\n" << date<< ": EvalConsensus done" << endl;
  cout << Date( ) << ": EvalConsensus done" << endl;
  log.close( );
  
}
