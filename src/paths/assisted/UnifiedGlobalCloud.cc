///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "VecUtilities.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"
#include "paths/KmerPathInterval.h"
#include "paths/assisted/KmerPathPerfectOverlaps.h"
#include "util/RunCommand.h"

void MapToGlobal( const vecKmerPath &unipaths,
		  const vec<tagged_rpint> &unipathsdb,
		  const vecKmerPath &gcloud,
		  vec<int> &to_gcid );

void TableReport( const String &descriptor,
		  const vecKmerPath &paths,
		  const vec<bool> &select,
		  vec< vec<String> > &table );

/**
 * UnifiedGlobalCloud
 *
 * Build a unified cloud of unipaths and reads, by selecting first
 * unipaths longer than MIN_ULEN (long unipaths), and then reads that
 * are not fully embedded in a long unipath. All operations carried
 * out in kmer space.
 *
 * NB: both fw and rc copies for each unipath and read are stored in
 *  the globlal cloud (unless palyndromic, in which case only one copy
 *  is stored).
 *
 * UNIPATHS: head name relative to RUN
 * READS: head name of reads, relative to RUN
 * HEAD_OUT: head name of output files, relative to RUN
 * MIN_ULEN: min length of unipaths to be selected (kmers size)
 */
int main( int argc, char *argv[] )
{
  RunTime( );
 
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_Int_OrDefault( K, 96 );
  CommandArgument_String_OrDefault( UNIPATHS, "all_reads" );
  CommandArgument_String_OrDefault( READS, "filled_reads_filt" );
  CommandArgument_String_OrDefault( HEAD_OUT, "global_cloud" );
  CommandArgument_Int_OrDefault( MIN_ULEN, 250 );
  EndCommandArguments;
  
  // Dir and file names.
  String strK = ToString( K );
  String run_dir = PRE + "/" + DATA + "/" + RUN;

  String unipaths_head = run_dir + "/" + UNIPATHS;
  String unipaths_file = unipaths_head + ".unipaths.k" + strK;
  String unipaths_rc_file = unipaths_head + ".unipaths_rc.k" + strK;
  String unipathsdb_file = unipaths_head + ".unipathsdb.k" + strK;
  String unibases_file = unipaths_head + ".unibases.k" + strK;
  
  String reads_head = run_dir + "/" + READS;
  String reads_paths_file = reads_head + ".paths.k" + strK;

  String out_head = run_dir + "/" + HEAD_OUT;
  String out_paths_file = out_head + ".paths.k" + strK;
  String out_rc_file = out_head + ".to_rc_paths.k" + strK;
  String out_bases_file = out_head + ".paths.k" + strK + ".fastb";
  String out_map_file = out_head + ".to_gcid";
  String log_file = out_head + ".UnifiedGlobalCloud.log";

  vec<String> needed;
  needed.push_back( unipaths_file );
  needed.push_back( unipathsdb_file );
  needed.push_back( reads_paths_file );
  needed.push_back( unipaths_rc_file );
  needed.push_back( unibases_file );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;

  // Load.
  cout << Date( ) << ": loading paths and unipaths" << endl;
  vecKmerPath reads_paths( reads_paths_file );
  vecKmerPath unipaths( unipaths_file );
  vecbvec unibases( unibases_file );
  vecKmerPath unipaths_rc( unipaths_rc_file );
  BREAD2( unipathsdb_file, vec<tagged_rpint>, unipathsdb );
  KmerBaseBroker kbb( K, unipaths, unipaths_rc, unipathsdb, unibases );

  // Selected unipaths (both fw and rc)
  cout << Date( ) << ": selecting objects in the global cloud" << endl;
  size_t n_sel = 0;
  vec<bool> uselect( unipaths.size( ), false );
  for (size_t ii=0; ii<unipaths.size( ); ii++) {
    if ( unipaths[ii].KmerCount( ) < MIN_ULEN ) continue;
    uselect[ii] = true;
    n_sel++;
  }
  
  // Selected reads.
  vec<bool> rselect( reads_paths.size( ), false );
  for (size_t ii=0; ii<reads_paths.size( ); ii++) {
    const KmerPathInterval &kpi = reads_paths[ii].Segment( 0 );
    const longlong kmer0 = kpi.Start( );
    vec<longlong> locs;
    Contains( unipathsdb, kmer0, locs );
    ForceAssert( locs.size( ) > size_t( 0 ) );
    ForceAssert( locs.size( ) < size_t( 3 ) );
    
    longlong loc = locs[0];
    if ( unipathsdb[ locs[0] ].Rc( ) ) {
      ForceAssert( locs.size( ) == size_t( 2 ) );
      ForceAssert( unipathsdb[ locs[1] ].Fw( ) );
      loc = locs[1];
    }
    const tagged_rpint &t_rpi = unipathsdb[loc];
    const int u_id = t_rpi.PathId( );
    const int u_pos = t_rpi.PathPos( );
    if ( IsReadEmbedded( reads_paths[ii], unipaths[u_id], 0, u_pos ) )
      continue;

    rselect[ii] = true;
    n_sel++;
  }
  
  // Build cloud paths. NB: both fw and rc copies for unipaths and reads.
  cout << Date( ) << ": adding rc copies, and removing duplicates" << endl;
  vecKmerPath gcloud;
  gcloud.reserve( n_sel );
  for (size_t ii=0; ii<uselect.size( ); ii++)
    if ( uselect[ii] )
      gcloud.push_back( unipaths[ii] );
  for (size_t ii=0; ii<rselect.size( ); ii++) {
    if ( ! rselect[ii] ) continue;
    KmerPath rc_path = reads_paths[ii];
    rc_path.Reverse( );
    gcloud.push_back( reads_paths[ii] );
    gcloud.push_back( rc_path );
  }
  
  // Sort unipaths, remove identical copies.
  sort( gcloud.begin( ), gcloud.end( ) );
  gcloud.erase( unique( gcloud.begin( ), gcloud.end( ) ), gcloud.end( ) );
  const int n_clouds = gcloud.size( );
  const String str_n_clouds = ToStringAddCommas( n_clouds );
  
  // Build rc map.
  vec<int> to_gcloud_rc;
  {
    to_gcloud_rc.resize( n_clouds, -1 );
    vecKmerPath::iterator it;
    for (int ii=0; ii<n_clouds; ii++) {
      KmerPath kp = gcloud[ii];
      kp.Reverse( );
      it = lower_bound( gcloud.begin( ), gcloud.end( ), kp );
      int rcii = distance( gcloud.begin( ), it );
      ForceAssert( rcii < n_clouds );
      ForceAssert( *it == kp );
      to_gcloud_rc[rcii] = ii;
    }
    for (int ii=0; ii<to_gcloud_rc.isize( ); ii++)
      ForceAssertEq( ii, to_gcloud_rc[ to_gcloud_rc[ii] ] );
  }
  
  // Generate map uid to gcid.
  vec<int> to_gcid;
  MapToGlobal( unipaths, unipathsdb, gcloud, to_gcid );
  
  // As vecbvec.
  vecbvec bases;
  bases.reserve( gcloud.size( ) );
  for (size_t ii=0; ii<gcloud.size( ); ii++)
    bases.push_back( kbb.Seq( gcloud[ii] ) );

  // Save.
  cout << Date( ) << ": saving " << str_n_clouds << " reads" << endl;
  BinaryWriter::writeFile( out_rc_file, to_gcloud_rc );
  gcloud.WriteAll( out_paths_file );
  bases.WriteAll( out_bases_file );
  WRITE( out_map_file, to_gcid );
			
  // Report overall stats.
  cout << "\nSELECTED UNIPATHS STATS (MIN_ULEN = " << MIN_ULEN << ")\n\n";
  vec< vec<String> > utable;
  TableReport( "unipaths", unipaths, uselect, utable );
  PrintTabular( cout, utable, 3, "rrrr" );
  
  cout << "\nSELECTED READS STATS\n\n";
  vec< vec<String> > rtable;
  TableReport( "reads", reads_paths, rselect, rtable );
  PrintTabular( cout, rtable, 3, "rrrr" );
  cout << endl;
  
  cout << "NB: the number of selected reads is not equal to the number of\n"
       << "    saved reads, since for every selected read both fw and rc\n"
       << "    copies are saved (after having removed duplicate copies)\n\n";
  // Done.
  cout << Date( ) << ": UnifiedGlobalCloud done" << endl;

}

/**
 * MapToGlobal
 *
 * Build a map uid to gcid (unipath id to id in global cloud).
 */
void MapToGlobal( const vecKmerPath &unipaths,
		  const vec<tagged_rpint> &unipathsdb,
		  const vecKmerPath &gcloud,
		  vec<int> &to_gcid )
{
  to_gcid.clear( );

  size_t n_unipaths = unipaths.size( );
  size_t n_gcloud = gcloud.size( );
  to_gcid.resize( n_unipaths, -1 );

  for (int gcid=0; gcid<(int)n_gcloud; gcid++) {
    const KmerPathInterval &kpi = gcloud[gcid].FirstSegment( );
    vec<longlong> locs;
    Contains( unipathsdb, kpi.Start( ), locs );

    bool uid_found = false;
    int uid = -1;
    for (int ii=0; ii<locs.isize( ); ii++) {
      const tagged_rpint &rpi = unipathsdb[ locs[ii] ];
      if ( rpi.Rc( ) ) continue;
      uid = rpi.ReadId( );
      if ( unipaths[uid] == gcloud[gcid] ) {
	uid_found = true;
	break;
      }
    }
    
    if ( ! uid_found ) continue;

    ForceAssertEq( to_gcid[uid], -1 );
    to_gcid[uid] = gcid;
  }
  
}

/**
 * TableReport
 *
 * Generate a printable table reporting "selected" stats.
 */
void TableReport( const String &descriptor,
		  const vecKmerPath &paths,
		  const vec<bool> &select,
		  vec< vec<String> > &table )
{
  table.clear( );

  size_t n_selected = 0;
  longlong klen_selected = 0;
  longlong klen_total = 0;
  for (size_t ii=0; ii<paths.size( ); ii++) {
    longlong klen = paths[ii].KmerCount( );
    klen_total += klen;
    if ( select[ii] ) {
      klen_selected += klen;
      n_selected++;
    }
  }
  
  vec<String> legend = MkVec( descriptor,
			      String( "selected" ),
			      String( "total" ),
			      String( "% selected" ) );
  table.push_back( legend );
  
  { // count
    double ratio = SafeQuotient( n_selected, longlong( paths.size( ) ) );
    vec<String> line = MkVec( String( "count" ),
			      ToStringAddCommas( n_selected ),
			      ToStringAddCommas( paths.size( ) ),
			      ToString( 100. * ratio, 1 ) );
    table.push_back( line );
  }

  { // kmer length
    double ratio = SafeQuotient( klen_selected, klen_total );
    vec<String> line = MkVec( String( "kmer length" ),
			      ToStringAddCommas( klen_selected ),
			      ToStringAddCommas( klen_total ),
			      ToString( 100. * ratio, 1 ) );
    table.push_back( line );
  }
  
}

