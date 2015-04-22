/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Basevector.h"
#include "PrettyPrintTable.h"
#include "ReadLocation.h"
#include "STLExtensions.h"
#include "feudal/BinaryStream.h"
#include "paths/KmerPath.h"
#include "paths/Unipath.h"
#include "paths/simulation/Placement.h"

/**
 * ShowUnilocsBrief
 *
 * Print in a compact way informations about unilocs on unipaths. Sends
 * output to either cout, or to a log file.
 *
 * INPUT (relative to <PRE>/<DATA>/<RUN>):
 *   <READS>.paths.k<K>
 *   <READS>.unipaths.k<K>
 *   <READS>.unipathsdb.k<K>
 *   <READS>.unilocs.<K>.<MAX_COPY_NUMBER>.1
 *
 * ARCHIVE: send output to a file stream rather than cout
 * FW_ONLY: do not print info for rc unipaths
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_Int( K );
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String_OrDefault( READS, "reads" );
  CommandArgument_Int_OrDefault( MAX_COPY_NUMBER, 10 );
  CommandArgument_Bool_OrDefault( ARCHIVE, True );
  CommandArgument_Bool_OrDefault( FW_ONLY, False );
  EndCommandArguments;
  
  String strK = ToString( K );
  String strMaxCn = ToString( MAX_COPY_NUMBER ) + ".1";

  String data_dir = PRE + "/" + DATA;
  String run_dir = data_dir + "/" + RUN;
  String base_reads = run_dir + "/" + READS;

  String paths_file =  base_reads + ".paths.k" + strK;
  String unipaths_file = base_reads + ".unipaths.k" + strK;
  String unipathsdb_file = base_reads + ".unipathsdb.k" + strK;
  String unilocs_file = base_reads + ".unilocs." + strK + "." + strMaxCn;
  String log_file = base_reads + ".ShowUnilocsBrief.out";

  // Load.
  cout << Date( ) << ": loading paths" << endl;
  vecKmerPath paths( paths_file );
  int nreads = paths.size( );
  
  cout << Date( ) << ": loading unipaths lengths" << endl;
  vec<int> to_rc;
  vec<int> up_lens;
  {
    vec<tagged_rpint> unipathsdb;
    BinaryReader::readFile(unipathsdb_file, &unipathsdb);
    vecKmerPath unipaths(unipaths_file);
    UnipathInvolution( unipaths, unipathsdb, to_rc );
    up_lens.resize( unipaths.size( ) );
    for (size_t ii=0; ii<unipaths.size( ); ii++)
      up_lens[ii] = unipaths[ii].KmerCount( ) + K - 1;
  }
  
  cout << Date( ) << ": loading and sorting unilocs" << endl;
  BREAD2( unilocs_file, vec<read_location_short>, locs );
  if ( ! is_sorted( locs.begin( ), locs.end( ) ) ) {
    cout << Date( ) << ": sorting unilocs" << endl;
    sort( locs.begin( ), locs.end( ) );
  }

  cout << Date( ) << ": generating unilocs maps\n" << endl;
  vec<int> flocs( up_lens.size( ), -1 );
  for (int ii=locs.isize( )-1; ii>=0; ii--)
    flocs[ locs[ii].Contig( ) ] = ii;
      
  // Put unipaths and their rc together.
  vec< pair<int,int> > upsets;
  upsets.reserve( to_rc.isize( ) / 2 );
  for (int ii=0; ii<to_rc.isize( ); ii++)
    if ( ii < to_rc[ii] )
      upsets.push_back( make_pair( ii, to_rc[ii] ) );
  
  // Out stream.
  ofstream archive_out;
  if ( ARCHIVE ) {
    archive_out.open( log_file.c_str( ) );
    PrintCommandPretty( archive_out );
    cout << " Sending output to " << log_file << "\n" << endl;
  }
  ostream &out = ARCHIVE ? archive_out : * (ostream *) &cout;

  // Loop over all unipaths (put together fw and rc).
  for (int set_id=0; set_id<upsets.isize( ); set_id++) {
    const pair<int,int> &upids = upsets[set_id];
    int n_unilocs = 0;
    if ( flocs[ upids.first ] > -1 ) {
      for (uint ii=flocs[ upids.first ]; ii<locs.size( ); ii++) {
	if ( locs[ii].Contig( ) != upids.first ) break;
	n_unilocs++;
      }
    }
    
    out << "UNIPATHS " << upids.first
	<< " and " << upids.second
	<< " : " << up_lens[upids.first]
	<< " bp, " << n_unilocs
	<< " uniloc" << ( n_unilocs > 1 ? "s" : "" )
	<< "\n";
    
    // Unilocs info as a table.
    vec< vec<String> > table;
    vec<String> aline;

    // Fw-Rc loop.
    for (int direction=0; direction<2; direction++) {
      if ( FW_ONLY && direction == 1 ) continue;

      const int unipath_id = ( direction == 0 ) ? upids.first : upids.second;
      const int rc_unipath_id = ( direction == 1 ) ? upids.first : upids.second;
      int floc = flocs[unipath_id];
      if ( floc < 0 ) continue;
      
      // Loop over all unilocs for this unipath.
      for (uint ii=floc; ii<locs.size( ); ii++) {
	const read_location_short &rloc = locs[ii];
	if ( rloc.Contig( ) != unipath_id ) break;

	int unipath_len = up_lens[unipath_id];
	int begin = rloc.StartOnContig( );
	int len = paths[rloc.ReadId( )].KmerCount( ) + K - 1;
	int end = begin + len;
	bool fw = rloc.Fw( );
	String str_orient = fw ? "fw" : "rc";

	int head = ( begin < 0 ) ? -begin : 0;
	int tail = ( end > unipath_len ) ? end - unipath_len : 0;
	int true_begin = ( begin < 0 ) ? 0 : begin;
	int true_end = ( end > unipath_len ) ? unipath_len : end;
	
	aline.clear( );

	String str_unipath
	  = " unipath_" + ToString( unipath_id )
	  + "." + ToString( ii - floc )
	  + "/" + ToString( n_unilocs - 1 )
	  + " (" + ToString( unipath_len )
	  + " bp)";
	aline.push_back( str_unipath );
	
	String str_read
	  = " r_" + ToString( rloc.ReadId( ) )
	  + "_" + str_orient
	  + " (" + ToString( len )
	  + " bp)";
	aline.push_back( str_read );

	String str_head = ( head > 0 ? " HEAD_" + ToString( head ) : "  " );
	aline.push_back( str_head );

	String str_interval
	  = "  [" + ToString( true_begin )
	  + "," + ToString( true_end )
	  + ") ";
	aline.push_back( str_interval );

	String str_tail = ( tail > 0 ? " TAIL_" + ToString( tail ) : "  ");
	aline.push_back( str_tail );
	
	table.push_back( aline );
      }
    }

    // Print current pair of unipaths.
    BeautifyAndPrintTable( table, out );
    out << "\n";
  }
  
  // Done.
  cout << Date( ) << ": done" << endl;
  out << flush;
  
}
