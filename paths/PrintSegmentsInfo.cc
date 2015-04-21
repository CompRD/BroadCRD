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
#include "graph/Digraph.h"
#include "paths/CPathMapper.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"

/**
 * PrintSegmentsInfo
 *
 * Print informations on which reads overlap the segments of the given
 * unipaths. The three {paths, paths_rc, pathsdb} must belong to the
 * same set of reads, while contigs.uinpaths are unipaths representing
 * contigs. Reads and contigs must have been "CommonPather"-ed for
 * PrintSegmentInfo to work.
 *
 * input (all relative to /PRE/DATA/RUN/SUBDIR):
 *   ../CONTIGS.UNIPATHS.kK
 *   ../READS.PATHS.kK
 *   ../READS.PATHS_rc.kK
 *   ../READS.PATHSDB.kK
 *
 * output: sent to cout
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String( SUBDIR );
  CommandArgument_Int( K );
  CommandArgument_String_OrDefault( CONTIGS, "contigs" );
  CommandArgument_String_OrDefault( READS, "reads" );
  CommandArgument_String_OrDefault( UNIPATHS, "unipaths" );
  CommandArgument_String_OrDefault( PATHS, "paths" );
  CommandArgument_String_OrDefault( PATHSDB, "pathsdb" );

  EndCommandArguments;
  
  // Dir and file names.
  String data_dir = PRE + "/" + DATA;
  String run_dir = data_dir + "/" + RUN;
  String sub_dir = run_dir + "/" + SUBDIR;
  
  String cg_base = sub_dir + "/" + CONTIGS;
  String r_base = sub_dir + "/" + READS; 
  String kK = ".k" + ToString( K );

  String unipaths_file = cg_base + "." + UNIPATHS + kK;
  
  String rpaths_file = r_base + "." + PATHS + kK;
  String rpaths_rc_file = r_base + "." + PATHS + "_rc" + kK;
  String rpathsdb_file = r_base + "." + PATHSDB + kK;
  
  // Load.
  cout << Date( ) << ": loading unipaths" << endl;
  vecKmerPath unipaths( unipaths_file ); 

  cout << Date( ) << ": loading reads' paths" << endl;
  vecKmerPath rpaths( rpaths_file ); 
  vecKmerPath rpaths_rc( rpaths_rc_file ); 

  cout << Date( ) << ": loading reads' pathsdb" << endl;
  BREAD2( rpathsdb_file, vec<tagged_rpint>, rpathsdb );
  
  cout << Date( ) << ": done loading\n" << endl;
  
  // A CPathMapper.
  CPathMapper mapper( &rpaths, &rpaths_rc, &rpathsdb );

  // Loop over all unipaths.
  for (int unip_id=0; unip_id<(int)unipaths.size( ); unip_id++) {
    const KmerPath &unipath = unipaths[unip_id];
    
    int tot_len = 0;
    
    // Loop over all segments in unipath.
    for (int seg_id=0; seg_id<unipath.NSegments( ); seg_id++) {
      const KmerPathInterval &kpint = unipath.Segment( seg_id );
      longlong len = kpint.Length( );

      vec<CKmerAlign> matches;
      mapper.Matches( kpint, matches, unip_id, tot_len );
      tot_len += len;

      cout << "u" << unip_id
	   << " s" << seg_id << "." << unipath.NSegments( ) - 1
	   << " len=" << len
	   << " nmatches=" << matches.size( ) << "\n";
      
      for (int ii=0; ii<(int)matches.size( ); ii++)
	matches[ii].PrintCoreInfo( cout );
      cout << "\n";
    }
  }

  // Done.
  cout << Date( ) << ": done" << endl;

}
