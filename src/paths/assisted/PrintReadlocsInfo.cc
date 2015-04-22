///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"
#include "VecUtilities.h"
#include "lookup/LookAlign.h"
#include "paths/ReadLoc.h"

/**
 * PrintReadlocsInfo
 *
 * How many readlocs on each contig.
 */
int main( int argc, char *argv[] )
{
  RunTime( );
 
  BeginCommandArguments;
  CommandArgument_String( RUN_DIR );
  CommandArgument_String_OrDefault( ASSEMBLY, "ideal/final" );
  CommandArgument_Bool_OrDefault( ARCHIVE, True );
  EndCommandArguments;
  
  String sub_dir = RUN_DIR + "/ASSEMBLIES/assisted";
  String locs_head = sub_dir + "/" + ASSEMBLY + ".contigs";
  String contigs_file = sub_dir + "/" + ASSEMBLY + ".contigs.fastb";
  String aligns_file = sub_dir + "/" + ASSEMBLY + ".Eval/aligns.qlt";
  String log_file = sub_dir + "/" + ASSEMBLY + ".readlocs.info";

  ofstream archive_out;
  if ( ARCHIVE ) {
    archive_out.open( log_file.c_str( ) );
    PrintCommandPretty( archive_out );
    cout << " Sending output to " << log_file << "\n" << endl;
  }
  ostream &out = ARCHIVE ? archive_out : * (ostream *) &cout;

  out << Date( ) << ": loading contigs" << endl;
  vecbvec contigs( contigs_file );
  size_t n_contigs = contigs.size( );

  out << Date( ) << ": loading look_aligns" << endl;
  vec<look_align> aligns;
  LoadLookAligns( aligns_file, aligns );
  vec<int> n_placs( n_contigs, 0 );
  for (size_t ii=0; ii<aligns.size( ); ii++)
    n_placs[aligns[ii].query_id] += 1;
  
  out << Date( ) << ": parsing locs from " << n_contigs << " contigs" << endl;
  vec<int> rl_counts( n_contigs, 0 );
  read_locs_on_disk locs_parser( locs_head, RUN_DIR );

  for (int ii=0; ii<(int)n_contigs; ii++) {
    vec<read_loc> locs;
    locs_parser.LoadContig( ii, locs );
    rl_counts[ii] = locs.size( );
  }
  
  vec< vec<String> > table;
  vec<String> line;

  line = MkVec( String( "cid" ),
		String( "clen" ),
		String( "mult" ),
		String( "nlocs" ),
		String( "dens%" ) );
  table.push_back( line );

  for (int cid=0; cid<(int)n_contigs; cid++) {
    int clen = contigs[cid].size( );
    int mult = n_placs[cid];
    int nlocs = rl_counts[cid];
    double density = ( clen > 1 ) ? 100. * (double)nlocs / (double)clen : .0;
    
    line = MkVec( ToString( cid ),
		  ToString( clen ),
		  ToString( mult ),
		  ToString( nlocs ),
		  ToString( density, 2 ) );
    table.push_back( line );
  }
  
  out << "\n"
      << "LEGEND\n"
      << "  cid:     contig id\n"
      << "  clen:    contig length\n"
      << "  mult:    multiplicity of contig on reference\n"
      << "  nlocs:   number of locs in contig\n"
      << "  dens%:   % density of locs ( = 100.0 * nlocs / clen )\n"
      << "\n";
  PrintTabular( out, table, 3, "rrrrr" );
  out << endl;

  out << Date( ) << ": done" << endl;

}

