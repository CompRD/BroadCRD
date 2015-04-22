/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "MainTools.h"

#include "ParseSet.h"
#include "Superb.h"
#include "VecUtilities.h"
#include "paths/ReadLoc.h"

/**
 * CountReadLocs
 *
 * Print lib-by-lib stats for a given file of read locations
 *
 * SUPERS: count locs on a selected set of supers (do all if empty)
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_String_OrDefault(ASSEMBLY, "linear_scaffolds0.clean.patched");
  CommandArgument_String_OrDefault( SUPERS, "" );
  EndCommandArguments;

  // File names.
  String data_dir = PRE + "/" + DATA;
  String run_dir = data_dir + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR; 
  
  String head = sub_dir + "/" + ASSEMBLY;
  String supers_file = head + ".superb";

  // Load.
  read_locs_on_disk locs_file( head, run_dir );
  
  cout << Date( ) << ": loading supers" << endl;
  vec<superb> supers;
  ReadSuperbs( supers_file, supers );

  vec<bool> todo;
  int n_selected = 0;
  if ( SUPERS == "" ) {
    n_selected = supers.size( );
    todo.resize( n_selected, true );
  }
  else {
    vec<int> sel;
    ParseIntSet( SUPERS, sel );
    n_selected = sel.size( );
    todo.resize( supers.size( ), false );
    for (int ii=0; ii<sel.isize( ); ii++)
      todo[ sel[ii] ] = true;
  }

  // Loop over all contigs.
  vec< vec<uint64_t> > count;   // count[class][lib_id]
  
  cout << Date( ) << ": parsing " << n_selected << " supers\n" << endl;
  #pragma omp parallel for
  for (int super_id=0; super_id<supers.isize( ); super_id++) {
    if ( ! todo[super_id] ) continue;
    const superb &sup = supers[super_id];
    for (int cpos=0; cpos<sup.Ntigs( ); cpos++) {
      vec<read_loc> locs;
      #pragma omp critical
      locs_file.LoadContig( sup.Tig( cpos ), locs );
      for (int ii=0; ii<locs.isize( ); ii++) {
	int read_class = locs[ii].ReadClass( );
	int lib_id = locs[ii].LibraryId( );
	if ( count.isize( ) <= read_class )
	  count.resize( read_class + 1 );
	if ( count[read_class].isize( ) <= lib_id )
	  count[read_class].resize( lib_id + 1 );
	count[read_class][lib_id]++;
      }
    }
  }
  
  // Print result.
  if ( SUPERS != "" ) {
    cout << "selected supers: " << SUPERS << "\n";
    int n_tigs = 0;
    longlong tot_clen = 0;
    longlong tot_slen = 0;
    for (int ii=0; ii<supers.isize( ); ii++) {
      if ( ! todo[ii] ) continue;
      n_tigs += supers[ii].Ntigs( );
      tot_slen += supers[ii].TrueLength( );
      tot_clen += supers[ii].ReducedLength( );
    }
    cout << "  n_contigs:       " << ToStringAddCommas( n_tigs ) << "\n"
	 << "  contig_length:   " << ToStringAddCommas( tot_clen ) << "\n"
	 << "  scaffold_length: " << ToStringAddCommas( tot_slen ) << "\n"
	 << endl;
  }
  for (int ii=0; ii<count.isize( ); ii++) {
    String type = "";
    if ( ii == 0 ) type = "frag";
    else if ( ii == 1 ) type = "jump";
    else type = "long_jump";
    for (int jj=0; jj<count[ii].isize( ); jj++)
      cout << type << "[" << jj << "]\t" << count[ii][jj] << "\n";
  }
  cout << endl;

  // Done
  cout << Date( ) << ": done" << endl;
  
}
