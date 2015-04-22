///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "PairsManager.h"
#include "Superb.h"
#include "Vec.h"
#include "Vec.h"
#include "VecUtilities.h"
#include "feudal/BinaryStream.h"
#include "paths/Alignlet.h"

/**
 * LibCoverage2
 *
 * Compute statistics on valid pairs for qltoutlet aligns. Code stolen
 * from PlotLibsOnRef and LibCoverage.
 *
 * REMARK! An insert is considered for coverage if its lengths is
 * within MAX_DEV deviations from its lib size.
 *
 * READS: it loads <READS>.{pairs,fastb} in RUN
 * ALIGNS: it loads <ALIGNS>.qltoutlet{,.index} in SUBDIR
 * ASSEMBLY: it loads <ASSEMBLY>.{contigs.fastb,superb} in SUBDIR
 */
 
int main(int argc, char *argv[])
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_String_OrDefault( READS, "scaffold_reads" );
  CommandArgument_String_OrDefault( ALIGNS, "scaffold_reads_filtered" );
  CommandArgument_String_OrDefault( ASSEMBLY, "linear_scaffolds0.clean.patched" );
  CommandArgument_Double_OrDefault( MAX_DEV, 5.0 );
  EndCommandArguments;

  // Dir and file names.
  String data_dir = PRE + "/" + DATA;
  String run_dir = data_dir + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR; 

  String pairs_file = run_dir + "/" + READS + ".pairs";

  String aligns_file = sub_dir + "/" + ALIGNS + ".qltoutlet";
  String index_file = sub_dir + "/" + ALIGNS + ".qltoutlet.index";
  String contigs_file = sub_dir + "/" + ASSEMBLY + ".contigs.fastb";
  String supers_file = sub_dir + "/" + ASSEMBLY + ".superb";

  // Load.
  cout << Date( ) << ": loading supers" << endl;
  vec<superb> supers;
  ReadSuperbs( supers_file, supers );

  cout << Date( ) << ": loading pairs" << endl;
  PairsManager pairs( pairs_file );
  int nlibs = pairs.nLibraries( );

  cout << Date( ) << ": loading aligns" << endl;
  vec<alignlet> aligns;
  BinaryReader::readFile( aligns_file, &aligns );
  
  cout << Date( ) << ": loading index" << endl;
  vec<int> index;
  BinaryReader::readFile( index_file, &index );
  
  // Various maps.
  int n_tigs  = MastervecFileObjectCount( contigs_file );
  int n_supers = supers.size( );
  vec<int> to_super( n_tigs, -1 );
  vec<int> start_on_super( n_tigs, 0 );
  vec<int> super_len( n_supers, 0 );
  for (int ii=0; ii<supers.isize( ); ii++) {
    super_len[ii] = supers[ii].TrueLength( );
    int start = 0;
    for (int jj=0; jj<supers[ii].Ntigs( ); jj++) {
      int tig = supers[ii].Tig( jj );
      to_super[tig] = ii;
      start_on_super[tig] = start;
      start += supers[ii].Len( jj );
      if ( jj < supers[ii].Ntigs( )-1 ) start += supers[ii].Gap( jj );
    }
  }
  
  // Valid windows of separations for various libraries (NB: separations).
  vec< pair<int,int> > valid_win;
  valid_win.reserve( nlibs );
  for (int ii=0; ii<nlibs; ii++) {
    int sep = pairs.getLibrarySep( ii );
    int dev = pairs.getLibrarySD( ii );
    int radius = int( MAX_DEV * double( dev ) );
    valid_win.push_back( make_pair( sep - radius, sep + radius ) );
  }
  
  // Remap contig-based aligns to super coordinates.
  cout << Date( ) << ": mapping aligns to super coordinates" << endl;
  for (size_t ii=0; ii<index.size( ); ii++) {
    if ( index[ii] < 0 ) continue;
    alignlet &al = aligns[ index[ii] ];
    int tig_id = al.TargetId( );
    int super_id = to_super[tig_id];
    if ( super_id < 0 ) continue;
    al.SetTargetId( super_id );
    al.SetTargetLength( super_len[super_id], al.Fw1( ) );
    al.Shift( start_on_super[tig_id] );
  }

  // Reads in input (size in PairsManager means number of pairs).
  vec<size_t> n_reads = pairs.getLibrarySizes( );
  for (int ii=0; ii<n_reads.isize( ); ii++)
    n_reads[ii] *= 2;
  
  // Counters (per library).
  vec<size_t> n_assembled( nlibs );
  vec<size_t> n_valid_pairs( nlibs );
  
  // Loop over all pairs.
  cout << Date( ) << ": parsing pairs" << endl;
  for (size_t pair_id=0; pair_id<pairs.nPairs( ); pair_id++) {
    int lib_id = pairs.libraryID( pair_id );
    int id1 = pairs.ID1( pair_id );
    int id2 = pairs.ID2( pair_id );
    if ( index[id1] >= 0 ) n_assembled[lib_id]++;
    if ( index[id2] >= 0 ) n_assembled[lib_id]++;
    if ( index[id1] < 0 || index[id2] < 0 ) continue;

    const alignlet &al1 = aligns[ index[id1] ];
    const alignlet &al2 = aligns[ index[id2] ];
    if ( al1.TargetId( ) != al2.TargetId( ) ) continue;
    if ( al1.Fw1( ) == al2.Fw1( ) ) continue;

    const alignlet &al_fw = al1.Fw1( ) ? al1 : al2;
    const alignlet &al_rc = al1.Fw1( ) ? al2 : al1;
    int sep = al_rc.pos2( ) - al_fw.Pos2( );
    if ( sep < valid_win[lib_id].first || sep > valid_win[lib_id].second )
      continue;
    
    n_valid_pairs[lib_id]++;
  }
  
  // Generate table.
  cout << Date( ) << ": generating table" << endl;
  
  vec< vec<String> > table;
  vec<String> line;
  {
    line = MkVec( String( "lib_name" ),
		  String( "lib_stats" ),
		  String( "n_reads" ),
		  String( "%_reads" ),
		  String( "n_pairs" ) );
    table.push_back( line );
  }
  
  for (int lib_id=0; lib_id<nlibs; lib_id++) {
    int sep = pairs.getLibrarySep( lib_id );
    int dev = pairs.getLibrarySD( lib_id );
    double ratio = double( n_assembled[lib_id] ) / double( n_reads[lib_id] );
    String str_name = pairs.getLibraryName( lib_id );
    String str_stats = ToString( sep ) + " +/- " + ToString( dev );
    String str_nreads = ToStringAddCommas( n_assembled[lib_id] );
    String str_pcreads = ToString( 100.0 * ratio, 1 );
    String str_npairs = ToStringAddCommas( n_valid_pairs[lib_id] );
    line = MkVec( str_name, str_stats, str_nreads, str_pcreads, str_npairs );
    table.push_back( line );
  }
  
  cout << "\nLEGEND\n"
       << "     n_reads:   number of reads assembled\n"
       << "     %_reads:   % of reads assembled\n"
       << "     n_pairs:   number of valid pairs\n"
       << "\n";
  PrintTabular( cout, table, 3, "lrrrr" );
  cout << endl;

  // Done.
  cout << Date( ) << ": done" << endl;
  
}
