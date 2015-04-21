///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"

#include "Basevector.h"
#include "Superb.h"
#include "Vec.h"
#include "VecUtilities.h"
#include "math/Functions.h"
#include "paths/reporting/ReftigUtils.h"
#include "util/RunCommand.h"

/**
 * AlignSeedsToAssembly
 *
 * Align seed unibases to an assembly.
 *
 * FORCE: if false do not use cached aligns, but regenerate them
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_Int( K );
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_String_OrDefault( READS, "all_reads" );
  CommandArgument_String_OrDefault( ASSEMBLY, "final" );
  CommandArgument_Bool_OrDefault( FORCE, False );
  EndCommandArguments;

  // Dir and file names.
  String data_dir = PRE + "/" + DATA;
  String run_dir = data_dir + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
  String out_dir = sub_dir + "/" + ASSEMBLY + ".seeds";

  String kK = ".k" + ToString( K );
  String unibases_file = run_dir + "/" + READS + ".unibases" + kK;

  String seeds_file = sub_dir + "/seeds.ids";
  String contigs_head = sub_dir + "/" + ASSEMBLY + ".contigs";
  String contigs_file = sub_dir + "/" + ASSEMBLY + ".contigs.fastb";
  String supers_file = sub_dir + "/" + ASSEMBLY + ".superb";
  String lookup_file = sub_dir + "/" + ASSEMBLY + ".contigs.lookup";

  String query_file = out_dir + "/seeds.fastb";
  String aligns_file = out_dir + "/seeds.qlt";
  String placem_file = out_dir + "/seeds.placements";
  String log_file = out_dir + "/main.log";
  
  Mkpath( out_dir );

  // Load.
  int n_unibases = MastervecFileObjectCount( unibases_file );

  cout << Date( ) << ": loading ids of seeds" << endl;
  READ( seeds_file, vec<int>, seeds );
  ForceAssert( is_sorted( seeds.begin( ), seeds.end( ) ) );

  cout << Date( ) << ": loading supers" << endl;
  vec<superb> supers;
  ReadSuperbs( supers_file, supers );
  
  // Generate maps.
  int ntigs = 0;
  for (size_t ii=0; ii<supers.size( ); ii++)
    ntigs += supers[ii].Ntigs( );
  
  vec<int> to_super( ntigs, -1 );
  vec<int> tig_len( ntigs, 0 );
  vec<int> start_on_super( ntigs, 0 );
  for (int ii=0; ii<supers.isize( ); ii++) {
    int cursor = 0;
    for (int jj=0; jj<supers[ii].Ntigs( ); jj++) {
      int tig = supers[ii].Tig( jj );
      to_super[tig] = ii;
      tig_len[tig] = supers[ii].Len( jj );
      start_on_super[tig] = cursor;
      cursor += supers[ii].Len( jj );
      if ( jj < supers[ii].Ntigs( )-1 ) cursor += supers[ii].Gap( jj );
    }
  }
  
  // Query fastb is almost empty, but ids are preserved.
  if ( FORCE || ! IsRegularFile( query_file ) ) {
    cout << Date( ) << ": loading seeds unibases" << endl;
    vecbvec query;
    query.SparseRead( unibases_file, seeds, 0 );
    
    cout << Date( ) << ": saving (sparse) query fastb" << endl;
    query.WriteAll( query_file );
  }
  
  // Generate lookup table.
  if ( FORCE || ! IsRegularFile( lookup_file ) ) {
    cout << Date( ) << ": generating lookup table" << endl;
    String theCommand
      = "MakeLookupTable SOURCE=" + contigs_file
      + " OUT_HEAD=" + contigs_head
      + " LOOKUP_ONLY=true";
    RunCommandWithLog( theCommand, "/dev/null" );
  }

  // Align seeds to assembly contigs.
  vec<look_align> hits;
  if ( FORCE || ! IsRegularFile( aligns_file ) ) {
    cout << Date( ) << ": aligning seeds to assembly" << endl;
    GetAlignsFast( K, query_file, lookup_file, aligns_file, hits,
		   ! FORCE, sub_dir + "/tmp" );
    
    cout << Date( ) << ": sorting and rewriting aligns" << endl;
    sort( hits.begin( ), hits.end( ) );
    WriteLookAligns( aligns_file, hits );
  }
  else {
    cout << Date( ) << ": loading aligns" << endl;
    LoadLookAligns( aligns_file, hits );    
  }
  
  // Map each contig to the seeds that align it (multiple seeds allowed).
  vec< vec<int> > tig2hits( ntigs );
  for (int ii=0; ii<hits.isize( ); ii++) {
    int tig = hits[ii].target_id;
    tig2hits[tig].push_back( ii );
  }
  
  // Map each seed to its alignments (multiple hits are allowed).
  vec< vec<int> > unibase2hits( n_unibases );
  for (int ii=0; ii<hits.isize( ); ii++) {
    int qid = hits[ii].query_id;
    unibase2hits[qid].push_back( ii );
  }

  // Generate info (main log).
  ofstream log( log_file.c_str( ) );
  PrintCommandPretty( log );
  for (int super_id=0; super_id<supers.isize( ); super_id++) {
    const superb &sup = supers[super_id];
    const String str_sid = ToString( super_id );
    const String str_tot = ToString( sup.Ntigs( ) );

    log << "SUPER " << super_id
	<< "\t(" << sup.TrueLength( )
	<< " bp, " << sup.Ntigs( )
	<< " contigs)\n";

    int cursor = 0;

    vec< vec<String> > table;
    vec<String> line;
    for (int cgpos=0; cgpos<sup.Ntigs( ); cgpos++) {
      line.clear( );
      int beg = cursor;
      int end = cursor + sup.Len( cgpos );
      String str_tig = ToString( sup.Tig( cgpos ) );
      String str_pos = ToString( cgpos );
      String str_win = "[" + ToString( beg ) + ", " + ToString( end ) + ")";
      String str_seed = "";
      if ( tig2hits[sup.Tig( cgpos )].size( ) > 0 ) {
	const vec<int> &hit_ids = tig2hits[sup.Tig( cgpos )];
	str_seed = "##### " + ToString( hit_ids.size( ) ) + " seed(s): ";
	for (int ii=0; ii<hit_ids.isize( ); ii++) {
	  str_seed += "u" + ToString( hits[ hit_ids[ii] ].query_id );
	  if ( ii<hit_ids.isize( )-1 ) str_seed += ", ";
	}
      }
      
      line.push_back( "  c" + str_tig );
      line.push_back( "s" + str_sid + "." + str_pos + "/" + str_tot );
      line.push_back( ToString( sup.Len( cgpos ) ) );
      line.push_back( str_win );
      line.push_back( str_seed );
      table.push_back( line );

      if ( cgpos < sup.Ntigs( ) - 1 ) {
	int next = end + sup.Gap( cgpos );

	line.clear( );
	line.push_back( "gap" );
	line.push_back( "" );
	line.push_back( ToString( sup.Gap( cgpos ) ) );
	line.push_back( "[" + ToString( end ) + ", " + ToString( next ) + ")" );
	line.push_back( "" );
	table.push_back( line );
      }

      cursor += sup.Len( cgpos );
      if ( cgpos < sup.Ntigs( )-1 ) cursor += sup.Gap( cgpos );
    }
    
    PrintTabular( log, table, 3, "rrrrl" );
    log << endl;
  }
  log.close( );

  // Placement details.
  ofstream plout( placem_file.c_str( ) );

  int n_aligned_uni = 0;
  int n_aligned_multi = 0;
  int n_unaligned = 0;
  for (int ii=0; ii<seeds.isize( ); ii++) {
    int n_hits = unibase2hits[ seeds[ii] ].size( );
    if ( n_hits < 1 ) n_unaligned++;
    else if ( n_hits == 1 ) n_aligned_uni++;
    else n_aligned_multi++;
  }

  plout << "\n"
	<< "PLACEMENT STATS:\n"
	<< "\n"
	<< "n_seeds:   " << seeds.size( ) << "\n"
	<< "unplaced:  " << n_unaligned << "\n"
	<< "unique:    " << n_aligned_uni << "\n"
	<< "multiple:  " << n_aligned_multi << "\n"
	<< endl;
  
  plout << "PLACEMENT INFO:\n\n";
  for (int tig=0; tig<ntigs; tig++) {
    if ( tig2hits[tig].size( ) < 1 ) continue;
    vec< triple<int,int,int> > placs;   // begin, end, unibase_id
    placs.reserve( tig2hits[tig].size( ) );
    for (int ii=0; ii<tig2hits[tig].isize( ); ii++) {
      const look_align &la = hits[ tig2hits[tig][ii] ];
      int beg = la.pos2( );
      int end = la.Pos2( );
      int unibase_id = la.query_id;
      placs.push_back( triple<int,int,int>( beg, end, unibase_id ) );
    }
    sort( placs.begin( ), placs.end( ) );
    for (int ii=0; ii<placs.isize( ); ii++)
      plout << "c" << tig
	    << "\tu" << placs[ii].third
	    << "\t[" << placs[ii].first
	    << ", " << placs[ii].second
	    << ")_" << tig_len[tig]
	    << "\n";
  }
  plout << endl;

  plout.close( );
  
  // Done.
  cout << Date( ) << ": done" << endl;
  
}

