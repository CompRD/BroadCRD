///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Basevector.h"
#include "FetchReadsAmb.h"
#include "VecUtilities.h"
#include "lookup/LookAlign.h"
#include "math/HoInterval.h"
#include "paths/ReadHBVs.h"
#include "paths/FlattenHyperFastavector.h"
#include "paths/reporting/ReftigUtils.h"

/**
 * MapNhoods
 *
 * Map local assembly to a reference (or to an assembly).
 *
 *
 * --- WARNING ---
 *
 *    If using FLATTEN_HYPER do not select more than one nhood at a
 *    time: there are static global variables in HyperFastavector that
 *    break ScaffoldComponent, when it runs on different hypers.
 *   
 * --- WARNING ---
 * 
 *
 * NHOODS: parsed with ParseIntSet (or do all)
 * OUT_DIR: relative to <SUBDIR>
 * REF_HEAD: relative to <DATA> (it needs <REF_HEAD>.lookup)
 * REF_GAPS: show gaps on reference (it needs <REF_HEAD>.fasta)
 * USE_CACHE: if true, load cached aligns (if found)
 * FLATTEN_HYPER: run FlattenHperFastavector and save output
 * GENERATE_DOT: if true, generate dot files (one per nhood)
 * VERBOSITY_DOT: used to label edges (if GENERATE_DOT=True, see GenerateDot)
 * ORIGIN: shift starts on target by this amount
 */
int main( int argc, char *argv[] )
{
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_Int( K );
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_String_OrDefault( NHOODS, "all" );
  CommandArgument_String_OrDefault( OUT_DIR, "nhoods.map" );
  CommandArgument_String_OrDefault( REF_HEAD, "genome" );
  CommandArgument_Bool_OrDefault( REF_GAPS, True );
  CommandArgument_Bool_OrDefault( USE_CACHE, False );
  CommandArgument_Bool_OrDefault( FLATTEN_HYPER, False );
  CommandArgument_Bool_OrDefault( GENERATE_DOT, False );
  CommandArgument_Int_OrDefault( VERBOSITY_DOT, 1 );
  CommandArgument_Int_OrDefault( ORIGIN, 0 );
  EndCommandArguments;

  // Dir and file names.
  String data_dir = PRE + "/" + DATA;
  String run_dir = data_dir + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
  String out_dir = sub_dir + "/" + OUT_DIR;
  
  String full_ref_head = data_dir + "/" + REF_HEAD;
  String seeds_ids_file = sub_dir + "/seeds.ids";
  String ref_lookup_file = full_ref_head + ".lookup";
  String ref_fasta_file = full_ref_head + ".fasta";
  
  String log_file = out_dir + "/main.log";

  Mkpath( out_dir );
  ofstream log( log_file.c_str( ) );
  PrintCommandPretty( log );
  cout << "Sending log to " << log_file << "\n" << endl;
  
  // Selected nhoods.
  int n_seeds = FirstLineOfFile( sub_dir + "/seeds.ids" ).Int( );
  vec<int> n_hoods;
  vec<bool> selected;  
  if ( NHOODS == "all" ) {
    n_hoods.resize( n_seeds, vec<int>::IDENTITY );
    selected.resize( n_seeds, True );
  }
  else {
    ParseIntSet( NHOODS, n_hoods );
    selected.resize( n_seeds, False );
    for (int ii=0; ii<n_hoods.isize( ); ii++)
      selected[ n_hoods[ii] ] = True;
  }

  // Eearly exit (see WARNING above).
  if ( FLATTEN_HYPER && n_hoods.size( ) > 1 ) {
    cout << "Can only select one nhood, if FLATTEN_HYPER=True.\n" << endl;
    return 1;
  }
  
  // Load.
  log << Date( ) << ": loading " << n_hoods.size( ) << " nhoods" << endl;
  vec<HyperBasevector> hypers;
  readHBVs( sub_dir, selected, &hypers );
  
  vecbitvector amb;
  if ( REF_GAPS ) {
    log << Date( ) << ": loading reference's amb bases" << endl;
    FetchReadsAmb( amb, ref_fasta_file, AMB_EQ_Nn );
  }
  vecbitvector *p_amb = REF_GAPS ? &amb : 0;
  
  // Loop over all selected nhoods.
  log << "\n";
  for (int seed_id=0; seed_id<n_seeds; seed_id++) {
    if ( ! selected[seed_id] ) continue;
    const HyperBasevector &the_hyper = hypers[seed_id];
    
    // Dir and file names for seed_id.
    String nhood_out_dir = out_dir + "/" + ToString( seed_id / 1000 );
    String tmp_dir = nhood_out_dir;
    Mkpath( nhood_out_dir );

    String nhood_head = nhood_out_dir + "/" + ToString( seed_id );
    String fastb_file = nhood_head + ".fastb";
    String aligns_file = nhood_head + ".aligns";
    String dot_file = nhood_head + ".dot";
    
    // Genearate (and save) nhood fastb file.
    if ( ( ! USE_CACHE ) || ( ! IsRegularFile( fastb_file ) ) ) {
      vecbvec bases;
      the_hyper.GenerateVecbasevector( bases );
      bases.WriteAll( fastb_file );
    }

    // Generate (or load) aligns.
    vec<look_align> aligns;
    GetAlignsFast( K, fastb_file, ref_lookup_file, aligns_file,
		   aligns, USE_CACHE, tmp_dir );

    // Generate reftigs.
    digraph agraph;
    vec< pair<int,ho_interval> > reftigs;
    HyperToReftigsCore( K, the_hyper, aligns, reftigs, &agraph );

    // Print dot file.
    if ( ( ! USE_CACHE ) || ( ! IsRegularFile( dot_file ) ) ) {
      GenerateDot( ORIGIN, nhood_head, the_hyper, agraph,
		   aligns, reftigs, VERBOSITY_DOT );
    }
    
    // Flatten hyper and save assembly output.
    if ( FLATTEN_HYPER ) {
      String base_out = nhood_head + ".flattened";
      String log_file = base_out + ".log";
      ofstream local_log( log_file.c_str( ) );
      HyperFastavector hfv( the_hyper );
      FlattenHyperFastavector( local_log, hfv, base_out );
    }

    // Print info.
    log << "n" << seed_id << "\n";
    PrintReftigs( log, K, ORIGIN, reftigs, p_amb );
  }
  
  // Done.
  log << Date( ) << ": done" << endl;
}
