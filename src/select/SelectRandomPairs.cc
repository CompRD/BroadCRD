///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "random/Shuffle.h"
#include "PairsManager.h"
#include "feudal/IncrementalWriter.h"
#include "feudal/BinaryStream.h"
#include <map>
#include "lookup/LookAlign.h"

/**
 * SelectRandomPairs
 *
 * Randomly select pairs of reads, up to a given amount. A
 * correspondence map is also saved, to trace reads back to their
 * original ids. 
 *
 * READS_IN: it loads <READS_IN>.{fastb,qualb,pairs,qltout}
 * READS_OUT: it saves <READS_OUT>.{fastb,qualb,pairs,select,qltout}
 * ARCHIVE: if True, send log to a separate file (defaults to False)
 * N_PAIRS: select these many pairs
 * MB_TOTAL: select pairs up to this total length (Mb)
 * FRAC: select this fraction of total length (1.0 selects all reads)
 * SEED: seed for the random number generator
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( READS_IN );
  CommandArgument_String( READS_OUT );
  CommandArgument_Bool_OrDefault( ARCHIVE, False );

  // Criteria for selection (only one can be chosen).
  CommandArgument_UnsignedInt_OrDefault( N_PAIRS, 0 );
  CommandArgument_Double_OrDefault( MB_TOTAL, -1.0 );
  CommandArgument_Double_OrDefault( FRAC, -1.0 );
  CommandArgument_Bool_OrDefault( USE_LOCS, True );

  // Seed for random selection.
  CommandArgument_Bool_OrDefault( SHUFFLE, True );
  CommandArgument_UnsignedInt_OrDefault( SEED, 666 );
  EndCommandArguments;

  // File names.
  String in_bases_file = READS_IN + ".fastb";
  String in_quals_file = READS_IN + ".qualb";
  String in_pairs_file = READS_IN + ".pairs";
  String in_names_file = READS_IN + ".names";
  String in_qlt_file   = READS_IN + ".qltout";

  String log_file = READS_OUT + ".SelectRandomPairs.log";
  String out_bases_file  = READS_OUT + ".fastb";
  String out_quals_file  = READS_OUT + ".qualb";
  String out_select_file = READS_OUT + ".select";
  String out_pairs_file  = READS_OUT + ".pairs";
  String out_names_file  = READS_OUT + ".names";
  String out_qlt_file    = READS_OUT + ".qltout";

  // Output dir.
  Mkpath( Dirname( READS_OUT ) );

  // Log stream.
  ofstream archive_out;
  if ( ARCHIVE ) archive_out.open( log_file.c_str( ) );
  ostream &log = ARCHIVE ? archive_out : * (ostream *) &cout;
  if ( ARCHIVE ) {
    cout << "See log file " << log_file << "\n" << endl;
    PrintCommandPretty( log );
  }
  
  // Check arguments.
  {
    if ( ! ((  (MB_TOTAL > 0) && !(FRAC > 0) && !(N_PAIRS > 0) ) ||
                ( !(MB_TOTAL > 0) &&  (FRAC > 0) && !(N_PAIRS > 0) ) ||
                ( !(MB_TOTAL > 0) && !(FRAC > 0) &&  (N_PAIRS > 0) )) )
    {    PRINT3( MB_TOTAL, FRAC, N_PAIRS );    }
    ForceAssert((  (MB_TOTAL > 0) && !(FRAC > 0) && !(N_PAIRS > 0) ) ||
                ( !(MB_TOTAL > 0) &&  (FRAC > 0) && !(N_PAIRS > 0) ) ||
                ( !(MB_TOTAL > 0) && !(FRAC > 0) &&  (N_PAIRS > 0) ));
  }
  
  // Selected read ids.
  size_t n_reads = MastervecFileObjectCount( in_bases_file );
  
  vec<size_t> select;
  vec<longlong> maps_to( n_reads, -1 );
  
  // Pairs and bases - scoped for memory.
  {
    log << Date( ) << ": loading pairing info" << endl;
    PairsManager pairs( in_pairs_file );
    size_t n_pairs = pairs.nPairs();
    
    // Load bases and pairs.
    log << Date( ) << ": loading bases" << endl;
    vecbvec bases( in_bases_file );
    
    
    // Shuffle.
    log << Date( ) << ": shuffling pairs ids" << endl;
    vec<uint64_t> shuffled;
    if (SHUFFLE)
    {    if ( n_pairs > 0 ) 
              Shuffle64( (uint64_t)n_pairs, shuffled, (uint64_t)SEED );    }
    else shuffled = vec<uint64_t>( n_pairs, vec<uint64_t>::IDENTITY );
  
    size_t n_keepers = 0;
    if (FRAC == 1.0 || N_PAIRS > 0) {
      n_keepers = (FRAC == 1.0 ? n_pairs : N_PAIRS);
      
      size_t tot_length = 0;
      for (size_t ii = 0; ii < n_keepers; ii++) 
        tot_length += (bases[pairs.ID1(shuffled[ii])].size() +
                       bases[pairs.ID2(shuffled[ii])].size());

      log << endl
	  << "total length of selected reads: " << tot_length << endl
	  << "pairs in input:                 " << n_pairs << endl
	  << "pairs selected:                 " << n_keepers << endl
	  << endl;
    } 
    else {
      // Determine total Mb to keep.
      size_t required_length = 0;
      if (MB_TOTAL > 0) required_length = size_t(MB_TOTAL * 1000000.0);
      if (FRAC     > 0) required_length = uint64_t(round(double(bases.sumSizes()) * FRAC));
      
      // Determine how many pairs we want to keep.
      size_t tot_length = 0;
      for (size_t ii=0; ii<n_pairs; ii++) {
        if ( tot_length >= required_length ) break;
        tot_length += bases[ pairs.ID1( shuffled[ii] ) ].size( );
        tot_length += bases[ pairs.ID2( shuffled[ii] ) ].size( );
        n_keepers = ii+1;
      }
      log << endl
	  << "total length of selected reads: " << tot_length << endl
	  << "required total length:          " << required_length << endl
	  << "pairs in input:                 " << n_pairs << endl
	  << "pairs selected:                 " << n_keepers << endl
	  << endl;
    }
    
    // Save pairs, and populate select list.
    log << Date( ) << ": saving pairs" << endl;

    PairsManager sel_pairs( n_keepers * 2 );
    for (size_t ii=0; ii<n_keepers; ii++) {

      size_t id1 = select.isize();
      select.push_back( pairs.ID1( shuffled[ii] ) );
      maps_to[ pairs.ID1( shuffled[ii] ) ] = id1;
      
      size_t id2 = select.isize( );
      select.push_back( pairs.ID2( shuffled[ii] ) );
      maps_to[ pairs.ID2( shuffled[ii] ) ] = id2;

      sel_pairs.addPair( id1, id2, pairs.sep( shuffled[ii] ), pairs.sd( shuffled[ii] ),
			 pairs.libraryName( shuffled[ii] ) );
    }
    sel_pairs.Write( out_pairs_file );

    log << Date( ) << ": saving correspondence map (select)" << endl;

    BinaryWriter::writeFile( out_select_file.c_str(), select );

    log << Date( ) << ": saving bases" << endl;
    IncrementalWriter<bvec> sel_bases(out_bases_file.c_str());
    for (size_t ii=0; ii < select.size(); ii++)
      sel_bases.add( bases[select[ii]] );
    sel_bases.close();

    

  }


  // Select corresponding quals.
  if (IsRegularFile(in_quals_file)) {
    log << Date( ) << ": loading quals" << endl;
    vecqvec quals(in_quals_file);
    log << Date( ) << ": saving quals" << endl;
    IncrementalWriter<qvec> sel_quals(out_quals_file.c_str());
    for (size_t ii=0; ii < select.size(); ii++)
      sel_quals.add( quals[select[ii]] );
    sel_quals.close();
  }



  // Select corresponding names.
  if (IsRegularFile(in_names_file)) {
    log << Date( ) << ": loading names" << endl;
    vecString names(in_names_file);
    log << Date( ) << ": saving names" << endl;
    IncrementalWriter<String> sel_names(out_names_file.c_str());
    for (size_t ii=0; ii < select.size(); ii++)
      sel_names.add( names[select[ii]] );
    sel_names.close();
  }



  // Select aligns (this could be improved: no need to load all aligns).
  if ( USE_LOCS && IsRegularFile(in_qlt_file) ){
    log << Date() << ": loading alignments" << endl;
    vec<look_align> aligns_in;
    LoadLookAligns( in_qlt_file, aligns_in );

    log << Date( ) << ": selecting and saving alignments" << endl;
    vec<look_align> aligns_out;
    aligns_out.reserve( aligns_in.size( ) );
    for( size_t ii=0; ii<aligns_in.size( ); ii++) {
      if ( maps_to[ aligns_in[ii].query_id ] < 0 ) continue;
      look_align al = aligns_in[ii];
      al.query_id = maps_to[ aligns_in[ii].query_id ];
      aligns_out.push_back( al );
    }
    sort( aligns_out.begin( ), aligns_out.end( ) );
    
    WriteLookAligns( out_qlt_file, aligns_out );
  }

  // Done.
  log << Date( ) << ": SelectRandomPairs done" << endl;
  
}
