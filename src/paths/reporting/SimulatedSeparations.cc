/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "LocsHandler.h"
#include "PairsHandler.h"

/**
 * SimulatedSeparations
 *
 * Find mean and stdev of the separations between pairs.
 *
 * INPUT (all files in <IN_DIR>):
 *   <READS>.fastb      reads bases
 *   <READS>.ref.locs   placement of reads onto reference
 *   <READS>.pairto     pairing info
 *   <GENOME>.fastb     reference genome
 *
 * IN_DIR: where input files are
 * READS: base name for reads
 * GENOME: base name for reference genome
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( IN_DIR );
  CommandArgument_String( READS );
  CommandArgument_String_OrDefault( GENOME, "genome" );
  EndCommandArguments;

  // File names.
  String reads_file = IN_DIR + "/" + READS + ".fastb";
  String pairs_file = IN_DIR + "/" + READS + ".pairto";
  String locs_file = IN_DIR + "/" + READS + ".ref.locs";
  String genome_file = IN_DIR + "/" + GENOME + ".fastb";

  // Load.
  int n_reads = MastervecFileObjectCount( reads_file );
  int n_genomes = MastervecFileObjectCount( genome_file );
  
  cout << Date( ) << ": loading pairs" << endl;
  phandler pairs( n_reads, pairs_file );
  unsigned int n_pairs = pairs.Pairs( ).size( );

  cout << Date( ) << ": loading locs" << endl;
  lhandler locs( n_reads, n_genomes );
  locs.LoadFromFile( locs_file );

  // Loop over all pairs.
  int dotter = 100000;
  cout << Date( ) << ": parsing "
       << n_pairs << " pairs (.="
       << dotter << " pairs):"
       << endl;

  vec<float> separations;
  separations.reserve( n_pairs );

  for (uint pair_id=0; pair_id<n_pairs; pair_id++) {
    if ( pair_id % dotter == 0 ) Dot( cout, pair_id / dotter );
    
    int id1 = pairs[pair_id].id1;
    int id2 = pairs[pair_id].id2;
    const read_location *loc1 = locs.GetPlacement( id1 );
    const read_location *loc2 = locs.GetPlacement( id2 );
    if ( !loc1 ) continue;
    if ( !loc2 ) continue;
    if ( loc1->Contig( ) != loc2->Contig( ) ) continue;
    if ( loc1->Orient( ) == loc2->Orient( ) ) continue;

    const read_location *loc_fw = loc1->Fw( ) ? loc1 : loc2;
    const read_location *loc_rc = loc1->Fw( ) ? loc2 : loc1;

    int sep = loc_rc->StartOnContig( ) - loc_fw->StopOnContig( ) - 1;
    separations.push_back( (float)sep );
  }
  cout << endl;

  // Report.
  NormalDistribution nd = SafeMeanStdev( separations );

  cout << "\n"
       << "pairs in input:  " << n_pairs << "\n"
       << "valid pairs:     " << separations.size( ) << "\n"
       << "separation:      ( " << ToString( nd.mu_,2 )
       << " +/- " << ToString( nd.sigma_, 2 ) << " )\n"
       << endl;
    
  // Done.
  cout << Date( ) << ": done" << endl;
}

