/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"

#include "Superb.h"
#include "Vec.h"
#include "math/Functions.h"
#include "paths/Alignlet.h"
#include "paths/RemoveDuplicateAligns.h"
#include "paths/ScaffoldsUtils.h"
#include "paths/reporting/AllContigLinksCore.h"
#include "feudal/BinaryStream.h"

/**
 * AllContigLinks
 *
 * Realign reads to contigs (if needed, alignments are cached), and
 * find all possible clusters of consistent links between contigs (a
 * pair of oriented contigs is allowd to own multiple clusters)
 *
 * MIN_LINKS: only accept clusters with >= MIN_LINKS
 * SCAFFOLDS: used to find the correct contigs file
 * REMOVE_DUPLICATE: if true, run RemoveDuplicateAligns
 * FORCE: do not used cached aligns
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
  CommandArgument_String_OrDefault( READS, "scaffold_reads" );
  CommandArgument_String_OrDefault( SCAFFOLDS, "final" );
  CommandArgument_Int_OrDefault( MIN_LINKS, 2 );
  CommandArgument_Bool_OrDefault( REMOVE_DUPLICATES, True );
  CommandArgument_Bool_OrDefault( FORCE, False );
  EndCommandArguments;

  // Dir and file names.
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
  String out_dir = sub_dir + "/" + SCAFFOLDS + ".Aligns";

  String reads_file = run_dir + "/" + READS + ".fastb";
  String pairs_file = run_dir + "/" + READS + ".pairs";
  
  String supers_file = sub_dir + "/" + SCAFFOLDS + ".superb";
  String contigs_file = sub_dir + "/" + SCAFFOLDS + ".contigs.fasta";
  String aligns_file = out_dir + "/" + READS + ".qltoutlet";
  String index_file = out_dir + "/" + READS + ".qltoutlet.index";
  String links_file = out_dir + "/" + READS + ".all_links.txt";
  String log_file = out_dir + "/" + READS + ".AllContigLinks.log";
  
  Mkpath( out_dir );
  
  // Output streams.
  ofstream log( log_file.c_str( ) );
  PrintCommandPretty( log );
  
  cout << "\nSending output to " << out_dir << "\n" << endl;

  // Load core data.
  log << Date( ) << ": loading pairing info" << endl;
  PairsManager pairs( pairs_file );
  
  log << Date( ) << ": loading contigs fasta" << endl;
  vec<fastavector> contigs;
  LoadFromFastaFile( contigs_file, contigs );
  
  // Print lib info.
  log << "\nLIBRARY STATISTICS:\n\n";
  for (size_t ii=0; ii<pairs.nLibraries( ); ii++)
    log << ii << "\t"
	<< pairs.getLibraryName( ii ) << "\t"
	<< pairs.getLibrarySep( ii ) << " +/- "
	<< pairs.getLibrarySD( ii ) << "\n";

  // Align reads, or load aligns.
  vec<alignlet> aligns;
  vec<int> index;
  if ( FORCE || ! IsRegularFile( aligns_file ) ) {
    AlignReadsToContigs( K, out_dir, reads_file, contigs, aligns, index, log );

    if ( REMOVE_DUPLICATES )
      RemoveDuplicateAligns( pairs, aligns, index, log );

    log << Date( ) << ": saving aligns" << endl;
    BinaryWriter::writeFile( aligns_file, aligns );
    BinaryWriter::writeFile( index_file, index );
  }
  else {
    log << Date( ) << ": loading aligns" << endl;
    BinaryReader::readFile( aligns_file, &aligns );
    BinaryReader::readFile( index_file, &index );
  }
  
  // Run AllContigLinksCore.
  ofstream out( links_file.c_str( ) );
  AllContigLinksCore( MIN_LINKS, contigs, aligns, index, pairs, out, log );
  out.close( );

  // Done.
  String str_done = Date( ) + ": done";
  cout << str_done << endl;
  log << str_done << endl;
  log.close( );

}

