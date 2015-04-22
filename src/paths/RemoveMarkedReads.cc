///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/* RemoveMarkedReads
 *
 * A simple module to remove reads from a dataset in accordance with a list of
 * reads to keep, e.g., in  reads.is_strong.  This module is designed to be a
 * replacement for RemoveSuspiciousReads2 in the RunAllPathsLG pipeline.
 *
 *
 * Input files required:
 *
 * <READS_IN>.fastb
 * <READS_IN>.pairs
 * <READS_IN>.is_strong
 *
 * Output files created:
 *
 * <READS_OUT>.fastb
 * <READS_OUT>.pairs
 * <READS_OUT>.lengths
 *
 *
 * All of these files are in /PRE/DATA/RUN.  Note that any files not listed
 * above will NOT be created or updated at READS_OUT, even if the requisite
 * READS_IN files exist.
 *
 *
 * Josh Burton
 * August 2009
 *
 ******************************************************************************/

#include "MainTools.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "PairsManager.h"
#include "util/ReadTracker.h"
#include "feudal/BinaryStream.h"

int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_String(READS_IN);
  CommandArgument_String(READS_OUT);
  CommandArgument_String_OrDefault_Doc(READS_TO_REMOVE, "",
    "Remove reads marked True,  keep reads marked False, in this vec<bool> file. "
    "Alternatively, use the argument READS_TO_KEEP.");
  CommandArgument_String_OrDefault_Doc(READS_TO_KEEP, "",
    "Remove reads marked False, keep reads marked True,  in this vec<bool> file. "
    "Alternatively, use the argument READS_TO_REMOVE.");
  CommandArgument_Bool_OrDefault_Doc(WRITE_QUALS, True,
    "remove read quality scores from associated qualb file" );
  CommandArgument_Bool_OrDefault_Doc(WRITE_LENGTHS, False,
    "generate read length file" );
  CommandArgument_Bool_OrDefault_Doc(SAVE_IN_READ_IDS, False,
    "for each OUT read_ID save the IN read_ID to a text file" );
  CommandArgument_Bool_OrDefault_Doc(TRACK_READS, True,
    "use ReadTracker to ouput a mapping of old reads to new" );
  
  EndCommandArguments;


  if (READS_TO_REMOVE != "" && READS_TO_KEEP != "")
    InputErr("Use either READS_TO_REMOVE or READS_TO_KEEP, not both.");
  if (READS_TO_REMOVE == "" && READS_TO_KEEP == "")
    InputErr("You must supply either READS_TO_REMOVE or READS_TO_KEEP.");

  
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String in_head  = run_dir + "/" + READS_IN;
  String out_head = run_dir + "/" + READS_OUT;
  
  
  /*****************************************************************************
   *
   *                                 LOAD FILES
   *
   ****************************************************************************/
  cout << Date( ) << ": Loading bases..." << endl;
  
  vecbasevector reads( in_head + ".fastb" );
  size_t n_reads_in = reads.size();
  
  cout << Date( ) << ": Loading quals..." << endl;
  
  vecqualvector quals;
  if (WRITE_QUALS) {
    quals.ReadAll(in_head + ".qualb");
    ForceAssertEq(reads.size(), quals.size());
  }

  cout << Date( ) << ": Loading pairing info..." << endl;

  PairsManager pairs( in_head + ".pairs" );
  size_t n_pairs_in = pairs.nPairs( );

  cout << Date( ) << ": Loading marked read info..." << endl;

  vec<Bool> reads_to_remove;
  if (READS_TO_KEEP != "") {
    BinaryReader::readFile( run_dir + "/" + READS_TO_KEEP, &reads_to_remove );
    for (size_t i = 0; i < reads_to_remove.size(); ++i)
      reads_to_remove[i] = (reads_to_remove[i] ? False : True); // flip vec<Bool>
  }
  else
    BinaryReader::readFile( run_dir + "/" + READS_TO_REMOVE, &reads_to_remove );

  // Input sanity checks.
  // If one of these asserts fails, the input files are inconsistent.
  PRINT2( n_reads_in, pairs.nReads() );
  ForceAssertEq( n_reads_in, pairs.nReads() );
  ForceAssertEq( n_reads_in, reads_to_remove.size( ) );
  
  
  /*****************************************************************************
   *
   *                      REMOVE MARKED READS FROM DATASET
   *
   ****************************************************************************/
  
  size_t n_reads_to_delete = Sum( reads_to_remove );
  cout << Date( ) << ": " << n_reads_to_delete << " reads out of " << n_reads_in << " are marked for deletion." << endl;
  
  // Shift the elements of reads and lane_index into their new locations.
  // Along the way, create a mapping from old read ID to new read ID.
  cout << Date( ) << ": Removing marked reads..." << endl;
  
  vec<size_t> in_read_IDs;
  size_t out_read_ID = 0;
  
  for (size_t in_read_ID = 0; in_read_ID < n_reads_in; ++in_read_ID) {
    if ( !reads_to_remove[in_read_ID] ) {

      if (SAVE_IN_READ_IDS || TRACK_READS)
        in_read_IDs.push_back(in_read_ID);

      if (in_read_ID != out_read_ID) {
        reads.SwapElements(out_read_ID, in_read_ID);
	if (WRITE_QUALS)
	  quals.SwapElements(out_read_ID, in_read_ID);
      }
      ++out_read_ID;
    }
  }
  
  
  // Tally the number of reads being removed from each part of each pair.
  size_t pairAReadsRemoved = 0;
  size_t pairBReadsRemoved = 0;
  size_t pairABReadsRemoved = 0;
  vec<Bool> removed(n_pairs_in, False);
  for ( size_t p = 0; p < n_pairs_in; ++p ) {
    bool remove1 = reads_to_remove[ pairs.ID1(p) ];
    bool remove2 = reads_to_remove[ pairs.ID2(p) ];
    
    if (  remove1 && !remove2 ) pairAReadsRemoved++;
    if ( !remove1 &&  remove2 ) pairBReadsRemoved++;
    if (  remove1 &&  remove2 ) pairABReadsRemoved++;
  }
  
  
  // Resize the data structures.
  size_t n_reads_out = out_read_ID;
  ForceAssertEq( n_reads_out, n_reads_in - n_reads_to_delete );
  
  reads.resize( n_reads_out );
  if (WRITE_QUALS)
    quals.resize( n_reads_out );
  pairs.removeReads( reads_to_remove, true );
  
  ForceAssertEq( n_reads_out, reads.size() );
  ForceAssertEq( n_reads_out, quals.size() );
  ForceAssertEq( n_reads_out, pairs.nReads() );
  
  
  // Add up the tallies.
  size_t n_pairs_out = pairs.nPairs();
  size_t n_pairs_removed = n_pairs_in - n_pairs_out;
  
  size_t pairedReadsRemoved = pairAReadsRemoved + pairBReadsRemoved + pairABReadsRemoved * 2;
  size_t pairedReadsRemaining = n_pairs_in * 2 - pairedReadsRemoved;
  
  size_t unpairedReadsRemoved = n_reads_to_delete - pairedReadsRemoved;
  size_t unpairedReadsRemaining = n_reads_out - pairedReadsRemaining;
  
  // Report to cout.
  cout << endl;
  cout << Date( ) << ": REPORT" << endl;
  cout << "\tRemoved " << unpairedReadsRemoved << " unpaired reads, " << unpairedReadsRemaining << " remaining (not including newly unpaired)." << endl;
  cout << "\tRemoved " << pairedReadsRemoved << " paired reads, " << pairedReadsRemaining << " remaining (includes orphaned paired reads)." << endl;
  cout << "\tRemoved " << n_pairs_removed << " pairs, " << n_pairs_out << " remaining." << endl;
  cout << "\tRemoved " << pairAReadsRemoved << " A reads of pair only." << endl;
  cout << "\tRemoved " << pairBReadsRemoved << " B reads of pair only." << endl;
  cout << "\tRemoved " << pairABReadsRemoved << " A and B reads of pair." << endl;
  cout << "\tUnpaired reads remaining: " << n_reads_out - n_pairs_out * 2 << endl;
  cout << "\tReads remaining: " << n_reads_out << endl;
  cout << endl;
  
  // Sanity checks: make sure the output is not degenerate.
  if ( n_reads_in > 0 && n_reads_out == 0 )
    FatalErr( "RemoveMarkedReads removed all reads!" );
  if ( n_pairs_removed > 0 && pairs.empty( ) )
    FatalErr( "RemoveMarkedReads removed all read pairings!" );
  
  
  /*****************************************************************************
   *
   *                                   OUTPUT
   *
   ****************************************************************************/
  cout << Date( ) << ": Writing output files..." << endl;
  
  reads.WriteAll( out_head + ".fastb" );
  if (WRITE_QUALS) quals.WriteAll( out_head + ".qualb" );
  
  if (SAVE_IN_READ_IDS) {
    ofstream outfs;
    outfs.open((out_head + ".in_read_IDs").c_str());
    outfs << "# out_read_ID in_read_ID" << endl;
    for (size_t out_read_ID = 0; out_read_ID < n_reads_out; out_read_ID++)
      outfs << out_read_ID << " " << in_read_IDs[out_read_ID] << endl;
    outfs.close();
  }

  if (TRACK_READS) {
    ReadTracker rt;
    rt.AddReadSet( in_head, reads_to_remove );
    rt.Dump(out_head);
  }


 //  // TEMP: For now, save the read pairs in the old format ( vec<read_pairing> )
//   // rather than as a PairsManager.
//   {
//     vec<read_pairing> pairings = pairs.convert_to_read_pairings( );
//     WritePairs( run_dir, pairings, n_reads_out, False, READS_OUT );
//   }
  pairs.Write( out_head + ".pairs" );
  
  // Write read length information.
  if (WRITE_LENGTHS) {
    vec<int> lengths(n_reads_out);
    for (size_t i = 0; i < n_reads_out; i++)
      lengths[i] = reads[i].size( );
    BinaryWriter::writeFile( out_head + ".lengths", lengths );
  }

  cout << Date( ) << ": Done!" << endl;
  return 0;
}
