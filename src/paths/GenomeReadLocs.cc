/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// ---------------------------------------------------------------------------------

// GenomeReadLocs.  Generate read locations on genome, as a vec<read_location>
//
// Files In:
//  reads.fastb
//  genome.fastb
//
// Files Out:
//  reads.ref.locs

#include "Basevector.h"
#include "MainTools.h"
#include "ReadLocation.h"
#include "lookup/LookAlign.h"
#include "lookup/PerfectLookup.h"

int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_Int_OrDefault_Doc(K, 12, "Kmer size of lookup table");
  CommandArgument_String_OrDefault_Doc(READS, "reads",
    "Reads that we need to find locations on the genome for.");
  CommandArgument_String_OrDefault(GENOME, "genome");
  CommandArgument_Bool_OrDefault_Doc(UNIQUE_ONLY, False,
    "If True, only use uniqely place reads." );
  CommandArgument_String_OrDefault_Doc(ALIGNS_IN, "",
    "Use this set of alignments rather than align internally.");
  EndCommandArguments;

  // Set up directories.

  String data_dir = PRE + "/" + DATA;
  String run_dir = PRE + "/" + DATA + "/" + RUN;

  cout << "Loading Read and Genome Data\n";

  vecbasevector reads(run_dir + "/" + READS + ".fastb");
  int nreads = reads.size();

  // Record some statistics

  int unique = 0, unplaced = 0, multiple = 0;
     
  vec<read_location> locs_on_ref;

  cout << "Aligning " << nreads << " reads to reference\n";

  // Perform Alignment

  size_t chunk_size;
  int chunks;

  if (ALIGNS_IN == "") {
    chunk_size = 10000000;    // 10,000,000 reads at a time
    chunks = (reads.size() + chunk_size - 1) / chunk_size;
  } else {
    chunks = 1;
  }
    
  cout << Date( ) << ": Aligning in " << chunks << " chunk(s)" << endl;
    
  for (int i = 0; i < chunks; ++i) {
    vec<look_align> lookAligns;
    vec<vec<int> > index;
 
    int chunkStart, chunkEnd, readsInChunk;

    if (ALIGNS_IN == "") {
      chunkStart =  i * chunk_size;
      chunkEnd =  Min( (i+1) * chunk_size, reads.size() ) - 1;
      readsInChunk  = chunkEnd - chunkStart + 1;

      index.clear_and_resize(readsInChunk);
      PerfectLookup(K, reads, data_dir + "/" + GENOME + ".lookup", lookAligns,
		  FW_OR_RC, chunkStart, chunkEnd);

      cout << Date( ) << ": Chunk " << i+1 << " : Found " << lookAligns.size() 
	   << " potential aligns " << endl;

      // Build alignment index
      
      for (int i = 0; i < lookAligns.isize(); i++ ) {
	int readId = lookAligns[i].query_id;
	index[readId - chunkStart].push_back(i);
      }


    } else {

      LoadLookAligns( run_dir + "/" + ALIGNS_IN, lookAligns, index, nreads );
      cout << Date( ) << ": Imported alignments : Found " << lookAligns.size() 
	   << " potential aligns " << endl;

      chunkStart = 0;
      chunkEnd = nreads - 1;
      readsInChunk = nreads;
    }


    // Create read locations for uniquely placed reads

    for (int i = 0; i < readsInChunk; i++ ) {
      if (index[i].size() == 0)
	unplaced++;
      else {
	if (!UNIQUE_ONLY || index[i].size() == 1) {
	  int readId = i + chunkStart;
	  int alignId = index[i][0];
	  look_align& la = lookAligns[alignId];
	  int contigId = la.target_id;
	  int startPos = la.StartOnTarget();
	  int readLength = la.query_length;
	  int contigLength = la.target_length;
	  Bool rc = la.rc1;
	  locs_on_ref.push_back( read_location( readId, readLength, contigId, startPos,
						rc, contigLength) );
	}
	if (index[i].size() == 1)
	  unique++;	  
	else
	  multiple++;	  
      } 
    }
  }

  PRINT3(unique, unplaced, multiple);

  vecbasevector genome(data_dir + "/" + GENOME + ".fastb");

  Sort(locs_on_ref);
  WriteLocs( run_dir + "/" + READS + ".ref.locs", locs_on_ref, 
       genome.size( ), nreads );    
  
}
