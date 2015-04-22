///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"

#include "Alignment.h"
#include "Basevector.h"
#include "PairsManager.h"
#include "Superb.h"
#include "Vec.h"
#include "paths/Alignlet.h"
#include "paths/ReadLoc.h"
#include "feudal/IncrementalWriter.h"

/**
 * ReOrderReads
 *
 * Reorder the reads so that records aligned with same contigs are 
 * close 
 * 
 * READS: it loads <READS>.{fastb,pairs}
 * UNIBASES: it loads <UNIBASES>.unibases.k<K>
 * OUTBASE: output is saved in <READS>.<OUTBASE>
 * MIN_LENGTH: min length to align (after stripping last k-1 bases)
 * SAMPLE_SIZE: how many pairs to align (randomly selected)
 * MAX_MUTATION_RATE: discard aligns with excessive mutation rate
 * MAX_INDELS: discard aligns with too many indels
 * SEED: used to seed the randomizer
 * FORCE: do not use cached data, even if found
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String(READLOCS_IN);
  CommandArgument_String(SCAFFOLDS_IN);
  CommandArgument_String_OrDefault(FRAG_READS,"");
  CommandArgument_String_OrDefault(JUMP_READS,"");
  CommandArgument_String_OrDefault(LONG_JUMP_READS,"");
  CommandArgument_String_OrDefault(READS_DIR,"."); //where to search for frag,jump,long_jump reads
  CommandArgument_String_OrDefault(SCAFFOLDS_DIR,"."); //where to search for scaffolds and readlocs
  CommandArgument_String_OrDefault(READS_OUT_DIR, "."); //output directory
  CommandArgument_Int_OrDefault(MAX_READS, 500*1000*1000); //maximum 500M reads a time
  EndCommandArguments;
  
  // load reads from RUN directory 
  // deal with three different libraries
  vec<String> reads_name(3,"");
  reads_name[0] = FRAG_READS;
  reads_name[1] = JUMP_READS;
  reads_name[2] = LONG_JUMP_READS;
  // check reads 
  vec<uint64_t> n_reads_all(3,0);
  for(int ipass=0;ipass<3;ipass++)
  {
    String reads_fastb_file = READS_DIR+ "/" + reads_name[ipass] + ".fastb";
    if ( reads_name[ipass] != "" && IsRegularFile( reads_fastb_file) )
    {
      n_reads_all[ipass] = MastervecFileObjectCount( reads_fastb_file );
      cout << "Reads file " << reads_fastb_file << endl;
      cout << "   n_reads= " << n_reads_all[ipass] << endl;
    }
  }
  uint64_t n_total_reads = n_reads_all[0] + n_reads_all[1] + n_reads_all[2];
  if( n_total_reads <=0 ){
    cout << "Error! Found no reads. ";
    exit(1);
  }

  // load scaffold and readloc files from SUB directory
  String read_locs_header =  SCAFFOLDS_DIR + "/" + READLOCS_IN;
  String assembly_header =  SCAFFOLDS_DIR + "/" + SCAFFOLDS_IN;
  // Load the contigs
  String contig_file = assembly_header + ".contigs.fastb";
  vecbasevector tigs( contig_file );
  uint32_t ntigs = tigs.size( );
  cout << "Load contigs file " << contig_file << endl;
  cout << "  ntigs " << ntigs << endl;

  // the mapping from old read_id to new read_id, initialized to the total reads
  vec<uint64_t> new_frag_read_ids(n_reads_all[0], n_reads_all[0]); 
  vec<uint64_t> new_jump_read_ids(n_reads_all[1], n_reads_all[1]); 
  vec<uint64_t> new_long_jump_read_ids(n_reads_all[2], n_reads_all[2]); 
  // new index counter
  vec<uint64_t> new_read_id(3,0); 

  // load readlocs  and create the mapping
  //read_locs_on_disk locs_file( read_locs_header, READS_DIR);
  read_locs_on_disk locs_file( read_locs_header , "/dev/null"); 
  cout << "Load locs file " << read_locs_header + ".readlocs" << endl;
  int dotter = 100;
  for(uint contig_id=0;contig_id<tigs.size();contig_id++){
    if ( contig_id % dotter == 0 ) Dot( cout, contig_id/ dotter );
    vec<read_loc> locs;
    locs_file.LoadContig(contig_id, locs );
    // sort all the reads based on the starting position in the contig
    vec<pair<int,uint64_t> > sort_pairs;
    sort_pairs.reserve(locs.size());
    for (uint64_t jj=0; jj<locs.size( ); jj++) {
      int pos2 = locs[jj].Start( );
      sort_pairs.push_back(make_pair(pos2,jj));
    }
    sort(sort_pairs.begin(), sort_pairs.end());
    for(vec<pair<int,uint64_t> >::iterator it= sort_pairs.begin();it!= sort_pairs.end();++it)
    {
      uint64_t sub_index = it->second;
      read_loc & loc = locs[sub_index];
      uint64_t read_id = loc.ReadId( );
      if ( loc.Frag() ) {
	ForceAssertLt(read_id, n_reads_all[0]);
	new_frag_read_ids[read_id] = new_read_id[0] ++;
      }
      if (loc.Jump() ) {
	ForceAssertLt(read_id, n_reads_all[1]);
	new_jump_read_ids[read_id] = new_read_id[1] ++;
      }
      if (loc.LongJump() ) {
	ForceAssertLt(read_id, n_reads_all[2]);
	new_long_jump_read_ids[read_id] = new_read_id[2] ++;
      }
    }
  }
  cout << endl;
  // unaligned reads are at the end
  for(uint64_t i=0;i<n_reads_all[0]; i++) 
    if ( new_frag_read_ids[i] == n_reads_all[0] ) new_frag_read_ids[i] = new_read_id[0] ++;
  for(uint64_t i=0;i<n_reads_all[1]; i++) 
    if ( new_jump_read_ids[i] == n_reads_all[1] ) new_jump_read_ids[i] = new_read_id[1] ++;
  for(uint64_t i=0;i<n_reads_all[2]; i++) 
    if ( new_long_jump_read_ids[i] == n_reads_all[2] ) new_long_jump_read_ids[i] =new_read_id[2] ++;
  for(int ip=0;ip<3;ip++)
    AssertEq(new_read_id[ip], n_reads_all[ip]); // The mapping should be complete

  for(int ipass=0;ipass<3;ipass++)
  {
    // wrappers
    uint64_t n_reads = n_reads_all[ipass];
    if (n_reads == 0 ) continue;
    String reads_fastb_file = READS_DIR+ "/" + reads_name[ipass] + ".fastb";
    String reads_pairs_file = READS_DIR + "/" + reads_name[ipass] + ".pairs";
    String reads_fastb_ro_file = READS_OUT_DIR + "/" + reads_name[ipass] + "_ro.fastb";
    String reads_pairs_ro_file = READS_OUT_DIR + "/" + reads_name[ipass] + "_ro.pairs";
    cout << "Processing reads file " << reads_fastb_file << endl;
    cout << "Writing reordered reads file " << reads_fastb_ro_file << endl;
    vec<uint64_t>&  new_read_ids = ( ipass == 0 ? new_frag_read_ids : ( ipass == 1 ? 
      new_jump_read_ids : new_long_jump_read_ids) );
    // Reorder the reads based on the map file
    // the IncrementalWriter for creating large files
    IncrementalWriter<bvec> reads_ro(reads_fastb_ro_file.c_str());
    std::vector<bool> flag_new(n_reads, False); // if a record in new library is copied properly
    vecbvec output_pool; // the output pool
    output_pool.reserve(MAX_READS);
    vecbvec reads_in; // the input buffer
    reads_in.reserve(MAX_READS);
    // read from the input file and fill the output buffer
    for(uint64_t start = 0; start < n_reads; start += MAX_READS )
    {
      uint64_t length  = start + MAX_READS >= n_reads ? n_reads - start: MAX_READS;
      uint64_t end = start + length;
      cout << "Output chunk " << start << "-" << end << endl;
      output_pool.resize(length);
      // now load the library by chunks
      for(uint64_t start0 = 0; start0 < n_reads; start0 += MAX_READS )
      {
	uint64_t end0 = start0 + MAX_READS;
	if ( end0 > n_reads ) end0 = n_reads;
	reads_in.ReadRange(reads_fastb_file, start0, end0);
	//cout << "    Load chunk " << start0 << "-" << end0 << endl;
	// cout << _reads_in.size() << endl;
	for(size_t i=0; i<reads_in.size();i++)
	{
	  uint64_t old_read_id = i + start0;
	  uint64_t new_read_id = new_read_ids[old_read_id];
	  if ( new_read_id >= start && new_read_id < end)
	  {
	    flag_new[new_read_id] = True;
	    output_pool[new_read_id - start] = reads_in[i];
	  }
	}
	reads_in.clear();
      }
      reads_ro.add(output_pool.begin(), output_pool.end());
    }
    cout << "Check completeness ..." << endl;
    if ( find(flag_new.begin(), flag_new.end(), False) != flag_new.end() ) 
    {
      cout << "Warning, incomplete new library " << endl;
      exit(1);
    }
    cout << "Close output file and exit " << endl;
    reads_ro.close();

    // now deal with the pairs
    PairsManager pairs(reads_pairs_file );
    for(uint64_t pair_id=0;pair_id < pairs.nPairs(); pair_id ++)
    {
      longlong id1_old = pairs.ID1(pair_id);
      longlong id2_old = pairs.ID2(pair_id);
      AssertGe(id1_old,0);
      AssertGe(id2_old,0);
      longlong id1_new = new_read_ids[id1_old];
      longlong id2_new = new_read_ids[id2_old];
      pairs.SetIDs(pair_id, id1_new, id2_new);
    }
    pairs.Write(reads_pairs_ro_file);
  }

  // Write the new readloc files, with new read_id
  //
  // load readlocs  and create the mapping
  //read_locs_on_disk locs_file( read_locs_header, READS_DIR);
  Cp(read_locs_header + ".readlocs", read_locs_header + "_ro.readlocs");
  read_locs_on_disk locs_file_ro( read_locs_header + "_ro", "/dev/null"); 
  cout << "Write new readlocs file " << read_locs_header + "_ro.readlocs" << endl;
  for(uint contig_id=0;contig_id<tigs.size();contig_id++){
    vec<read_loc> locs;
    locs_file.LoadContig(contig_id, locs );
    // sort all the reads based on the starting position in the contig
    vec<pair<int,uint64_t> > sort_pairs;
    sort_pairs.reserve(locs.size());
    for (uint64_t jj=0; jj<locs.size( ); jj++) {
      int pos2 = locs[jj].Start( );
      sort_pairs.push_back(make_pair(pos2,jj));
    }
    sort(sort_pairs.begin(), sort_pairs.end());
    for(vec<pair<int,uint64_t> >::iterator it= sort_pairs.begin();it!= sort_pairs.end();++it)
    {
      uint64_t sub_index = it->second;
      read_loc & loc = locs[sub_index];
      uint64_t read_id = loc.ReadId( );
      uint64_t p_read_id = loc.PartnerReadId( );
      if ( loc.Frag() ) {
	ForceAssertLt(read_id, n_reads_all[0]);
	ForceAssertLt(p_read_id, n_reads_all[0]);
	loc.SetReadId(new_frag_read_ids[read_id]);
	loc.SetPartnerReadId(new_frag_read_ids[p_read_id]);
      }
      if (loc.Jump() ) {
	ForceAssertLt(read_id, n_reads_all[1]);
	ForceAssertLt(p_read_id, n_reads_all[1]);
	loc.SetReadId(new_jump_read_ids[read_id]);
	loc.SetPartnerReadId(new_jump_read_ids[p_read_id]);
      }
      if (loc.LongJump() ) {
	ForceAssertLt(read_id, n_reads_all[2]);
	ForceAssertLt(p_read_id, n_reads_all[2]);
	loc.SetReadId(new_long_jump_read_ids[read_id]);
	loc.SetPartnerReadId(new_long_jump_read_ids[p_read_id]);
      }
    }
    locs_file.WriteContig( read_locs_header + "_ro.readlocs", contig_id, locs );
  }
  cout << endl;

}
