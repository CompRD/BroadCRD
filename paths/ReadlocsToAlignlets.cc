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

#include "Alignment.h"
#include "Basevector.h"
#include "PairsManager.h"
#include "Superb.h"
#include "Vec.h"
#include "paths/Alignlet.h"
#include "paths/ReadLoc.h"
#include "feudal/BinaryStream.h"

/**
 * ReadlocsToAlignlets
 *
 * Convert a set of readlocs to alignlets.
 * Default input:
 *  SUB/initial_scaffolds.contig.fastb
 *  SUB/initial_scaffolds.readlocs
 *  RUN/jump_reads_ec.fastb 
 *  RUN/long_jump_reads_ec.fastb
 *  (note that {long,}jump_reads_ec.{fastb,pairs} combines to scaffold_reads)
 * Default Output:
 *   SUB/scaffold_reads.qltoutlet
 *   SUB/scaffold_reads.qltoutlet.index
 */ 
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( FULL_RUN ); //search for jump,long_jump reads
  CommandArgument_String_OrDefault(SUB, "test");
  CommandArgument_String_OrDefault(READLOCS_IN, "initial_scaffolds");
  CommandArgument_String_OrDefault(SCAFFOLDS_IN, "initial_scaffolds");
  CommandArgument_String_OrDefault(ALIGNS_OUT, "scaffold_reads");
  CommandArgument_String_OrDefault(JUMP_READS, "jump_reads_ec");
  CommandArgument_String_OrDefault(LONG_JUMP_READS, "long_jump_reads_ec");
  EndCommandArguments;
  // load reads from RUN directory 
  String jump_reads_head = FULL_RUN + "/" + JUMP_READS;
  String jump_reads_fastb_file = jump_reads_head + ".fastb";
  String long_jump_reads_head = FULL_RUN + "/" + LONG_JUMP_READS;
  String long_jump_reads_fastb_file = long_jump_reads_head + ".fastb";
  String jump_pairs_file = jump_reads_head + ".pairs";
  String long_jump_pairs_file = long_jump_reads_head + ".pairs";
  // load scaffold and readloc files from SUB directory
  String FULL_SUB = FULL_RUN + "/ASSEMBLIES/" + SUB;
  String assembly_header =  FULL_SUB + "/" + SCAFFOLDS_IN;
  String read_locs_header =  FULL_SUB + "/" + READLOCS_IN;
  String aligns_out_file = FULL_SUB + "/" + ALIGNS_OUT +".qltoutlet";
  String index_out_file =  aligns_out_file + ".index";
  // load reads 
  uint64_t n_jump_reads = 0;
  if ( JUMP_READS != "" && IsRegularFile( jump_reads_fastb_file) )
  {
    n_jump_reads = MastervecFileObjectCount( jump_reads_fastb_file );
    cout << "Load jump_read file " << jump_reads_fastb_file << endl;
    cout << "   n_jump_reads= " << n_jump_reads << endl;
  }
  uint64_t n_long_jump_reads = 0;
  if ( LONG_JUMP_READS != "" &&  IsRegularFile( long_jump_reads_fastb_file) )
  {
    n_long_jump_reads = MastervecFileObjectCount( long_jump_reads_fastb_file );
    cout << "Load long_jump_read file " << long_jump_reads_fastb_file << endl;
    cout << "   n_long_jump_reads= " << n_long_jump_reads << endl;
  }
  uint64_t n_scaffold_reads = n_jump_reads + n_long_jump_reads;
  if( n_scaffold_reads <=0 ){
    cout << "Error! Found no reads. ";
    exit(1);
  }
  // load readlocs 
  read_locs_on_disk locs_file( read_locs_header, FULL_RUN );
  cout << "Load locs file " << read_locs_header + ".readlocs" << endl;

  // Load the contigs
  String contig_file = assembly_header + ".contigs.fastb";
  vecbasevector tigs( contig_file );
  uint32_t ntigs = tigs.size( );
  cout << "Load contigs file " << contig_file << endl;

  int dotter = 100;
  cout << Date( ) << ": " << ntigs << " contigs (.=" << dotter
    << " contigs):" << endl;
  vec<alignlet> aligns;
  vec<int> index( n_scaffold_reads, -1 );
  aligns.reserve( n_scaffold_reads);
  uint64_t n_long_jump_aligns = 0;
  uint64_t n_total_locs_scaffold = 0;
  for(uint contig_id=0;contig_id<tigs.size();contig_id++){
    if ( contig_id % dotter == 0 ) Dot( cout, contig_id/ dotter );
    vec<read_loc> locs;
    locs_file.LoadContig(contig_id, locs );
    for (int jj=0; jj<locs.isize( ); jj++) {
      n_total_locs_scaffold++;
      int read_id = locs[jj].ReadId( );
      int pos2 = locs[jj].Start( );
      int Pos2 = locs[jj].Stop( );
      uint target_id = locs[jj].ContigId( ); 
      int target_len = tigs[contig_id].isize();
      Bool fw1 = locs[jj].Fw( );
      ForceAssertEq( target_id, contig_id );
      // the read_id of long_jump_reads are numbered after jump_reads
      if ( locs[jj].LongJump() )
      {
	n_long_jump_aligns++;
	read_id += n_jump_reads;
      }
      if ( ! locs[jj].Jump( ) &&  ! locs[jj].LongJump( )) continue;
      // take only good aligns
      if ( pos2 < 0 || Pos2 > target_len ) { continue; }
      //if (locs[jj].ReadLength( ) < 96 ) { continue; }
      index[read_id] = aligns.isize( );
      aligns.push_back( alignlet( pos2, Pos2, target_id, target_len, fw1 ) );
    }
  }


  cout << endl;
  cout << "Total locs aligned with the scaffold: " << n_total_locs_scaffold << endl;
  cout << "Find " << n_long_jump_aligns << " long jump reads aligns " << endl;

  cout << Date( ) << ": saving " << aligns.size( ) << " aligns" << endl;
  BinaryWriter::writeFile( aligns_out_file, aligns );
  BinaryWriter::writeFile( index_out_file, index );

}
