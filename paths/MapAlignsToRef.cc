///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"
#include "Vec.h"
#include "feudal/BinaryStream.h"
#include "lookup/LookAlign.h"
#include "paths/Alignlet.h"
#include "paths/MapAlignletToRef.h"

/**
 * MapAlignsToRef
 *
 * Given alignments of reads onto contigs, and of contigs onto a
 * reference, generate the set of implied alignments of reads onto
 * reference.
 *
 * This code is just a wrapper-example for paths/MapAlignletToRef.
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_String_OrDefault( ASSEMBLY, "initial_scaffolds" );
  CommandArgument_String_OrDefault( READS, "scaffold_reads" );
  CommandArgument_String_OrDefault( ALIGNS, "scaffold_reads_filtered" );
  
  // Dir and file names.
  String data_dir = PRE + "/" + DATA;
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
  String eval_dir = sub_dir + "/" + ASSEMBLY + ".Eval";

  String genome_fastb_file = data_dir + "/genome.fastb";

  String reads_fastb_file = run_dir + "/" + READS + ".fastb";

  String contigs_fastb_file = sub_dir + "/" + ASSEMBLY + ".contigs.fastb";
  String reads_aligns_file = sub_dir + "/" + ALIGNS + ".qltoutlet";  
  String reads_index_file = sub_dir + "/" + ALIGNS + ".qltoutlet.index";
  
  String contigs_aligns_file = eval_dir + "/aligns.qlt";
  
  // Leave early if reference is not found.
  if ( ! IsRegularFile( genome_fastb_file ) ) {
    cout << "Reference's fastb not found.\n" << endl;
    return 0;
  }

  if ( ! IsRegularFile( contigs_aligns_file ) ) {
    cout << "Aligns of contigs on reference missing. Run EvalScaffolds.\n";
    return 0;
  }

  // Load.
  int n_contigs = MastervecFileObjectCount( contigs_fastb_file );
  
  cout << Date( ) << ": loading aligns of contigs onto reference" << endl;
  vec<look_align_plus> cg_hits;
  LoadLookAlignPlus( contigs_aligns_file, cg_hits );
  
  cout << Date( ) << ": loading aligns of reads onto contigs" << endl;
  vec<alignlet> aligns;
  BinaryReader::readFile( reads_aligns_file, &aligns );
  
  cout << Date( ) << ": loading index of reads' aligns" << endl;
  vec<int> index;
  BinaryReader::readFile( reads_index_file, &index );
  
  cout << Date( ) << ": loading reads bases" << endl;
  vecbvec reads( reads_fastb_file );

  cout << Date( ) << ": loading contigs bases" << endl;
  vecbvec contigs( contigs_fastb_file );

  cout << Date( ) << ": loading reference bases" << endl;
  vecbvec genome( genome_fastb_file );

  cout << Date( ) << ": done loading\n" << endl;
  
  // Maps.
  vec< vec<int> > to_hits( n_contigs );   // contig_id -> ids in cg_hits
  for (int ii=0; ii<cg_hits.isize( ); ii++)
    to_hits[ cg_hits[ii].query_id ].push_back( ii );
  
  // Loop over all read alignments.
  for (size_t read_id=0; read_id<index.size( ); read_id++) {
    if ( index[read_id] < 0 ) continue;
    
    vec<alignlet> result;
    vec<int> result_hit_ids;
    const alignlet &al = aligns[ index[read_id] ]; 
    MapAlignletToRef( al, cg_hits, to_hits, result, &result_hit_ids );
    
    cout << "ALIGNS of r" << read_id << "\n";
    for (int ii=0; ii<result.isize( ); ii++) {
      bool contig_rc = cg_hits[ result_hit_ids[ii] ].Rc1( );
      bool read_rc = ! al.Fw1( );
      const bvec &theGenome = genome[result[ii].TargetId( )];
      const bvec &theContig = contigs[al.TargetId( )];
      const bvec &theRead = reads[read_id];
      int len = theRead.size( );

      cout << "c" << ( contig_rc ? "-" : "+" ) << ", "
	   << "r" << ( read_rc ? "-" : "+" ) << "\n";
      theGenome.PrintBases( cout, result[ii].pos2( ), len, false, 120 );
      theContig.PrintBases( cout, al.pos2( ), len, contig_rc, 120 );
      theRead.PrintBases( cout, 0, len, contig_rc != read_rc, 120 );
      
      cout << "\n";
    }
  }
    
  // Done.
  cout << Date( ) << ": done" << endl;

}

