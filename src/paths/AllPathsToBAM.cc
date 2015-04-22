///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// AllPathsToBAM.  Convert assembled sequence to BAM format for use in external
// processing and visualizers.  MakeReadLocs must have been run first.  This is
// automatic in the RunAllPathsLG pipeline if PATCH_SCAFFOLDS=True.
//
// The bandwidth entries in the read locations have to be set.  If we wanted to
// work around this we could pass a default bandwidth to the alignment code.

#include "Basevector.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "Superb.h"
#include "paths/ReadLoc.h"
#include "util/RunCommand.h"
#include <sys/wait.h>

#include <omp.h>
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

int main(int argc, char *argv[]) {
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_String_Doc(SAMPLE, "Project G number");
  CommandArgument_String_OrDefault(SUBDIR, "test");
  CommandArgument_String_OrDefault_Doc(OUT_HEAD, "",
       "If unspecified, defaults to:\n"
       "<PRE>/<DATA>/<RUN>/ASSEMBLIES/<SUBDIR>/<ASSEMBLY>.[contig|super]s.bam\n"
       "Otherwise, will use complete path <OUT_HEAD>.bam");
  CommandArgument_String_OrDefault(ASSEMBLY, "linear_scaffolds0.clean.patched");
  CommandArgument_String_OrDefault_Doc(TIGS, "", 
       "if unspecified, process all contigs;"
       " otherwise it is one of the following: \n"
       "(a) a list of contig ids (in ParseIntSet format) or \n"
       "(b) the letter s followed by a list of scaffolds or \n"
       "(c) s<scaffold id>.<list of indices of contigs in the scaffold");
  CommandArgument_Bool_OrDefault_Doc(ORIG_READ_NAMES, False,
       "Tracing the original read names uses up a lot of overhead.\n"
       "When False, the names will have the format:\n\t"
       "frag|jump|long.<library>.<id>\n"
       "IMPLEMENTATION PENDING CAPABILITY TO STORE ORIGINAL READ NAMES");
  CommandArgument_Bool_OrDefault_Doc(SINGLE_SCAFFOLD_BAMS, False,
       "Create a separate BAM for each scaffold\n"
       "NOT YET IMPLEMENTED");
  CommandArgument_Bool_OrDefault_Doc(USE_SCAFFOLD_REFS, False,
       "Use scaffold alignments instead of contig alignments\n"
       "NOT YET IMPLEMENTED");
  CommandArgument_UnsignedInt_OrDefault_Doc(GAP_FLOOR, 100, 
       "Gaps below this value are extended to this size\n"
       "NOT YET IMPLEMENTED");
  CommandArgument_String_OrDefault_Doc(SORT_MEM, "8G", 
       "-m argument for samtools sort\n");
  CommandArgument_String_OrDefault_Doc(READLOCS_PREFIX, "",
       "if specified, file extension is .READLOCS_PREFIX.readlocs "
       "instead of .readlocs");
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
       "Number of threads to use (use all available processors if set to 0)");
  EndCommandArguments;

  // Thread control
   
  NUM_THREADS = configNumThreads(NUM_THREADS);
  omp_set_num_threads( NUM_THREADS );

  cout << "Using " << NUM_THREADS << " threads." << endl;

  // Define directories.

  String data_dir = PRE + "/" + DATA, run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;

  // Quietly disallow activation of incomplete implementation(s)

  if ( ORIG_READ_NAMES ) ORIG_READ_NAMES = False;
  if ( SINGLE_SCAFFOLD_BAMS ) SINGLE_SCAFFOLD_BAMS = False;
  if ( USE_SCAFFOLD_REFS ) USE_SCAFFOLD_REFS = False;

  // Fork off indexing of fasta
  String fasta_file = sub_dir + "/" + ASSEMBLY +
    ( USE_SCAFFOLD_REFS ? ".super" : ".contig" ) + "s.fasta";
  pid_t pid = Fork( "samtools faidx " + fasta_file );

  // Parse TIGS.

  vec<superb> scaffolds;
  vec<int> tigs, supers;
  if ( TIGS == "" ) {
    int n_tigs = MastervecFileObjectCount( sub_dir + "/" + ASSEMBLY
					   + ".contigs.fastb" );
    for ( int j = 0; j < n_tigs; j++ )
      tigs.push_back(j);
  } else if (TIGS.Contains("s", 0)) {
    TIGS = TIGS.After("s");
    ReadSuperbs(sub_dir + "/" + ASSEMBLY + ".superb", scaffolds);
    if (TIGS.Contains(".")) {
      if ( USE_SCAFFOLD_REFS ) USE_SCAFFOLD_REFS = False;
      int scaffold = TIGS.Before(".").Int();
      ForceAssertLt(scaffold, scaffolds.isize());
      vec<int> spos;
      ParseIntSet(TIGS.After("."), spos);
      for (int j = 0; j < spos.isize(); j++)
	tigs.push_back(scaffolds[scaffold].Tig(spos[j]));
    } else {
      ParseIntSet(TIGS, supers);
      UniqueSort( supers );
      for (int i = 0; i < supers.isize(); i++) {
	int scaffold = supers[i];
	ForceAssertLt(scaffold, scaffolds.isize());
	for (int j = 0; j < scaffolds[scaffold].Ntigs(); j++)
	  tigs.push_back(scaffolds[scaffold].Tig(j));
      }
    }
    UniqueSort( tigs );
  } else {
    if ( SINGLE_SCAFFOLD_BAMS ) SINGLE_SCAFFOLD_BAMS = False;
    if ( USE_SCAFFOLD_REFS ) USE_SCAFFOLD_REFS = False;
    ParseIntSet(TIGS, tigs);
    UniqueSort( tigs );
  }

  // Load scaffolds and SAM Reference sequence dictionary.

  vec<int> to_super, to_super_pos;
  vec< vec<int> > to_super_offset;
  vec<String> sam_header(1, "@HD\tVN:1.3\tSO:coordinate" );
  if ( USE_SCAFFOLD_REFS ) {
    int ntigs 
      = MastervecFileObjectCount( sub_dir + "/" + ASSEMBLY + ".contigs.fastb" );
    int header_idx = 0;
    to_super.resize( ntigs, -1 );
    to_super_pos.resize( ntigs, -1 );
    to_super_offset.resize( scaffolds.isize( ) );
    if ( supers.nonempty( ) ) {
      sam_header.resize( supers.isize( ) + 1 );
    } else {
      sam_header.resize( scaffolds.isize( ) + 1 );
    }
    for ( int i = 0; i < scaffolds.isize( ); i++ ) {
      int n = scaffolds[i].Ntigs( );
      int cur_offset = 0;
      to_super_offset[i].resize( n, -1 );
      for ( int j = 0; j < n; j++ ) {
	to_super[ scaffolds[i].Tig(j) ] = i;
	to_super_pos[ scaffolds[i].Tig(j) ] = j;
	to_super_offset[i][j] = cur_offset;
	int gap = ( j < n - 1 ) ? ( ( scaffolds[i].Gap(j) < (int)GAP_FLOOR ) ?
				    GAP_FLOOR : scaffolds[i].Gap(j) ) : 0;
	cur_offset += gap + scaffolds[i].Len(j);
      }
      if ( supers.nonempty( ) && BinMember(supers, i) ) {
	header_idx++;
	sam_header[header_idx] = "@SQ\tSN:scaffold_" + ToString( i ) + "\tLN:" +
	  ToString( cur_offset );
      } else {
	sam_header[i + 1] = "@SQ\tSN:scaffold_" + ToString( i ) + "\tLN:" +
	  ToString( cur_offset );
      }
    }
  } else {
    String contigs_file = sub_dir + "/" + ASSEMBLY + ".contigs.fastb";
    typedef VirtualMasterVec<basevector> VVecBVec;
    VVecBVec contigs( (contigs_file).c_str( ) );
    VVecBVec::const_iterator contigItr = contigs.begin( );
    int prev_tig = 0;
    sam_header.resize( tigs.isize( ) + 1 );
    for ( int i = 0; i < tigs.isize( ); i++ ) {
      contigItr = contigItr + ( tigs[i] - prev_tig );
      sam_header[i + 1] = "@SQ\tSN:contig_" + ToString( tigs[i] ) + "\tLN:" +
	ToString( (uint)(contigItr->size( )) );
      prev_tig = tigs[i];
    }
  }

  String head = sub_dir + "/" + ASSEMBLY;
  String sam_file =
    ( OUT_HEAD != "" ) ? ( OUT_HEAD + ".sam" ) :
    ( head + ( ( USE_SCAFFOLD_REFS ) ? ".super" : ".contig" ) + "s.sam" );
  ofstream sam;
  OpenOfstream( sam, sam_file );

  // Print Header Line and Sequence Dictionary
  for ( int i = 0; i < sam_header.isize( ); i++ )
    sam << sam_header[i] << "\n";

  // Process Read Groups
  if ( READLOCS_PREFIX != "" ) head += "." + READLOCS_PREFIX;
  read_locs_on_disk locs_file( head, run_dir );
  if (true) { // Force garbage collection
    vec<Bool> lib_types( 3, True );
    vec<String> pairs_files( 3, "" );
    for ( int pass = 0; pass < 3; pass++ ) {
      String file_name = run_dir + "/" +
	String( pass == 0 ? "frag" : ( pass == 1 ? "jump" : "long_jump" ) ) +
	"_reads_filt";
      file_name += String( IsRegularFile( file_name + "_cpd.pairs" ) ?
			   "_cpd" : "" ) + ".pairs";
      
      if ( ! IsRegularFile( file_name ) ) {
	cout << "Unable to find " << String( pass == 0 ? "frag" :
					     ( pass == 1 ? "jump" :
					       "long_jump" ) )
	     << " pairs file" << endl;
	lib_types[pass] = False;
      } else {
	pairs_files[pass] = file_name;
      }
    }

    vec< vec<int> > lib_inserts( 3 );
    vec< vec<int> > lib_seps( 3 );
    vec< vec<String> > read_groups( 3 );
    int lib_count = 0;

    for ( int pass = 0; pass < 3; pass++ ) {
      if ( lib_types[pass] ) {
	longlong nreads;
	vec<String> lib_names;
	vec<int> lib_sd;
	// These names match what is given by read_loc::ReadClassName( ), if
	// that ever changes, this will break.
	String type_name = pass == 0 ? "frag" : ( pass == 1 ? "jump" : "long" );
	ReadPairsManagerLibInfo( pairs_files[pass], nreads, lib_names,
				 lib_seps[pass], lib_sd );
	lib_count += lib_names.isize( );
	lib_inserts[pass].resize( lib_seps[pass].isize( ), -1 );
	read_groups[pass].resize( lib_names.isize( ) );
	for ( int i = 0; i < lib_names.isize( ); i++ ) {
	  read_groups[pass][i] = "@RG\tID:" + type_name + "." + lib_names[i] +
	    "\tPL:ILLUMINA\tSM:" + SAMPLE;
	}
      }
    }

    // Determine Insert/Fragment sizes
    for ( int64_t id = 0; id < locs_file.getLocsSize() && lib_count > 0;
	  id++ ) {
      read_loc rl;
      locs_file.LoadOne( id, rl );
      // The first should always evaluate True, but just in case...
      if ( lib_types[rl.ReadClass( )] &&
	   lib_inserts[rl.ReadClass( )][rl.LibraryId( )] < 0 ) {
	lib_count--;
	lib_inserts[rl.ReadClass( )][rl.LibraryId( )] =
	  rl.ReadLength( ) * 2 + lib_seps[rl.ReadClass( )][rl.LibraryId( )];
      }
    }
	
    // Print the Read Groups
    for ( int pass = 0; pass < 3; pass++ ) {
      if ( lib_types[pass] ) {
	for ( int i = 0; i < read_groups[pass].isize( ); i++ ) {
	  // Skip unassembled Read Groups
	  if ( lib_inserts[pass][i] != -1 ) {
	    sam << read_groups[pass][i] << "\tPI:" << lib_inserts[pass][i]
		<< "\n";
	  }
	}
      }
    }
    sam.flush( );
  }


  // Process contigs.

  cout << Date( ) << ": parsing " << tigs.isize( ) << " contigs" << endl;
  for ( int it = 0; it < tigs.isize( ); /*it++*/ ) {
    vec< vec<read_loc> > locs( NUM_THREADS );
    uint read_ct = 0;
    //Dot( cout, it );
    //int TIG = tigs[it];
    //vec<read_loc> locs;
    //locs_file.LoadContig( TIG, locs );
    while ( read_ct < 5000000 && it < tigs.isize( ) ) {
      Dot( cout, it );
      vec<read_loc> all_locs;
      int TIG = tigs[it];
      locs_file.LoadContig( TIG, all_locs );
      read_ct += (uint)( all_locs.size( ) );
      it++;
      for ( int j = 0; j < all_locs.isize( ); j++ ) {
	locs[j % NUM_THREADS].push_back( all_locs[j] );
      }
    }
    #pragma omp parallel
    {
      int id = omp_get_thread_num();
      ostringstream sam_aligns;
      PrintSAMAligns( sam_aligns, locs[id], to_super, to_super_pos,
		      to_super_offset, run_dir, sub_dir, ASSEMBLY );
      #pragma omp critical
      {
	sam << sam_aligns.str( );
      }
    }
    sam.flush( );
  }
  sam.close( );
  cout << endl;
  waitpid( pid, NULL, 0 );
  String bam_file = sam_file.Before( "sam" ) + "bam";
  String sorted_prefix = sam_file.Before( "sam" ) + "sorted";
  cout << "Generating BAM: " << bam_file << endl;
  String theCommand = "samtools calmd -S -b " + sam_file + " " + fasta_file +
    " > " + bam_file;
  RunCommand( theCommand );
  cout << "Removing intermediate SAM" << endl;
  Remove( sam_file );
  cout << "Sorting BAM" << endl;
  theCommand = "samtools sort -m " + ToString( SORT_MEM.Int( ) ) + " " +
    bam_file + " " + sorted_prefix;
  RunCommand( theCommand );
  Mv( sorted_prefix + ".bam", bam_file );
  cout << "Indexing BAM" << endl;
  theCommand = "samtools index " + bam_file;
  RunCommand( theCommand );
  cout << Date( ) << ": done" << endl;
}
