///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeReadLocs.  Generate a file of read locations.  Experimental.  Does not
// include long jumps.

#include "Basevector.h"
#include "FetchReads.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "Vec.h"
#include "feudal/BinaryStream.h"
#include "paths/ReadLoc.h"
#include "paths/UnipathFixerTools.h"
#include "util/SearchFastb2Core.h"

Bool cmp_pos( const triple<int64_t,int,int>& x1, const triple<int64_t,int,int>& x2 )
{    if ( x1.second < x2.second ) return True;
     if ( x1.second > x2.second ) return False;
     if ( x1.third < x2.third ) return True;
     if ( x1.third > x2.third ) return False;
     return x1.first < x2.first;    }

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(SUBDIR, "test");
     CommandArgument_Int(K);
     CommandArgument_String_OrDefault(HEAD, "extended40.shaved");
     CommandArgument_String_OrDefault(SCAFFOLDS_IN, "linear_scaffolds0.patched");
     CommandArgument_String_OrDefault(UNIPATH_PATCH_DIR, "unipath_patch");
     CommandArgument_String_OrDefault(JUMP_READS, "jump_reads_filt_cpd");
     CommandArgument_String_OrDefault(FRAG_READS, "frag_reads_filt_cpd");
     EndCommandArguments;

     // Define directories, etc.

     String data_dir = PRE + "/" + DATA, run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
     cout << Date( ) << ": " << run_dir << endl;

     // Note some files that are needed.

     String unibases_file = run_dir + "/" + HEAD + ".unibases.k" + ToString(K);
     String reads_file = run_dir + "/" + FRAG_READS + ".fastb";
     String jreads_file = run_dir + "/" + JUMP_READS + ".fastb";

     // Get total number of contigs.

     size_t ntigs = MastervecFileObjectCount( 
          sub_dir + "/" + SCAFFOLDS_IN + ".contigs.vecfasta" );

     // Compute fragment read lengths.

     vec<uint16_t> read_len, jread_len;
     cout << Date( ) << ": computing fragment read lengths" << endl;
     {    vecbasevector reads(reads_file);
          read_len.resize( reads.size( ) );
          for ( size_t i = 0; i < reads.size( ); i++ )
          {    ForceAssertLt( reads[i].size( ), 65536u );
               read_len[i] = reads[i].size( );    }    }

     // Map the unibases to the contigs.  There is a potential problem with the way 
     // we're doing this: we require that each unibase maps perfectly.

     cout << Date( ) << ": mapping unibases to contigs" << endl;
     vec< triple<int64_t,int64_t,int> > UALIGNS;
     const int max_placements = 1;
     SearchFastb2( unibases_file, sub_dir + "/" + SCAFFOLDS_IN + ".contigs.fastb", 
          K, &UALIGNS, 0, max_placements );

     // Load the unibases.

     vecbasevector unibases( unibases_file );

     // Index the unibase alignments.

     cout << Date( ) << ": indexing unibase alignments" << endl;
     vec<size_t> U_START( unibases.size( ) + 1 );
     size_t POS = 0;
     for ( int64_t u = 0; u <= (int64_t) unibases.size( ); u++ )
     {    while( POS < UALIGNS.size( ) && UALIGNS[POS].first < u ) ++POS;
          U_START[u] = POS;    }

     // Map fragment reads to the contigs.  Note that we assume reads are fw on the
     // unibases.

     vec< triple<int64_t,int,int> > RALIGNS, JRALIGNS; // (rid, tig, pos)
     cout << Date( ) << ": loading segments" << endl;
     {    vec<segalign> SEGS; 
          String SEGS_file = run_dir + "/" + UNIPATH_PATCH_DIR + "/UnipathPatcher.SEGS";
          BinaryReader::readFile( SEGS_file, &SEGS );
          cout << Date( ) << ": start process of mapping reads to contigs" << endl;
          for ( size_t i = 0; i < SEGS.size( ); i++ )
          {    const segalign& a = SEGS[i];
               int64_t rid = a.rid; 
               int u = a.u, rpos = a.rpos, upos = a.upos;
               if ( U_START[u] < U_START[u+1] )
               {    int tig = UALIGNS[ U_START[u] ].second;
                    int tpos = UALIGNS[ U_START[u] ].third;
                    int read_start_on_tig;
                    if ( tpos >= 0 ) read_start_on_tig = tpos + upos - rpos;
                    else
                    {    read_start_on_tig = -tpos-1 + unibases[u].isize( ) - upos
                              - read_len[rid] + rpos;    }
                    read_start_on_tig 
                         = Max( 0, read_start_on_tig ); // don't like!
                    if ( tpos >= 0 ) RALIGNS.push(rid, tig, read_start_on_tig);
                    else RALIGNS.push(rid, tig, -read_start_on_tig-1);    }    }    }
     ParallelUniqueSort(RALIGNS);
     cout << Date( ) << ": found " << RALIGNS.size( )
          << " alignments of fragment reads" << endl;

     // Build read locations for fragment reads.

     vec<read_loc> raligns;
     int64_t raligns_count = 0;
     int64_t unpaired_frag_count = 0, total_frag_count = 0;
     cout << Date( ) << ": generating fragment read locations" << endl;
     {    ParallelSort( RALIGNS, cmp_pos );
          PairsManager pairs( run_dir + "/" + FRAG_READS + ".pairs" );
          raligns.resize( RALIGNS.size( ) );
          for ( size_t i = 0; i < raligns.size( ); i++ )
          {    int pos = RALIGNS[i].third;
               if ( pos < 0 ) pos = -pos-1;
               int64_t rid = RALIGNS[i].first;
               int64_t pid = pairs.getPairID(rid);
               total_frag_count++;
               if ( pid < 0 )
               {    unpaired_frag_count++;
                    continue;    }
               raligns[raligns_count++] = read_loc( rid, RALIGNS[i].second, 
                    pos, RALIGNS[i].third >= 0, 0, pairs.libraryID(pid),
                    read_len[rid] );    }    }
     cout << Date( ) << ": " << PERCENT_RATIO( 3, unpaired_frag_count,
          total_frag_count ) << " of fragment reads are unpaired" << endl;
     Destroy(RALIGNS), Destroy(read_len);

     // Compute jump read lengths.

     cout << Date( ) << ": computing jump read lengths" << endl;
     {    vecbasevector jreads(jreads_file);
          jread_len.resize( jreads.size( ) );
          for ( size_t i = 0; i < jreads.size( ); i++ )
          {    ForceAssertLt( jreads[i].size( ), 65536u );
               jread_len[i] = jreads[i].size( );    }    }

     // Map jump reads to the contigs.  Note that we assume reads are fw on the
     // unibases.

     {    vec<segalign> JSEGS;
          String JSEGS_file = run_dir + "/" + UNIPATH_PATCH_DIR + "/UnipathPatcher.JSEGS";
          BinaryReader::readFile( JSEGS_file, &JSEGS );
          cout << Date( ) << ": mapping jump reads" << endl;
          for ( size_t i = 0; i < JSEGS.size( ); i++ )
          {    const segalign& a = JSEGS[i];
               int64_t rid = a.rid; 
               int u = a.u, rpos = a.rpos, upos = a.upos;
               if ( U_START[u] < U_START[u+1] )
               {    int tig = UALIGNS[ U_START[u] ].second;
                    int tpos = UALIGNS[ U_START[u] ].third;
                    int read_start_on_tig;
                    if ( tpos >= 0 ) read_start_on_tig = tpos + upos - rpos;
                    else
                    {    read_start_on_tig = -tpos-1 + unibases[u].isize( ) - upos
                              - jread_len[rid] + rpos;    }
                    read_start_on_tig = Max( 0, read_start_on_tig ); // don't like!
                    if ( tpos >= 0 ) JRALIGNS.push(rid, tig, read_start_on_tig);
                    else JRALIGNS.push(rid, tig, -read_start_on_tig-1);   }   }   }
     ParallelUniqueSort(JRALIGNS);
     Destroy(unibases), Destroy(U_START), Destroy(UALIGNS);

     // Build read locations for jump reads.

     int64_t unpaired_jump_count = 0, total_jump_count = 0;
     cout << Date( ) << ": generating jump read locations" << endl;
     {    ParallelSort( JRALIGNS, cmp_pos );
          PairsManager jpairs( run_dir + "/" + JUMP_READS + ".pairs" );
          size_t frag_count = raligns.size( );
          raligns.resize( frag_count + JRALIGNS.size( ) );
          for ( size_t i = 0; i < JRALIGNS.size( ); i++ )
          {    int pos = JRALIGNS[i].third;
               if ( pos < 0 ) pos = -pos-1;
               int64_t rid = JRALIGNS[i].first;
               int64_t pid = jpairs.getPairID(rid);
               total_jump_count++;
               if ( pid < 0 )
               {    unpaired_jump_count++;
                    continue;    }
               raligns[raligns_count++] = read_loc( rid, JRALIGNS[i].second, 
                    pos, JRALIGNS[i].third < 0, 1,
                    jpairs.libraryID(pid), jread_len[rid] );    }    }
     cout << Date( ) << ": " << PERCENT_RATIO( 3, unpaired_jump_count,
          total_jump_count ) << " of jump reads are unpaired" << endl;
     Destroy(JRALIGNS), Destroy(jread_len);
     raligns.resize(raligns_count);

     // Sort and write read locations.

     ParallelSort(raligns);
     WriteReadLocs( sub_dir + "/" + SCAFFOLDS_IN, ntigs, raligns );    }
