/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/// MarkUnique.  For a given genome, mark the starting positions of N-base perfect
/// reads that can be aligned uniquely to it.
///
/// If a read has a false placement having <= D mismatches, the read is counted
/// as not aligning uniquely.  Only alignments subsuming a K-base perfect match
/// are seen.
///
/// If ALT=True, we require that D = 0 and that N is divisible by 4, and use an
/// alternate method, which is faster but requires more memory.  Ambiguous bases are
/// treated incorrectly by this method.
///
/// Files GENOME.fastb and GENOME.lookup must be provided, as well as GENOME.fastamb
/// if there are ambiguous bases.
///
/// If DIFFERENT_LOOKUP is specified, use it instead of GENOME.lookup.
///
/// OUT = file to write answer to, as vecbitvector.
///
/// WORKDIR = work directory, must be specified (unless ALT=True)
///
/// Parallelized for blade farm.  If there are no failures on the farm, you can run 
/// the code from beginning to end by launching a single job - without specifying
/// PHASE or BATCH_ID.  Otherwise you may have to intervene by launching PHASE 2
/// jobs (with BATCH_ID specified), and conceivably by launching a PHASE 3 job to
/// finish up, if the original job has died.

#include "Basevector.h"
#include "Bitvector.h"
#include "FeudalMimic.h"
#include "MainTools.h"
#include "lookup/ImperfectLookup.h"
#include "lookup/LookAlign.h"
#include "math/HoInterval.h"
#include "paths/ReadsToPathsCoreX.h"
#include "feudal/BinaryStream.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(GENOME, "reads of length N from every position in "
          "GENOME.fastb will be looked up in GENOME.lookup and marked if they are "
          "unique");
     CommandArgument_String_OrDefault_Doc(DIFFERENT_LOOKUP, "",
          "if specified, overrides GENOME.lookup, reads are still generated from "
          "GENOME.fastb");
     CommandArgument_Int_Doc(N, "generate reads of length N");
     CommandArgument_Int_Doc(D,"Defines uniqueness criterion: perfect alignment is "
          "unique if second best has > D mismatches");
     CommandArgument_Bool_OrDefault(ALT, False);
     CommandArgument_Int_OrDefault(K, 12);
     CommandArgument_Int_OrDefault(PHASE, 1);
     CommandArgument_Int_OrDefault(READS_PER_BATCH, 1000000);
     CommandArgument_Int_OrDefault(BATCH_ID, -1);
     CommandArgument_String_OrDefault(WORKDIR, "");
     CommandArgument_String(OUT);
     CommandArgument_String_OrDefault(QUEUE, "");
     CommandArgument_Bool_OrDefault_Doc(SKIPDONE, False, "Skip directories with marks.done.  To only be used if started PHASE 2.");
     EndCommandArguments;

     if (ALT)
     {    ForceAssert( D == 0 && N % 4 == 0 );
          vecbasevector genome( GENOME + ".fastb" );
          vecKmerPath paths, pathsrc;
          vec<big_tagged_rpint> pathsdb;
          ReadsToPathsCoreY( genome, N, paths, pathsrc, pathsdb );
          vecbitvector cov;
          Mimic( genome, cov );
          vec< vec<int> > starts( paths.size( ) );
          for ( size_t i = 0; i < paths.size( ); i++ )
          {    int S = 0;
               for ( int j = 0; j < paths[i].NSegments( ); j++ )
               {    starts[i].push_back(S);
                    S += paths[i].Length(j);    }    }
          longlong START = 0;
          for ( int e = 0; e < pathsdb.isize( ); e++ )
          {    const big_tagged_rpint& E = pathsdb[e];
               if ( E.Rc( ) ) continue;
               longlong Estart = E.Start( ), Estop = E.Stop( );
               int t = E.ReadId( );
               int S = starts[t][ E.PathPos( ) ];
               if ( e > 0 && START > Estart )
               {    S += START - Estart;
                    Estart = START;    }
               if ( Estart > Estop ) continue;
               KmerPathInterval I( Estart, Estop );
               static vec<longlong> con;
               Contains( pathsdb, I, con );
               static vec<ho_interval> C;
               C.clear( );
               for ( int x = 0; x < con.isize( ); x++ )
               {    const big_tagged_rpint& t = pathsdb[ con[x] ];
                    int start = t.Start( ) - I.Start( );
                    int stop = t.Stop( ) + 1 - I.Start( );
                    if ( start < 0 ) start = 0;
                    if ( stop > I.Length( ) ) stop = I.Length( );
                    C.push( start, stop );    }
               static vec< pair<ho_interval,int> > C2;
               CondenseIntervals( I.Length( ), C, C2 );
               for ( int x = 0; x < C2.isize( ); x++ )
               {    if ( C2[x].second == 1 )
                    {    for ( int u = C2[x].first.Start( );
                              u < C2[x].first.Stop( ); u++ )
                         {    cov[t].Set( S + u, True );    }    }    }
               START = Estop + 1;    }
          cov.WriteAll(OUT);    
          exit(0);    }

     ForceAssert( WORKDIR != "" );

  if ( PHASE == 1 ) {    
    ForceAssertEq( BATCH_ID, -1 );
    // We do not check when SKIPDONE is true 
    // since we want to go through it again.
    if ( !SKIPDONE && IsDirectory(WORKDIR) ) {    
      cout << "WORKDIR=" << WORKDIR << " exists.\n";
      cout << "To avoid accidental overwrites, you must remove "
	   << "it yourself, then call this program again." << endl;
      exit(1);    }
    if (!SKIPDONE)
        Mkdir777(WORKDIR);
    vecbasevector genome( GENOME + ".fastb" );
    vecbitvector amb;
    if ( IsRegularFile( GENOME + ".fastamb" ) ) 
      amb.ReadAll( GENOME + ".fastamb" );
    int batch = 1;
    vecbasevector reads;
    for ( size_t i = 0; i < genome.size( ); i++ ) {
      for ( int j = 0; j <= genome[i].isize( ) - N; j++ ) {    
	if ( amb.size( ) > 0 ) {    
	  int k;
	  for ( k = 0; k < N; k++ )
	    if ( amb[i][j+k] ) break;
	  if ( k < N ) continue;    }
	static basevector b;
	b.SetToSubOf( genome[i], j, N );
	reads.push_back_reserve(b);    
	if ( reads.size( ) == static_cast<size_t>(READS_PER_BATCH) ) {
	  String dir = WORKDIR + "/batch_" + ToString(batch++);
          if ( !IsDirectory(dir) )
	     Mkdir777(dir);
          if (!SKIPDONE)
	     reads.WriteAll( dir + "/reads.fastb" );
	  reads.clear( );    
	}    
      }    
    }
    if ( !reads.empty( ) ) {    
      String dir = WORKDIR + "/batch_" + ToString(batch++);
      if ( !IsDirectory(dir) )
         Mkdir777(dir);
      if (!SKIPDONE)
         reads.WriteAll( dir + "/reads.fastb" );    
    }
    for ( int b = 1; b < batch; b++ ) {    
      String dir = WORKDIR + "/batch_" + ToString(b);
      String queue = ( QUEUE == "" ? String("") : "-q " + QUEUE + " " );
      String difflook = ( DIFFERENT_LOOKUP == "" ? String("") 
          : "DIFFERENT_LOOKUP=" + DIFFERENT_LOOKUP );
      String markstemp = "";
      if (SKIPDONE) markstemp = " SKIPDONE=True ";
      SystemSucceed( "bsub -o " + dir + "/MarkUnique.bsub.out " + queue + "-e " 
		     + dir + "/MarkUnique.bsub.err " + " \"MarkUnique GENOME=" 
		     + GENOME + " N=" + ToString(N) + " D=" + ToString(D) + " K=" 
		     + ToString(K) + " PHASE=2 READS_PER_BATCH=" 
		     + ToString(READS_PER_BATCH) + " BATCH_ID=" + ToString(b) 
                     + " " + difflook + markstemp
		     + " WORKDIR=" + WORKDIR + " OUT=" + OUT + " > " 
		     + dir + "/MarkUnique.out\"" );    
    }
    PHASE = 3;    
  }

  if ( PHASE == 2 ) {    
    ForceAssertGe( BATCH_ID, 0 );
    String bdir = WORKDIR + "/batch_" + ToString(BATCH_ID);
    if (!( (SKIPDONE == True) && (IsRegularFile( bdir + "/marks.done" )) )) {
       vecbasevector reads( bdir + "/reads.fastb" );
       vec<look_align> aligns;
       vec<int> min_errs;
       int max_errors = 0;
       vec<Bool> unique_aligned;
       String lookup 
            = ( DIFFERENT_LOOKUP == "" ? GENOME + ".lookup" : DIFFERENT_LOOKUP );
       ImperfectLookup( K, reads, lookup, aligns, min_errs, FW_OR_RC,
		        D, max_errors, &unique_aligned );
       BinaryWriter::writeFile( bdir + "/marks", unique_aligned );
       Echo( "done", bdir + "/marks.done" );    
    }
  }
          
  if ( PHASE == 3 ) {    
    ForceAssertEq( BATCH_ID, -1 );
    vecbasevector genome( GENOME + ".fastb" );
    vecbitvector amb, cov;
    if ( IsRegularFile( GENOME + ".fastamb" ) ) 
      amb.ReadAll( GENOME + ".fastamb" );
    Mimic( genome, cov );
    int batch = 1, nreads = 0;
    vec<Bool> unique_aligned;
    String bdir = WORKDIR + "/batch_" + ToString(batch);
    while( !IsRegularFile( bdir + "/marks.done" ) ) sleep(10);
    BinaryReader::readFile( bdir + "/marks", &unique_aligned );
    for ( size_t i = 0; i < genome.size( ); i++ ) {
      for ( int j = 0; j <= genome[i].isize( ) - N; j++ ) {    
	if ( amb.size( ) > 0 ) {    
	  int k;
	  for ( k = 0; k < N; k++ )
	    if ( amb[i][j+k] ) break;
	  if ( k < N ) continue;    
	}
	cov[i].Set( j, unique_aligned[nreads++] );
	if ( nreads == READS_PER_BATCH ) {    
	  Dot(cout, batch);
	  ++batch;
	  bdir = WORKDIR + "/batch_" + ToString(batch);
	  while( !IsRegularFile( bdir + "/marks.done" ) ) sleep(10); 
	  BinaryReader::readFile( bdir + "/marks", &unique_aligned );
	  nreads = 0;    
	}    
      }    
    }
    cov.WriteAll(OUT);    
  }
  cout << endl;
}
