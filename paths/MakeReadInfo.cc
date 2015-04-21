///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeReadInfo.  Generate a file showing information about reads, organized
// by position on the genome.  This only works with simulations or assemblies
// from reads picked out of a BAM.

#include "FastIfstream.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "ReadLocationLG.h"
#include "lookup/LookAlign.h"
#include "util/ReadTracker.h"
#include "feudal/BinaryStream.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_Bool_OrDefault(FIX_ORIGIN, True);
     EndCommandArguments;

     // Define directories.

     String data_dir = PRE + "/" + DATA;
     String run_dir = PRE + "/" + DATA + "/" + RUN;

     // Load read locations.  Two test cases at present.

     vec<ReadLocationLG> locs;
     if ( IsRegularFile( data_dir + "/frag_reads_orig.ref.locs" ) )
          BinaryReader::readFile( data_dir + "/frag_reads_orig.ref.locs", &locs );
     else if ( IsRegularFile( data_dir + "/../region" )
          && IsRegularFile( data_dir + "/frag_reads_orig.qltout" ) )
     {    
          // Define region, then convert alignments to read locations.  The handling
          // of tig_id and origin is not really right, and likely to break.

          String region = FirstLineOfFile( data_dir + "/../region" );
          int tig_id = region.Before( ":" ).Int( );
          int origin = region.Between( ":", "-" ).Int( );
          if ( !FIX_ORIGIN ) origin = 0;
          vec<look_align> aligns;
          LoadLookAligns( data_dir + "/frag_reads_orig.qltout", aligns );
          for ( int i = 0; i < aligns.isize( ); i++ )
          {    const look_align& la = aligns[i];
               locs.push( la.query_id, la.target_id - tig_id, la.pos2( ) - origin,
                    ( la.Fw1( ) ? ORIENT_FW : ORIENT_RC ) );    }    }
     else
     {    cout << "The files needed to place reads on the reference are "
               << "unavailable." << endl;
          cout << "Abort." << endl;
          return 1;    }
     Sort(locs);

     // Load read pairs.

     PairsManager pairs( data_dir + "/frag_reads_orig.pairs" );
     PairsManager pairs_corr( run_dir + "/frag_reads_corr.pairs" );

     // Get number of reads.

     size_t nreads = MastervecFileObjectCount( data_dir + "/frag_reads_orig.fastb" );
     size_t nreads_corr 
          = MastervecFileObjectCount( run_dir + "/frag_reads_corr.fastb" );

     // Generate map of reads to "corr" reads.

     vec<int> id_corr( nreads, -1 ); 
     {    ReadTracker rt;
          rt.Load( run_dir + "/frag_reads_corr" );
          for ( size_t r = 0; r < nreads_corr; r++ )
               id_corr[ rt.GetReadIndex(r) ] = r;    }

     // Show which reads are corrected.

     vec<Bool> corrected( nreads, False );
     {    for ( size_t i = 0; i < nreads; i++ )
               corrected[i] = id_corr[i] >= 0;    }

     // Index read locations.

     vec<int> locs_index( nreads, -1 );
     {    for ( int i = 0; i < locs.isize( ); i++ )
          {    const ReadLocationLG& l = locs[i];
               locs_index[ l.ReadId( ) ] = i;    }    }

     // Figure out which reads went into filled pairs.

     vec<int> filled( nreads_corr, -1 );
     {    String line;
          fast_ifstream fin( run_dir + "/filled_reads.info" );
          while(1)
          {    getline( fin, line );
               if ( fin.fail( ) ) break;
               int pid = line.Before( " " ).Int( ); // .corr pair id
               int frag_id = line.Between( " ", " " ).Int( ); // filled frag id
               int cid1 = pairs_corr.ID1(pid), cid2 = pairs_corr.ID2(pid);
               filled[cid1] = filled[cid2] = frag_id;    }    }

     // Print read info.
     
     Ofstream( out, run_dir + "/frag_reads_orig.info" );
     for ( int i = 0; i < locs.isize( ); i++ )
     {    const ReadLocationLG& l1 = locs[i];
          int id1 = l1.ReadId( );
          if ( !l1.Fw( ) ) continue;
          int id2 = pairs.getPartnerID(id1);
          out << l1.ReadId( ) << ( l1.Fw( ) ? "fw" : "rc" ) << " at " 
               << l1.Contig( ) << "." << l1.Start( ) << "-" 
               << l1.Start( ) + 101 << " --> ";
          if ( id2 >= 0 ) out << id2;
          else out << "no partner";
          if ( id2 >= 0 && locs_index[id2] >= 0 )
          {    const ReadLocationLG& l2 = locs[ locs_index[id2] ];
               out << ( l2.Fw( ) ? "fw" : "rc" )
                    << " at " << l2.Contig( ) << "." << l2.Start( ) 
                    << "-" << l2.Start( ) + 101;    }
          out << " " << ( corrected[id1] ? "C" : "U" );
          if ( id2 >= 0 ) out << "/" << ( corrected[id2] ? "C" : "U" );
          if ( id2 >= 0 && corrected[id1] && corrected[id2] 
               && filled[ id_corr[id1] ] >= 0 ) 
          {    out << " filled/" << filled[ id_corr[id1] ];    }
          out << "\n";    }    }
