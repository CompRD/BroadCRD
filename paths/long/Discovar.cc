///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "TokenizeString.h"
#include "efasta/EfastaTools.h"
#include "paths/long/DataSpec.h"
#include "paths/long/DiscovarTools.h"
#include "paths/long/Heuristics.h"
#include "paths/long/ImproveLongHyper.h"
#include "paths/long/LoadCorrectCore.h"
#include "paths/long/Logging.h"
#include "paths/long/LongHyper.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/PairInfo.h"
#include "paths/long/RefTrace.h"
#include "paths/long/ReadOriginTracker.h"
#include "paths/long/SupportedHyperBasevector.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(READS,
          "Comma-separated list of one or more bam files, each ending in .bam.\n"
          "Alternatively, this may have the form @fn, where fn is a file\n"
          "containing a list of bam file names, one per line.");
     CommandArgument_String_Doc(REGIONS,
          "Regions to be extracted from bam files: a comma-separated list of one "
          "or more region specifications chr:start-stop, where chr is a chromosome "
          "name (consistent with usage in the bam files), and start-stop defines "
          "a range of bases on chr (zero based).  If REGIONS = all, "
          "bam files will be used in their entirety.");
     CommandArgument_String_Doc(TMP, "Directory to put temporary files in.");
     CommandArgument_String_Doc(OUT_HEAD, "Full path prefix for output files.");
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0,
          "Number of threads to use (use all available processors if set to 0).");
     CommandArgument_String_OrDefault_Doc(REFERENCE, "",
          "FASTA file containing reference - used for variant calling.");
     CommandArgument_Bool_OrDefault_Doc(STATUS_LOGGING, False,
          "if set to True, generate cryptic logging that reports on the status of "
          "intermediate calculations");
     CommandArgument_Bool_OrDefault(USE_OLD_LRP_METHOD, True);
     CommandArgument_Bool_OrDefault_Doc(DRY_RUN, False,
          "Set to True for a dry run to check input parameters.");
     CommandArgument_LongLong_OrDefault_Doc(MAX_MEMORY_GB,0,
          "Try not to use more than this amount of memory.");
     EndCommandArguments;

     SetThreads(NUM_THREADS);
     SetMaxMemory(MAX_MEMORY_GB<<30);

     // Note time and check requirements.

     double clock = WallClockTime( );
     DiscovarTools::CheckDiscovarSystemRequirements( );

     // Check READS and REGIONS arguments.

     DiscovarTools::CheckDiscovarRegionsInput(REGIONS);
     if ( REGIONS == "all" ) REGIONS = "";
     vec<String> bams;
     DiscovarTools::CheckDiscovarReads(READS);
     if ( READS.Contains( "," ) && !READS.Contains( "{", 0 ) )
          READS = "{" + READS + "}";
     ParseStringSet( READS, bams );
     DiscovarTools::CheckDiscovarBams(bams);
     DiscovarTools::CheckDiscovarRegions(REGIONS);
     DiscovarTools::CheckDiscovarOutHead(OUT_HEAD);
     DiscovarTools::CheckDiscovarTmp(TMP);
     DiscovarTools::TestDiscovarRegionsBamsCompatibility( REGIONS, bams );

     // Handle TMP. Clear out old frag_reads_orig files.

     if ( !IsDirectory(TMP) )  Mkdir777(TMP);
     else  System("rm -f " + TMP + "/frag_reads_orig.*");
 
     // Set up logging etc.

     ref_data ref;
     long_logging logc( "", "" );
     logc.STATUS_LOGGING = STATUS_LOGGING;
     long_data_spec spec( "" );

     DiscovarTools::CheckReferenceInput(REFERENCE,OUT_HEAD);
     DiscovarTools::DiscovarRefTraceControl ref_trace_control(REFERENCE,REGIONS,logc,OUT_HEAD+".final.variant");
     String sHeurControl="";

     if( ref_trace_control.IsRefTraceOn() ){
         ref_trace_control.CheckConsistency();
         if (logc.STATUS_LOGGING)
         {    std::cout << Date() << ": converting sequences from " << REFERENCE 
                  << " to internal data structure" << std::endl;    }
         //note that BuildRefData was never called to initialized 'ref's data before variant calls are added
         BuildRefData_Discovar( ref, ref_trace_control ); //note that TestDiscovarRegionsBamsCompatibility must be run before this line
         logc.REFTRACE="True";
         logc.REFTRACE_VARIANTS=True;
         ForceAssert(sHeurControl==""); // change this nicely if one wants other things for heur control
         sHeurControl = "{}";
         ref_trace_control.getRefHead() = REFERENCE.Before(".fasta");
         ref_trace_control.setReadHead(TMP+"/frag_reads_orig");
     }

     vec<ref_loc> readlocs;
     long_logging_control log_control( ref, &readlocs, OUT_HEAD, "" );
     long_heuristics heur( sHeurControl );

     // Load and correct reads.

     vec<String> heads;
     for ( int i = 0; i < bams.isize( ); i++ )
     {    String regions = REGIONS;
          regions.GlobalReplaceBy( ",", " " );
          String getsam = "samtools view -h " + bams[i] + " " + regions;
          SamIAm( i, getsam, TMP );
          heads.push_back( ToString(i) );    }

     if(DRY_RUN){
         std::cout << "\nDry run has been completed successfully.\n" << std::endl;
         PrintTail( command, "", "", clock, logc );
         return 0;
     }

     MergeReadSets( heads, TMP, logc );
     SelectRandom( TMP, 1, logc, spec );
     VecEFasta corrected;
     vec<pairing_info> cpartner;
     vec<int> cid;
     vecbasevector creads;
     CorrectionSuite( TMP, heur, logc, log_control, creads, corrected, cid, 
          cpartner, NUM_THREADS, "", clock, USE_OLD_LRP_METHOD );
     vec<Bool> to_delete( creads.size( ), False );
     DefinePairingInfo( TMP, creads, to_delete, cid, corrected, cpartner, logc );

     // Build and simplify graph.

     SupportedHyperBasevector shb;
     if ( !LongHyper( READS, corrected, cpartner, shb, heur,
          log_control, logc, TMP, USE_OLD_LRP_METHOD ) )
     {    DiscovarTools::ExitPathsEmpty( );    }
     ReportAssemblyStats(shb);
     ImproveLongHyper( "", "", shb, TMP, READS, spec, heur, NUM_THREADS, 
          log_control, logc, 2, USE_OLD_LRP_METHOD );

     // Generate output files.

     shb.DumpFiles( OUT_HEAD+".final", log_control, logc );
     Ofstream( fout, OUT_HEAD + ".final.fasta0" );
     vec<int> to_left, to_right;
     shb.ToLeft(to_left), shb.ToRight(to_right);
     for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
     {    int v = to_left[e], w = to_right[e];
          fout << ">edge_" << e << " " << v << ":" << w
               << " bases=" << shb.EdgeLengthBases(e) - ( shb.K( ) - 1 ) << "\n";
          basevector b = shb.EdgeObject(e);
          if ( shb.From(w).nonempty( ) ) b.resize( b.isize( ) - ( shb.K( ) - 1 ) );
          b.Print(fout);    }


     //call ref trace if needed
     if( ref_trace_control.IsRefTraceOn() ){
         if (logc.STATUS_LOGGING)
              std::cout << Date() << ": Performing reference tracing" << std::endl;
         ref_trace_control.ReadySampleLookup();
         ReadOriginTracker read_tracker(ref_trace_control);
         Ofstream( callsOut, OUT_HEAD + ".final.calls" );
         RefTraceAndCallVaraint( ref, shb, shb.Inv( ),
                   logc.verb[ "REFTRACE" ], logc,
                   cout, callsOut, RefTraceHeuristics(), "", False,
                   ref_trace_control,
                   &read_tracker );
         callsOut.close();
     }

     // Print a brief report.

     cout << "\n=================================================================="
          << "==================\n\n";
     cout << "DISCOVAR SUMMARY STATS\n\n";
     cout << shb.NComponents( ) << " components" << endl;
     cout << shb.EdgeObjectCount( ) << " edges" << endl;
     int64_t total_kmers = 0;
     for ( int i = 0; i < shb.EdgeObjectCount( ); i++ )
          total_kmers += shb.EdgeLengthKmers(i);
     cout << total_kmers << " kmers" << endl;
     PrintTail( command, "", "", clock, logc );    }
