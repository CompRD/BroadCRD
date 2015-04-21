///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// LongProto is the internal version of Discovar.  Both LongProto.cc and Discovar.cc
// should be as short and synchronous as possible.  We will have a nightly test
// to ensure that the assemblies they generate are identical.
//
// There is a google doc that describes the LongProto approach and targets for
// further improvement.  We run a nightly dexter 'longread' test to assay 
// performance.

// MakeDepend: private

#include "Basevector.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "ParseSet.h"
#include "paths/HyperEfasta.h"
#include "paths/long/CleanEfasta.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/DataSpec.h"
#include "paths/long/EvalByReads.h"
#include "paths/long/EvalCorrected.h"
#include "paths/long/Heuristics.h"
#include "paths/long/ImproveLongHyper.h"
#include "paths/long/LoadAndCorrect.h"
#include "paths/long/LongHyper.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/MakeSimReads.h"
#include "paths/long/PairInfo.h"
#include "paths/long/Simulation.h"
#include "paths/long/ValidateRef.h"
#include "paths/long/Variants.h"
#include "paths/long/DiscovarTools.h"
#include "paths/long/RefTrace.h"
#include "paths/long/ultra/MakeCorrectedReads.h"
#include "fastg/FastgGraph.h"
#include "paths/long/EvalAssembly.h"

int main(int argc, char *argv[])
{
     RunTime( );
     double clock = WallClockTime( );

     BeginCommandArguments;

     // CONTROL OVER REGION:

     CommandArgument_String_OrDefault_Doc(SAMPLE, "", "see local/SAMPLE");
     CommandArgument_String_OrDefault_Doc(X, "", "see local/X");
     CommandArgument_Double_OrDefault_Doc(GENOME_SUB_PERCENT, 0.0,
          "if specified, create a diploid genome by introducing substitutions "
          "at the specified rate; only makes sense for simulated reads");
     CommandArgument_String_OrDefault_Doc(DATA_SPEC, "",
          "list of one or more data specification parameters, in the form arg=val "
          "or {arg1=val1,...,argn=valn} -- see DataSpec.h for list of arguments");

     // CONTROL OVER READS:

     CommandArgument_String_OrDefault_Doc(READS, "",
          "if specified, use these reads (fastb) instead of making simulated "
          "reads; comma-separated list of one or more BAM files is also allowed "
          "(use of @file reads a list of .bam files); "
          " list of fastq files are also allowed: if 1 fastq file is supplied, read 2i and 2i+1, are treated as pairs. "
          " if an even number of file is supplied, the reads from file 2i and 2i+1 are treated as pairs. "
          "in the special case where READS=#picard, use a predefined set "
          "of real Illumina reads and error correct them using the ALLPATHS-LG "
          "process -- such data are available only when SAMPLE=rhody, plasmo, "
          "human, hpool1, hpool2, hpool3, human.hpool2 or human.hpool3, and one "
          "must specify the TMP argument; usable only by development team");
     CommandArgument_String_OrDefault_Doc(DATASET, "", "comma-separated list; for "
          "certain SAMPLEs, this may be specified to further define the dataset");
     CommandArgument_Bool_OrDefault_Doc(HAVE_PICARD, False, 
          "if specified and using READS=#picard, assume that Illumina reads have "
          "already been loaded, initially error corrected, and written to temp "
          "files");
     CommandArgument_String_OrDefault_Doc(RID, "",
          "if specified, only correct the reads described by this list "
          " (in ParseIntSet format)");
     CommandArgument_String_OrDefault_Doc(SIM, "",
          "list of one or more simulation arguments, in the form arg=val or "
          "{arg1=val1,...,argn=valn} -- see Simulation.h for list of arguments");

     // ALGORITHMIC OPTIONS:

     CommandArgument_String_OrDefault_Doc(HEURISTICS, "",
          "list of one or more algorithmic heuristics, in the form arg=val or "
          "{arg1=val1,...,argn=valn} -- see Heuristics.h for list of arguments");

     // LOGGING OPTIONS:

     CommandArgument_String_OrDefault_Doc(VERB, "",
          "to control verbosity level of particular modules, may set to x1 or "
          "{x1,...,xn} where each xi is of the form flag:value (see doc/verb_flags "
          "for list of allowed flags) and value is a nonnegative integer, however "
          "setting to zero would do nothing");
     CommandArgument_String_OrDefault_Doc(LOGGING, "",
          "list of one or more logging options, in the form arg=val or "
          "{arg1=val1,...,argn=valn} -- see Logging.h for list of arguments");
     CommandArgument_String_OrDefault_Doc(IN_GENOME, "", 
          "fasta or fastb or hbv file for genome to use as reference sequence");
     CommandArgument_Bool_OrDefault_Doc(KEEP_NAMES, False, "keep read names");

     // SPECIAL OPTIONS FOR RUNNING ONLY PART OF THE CODE (SEE ALSO OUTPUT OPTIONS):

     CommandArgument_String_OrDefault_Doc(IN_EFASTA_READS, "", 
          "optional input file for corrected reads; can be set to\n"
          "#dexter:sample:date where sample is a case-insensitive dexter longread "
          "test name and date is the date of a dexter run; this will also port "
          "values for region definition variables; "
          "usable only by development team");
     CommandArgument_String_OrDefault_Doc(IN_SHBV, "",
          "optional input file for SupportedHyperBasevector; can be set to\n"
          "#dexter:sample:date as in IN_EFASTA_READS; if of the form *.n.shbv "
          "will start graph simplification at step n+1.");
     CommandArgument_Int_OrDefault_Doc(START_STEP, 1, "if set, force start step "
          "for graph simplification");
     CommandArgument_String_OrDefault_Doc(IN_SHBV_FINAL, "",
          "optional input file for final SupportedHyperBasevector; can be set to\n"
          "#dexter:sample:date -- see IN_EFASTA_READS");
     CommandArgument_String_OrDefault_Doc(EXIT, "", "can be set to CORRECT, "
          "IMPROVE, LOAD, LONG_HYPER, MAIN, NOMINAL_COV, OUT_EFASTA_READS, "
          "OUT_GENOME, OUT_SIM_READS or REFTRACE");

     // OUTPUT OPTIONS:

     CommandArgument_String_OrDefault_Doc(OUT_GENOME, "", 
          "optional fastb output file for genome");
     CommandArgument_String_OrDefault_Doc(OUT_SIM_READS, "",
          "optional fastb output file for simulated reads");
     CommandArgument_String_OrDefault_Doc(OUT_EFASTA_READS, "", 
          "optional output file for corrected reads");
     CommandArgument_String_OrDefault_Doc(OUT_INT_HEAD, "", 
          "generate intermediate assembly files "
          "OUT_INT_HEAD.n.{shbv,dot,support,fasta} for n = 1,2...,final");
     CommandArgument_String_OrDefault_Doc(OUT_HEAD, "", 
          "output assembly files: OUT_HEAD.efasta, OUT_HEAD.hbv and OUT_HEAD.dot");
     
     // OTHER:
     CommandArgument_Bool_OrDefault(USE_OLD_LRP_METHOD,True);

     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0,
          "Number of threads to use (use all available processors if set to 0)");
     CommandArgument_String_OrDefault_Doc(TMP, "", 
          "scratch directory, only used when READS=#picard");

     CommandArgument_LongLong_OrDefault_Doc(MAX_MEMORY_GB,0,
          "Try not to use more than this amount of memory.");
     EndCommandArguments;

     // Initial argument handling.
     
     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads( NUM_THREADS );
     SetMaxMemory(MAX_MEMORY_GB<<30);
     DiscovarTools::CheckDiscovarSystemRequirements( );
     long_logging logc( LOGGING, VERB );
     long_data_spec spec( DATA_SPEC );

     if(READS.Contains("@")){
         vec<String> bams;
         ParseStringSet( READS, bams );
         for(const auto& entry: bams){ForceAssert(entry.Contains(".bam",String::npos));}
         if( bams.size()==0){
             std::cout<<"Warning: zero entries parsed from " << READS << std::endl;
             READS="";
         }
         else{
             READS=bams[0];
             for(size_t ii=1;ii<bams.size();++ii){ READS+=","+bams[ii]; }
         }
     }
     if ( GENOME_SUB_PERCENT > 0.0 && READS != "" )
          FAIL_MSG( "If you specify READS, you can't use GENOME_SUB_PERCENT." );
     if ( READS == "#picard" && TMP == "" )
          FAIL_MSG( "If READS=#picard, TMP must be specified." );
     DexterizeAll( IN_SHBV, IN_SHBV_FINAL, IN_EFASTA_READS, START_STEP, SAMPLE, X, 
          GENOME_SUB_PERCENT, TMP, READS );
     if ( spec.HUMAN_CONTROLS != "" ) logc.REFTRACE = "True";
     logc.REFTRACE_VARIANTS_SUMMARY=True;
     TestSampleAndReads( SAMPLE, READS, logc );

     // Create genome.

     ref_data ref;
     String sVOut;
     if( OUT_INT_HEAD!="") { sVOut = OUT_INT_HEAD+".final.variant"; }
     else if( OUT_HEAD!=""){ sVOut = OUT_HEAD+".variant"; }
     else                  { sVOut = "variant"; }
     RefTraceControl RTCtrl(sVOut);
     if( READS=="#picard" || READS.Contains(".bam",-1)) RTCtrl.setReadHead(TMP+"/frag_reads_orig");

     String X_actual;
     if ( SAMPLE != "unknown" && !logc.TREAT_AS_UNKNOWN )
     {    BuildRefData( IN_GENOME, SAMPLE, X, spec.HUMAN_CONTROLS, X_actual, ref,
               GENOME_SUB_PERCENT, logc, RTCtrl );
          if ( OUT_GENOME != "" ) ref.G.WriteAll(OUT_GENOME);
          if ( EXIT == "OUT_GENOME" ) Done(clock);    }

     // Define data structure for read locations.  This is defined for simulated
     // reads and partially defined for reads that come from a bam file, provided
     // that KEEP_LOCS = True.

     vec<ref_loc> readlocs;

     // Set up logging control, heuristics, and simulation parameters.

     long_logging_control log_control( ref, &readlocs, OUT_INT_HEAD, VERB );
     long_heuristics heur( HEURISTICS );
     if ( SAMPLE == "rhody" || SAMPLE == "tb" ) heur.INJECT_REF = False;
     long_sim sim( SIM );

     // Load PacBio data if needed.

     if ( heur.PACBIO_PATCH && ( SAMPLE == "hpool2" || SAMPLE == "hpool3" ) )
     {    vecbasevector pb;
          LoadPacBio( X, TMP, pb, log_control, logc, spec.PACBIO_MORE );
          pb.WriteAll( TMP + "/pb.fastb" );    }
     if ( EXIT == "PB" ) Done(clock);

     // Load corrected reads if provided.  Load and error correct Illumina reads 
     // if requested.

     VecEFasta corrected;
     vec<pairing_info> cpartner;
     vec<int> cid;
     Bool have_corrected = ( IN_EFASTA_READS != "" 
          || READS == "#picard" || READS.Contains( ".bam", -1 )
          || READS.Contains( ".fastq", -1 ) );
     REPORT_TIME( clock, "used in initial setup" );
     if ( IN_EFASTA_READS != "" )
          LoadEfastaLong( IN_EFASTA_READS, corrected, cid, cpartner, logc );
     else if ( ( READS == "#picard" || READS.Contains( ".bam", -1 ) || READS.Contains( ".fastq", -1 ) )
          && IN_SHBV == "" && IN_SHBV_FINAL == "" )
     {    LoadAndCorrectIllumina( SAMPLE, READS, DATASET, X, TMP, spec, corrected, 
               cid, cpartner, HAVE_PICARD, heur, log_control, logc, NUM_THREADS, 
               EXIT, clock, USE_OLD_LRP_METHOD, KEEP_NAMES );    }
     if ( EXIT == "CORRECT" ) Done(clock);

     // Create reads and correct them.

     vecbasevector reads;
     if ( !have_corrected && IN_SHBV == "" && IN_SHBV_FINAL == "" )
     {    
          // Generate simulated reads.
     
          double sclock = WallClockTime( );
          if ( READS != "" ) reads.ReadAll(READS);
          else
          {    cout << Date( ) << ": making SIMULATED reads" << endl;
               MakeSimReads( sim, ref, reads, readlocs );
               if ( OUT_SIM_READS != "" ) reads.WriteAll(OUT_SIM_READS);    
               if ( EXIT == "OUT_SIM_READS" ) Done(clock);
               if (logc.DUMP_SIM_LOCS) DumpSimLocs( reads, readlocs );
               FilterSimReads( sim, ref.G, readlocs, RID );    }
          int N = reads.size( );
          DATE_MSG( "there are " << ToStringAddCommas(N) << " reads" );
          if ( RID.Contains( "random:", 0 ) )
          {    int n = RID.After( "random:" ).Int( );
               if ( n > N )
               {    FAIL_MSG("\nYou've requested " << n << " random reads, but "
                              "there are only " << N << " in total.");    }    }
          REPORT_TIME( sclock, "used getting initial reads" );

          // Make corrected reads.  For now they are declared unpaired.

          MakeCorrectedReads( sim, reads, RID, NUM_THREADS, heur, 
               log_control, logc, ref, corrected, cid, TMP );
          cpartner.resize( corrected.size( ), pairing_info( 0, -1, -1 ) );    }

     // Write and evaluate corrected reads.  

     if ( IN_SHBV == "" && IN_SHBV_FINAL == "" )
     {    if ( OUT_EFASTA_READS != "" )
          {    Ofstream( out, OUT_EFASTA_READS )
               for ( size_t i = 0; i < corrected.size( ); i++ )
                    corrected[i].Print( out, "corrected_" + ToString(cid[i]) );    }
          if ( EXIT == "OUT_EFASTA_READS" ) Done(clock);
          if (logc.EVAL_UNCORRECTED)
          {    VecEFasta uncorrected;
               uncorrected.reserve(cid.size());
               for ( int i = 0; i < cid.isize( ); i++ )
                    uncorrected.push_back( efasta( reads[ cid[i] ] ) );
               cout << "\nEvaluation of uncorrected reads:\n\n";
               EvalCorrected( uncorrected, cid, ref, log_control, logc );    }
          PrintCorrectedReadStats(corrected);
          if (logc.EVAL_CORRECTED)
          {    cout << "\nEvaluation of corrected reads:\n\n";
               EvalCorrected( corrected, cid, ref, log_control, logc );    }    }
     if ( EXIT == "MAIN" ) Done(clock);

     // Load SupportedHyperBasevector, create and improve HyperBasevector assembly.

     SupportedHyperBasevector shb;
     if ( IN_SHBV != "" ) BinaryReader::readFile( IN_SHBV, &shb );
     if ( IN_SHBV_FINAL != "" ) BinaryReader::readFile( IN_SHBV_FINAL, &shb );
     if ( IN_SHBV == "" && IN_SHBV_FINAL == "" ) 
     {    if ( !LongHyper( READS, corrected, cpartner, shb, heur,
               log_control, logc, TMP, USE_OLD_LRP_METHOD ) )
          {    DiscovarTools::ExitPathsEmpty( );    }    }
     ReportAssemblyStats(shb);
     if ( EXIT == "LONG_HYPER" ) Done(clock);
     if ( IN_SHBV_FINAL == "" )
     {    ImproveLongHyper( SAMPLE, X, shb, TMP, READS, spec, heur, NUM_THREADS, 
               log_control, logc, StartStep(IN_SHBV, START_STEP),
                                           USE_OLD_LRP_METHOD );    }
     if ( EXIT == "IMPROVE" ) Done(clock);

     // Detect variants, trace edges, analyze assembly.

     vec<VariantSignature> snp_bubbles;
     if ( READS == "#picard" )
     {    vecbasevector bases( TMP + "/frag_reads_orig.fastb" );
          vecqualvector quals( TMP + "/frag_reads_orig.qualb" );
          if ( heur.DETECT_VARIANTS && OUT_INT_HEAD != "" ) 
          {    vec<int> variant_edge_list;
               shb.DetectVarients( bases, quals, logc, &snp_bubbles, 
                   &variant_edge_list );
               VariantsMarkedDot( OUT_INT_HEAD + ".final.marked", shb, 
                   variant_edge_list, logc );    }
          TraceEdges( shb, logc.TRACE_EDGES, bases, quals );    }
     if (logc.ANALYZE) AnalyzeAssembly( shb, ref.G, ref.LG, ref.Glocs );

     // Form efasta assembly.  Hide bubbles.  Identify reverse complement edges
     // to be ignored.  Identify hidden edges.  Mark variants back.

     HyperEfasta he(shb);
     vec<int> inv = shb.Inv( );
     vec<Bool> hide;
     MakeEfastaAssembly( he, inv, hide, snp_bubbles, heur, logc );

     // Output stuff, trace reference, validate.
     if ( log_control.OUT_INT_HEAD != "" )           
          shb.DumpFiles( log_control.OUT_INT_HEAD + ".final", log_control, logc );

     AssessAssembly( SAMPLE, shb, he, hide, TMP, ref,
          spec.HUMAN_CONTROLS, logc, NUM_THREADS, RTCtrl );
     if ( EXIT == "REFTRACE" ) Done(clock);

     if (logc.EVAL_ASSEMBLY) 
     {    const int iSWAk=0; // set to 501 to use hash-based aligner to avoid n^2 blow up
          AssemblyError e = EvalAssembly( ref, shb, shb.Inv(), cout, iSWAk, 0 );
          vec<int> gaps = e.GetGaps( ), endgaps = e.GetLeftGaps( );
          endgaps.append( e.GetRightGaps( ) );
          vec<int> indels = e.GetIndels( );
          int subs = e.GetSubs( );
          ReverseSort(gaps), ReverseSort(endgaps), ReverseSort(indels);
          cout << "SUMMARY GAPS: " << printSeq(gaps) << endl;
          cout << "SUMMARY ENDGAPS: " << printSeq(endgaps) << endl;
          cout << "SUMMARY INDELS: " << printSeq(indels) << endl;
          cout << "SUMMARY SUBS: " << subs << endl;    }

     if ( OUT_HEAD != "" ){
         if( logc.MAKE_FASTG) fastg::WriteFastg(OUT_HEAD+".fastg",he);
         DumpEfastaAssembly( shb, he, inv, hide, OUT_HEAD, logc );
     }

     // Done.

     PrintTail( command, X, X_actual, clock, logc );
     if ( logc.VALIDATE >= 0 ) ValidateRef( shb, TMP, log_control, logc );    }
