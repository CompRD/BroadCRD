///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// =================================================================================
//
// INTRODUCTION.
//
// SomaticCall is a program to detect somatic mutations using sequence data
// from tumor and matched normal samples.  It only looks for substitutions, and it
// is designed to be highly stringent, so as to achieve a low false positive rate.
//
// It takes as input aligned reads in BAM format.  The following are required:
// (1) quality scores have been recalibrated;
// (2) duplicates have been marked.
// These requirements are not tested.  Details of implementation of (1) and (2), as
// well as details of the alignment algorithm may significantly affect the results
// of this program.
//
// The final output of the program is a file OUT/mutation_reports.  The first
// column of this file provides zero-based coordinates for putative mutations.
// The other columns provide statistics about the event.
//
// The algorithm is outlined below and implemented in this file, and the files
// SomaticCallSlave.cc, SomaticCallTools.cc, FindMoreReadsSlave.cc.
//
// =================================================================================
//
// SYSTEM REQUIREMENTS.
//
// (1) gcc 4.3.3
// (2) GMP with C++ support
// (3) samtools must be in your path
// (4) various executables from this codebase including SomaticCallSlave need
// to be in your path.  See required_execs, below, or just "make SomaticCall" and
// put the bin directory from that make in your path.
//
// =================================================================================
//
// DESCRIPTION OF THE ALGORITHM.
//
// (Note that this isn't quite right anymore: the last file created is
// mutation_reports4.)
//
// (A) Handling of alignments.  They are discarded if they are marked as duplicates 
// or assigned mapping quality score 0.  The mapping quality score is not otherwise 
// used.
//
// (B) Core statistical computation.  We describe a test that is repeated several
// times in the process.  Given aligned reads for the tumor and normal, the code 
// SomaticMutation in SomaticCallTools.cc tests for a likely somatic mutation 
// according to the following criteria: 
// (i) either an adjusted quality score sum in the tumor for the mutant base must be 
//     at least 100 or the LOD score for mutant:ref + mutant:mutant vs ref:ref must 
//     be at least 6.3;
// (ii) the quality score sum for the mutant base in the normal must be < 50 and the
//      LOD score for ref:ref vs mutant:ref + mutant:mutant must be at least 2.3.
// - This test is applied several times as the alignments are modified in the course
// of the program.  Once a putative mutation fails the test, it is not considered
// again.
// - There is also a pretest which checks that the adjusted quality score sum in the
// tumor is at least 60, but this is intended only as an optimization.
// - Note that the adjusted quality score sum test of (i) does not actually affect
// the final output, but is included for experimental purposes.
//
// (C) Overall steps in the process.
//
// (i) First, dbSNP position are excluded, and a list of putative somatic mutations 
// is obtained by applying SomaticMutation (yielding mutation_reports0).
//
// (ii) Then we require that a graph assembly of the locus has no cycles and that 
// all paths across it have the same length (yielding mutation_reports1).
//
// (iii) All reads supporting the mutant allele in the tumor are then realigned at 
// high stringency.  We seed on 12-mers of frequency up to 1000, allow up to 4 
// errors (including indels of size up to 2), and require that the next best 
// placement has // more than 6 more errors.  Then we call SomaticMutation again
// (yielding mutation_reports2).
//
// (iv) Next we eliminate mutations for which the alternative base is supported by 
// only one read start position (yielding mutation_reports3).
//
// (v) Then we search the entire set of reads from the normal sample with the goal 
// of finding more reads at each putative mutation.  
//    (a) To do this we first examine the given tumor alignments to find 21-mers 
//        from the tumor data that are centered at the mutation.  Find the most 
//        frequent 3 such 21-mers.  If there are no such 21-mers, reject the locus.
//        (At this point we create mutation_reports4.)
//    (b) For each of the two 20-mers in the above 21-mers, find every instance in 
//        the normal reads.
//    (c) For each such instance, align the read to the reference at the locus of 
//        the somatic mutation, allowing up to 6 mismatches and indels.
//    (d) Align the read to the whole genome using 12-mers, not allowing indels.  If
//        an alignment is found that is nearly as good or better than the local 
//        placement, discard the read.
//    (e) Return the local placements of the surviving reads.
//    (f) Add in the additional normal reads and recall mutations.  
// This process yields mutation_reports5.
//
// (vi) Finally we filter for tumor score >= 6.3, normal_score 2.3, yielding 
// mutation_reports5f.  We also create a link mutation_reports that points to this.
// Each entry in these files has a reference to report_dir, a directory that contains
// files t.visual and n.visual, that can be viewed with less -r to show the reads
// that support a given mutation.
//
// =================================================================================
//
// NOTES ON PARALLELIZATION.
//
// - In one test on a small data set (human hybrid selection), the optimum value for
//   MAX_THREADS was found to be about 10.  Adding more processors increased wall
//   clock time for the run.  As compared to MAX_THREADS=1, MAX_THREADS=10 was
//   about 4 times faster.
//
// - For human WGS, for 16-processor machine, recommend MAX_THREADS=16.
//
// - The optimal value for MAX_THREADS may depend on the type of filesystem used 
//   and system load at a particular time.
//
// - MAX_THREADS is applied at three different points in the code, and these three
//   points may have different optima.
//
// =================================================================================
//
// HOW TO GENERATE REFERENCE FILE FORMATS.
//
// - To generate the REFHEAD.fastb and REFHEAD.lookup files, run
//   MakeLookupTable SOURCE=g OUT_HEAD=REFHEAD
// where
//   * g is the fasta file for the genome reference
//   * REFHEAD is the output "head"; will create REFHEAD.fastb and REFHEAD.lookup.
//
// =================================================================================
//
// KNOWN PROBLEMS.
//
// 1. Naming and numbering of chromosomes in report files is tailored to human
// and actually depends on the order in the reference file.
//
// =================================================================================

#include <map>
#include <strstream>

#include "Basevector.h"
#include "Bitvector.h"
#include "FastIfstream.h"
#include "FeudalMimic.h"
#include "MainTools.h"
#include "ParseSet.h"
#include "Qualvector.h"
#include "Vec.h"
#include "VecUtilities.h"
#include "cancer/SomaticCallTools.h"
#include "lookup/AlignCollector.h"
#include "lookup/LookAlign.h"
#include "lookup/ShowAlignmentPileTools.h"
#include "math/Functions.h"
#include "system/Parallel.h"
#include "feudal/BinaryStream.h"

// MakeDepend: dependency MakeLookupTable
// MakeDepend: dependency ShortQueryLookup
// MakeDepend: dependency SomaticCallSlave

void BuildQuals( const vec< pair<int,int> >& mutations, const vecbasevector& tbases,
     const vecqualvector& tquals, const vec<look_align>& taligns,
     const vecbasevector& nbases, const vecqualvector& nquals,
     const vec<look_align>& naligns, 
     vec< vec<unsigned char> >& TQ, vec< vec<unsigned char> >& NQ,
     vec< vec<unsigned char> >& TQ_rc, vec< vec<unsigned char> >& NQ_rc )
{    int nmut = mutations.size( );
     TQ.clear_and_resize(4*nmut), NQ.clear_and_resize(4*nmut);
     TQ_rc.clear_and_resize(4*nmut), NQ_rc.clear_and_resize(4*nmut);
     cout << Date( ) << ": start building Q score vectors" << endl;
     for ( int pass = 1; pass <= 2; pass++ )
     {    vec< vec<unsigned char> >& Q = ( pass == 1 ? TQ : NQ );
          vec< vec<unsigned char> >& Q_rc = ( pass == 1 ? TQ_rc : NQ_rc );
          const vec<look_align>& aligns = ( pass == 1 ? taligns : naligns );
          const vecbasevector& bases = ( pass == 1 ? tbases : nbases );
          const vecqualvector& quals = ( pass == 1 ? tquals : nquals );
          for ( int i = 0; i < aligns.isize( ); i++ )
          {    const look_align& la = aligns[i];
               int p1 = la.pos1( ), p2 = la.pos2( );
               basevector rd1 = bases[ la.query_id ];
               qualvector q1 = quals[ la.query_id ];
               if ( la.Rc1( ) )
               {    rd1.ReverseComplement( );
                    q1.ReverseMe( );    }
               for ( int j = 0; j < la.a.Nblocks( ); j++ )
               {    if ( la.a.Gaps(j) > 0 ) p2 += la.a.Gaps(j);
                    if ( la.a.Gaps(j) < 0 ) p1 -= la.a.Gaps(j);
                    for ( int x = 0; x < la.a.Lengths(j); x++ ) 
                    {    int mpos = BinPosition( 
                              mutations, make_pair( la.target_id, p2 ) );
                         if ( mpos >= 0 ) 
                         {    Q[ 4*mpos + rd1[p1] ].push_back( q1[p1] );
                              Q_rc[ 4*mpos + rd1[p1] ].push_back( la.Rc1( ) );    }
                         ++p1; ++p2;    }    }    }    }
     return;    }

int main(int argc, char **argv) 
{
     RunTime( );

     BeginCommandArguments;
     CommandDoc( "Detect somatic mutations using sequence data from tumor and "
          "matched normal samples.  See SomaticCall.cc for detailed "
          "documentation." );
     CommandArgument_String_Doc(OUT, 
          "output directory, for exclusive use of this process");
     CommandArgument_String_Doc(TUMOR_BAM, ".bam file containing the de-dupped "
          "aligned reads for the tumor; index .bam.bai must exist too");
     CommandArgument_String_Doc(NORMAL_BAM, ".bam file containing the de-dupped "
          "aligned reads for the normal; index .bam.bai must exist too");
     CommandArgument_String_Doc(REFHEAD, "reference sequence head; "
          "there should exist files REFHEAD.fasta, REFHEAD.fastb, and "
          "REFHEAD.lookup; see documentation in SomaticCall.cc for how to "
          "generate the latter two files");
     CommandArgument_String_Doc(DBSNP, "name of file containing one line for each "
          "known SNP, where each line has two fields: the chromosome name, and the "
          "one-based position of the SNP on the chromosome");
     CommandArgument_String_OrDefault_Doc(CHRS,
          "chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y}",
          "list of chromosomes to search for somatic mutations; the default list "
          "is appropriate for a human male; note exclusion of mitochondrial genome "
          "and 'random' chromosomes");
     CommandArgument_Int_Doc(MAX_THREADS, "maximum number of "
          "parallel threads; see notes in SomaticCall.cc");
     CommandArgument_Int_OrDefault_Doc(MAX_ERRS_TUMOR, 4,
          "exclude tumor alignments having more than this number of substitutions");
     CommandArgument_Int_OrDefault_Doc(MAX_ERRS_NORMAL, 10,
          "exclude normal alignments having more than this number of substitutions");
     CommandArgument_Int_OrDefault_Doc(MAX_QUAL_DIFF, 150,
          "exclude alignments for which the sum of read quality scores at "
          "substitutions exceeds this value");
     CommandArgument_Int_OrDefault_Doc(LEN, 12000000, "chunk size for first part "
          "of code - only of interest for optimizing run time; optimized for "
          "hybrid selection; recommend 3M for whole genomes");
     CommandArgument_Int_OrDefault_Doc(PHASE, 1, "phase to start processing, "
          "either 1 or 2 or 3 or 4 - only of interest for experimentation or "
          "debugging; set PHASE=3 to restart at the point where the "
          "FindMoreReadsSlaves are run");
     CommandArgument_Int_OrDefault(MIN_MUTANT_SUM_PRETEST, 60);
     CommandArgument_Int_OrDefault(MIN_MUTANT_SUM, 100);
     CommandArgument_Double_OrDefault(TUMOR_THRESHOLD, 10.0);
     CommandArgument_Double_OrDefault(NORMAL_THRESHOLD, 4.0);
     CommandArgument_Bool_OrDefault_Doc(DUMP_MAPPING_SCORES, False,
          "show mapping scores of all tumor alignments that are kept, and those "
          "that are discarded by realignment step");
     CommandArgument_String_OrDefault_Doc(FAKE_MUTATIONS, "",
          "ignore real mutations but introduce this file of fake mutations; "
          "format is one mutation per line "
          "with four fields: (1) zero-based chromosome index; (2) zero-based "
          "position; (3) mutated base; (4) frequency");
     CommandArgument_Bool_OrDefault_Doc(REALIGN, False,
          "realign reads to make sure they're in the right place");
     EndCommandArguments;

     // Start!

     cout << Date( ) << ": SomaticCall starting" << endl;
     double clock = WallClockTime( );
     Mkdir777(OUT); 
     {    Ofstream( out, OUT + "/command" );
          command.PrintTheCommandPretty( out );    }

     // Sanity checks.

     vec<String> required_execs, required_files;
     required_execs.push_back( "samtools", "SomaticCallSlave", "ShortQueryLookup" );
     for ( int i = 0; i < required_execs.isize( ); i++ )
     {    if ( LineOfOutput( "which " + required_execs[i] ).Contains( "not found" ) )
          {    cout << "You need to have " + required_execs[i] 
                    + " in your path." << endl;
               exit(1);    }    }
     required_files.push_back( TUMOR_BAM, TUMOR_BAM + ".bai", NORMAL_BAM, DBSNP,
          NORMAL_BAM + ".bai", REFHEAD + ".fastb", REFHEAD + ".lookup" );
     for ( int i = 0; i < required_files.isize( ); i++ )
     {    if ( !IsRegularFile( required_files[i] ) )
          {    cout << "Can't find file " << required_files[i] << "." << endl;
               exit(1);    }    }
     if ( IsOlder( TUMOR_BAM + ".bai", TUMOR_BAM ) )
     {    cout << "Index for TUMOR_BAM predates it." << endl;
          exit(1);    }
     if ( IsOlder( NORMAL_BAM + ".bai", NORMAL_BAM ) )
     {    cout << "Index for NORMAL_BAM predates it." << endl;
          exit(1);    }
     if ( MAX_THREADS <= 0 ) 
     {    cout << "Illegal value for MAX_THREADS." << endl;
          exit(1);    }

     // Load reference.

     vecbasevector ref( REFHEAD + ".fastb" );
     vec<String> refdict, chrs;
     String line, chr, REFLIST = OUT + "/ref_list";
     {    Ofstream( rout, REFLIST );
          fast_ifstream rin( REFHEAD + ".fasta" );
          int count = 0;
          while(1)
          {    Bool match = getline_if_match( rin, line, ">" );
               if ( rin.fail( ) ) break;
               if ( !match ) continue;
               line = line.After( ">" );
               if ( line.Contains( " " ) ) line = line.Before( " " );
               rout << line << " " << ref[count++].size( ) << "\n";
               refdict.push_back(line);    }    }
     ParseStringSet( CHRS, chrs );
     vec<int> chrsi;
     for ( int i = 0; i < chrs.isize( ); i++ )
     {    if ( !Member( refdict, chrs[i] ) )
          {    cout << "Your argument for CHRS includes " << chrs[i]
                    << ", which is not present in REFHEAD.fasta." << endl;
               cout << "This will need to be fixed either by changing CHRS "
                    << "or REFHEAD.fasta." << endl;
               cout << "Abort." << endl;
               CRD::exit(1);    }
          chrsi.push_back( Position( refdict, chrs[i] ) );    }
     UniqueSort(chrsi);

     // Load DBSNP.

     vecbitvector dbsnp;
     Mimic( ref, dbsnp );
     fast_ifstream din(DBSNP);
     String dchr;
     int snps = 0;
     while(1)
     {    getline( din, line );
          if ( din.fail( ) ) break;
          line.GlobalReplaceBy( "\t", " " );
          if ( !line.Contains( " " ) || !line.After( " " ).IsInt( ) )
          {    cout << "Your DBSNP file is in the wrong format." << endl;
               cout << "It contains the following line: \"" << line << "\"" << endl;
               cout << "Each line should have two fields, separated by a single "
                    << "blank." << endl;
               cout << "Please see the manual." << endl;
               CRD::exit(1);    }
          istrstream iline( line.c_str( ) );
          int pos;
          iline >> dchr >> pos;
          int chr = Position( refdict, dchr );
          if ( chr < 0 ) continue;
          snps++;
          dbsnp[chr].Set( pos-1, True );    }
     cout << Date( ) << ": " << snps << " SNPs marked" << endl;
     DBSNP = OUT + "/dbsnp.vecbitvector";
     dbsnp.WriteAll(DBSNP);

     // Start phase 1.

     if ( PHASE == 1 )
     {    Mkdir777( OUT + "/SomaticCall" );
          vec<String> commands;
          for ( int i = 0; i < chrsi.isize( ); i++ )
          {    int C = chrsi[i];
               for ( int start = 0; start < ref[C].isize( ); start += LEN )
               {    int len = Min( LEN, ref[C].isize( ) - start );
                    String ns = ToString( commands.size( ) );
                    Mkdir777( OUT + "/SomaticCall/" + ns );
                    commands.push_back( "SomaticCallSlave" 
                         + ARG(REF, REFHEAD + ".fastb") + ARGC(REFLIST)
                         + ARGC(TUMOR_BAM) + ARGC(NORMAL_BAM)
                         + ARGC(DBSNP) + ARGC(C) + ARG(START, start) 
                         + ARG(LEN, len) + ARGC(MAX_ERRS_TUMOR) 
                         + ARGC(MAX_ERRS_NORMAL) + ARGC(MAX_QUAL_DIFF)
                         + ARGC(MIN_MUTANT_SUM_PRETEST) + ARGC(MIN_MUTANT_SUM)
                         + ARGC(TUMOR_THRESHOLD) + ARGC(NORMAL_THRESHOLD)
                         + ARGC(FAKE_MUTATIONS)
                         + ARG(OUT_HEAD, OUT + "/SomaticCall/" + ns + "/")
                         + " >& " + OUT + "/SomaticCall/" + ns + "/out" );    }    }
          Echo( ToString( commands.size( ) ), OUT + "/SomaticCall.ncommands" );
          cout << Date( ) << ": running " << commands.size( )
               << " instances of SomaticCallSlave" << endl;
          RunCommandsInParallel( commands, Max( 1, MAX_THREADS-1 ) );    }

     // Load tumor and normal reads that cover a putative mutation.

     cout << Date( ) << ": start phase 2" << endl;
     int ncommands = FirstLineOfFile( OUT + "/SomaticCall.ncommands" ).Int( );
     vecbasevector tbases, tbases_hot, nbases;
     vecqualvector tquals, nquals;
     vec<look_align> taligns, naligns;
     vec<int> tmappingscore, nmappingscore;
     vec<Bool> thot;
     vec<unsigned char> altbase;
     vec< pair<int,int> > mutations;
     vec<String> mutation_reports0, mutation_reports1;
     for ( int i = 0; i < ncommands; i++ )
     {    String head = OUT + "/SomaticCall/" + ToString(i) + "/";
          if ( !IsRegularFile( head + "mutations" ) ) continue;
          thot.appendFromBinaryFile(head + "thot");
          altbase.appendFromBinaryFile(head + "altbase");
          vec<look_align> talignsi, nalignsi;
          LoadLookAligns( head + "tqltout", talignsi, True );
          LoadLookAligns( head + "nqltout", nalignsi, True );
          for ( int j = 0; j < talignsi.isize( ); j++ )
               talignsi[j].query_id += tbases.size( );
          for ( int j = 0; j < nalignsi.isize( ); j++ )
               nalignsi[j].query_id += nbases.size( );
          taligns.append(talignsi), naligns.append(nalignsi);
          nbases.ReadAll( head + "nbases", True );
          tbases.ReadAll( head + "tbases", True );
          tquals.ReadAll( head + "tquals", True );
          nquals.ReadAll( head + "nquals", True );
          ReadStrings( head + "mutation_reports0", mutation_reports0 );
          ReadStrings( head + "mutation_reports1", mutation_reports1 );
          mutations.appendFromBinaryFile(head + "mutations");
          tmappingscore.appendFromBinaryFile(head + "tmappingscore");
          nmappingscore.appendFromBinaryFile(head + "nmappingscore");    }
     if ( PHASE < 3 )
     {    BinaryWriter::writeFile( OUT + "/mutations", mutations );
          WriteStrings( OUT + "/mutation_reports0", mutation_reports0 );
          WriteStrings( OUT + "/mutation_reports1", mutation_reports1 );    }

     // Realign the reads that support a mutation.

     vec<int> thot_ids;
     for ( size_t i = 0; i < tbases.size( ); i++ )
     {    if ( thot[i] ) 
          {    tbases_hot.push_back_reserve( tbases[i] );
               thot_ids.push_back(i);    }    }
     const size_t reads_per_batch = 500;
     size_t N = tbases_hot.size( );
     if ( PHASE < 3 && REALIGN )
     {    vec<String> commands2;
          Mkdir777( OUT + "/SomaticCall2" );
          for ( size_t i = 0; i < tbases_hot.size( ); i += reads_per_batch )
          {    Mkdir777( OUT + "/SomaticCall2/" + ToString(i/reads_per_batch) );
               String head2 = OUT + "/SomaticCall2/" + ToString(i/reads_per_batch) 
                    + "/reads";
               vecbasevector batch;
               batch.Append( tbases_hot, i, Min( i + reads_per_batch, N ) );
               batch.WriteAll( head2 + ".fastb" );
               commands2.push_back( "ShortQueryLookup" + ARG(SEQS, head2 + ".fastb")
                    + ARG(L, REFHEAD + ".lookup")
                    + ARG(OUT_SUFFIX, ".qltout") + ARG(MAX_FREQ, 1000)
                    + ARG(MAX_ERRS, 10) + ARG(ERR_DIFF, 6) + ARG(UNIQUE, True) 
                    + " >& " + head2 + ".log" );    }
          cout << Date( ) << ": running " << commands2.size( ) 
               << " instances of ShortQueryLookup" << endl;
          RunCommandsInParallel( commands2, Max( 1, MAX_THREADS-1 ) );    }

     // Boring bookkeeping to load the unique alignments, then discard the 
     // non-unique ones.  This code touches tbases, tquals, taligns, and nothing 
     // else.

     if (REALIGN)
     {    vec<look_align> new_taligns;
          size_t n = tbases.size( );
          vec<Bool> tdiscard( n, False ), treplace( n, False );
          for ( size_t i = 0; i < tbases_hot.size( ); i += reads_per_batch )
          {    String head2 = OUT + "/SomaticCall2/" + ToString(i/reads_per_batch)
                    + "/reads";
               Remove( head2 + ".qltout" ), Remove( head2 + ".minAlignErrors.txt" );
               vec<look_align> aligns;
               LoadLookAligns( head2 + ".qltout.unique", aligns );
               vec<Bool> aligned( Min( reads_per_batch, N-i ), False );
               for ( int j = 0; j < aligns.isize( ); j++ )
               {    aligned[ aligns[j].query_id ] = True;
                    aligns[j].query_id = thot_ids[ i + aligns[j].query_id ];
                    treplace[ aligns[j].query_id ] = True;    }
               new_taligns.append(aligns);
               for ( int j = 0; j < aligned.isize( ); j++ )
                    if ( !aligned[j] ) tdiscard[ thot_ids[i+j] ] = True;    }
          vec<Bool> delete_taligns( taligns.size( ), False );
          if (DUMP_MAPPING_SCORES) cout << "\n";
          for ( int i = 0; i < taligns.isize( ); i++ )
          {    int id = taligns[i].query_id;               
               if (DUMP_MAPPING_SCORES)
               {    cout << "tumor alignment " << i << " mapping score " 
                         << tmappingscore[i];
                    if ( tdiscard[id] ) cout << " discarding\n";
                    else cout << " keeping\n";    }
               if ( treplace[id] || tdiscard[id] ) delete_taligns[i] = True;    }
          if (DUMP_MAPPING_SCORES) cout << "\n";
          EraseIf(taligns, delete_taligns);
          taligns.append(new_taligns);
          vec<int> to_new_ids( tbases.size( ), -1 );
          int count = 0;
          for ( size_t i = 0; i < n; i++ )
               if ( !tdiscard[i] ) to_new_ids[i] = count++;
          for ( size_t i = 0; i < taligns.size( ); i++ )
               taligns[i].query_id = to_new_ids[ taligns[i].query_id ];
          tbases.EraseIf(tdiscard), tquals.EraseIf(tdiscard);    }

     // Recall mutations.

     vec< vec<unsigned char> > TQ, NQ, TQ_rc, NQ_rc;
     BuildQuals( mutations, tbases, tquals, taligns, nbases, nquals,
          naligns, TQ, NQ, TQ_rc, NQ_rc );
     int nmut = mutations.size( );
     vec<Bool> rejected( nmut, False );
     vec<String> mutation_reports2 = mutation_reports1;
     for ( int i = 0; i < nmut; i++ )
     {    int chr = mutations[i].first, pos = mutations[i].second;
          unsigned char refbase = ref[chr][pos];
          String info;
          if ( !SomaticMutation( chr, pos, ref, &TQ[4*i], &NQ[4*i], &TQ_rc[4*i], 
               &NQ_rc[4*i], refbase, MIN_MUTANT_SUM_PRETEST, MIN_MUTANT_SUM, 
               TUMOR_THRESHOLD, NORMAL_THRESHOLD, altbase[i], info ) )
          {    rejected[i] = True;    }    
          else mutation_reports2[i] = info;    }
     EraseIf( mutation_reports2, rejected );
     if ( PHASE < 3 ) WriteStrings( OUT + "/mutation_reports2", mutation_reports2 );
          
     // Build index to tumor alignments.

     vec< vec<int> > tindex(nmut);
     for ( int i = 0; i < taligns.isize( ); i++ )
     {    const look_align& la = taligns[i];
          int p2 = la.pos2( ), P2 = la.Pos2( ), t = la.target_id;
          int low = LowerBound( mutations, make_pair( t, p2 ) );
          int high = LowerBound( mutations, make_pair( t, P2 ) );
          for ( int j = low; j < high; j++ )
               tindex[j].push_back(i);    }

     // Find mutations supported by only one read start position.  This should be 
     // handled differently, by eliminating "unpaired duplicates" in the upstream 
     // process.

     for ( int i = 0; i < nmut; i++ )
     {    vec< vec<int> > starts(4);
          int tig = mutations[i].first, pos = mutations[i].second;
          for ( int j = 0; j < tindex[i].isize( ); j++ )
          {    const look_align& la = taligns[ tindex[i][j] ];
               basevector rd = tbases[ la.query_id ];
               if ( la.Rc1( ) ) rd.ReverseComplement( );
               int p1 = la.a.pos1( ), p2 = la.a.pos2( );
               if ( p2 == pos ) starts[ rd[p1] ].push_back( la.Fw1( ) ? p1 : -p1-1 );
               else
               {    for ( int j = 0; j < la.a.Nblocks( ); j++ ) 
                    {    if ( la.a.Gaps(j) > 0 )  
                         {    p2 += la.a.Gaps(j);    
                              if ( p2 == pos ) 
                              {    starts[ rd[p1] ].push_back( 
                                        la.Fw1( ) ? p1 : -p1-1 );
                                   break;    }    }
                         if ( la.a.Gaps(j) < 0 ) p1 -= la.a.Gaps(j);    
                         for ( int x = 0; x < la.a.Lengths(j); x++ ) 
                         {     ++p1; ++p2;    
                              if ( x < la.a.Lengths(j) - 1 && p2 == pos ) 
                              {    starts[ rd[p1] ].push_back( 
                                        la.Fw1( ) ? p1 : -p1-1 );
                                   break;    }    }
                         if ( p2 == pos ) break;     }    }    }
          Bool supported = False;
          for ( int j = 0; j < 4; j++ )
          {    if ( j != ref[tig][pos] )
               {    UniqueSort( starts[j] );
                    if ( starts[j].size( ) > 1 ) supported = True;    }    }
          if ( !supported ) rejected[i] = True;    }
     vec<String> mutation_reports3 = mutation_reports1;
     EraseIf( mutation_reports3, rejected );
     if ( PHASE < 3 ) WriteStrings( OUT + "/mutation_reports3", mutation_reports3 );

     // Perform main compute to find more reads.

     vec<vecbasevector> reads_by_mutation(nmut);
     vec<vecqualvector> quals_by_mutation(nmut);
     vec< vec<look_align> > aligns_by_mutation(nmut);
     String OUTHEAD = OUT + "/FindMoreReads";
      
     // Find the 21-mers that characterize each mutation.

     cout << Date( ) << ": logging local contexts to OUT/local_context" << endl;
     Ofstream( context, OUT + "/local_context" );
     vec< vec<basevector> > m21(nmut);
     vec<int> m21_mutation_index;
     for ( int i = 0; i < nmut; i++ )
     {    if ( rejected[i] ) continue;
          int tig = mutations[i].first, pos = mutations[i].second;
          basevector R;
          R.SetToSubOf( ref[tig], pos - 10, 21 );
          m21[i].push_back(R);
          m21_mutation_index.push_back(i);
          vec<basevector> rs;
          for ( int j = 0; j < tindex[i].isize( ); j++ )
          {    const look_align& la = taligns[ tindex[i][j] ];
               basevector rd = tbases[ la.query_id ];
               if ( la.Rc1( ) ) rd.ReverseComplement( );
               int p1 = la.a.pos1( ), p2 = la.a.pos2( );
               int readpos = -1;
               if ( p2 == pos && rd[p1] != ref[tig][p2] ) readpos = p1;
               for ( int j = 0; j < la.a.Nblocks( ); j++ ) 
               {    if ( la.a.Gaps(j) > 0 )  
                    {    p2 += la.a.Gaps(j);    
                         if ( p2 == pos && rd[p1] != ref[tig][p2] ) 
                              readpos = p1;    }
                    if ( la.a.Gaps(j) < 0 ) p1 -= la.a.Gaps(j);    
                    for ( int x = 0; x < la.a.Lengths(j); x++ ) 
                    {     ++p1; ++p2;    
                         if ( p2 == pos && rd[p1] != ref[tig][p2] ) 
                              readpos = p1;    }    }
               if ( readpos >= 0 && readpos >= 10 && readpos+10 < rd.isize( ) )
               {    basevector r;
                    r.SetToSubOf( rd, readpos - 10, 21 );
                    if ( r == R ) continue;
                    rs.push_back(r);    }     }
          context << "\nmutation " << i << " = " << tig << "." << pos << endl;
          context << "                *" << endl;
          context << "ref = " << R.ToString( ) << endl;
          Sort(rs);
          vec<basevector> rs_red;
          vec<int> rs_count;
          for ( int j = 0; j < rs.isize( ); j++ )
          {    int k = rs.NextDiff(j);
               rs_red.push_back( rs[j] );
               rs_count.push_back(k-j);    
               j = k - 1;    }
          ReverseSortSync( rs_count, rs_red );
          for ( int j = 0; j < Min( 3, rs_red.isize( ) ); j++ )
          {    context << "mut = " << rs_red[j].ToString( ) 
                    << " [" << rs_count[j] << "]" << endl;    
               m21[i].push_back( rs_red[j] );    
               m21_mutation_index.push_back(i);    }    
          if ( rs.empty( ) )
          {    context << "NO ALTERNATES FOUND!\n";
               rejected[i] = True;    }    }
     vec<String> mutation_reports4 = mutation_reports1;
     EraseIf( mutation_reports4, rejected );
     if ( PHASE < 3 ) WriteStrings( OUT + "/mutation_reports4", mutation_reports4 );

     // Define a directory for each mutation.  This is complicated by not wanting
     // to have directories with too many entries.

     Mkdir777( OUT + "/MutationsDir" );
     int nmutf = mutation_reports4.size( );
     vec<String> dirs(nmutf);
     if ( nmutf <= 100 )
     {    for ( int i = 0; i < nmutf; i++ )
          {    dirs[i] = OUT + "/MutationsDir/" + ToString(i);
               Mkdir777( dirs[i] );    }    }
     else
     {    int nsubs = int( round( sqrt(nmutf) ) );
          for ( int x = 0; x < nmutf; x++ )
          {    int i = x / nsubs, j = x % nsubs;
               dirs[x] = OUT + "/MutationsDir/" + ToString(i) + "/" + ToString(j);
               Mkdir777( OUT + "/MutationsDir/" + ToString(i) );
               Mkdir777( dirs[x] );    }    }

     // Create files of supporting data for each mutation.

     for ( int i = 0; i < nmutf; i++ )
     {    String loc = mutation_reports4[i].After( 
               "reference_fasta_record_" ).Before( " " );
          int chr = loc.Before( ":" ).Int( ), pos = loc.After( ":" ).Int( );
          for ( int pass = 1; pass <= 2; pass++ )
          {    vecbasevector basesx;
               vecqualvector qualsx;
               vec<look_align> alignsx;
               const vecbasevector& bases = ( pass == 1 ? tbases : nbases );
               const vecqualvector& quals = ( pass == 1 ? tquals : nquals );
               const vec<look_align>& aligns = ( pass == 1 ? taligns : naligns );
               int count = 0;
               for ( int j = 0; j < aligns.isize( ); j++ )
               {    look_align la = aligns[j];
                    if ( la.target_id != chr ) continue;
                    if ( la.Pos2( ) <= pos || la.pos2( ) > pos ) continue;
                    basesx.push_back_reserve( bases[la.query_id] );
                    qualsx.push_back_reserve( quals[la.query_id] );
                    la.query_id = count++;
                    alignsx.push_back(la);    }
               String sample = ( pass == 1 ? "t" : "n" );
               basesx.WriteAll( dirs[i] + "/" + sample + ".fastb" );
               qualsx.WriteAll( dirs[i] + "/" + sample + ".qualb" );
               Ofstream( aout, dirs[i] + "/" + sample + ".qltout" );
               for ( int j = 0; j < alignsx.isize( ); j++ )
                    alignsx[j].PrintParseable(aout);    }    }

     // Create visual piles for each mutation.  This is very inefficient as it
     // reads the reference twice for each mutation.

     for ( int i = 0; i < nmutf; i++ )
     {    String loc = mutation_reports4[i].After( 
               "reference_fasta_record_" ).Before( " " );
          int chr = loc.Before( ":" ).Int( ), pos = loc.After( ":" ).Int( );
          for ( int pass = 1; pass <= 2; pass++ )
          {    String sample = ( pass == 1 ? "t" : "n" );
               int read_length = 0;
               vecbasevector reads( dirs[i] + "/" + sample + ".fastb" );
               for ( size_t j = 0; j < reads.size( ); j++ )
                    read_length = Max( read_length, reads[j].isize( ) );

               // Code copied from lookup/ShowAlignmentPile.cc.

               String HEAD = dirs[i] + "/" + sample;
               String COORDINATES = ToString(chr) + ":" + ToString(pos);
               Ofstream( fout, OUT + "/pile" );
	       vecbitvector basesUsed;
	       vec< pair<bool,unsigned int> > specs; 
	       MaxErrDiffAlignCollector alignments 
                    = LoadAlignments( HEAD + ".qltout", reads.size(), 1000000000);
	       vec<GenomeCoordinateRange> coordinates 
                    = LoadCoordinates(COORDINATES, read_length, specs);
	       map< GenomeCoordinateRange, vec<look_align> > coverage;
	       for (unsigned int a = 0; a < alignments.size(); a++) 
               {    for (MaxErrDiffAlignCollector::Iterator qltit 
                         = alignments.Begin(a); qltit != alignments.End(a); qltit++) 
                    {    for (unsigned int c = 0; c < coordinates.size(); c++) 
                         {    if (qltit->target_id != 
                                   static_cast<int>(coordinates[c]._contig) ) 
                              {    continue;    }
		              if ( qltit->Pos2() <= static_cast<int>
                                   (coordinates[c]._start) || qltit->pos2() >= 
                                   static_cast<int>(coordinates[c]._end) ) 
                              {    continue;    }
		              if ( specs[c].first && ( static_cast<unsigned int>
                                   (qltit->Pos2()) <= specs[c].second || 
                                   static_cast<unsigned int> (qltit->pos2()) 
                                   > specs[c].second ) ) 
                              {    continue;    }
		              coverage[coordinates[c]].push_back(
                                   *qltit);    }    }    }
	       for (map< GenomeCoordinateRange, vec<look_align> >::iterator cit 
                    = coverage.begin(); cit != coverage.end(); cit++) 
               {    unsigned int length = cit->first._end - cit->first._start;
                    if (cit->first._start + length > ref[cit->first._contig].size()) 
                    {    length = ref[cit->first._contig].size() 
                              - cit->first._start;    }
		    basevector refsection(ref[cit->first._contig], 
                         cit->first._start, length);
		    PileOfAlignments pile(cit->first, refsection);
		    vec<look_align> vla = cit->second;
		    for (unsigned int i = 0; i < vla.size(); i++) 
                         pile.AddAlignment(vla[i], reads[vla[i].query_id]);
		    pile.Print(fout, 0, basesUsed, "dot");    }
	       fout.close();

               // Parse the "ShowAlignmentPile" output just generated.

               Ofstream( out, dirs[i] + "/" + sample + ".visual" );
               String line0;
               fast_ifstream inx( OUT + "/pile" );
               Bool first = True;
               int pluses = 0, saved_pluses = 0;
               while(1)
               {    getline( inx, line );
                    if ( inx.fail( ) ) break;
                    if ( line0 == "" ) line0 = line;
                    else first = False;
                    Bool saw_letter = False;
                    for ( int j = 0; j < line.isize( ); j++ )
                    {    if ( isspace( line[j] ) )
                         {    if (saw_letter) break;    }
                         else saw_letter = True;
                         if ( first && line[j] == '+' ) ++pluses;
                         if ( j != read_length + pluses )
                         {    if ( j > 2*read_length || line0[j] == line[j] ) 
                                   out << line[j];
                              else out << START_CYAN << line[j] << END_ESCAPE;    }
                         else
                         {    first = False;
                              if ( line[j] == 'A' )
                                   out << START_MAGENTA << "A" << END_ESCAPE;
                              else if ( line[j] == 'C' )
                                   out << START_GREEN << "C" << END_ESCAPE;
                              else if ( line[j] == 'G' )
                                   out << START_BLUE << "G" << END_ESCAPE;
                              else if ( line[j] == 'T' )
                                   out << START_RED << "T" << END_ESCAPE;
                              else out << line[j];    }    }
                    out << "\n";    }    
               Remove( OUT + "/pile" );    }    }

     // Add directory entries to final reports.

     if ( mutation_reports4.nonempty( ) )
     {    String s = mutation_reports4[0];
          int lastcol = s.RevBefore( " " ).RevAfter( " " ).Before( "." ).Int( );
          for ( int i = 0; i < nmutf; i++ )
          {    mutation_reports4[i] +=
                    " " + ToString(lastcol+2) + ".report_dir " + dirs[i];    }    }

     // Write final reports.

     Cp( OUT + "/mutation_reports4", OUT + "/mutation_reports4.backup" );
     vec<String> mutation_reports4_export;
     for ( int i = 0; i < mutation_reports4.isize( ); i++ )
     {    String tail = mutation_reports4[i].After( "reference_fasta_record_" );
          int record_id = tail.Before( ":" ).Int( );
          mutation_reports4_export.push_back( ToString( refdict[record_id] ) 
               + ":" + mutation_reports4[i].After( ":" ) );    }
     WriteStrings( OUT + "/mutation_reports4", mutation_reports4_export );
     SymlinkForce( OUT + "/mutation_reports4", OUT + "/mutation_reports" );
     cout << Date( ) << ": SomaticCall finished, elapsed time = "
          << TimeSince(clock) << endl;    }
