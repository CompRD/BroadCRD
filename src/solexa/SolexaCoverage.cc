/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/// SolexaCoverage.  Read in Solexa data, as determined by HEAD, which can be a
/// list.  Then generate information about coverage of the reference:
///
/// 1. Base-by-base display, along the reference, at each position showing a dot
/// if a read has a base that agrees with the reference there, A, C, G, or T if
/// an incorrect base appears in the read, and D or I if the read is deleted or
/// inserted relative to the reference at that position.
///
/// 2. Coverage by GC content.
///
/// Output to HEAD.SolexaCoverage[.OUT_SUFFIX] if HEAD is a singleton and OUT is
/// not specified.  Otherwise output to OUT.OUT_SUFFIX.
///
/// SHOW_INTENSITY: if True, print average intensity (K) of all bases matching
/// a given reference base.
///
/// MAX_ERRORS: if set, do not use alignments having more than this number of errors.
///
/// LEFT_TRIM: if set, declare that reads have been trimmed on left, infer trim
/// amounts from the read length and the run read length.
///
/// Whatever the output file X is, also generate another file X.binary which is a
/// VecDumbcallVec showing the base calls along the reference.
/// TSV_ONLY: if true have TSV output: contig \t pos \t base \t base cover \t
///           A cover \t C cover \t G cover \t T cover \t
///           DIFFERENT (if DIFF options satisfied) else space \n
///    The cover value is 0 for the base.  For example if the reference base
///    is covered 10 times and is a C, we may see 4 0 20 6, where the C column
///    is 0.  Of course, the value is actually shown in the base cover column.
///    Above default is False (original behavior).

#include "Basevector.h"
#include "FetchReads.h"
#include "FeudalMimic.h"
#include "Floatvector.h"
#include "math/Functions.h"
#include "Intvector.h"
#include "MainTools.h"
#include "PackAlign.h"
#include "Qualvector.h"
#include "lookup/LookAlign.h"
#include "math/HoInterval.h"
#include "paths/KmerPath.h"
#include "paths/ReadsToPathsCoreX.h"
#include "polymorphism/DumbCall.h"
#include "solexa/FourBase.h"
#include "solexa/SolexaMetrics.h"
#include "solexa/SolexaTools.h"
#include "dna/Bases.h"
#include "random/Random.h"

String globalPipelineDir;


bool dirSortOp(const String& x, const String& y) {
        String FCx = x.RevAfter("/");
        String FCy = y.RevAfter("/");

        if (FCx != FCy)
                return(FCx < FCy);
        else {
                String gPD = globalPipelineDir + "/";
                String Datex1 = x.After(gPD);
                String Datey1 = y.After(gPD);
                String Datex2 = Datex1.Before("_");
                String Datey2 = Datey1.Before("_");
                return(Datex2 < Datey2);
        }
}



int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(HEAD);
     CommandArgument_Bool_OrDefault(TSV_ONLY, False);
     CommandArgument_String_OrDefault(QLTOUT, "qltout");
     CommandArgument_String_OrDefault(QUAL_PREFIX, "");
     CommandArgument_String_OrDefault(REF, "");
     CommandArgument_Double_OrDefault(MIN_RATIO, 1.0);
     CommandArgument_Double_OrDefault(MIN_RATIO_10, 1.0);
     CommandArgument_Double_OrDefault(MIN_INTENSITY, 0.0);
     CommandArgument_Double_OrDefault(MIN_CYCLE_QUAL, 0.0);
     CommandArgument_Double_OrDefault(MIN_BASE_QUAL, 0.0);
     CommandArgument_Int_OrDefault(GC_WINDOW, 0);
     CommandArgument_Int_OrDefault(MAX_ERRORS, 0);
     CommandArgument_String_OrDefault(OUT, "");
     CommandArgument_String_OrDefault(OUT_SUFFIX, "");
     CommandArgument_String_OrDefault(HEAD_SUFFIX, "");
     CommandArgument_Bool_OrDefault(PAIRED, False);
     CommandArgument_String_OrDefault_Doc(PAIRED_HEADS_QLTOUT, "",
          "Very kludgy way of overriding naming of qltout files in PAIRED mode.");
     CommandArgument_Double_OrDefault_Doc(PAIRED_WIDTH, 4.0, "Number of stddev away from mean is still considered an acceptable pair.");
     CommandArgument_Bool_OrDefault_Doc(PAIRED_UNIQUE, False, "Both ends of the pair must align uniquely.");
     CommandArgument_Bool_OrDefault_Doc(PAIRED_UNIQUE_ONE_END, False, "One end of the pair must align uniquely.");
     CommandArgument_Bool_OrDefault_Doc(PAIRED_USE_MEDIAN, True, "Use pairs with distance closest to median if true, else if false use average.");
     CommandArgument_Bool_OrDefault(SHOW_INTENSITY, False);
     CommandArgument_Bool_OrDefault(COVERAGE_HISTOGRAM, False);
     CommandArgument_Bool_OrDefault(UNIQUE_ONLY, False);
     CommandArgument_Bool_OrDefault(COVERABILITY_ANALYSIS, False);
     CommandArgument_Bool_OrDefault(LEFT_TRIM, False);
     CommandArgument_Bool_OrDefault_Doc(NO_DUMB_CALLS, False, "If True, does not output the binary file.");
     CommandArgument_Int_OrDefault(TRUNC, 0);
     CommandArgument_Int_OrDefault(DIFF_MIN_COVER, 5);
     CommandArgument_Double_OrDefault(DIFF_MIN_FRAC, 0.1);
     CommandArgument_String_OrDefault_Doc(TRACE, "",
          "Find reads that support a particular base.  Format is "
          "tig.pos.base, e.g. 0.456.C.");
     CommandArgument_String_OrDefault(ALIGNS_FILE, "");
     CommandArgument_String_OrDefault(ALIGNS_DIRECTION, "");
     CommandArgument_String_OrDefault(PIPELINE, Getenv("SOLEXA_PIPELINE_DIR"));
     CommandArgument_Bool_OrDefault_Doc(SKIP_INDELS, False, "If true, exclude any reads containing indels.");
     CommandArgument_String_OrDefault_Doc(COV_BASES_USED, "", "If supplied, will create a vecbitvector of this name with the same dimensions as the READS indicating whether each base was included in coverage");
     EndCommandArguments;

     Bool loadIntensities, dummyCall;
     if (MIN_INTENSITY > 0.0) {
        loadIntensities = True;
        dummyCall = False;
     } else {
        loadIntensities = False;
        dummyCall = True;
     }

     globalPipelineDir = PIPELINE;
     int gPDlen = globalPipelineDir.isize();
     if (gPDlen > 0) {
        char lastchar = globalPipelineDir[gPDlen-1];
        if (lastchar == '/')
           globalPipelineDir.resize(gPDlen-1);
        lastchar = globalPipelineDir[globalPipelineDir.size()-1];
     }

     // Find HEADS, setup output.

     vec<String> HEADS;
     ExpandHead( HEAD, HEADS, TRUNC, PIPELINE, PAIRED );
     // Using dirSortOp, directories arre sorted by FC.lane, and within
     // each FC.lane, are sorted by date.
     sort(HEADS.begin(), HEADS.end(), dirSortOp);
     vec<String> paired_heads_qltout;
     if ( PAIRED_HEADS_QLTOUT != "" )
          ParseStringSet( PAIRED_HEADS_QLTOUT, paired_heads_qltout );

     if ( !HEADS.solo( ) ) ForceAssert( OUT != "" );
     String outfile = ( OUT == "" ? HEADS[0] + ".SolexaCoverage" : OUT );
     if ( OUT_SUFFIX != "" ) outfile += "." + OUT_SUFFIX;
     cout << "Writing to " << outfile << "." << endl;
     Ofstream( out, outfile );
     command.PrintTheCommandPretty( out );

     // Parse TRACE argument.

     int trace_tig = -1, trace_pos = -1, trace_base = -1;
     if ( TRACE != "" )
     {    trace_tig = TRACE.Before( "." ).Int( );
          trace_pos = TRACE.Between( ".", "." ).Int( );
          trace_base = as_char( TRACE.After( "." ).After( "." )[0] );    }

     // Load reference.

     int align_direction = 0 ;
     if ( ! ALIGNS_DIRECTION.empty() ) {
       if ( ALIGNS_DIRECTION == "FW" ) {
	 align_direction = 1;
       } else {
	 if ( ALIGNS_DIRECTION == "RC" ) {
	   align_direction = 2;
	 } else {
	   cout << "Illegal value passed for ALIGNS_DIRECTION argument" << endl;
	 }
       }
     }
     vecbasevector ref;
     if ( REF == "" ) RequireSameRef(HEADS);
     FetchReads( ref, 0, REF == "" ? HEADS[0] + ".reference.fasta" : REF );

     // Set up to track coverage of reference.

     longlong total_bases = 0;
     vec<longlong> cov_this_tig( ref.size( ), 0 );
     vec<VecIntVec> calls(6);
     vec< VecFloatVec > Icalls(4);
     for ( int i = 0; i < 6; i++ )
          Mimic( ref, calls[i] );
     for ( int i = 0; i < 4; i++ )
          Mimic( ref, Icalls[i] );
     vec<int> cov(1000000, 0);
     VecDumbcallVec dumbcalls;
     Mimic( ref, dumbcalls );



     vecbitvector pairalignOK;
     vector< vector<int> > pairWhichAlignA(HEADS.isize( ));
     vector< vector<int> > pairWhichAlignB(HEADS.isize( ));
     if (PAIRED) {
             pairalignOK.clear();
             int numReads = 0;
             for ( int ih = 0; ih < HEADS.isize( ); ih++ ) {
                String temphead = HEADS[ih];
                if ( HEAD_SUFFIX != "" ) temphead += "." + HEAD_SUFFIX;
                vecbasevector tempbases(temphead + ".fastb");
                numReads += tempbases.size();
             }
             pairalignOK.Reserve(numReads + HEADS.isize(), numReads);
             pairWhichAlignA.resize( HEADS.isize() );
             pairWhichAlignB.resize( HEADS.isize() );
             static bitvector b;
             for ( int ih = 0; ih < HEADS.isize( ); ih +=2 ) {
                String HEADA = HEADS[ih];
                String HEADB = HEADS[ih+1];
                String HEADPLUSA = HEADA;
                if ( HEAD_SUFFIX != "" ) HEADPLUSA += "." + HEAD_SUFFIX;
                String HEADPLUSB = HEADB;
                if ( HEAD_SUFFIX != "" ) HEADPLUSB += "." + HEAD_SUFFIX;
                vecbasevector basesA, basesrcA;
                VecFourBaseVec IA;
                vec<look_align> alignsA;
                vec< vec<int> > aligns_indexA;
                LoadSolexaData( HEADPLUSA, basesA, basesrcA, IA, loadIntensities );
                int nreadsA = basesA.size( );
                vec<float> minq10A;
                GetMinQ10( IA, nreadsA, minq10A, dummyCall );
                String aligns_fileA;
                // ALIGNS_FILE makes no sense in paired mode.
                aligns_fileA = HEADPLUSA + "." + QLTOUT;
                if ( PAIRED_HEADS_QLTOUT != "" )
                     aligns_fileA = paired_heads_qltout[ih] + "." + QLTOUT;
                LoadLookAligns( aligns_fileA, alignsA, aligns_indexA, nreadsA );
                solexa_metric_db db( HEADA + ".metrics" );
                vecbasevector basesB, basesrcB;
                VecFourBaseVec IB;
                vec<look_align> alignsB;
                vec< vec<int> > aligns_indexB;
                LoadSolexaData( HEADPLUSB, basesB, basesrcB, IB, loadIntensities );
                int nreadsB = basesB.size( );
                vec<float> minq10B;
                GetMinQ10( IB, nreadsB, minq10B, dummyCall );
                String aligns_fileB;
                // ALIGNS_FILE makes no sense in paired mode.
                aligns_fileB = HEADPLUSB + "." + QLTOUT;
                if ( PAIRED_HEADS_QLTOUT != "" )
                     aligns_fileB = paired_heads_qltout[ih+1] + "." + QLTOUT;
                LoadLookAligns( aligns_fileB, alignsB, aligns_indexB, nreadsB );
                // We do not need the second db - should be same as first.
                // solexa_metric_db dbB( HEADB + ".metrics" );

                AssertEq(basesA.size(),basesB.size());

                b.clear().resize( basesA.size( ) );
                pairalignOK.push_back(b);
                b.clear().resize( basesB.size( ) );
                pairalignOK.push_back(b);

                pairWhichAlignA[ih].resize(basesA.size( ));
                pairWhichAlignB[ih+1].resize(basesB.size( ));

                String tempstr;
                db.GetValue("paired_read_average_fragment_length", tempstr);
                double avgfrag = (double)tempstr.Double();
                db.GetValue("paired_read_median_fragment_length", tempstr);
                double medfrag = (double)tempstr.Double();
                //  double maxfrag = (avgfrag < medfrag) ? medfrag : avgfrag;
                db.GetValue("paired_read_std_dev_fragment_length", tempstr);
                double stddevfrag = (double)tempstr.Double();
                db.GetValue("paired_read_avg_dev_fragment_length", tempstr);
                double avgdevfrag = (double)tempstr.Double();
                // double maxdev = (stddevfrag < avgdevfrag) ? avgdevfrag : stddevfrag;
                double theAvg = 0, maxdist = 0, mindist = 0;
                if (PAIRED_USE_MEDIAN) {
                        theAvg = medfrag;
                        maxdist = medfrag + PAIRED_WIDTH*avgdevfrag;
                        mindist = medfrag - PAIRED_WIDTH*avgdevfrag;
                } else {
                        theAvg = avgfrag;
                        maxdist = avgfrag + PAIRED_WIDTH*stddevfrag;
                        mindist = avgfrag - PAIRED_WIDTH*stddevfrag;
                }

                // We do not do if(mindist < 0)mindist=0, since
                // David says overlapping reads are OK.


                for (size_t id = 0; id < basesA.size(); id++) {
                   if ( minq10A[id] < MIN_RATIO_10 ) continue;
                   if ( minq10B[id] < MIN_RATIO_10 ) continue;

                   int laAs = aligns_indexA[id].isize();
                   int laBs = aligns_indexB[id].isize();

                   // Obviously...
                   if (laAs == 0 || laBs == 0) continue;

                   if (PAIRED_UNIQUE && (laAs != 1 || laBs != 1)) continue;
                   if (PAIRED_UNIQUE_ONE_END && (laAs != 1 && laBs != 1)) continue;

                   int min_errors = 1000000000;
                   vec<int> at_minA, at_minB;
                   for ( int j = 0; j < laAs; j++ )
                   {    const look_align& la = alignsA[ aligns_indexA[id][j] ];
                        min_errors = Min( min_errors, la.Errors( ) );    }
                   for ( int j = 0; j < laAs; j++ )
                   {    const look_align& la = alignsA[ aligns_indexA[id][j] ];
                        if (la.Errors( ) == min_errors) at_minA.push_back(j); }
                   min_errors = 1000000000;
                   for ( int j = 0; j < laBs; j++ )
                   {    const look_align& la = alignsB[ aligns_indexB[id][j] ];
                        min_errors = Min( min_errors, la.Errors( ) );    }
                   for ( int j = 0; j < laBs; j++ )
                   {    const look_align& la = alignsB[ aligns_indexB[id][j] ];
                        if (la.Errors( ) == min_errors) at_minB.push_back(j); }

                   int amAs = at_minA.isize();
                   int amBs = at_minB.isize();

                   if (amAs == 0 || amBs == 0) continue;

                   // Below: find the closest pair to average distance.
                   double best_dist = 1000000000.0;
                   int d = 1000000000;
                   int best_indexA = -1, best_indexB = -1;
                   for (int i = 0; i < amAs; ++i) {
                      for (int j = 0; j < amBs; ++j) {
                         const look_align& laA = alignsA[ aligns_indexA[id][at_minA[i] ] ];
                         const look_align& laB = alignsB[ aligns_indexB[id][at_minB[j] ] ];
                         int dist = 1000000000;
                         if (laA.target_id == laB.target_id &&
                             laA.rc1 != laB.rc1) {
                            ho_interval exA, exB;
                            exA = laA.Extent2();
                            exB = laB.Extent2();
                            dist = Distance(exA, exB);
                         }
                         double temp_dist = fabs(theAvg - dist);
                         if (temp_dist < best_dist) {
                            d = dist;
                            best_dist = temp_dist;
                            best_indexA = i;
                            best_indexB = j;
                         }
                      }
                   }

                   if (d >= mindist && d <= maxdist) {
                     pairalignOK[ih].Set(id, 1);
                     pairalignOK[ih+1].Set(id, 1);
                     pairWhichAlignA[ih][id] = at_minA[best_indexA];
                     pairWhichAlignB[ih+1][id] = at_minB[best_indexB];
                   } else {
                     // Below will give an error if ever used.
                     pairWhichAlignA[ih][id] = -1;
                     pairWhichAlignB[ih+1][id] = -1;
                   }
                }
             }
     }




     // Go through the heads.  Compute reference coverage.

     for ( int ih = 0; ih < HEADS.isize( ); ih++ )
     {    String HEAD = HEADS[ih];
          String HEADPLUS = HEAD;
          if ( HEAD_SUFFIX != "" ) HEADPLUS += "." + HEAD_SUFFIX;

          // Set up data structures.

          cout << "loading data from " << HEAD << endl;
          double clock = WallClockTime( );
          vecbasevector bases, basesrc;

          VecFourBaseVec I;
          vec<look_align> aligns;
          vec< vec<int> > aligns_index;
          LoadSolexaData( HEADPLUS, bases, basesrc, I, loadIntensities );

          // if requested build a vecbitvector the size of the READS
          // that indicates if a base in a read was used for coverage
          vecbitvector basesUsed;
          if (COV_BASES_USED != "") {
              Mimic(bases, basesUsed);
          }

          vecqualvector quals;
          if ( MIN_BASE_QUAL > 0.0 )
               quals.ReadAll( HEAD + QUAL_PREFIX + ".qualb" );
          String aligns_file;
            if (ALIGNS_FILE != "") { aligns_file = ALIGNS_FILE; }
            else                   { aligns_file = HEADPLUS + "." + QLTOUT; }
          if ( PAIRED_HEADS_QLTOUT != "" )
               aligns_file = paired_heads_qltout[ih] + "." + QLTOUT;
          int nreads = bases.size( );
          LoadLookAligns( aligns_file, aligns, aligns_index, nreads );
          solexa_metric_db db( HEAD + ".metrics" );

          vec< vec<double> > cycle_qual;
          if (db.Defined("cycle_qual"))
          {
            cycle_qual = db.ValueVecVecDouble( "cycle_qual" );
          }
          else
          {
            cycle_qual.resize(db.ValueInt("read_length"));
            for (int i = 0; i < db.ValueInt("read_length"); i++)
            {
               cycle_qual[i].resize(3);
               cycle_qual[i][0] = 100.0;
               cycle_qual[i][1] = 100.0;
               cycle_qual[i][2] = 100.0;
            }
          }

          int read_length = 0;
          if (LEFT_TRIM) read_length = db.ValueInt( "read_length" );
          cout << "done, time used = " << TimeSince(clock) << endl;

          // For each read, find the minimum quality of its first 10 bases.

          vec<float> minq10;
          GetMinQ10( I, nreads, minq10, dummyCall );

          // Look at coverage of reference.

          for ( int i = 0; i < nreads; i++ )
               total_bases += bases[i].size( );
          for ( int id = 0; id < nreads; id++ ) {
	       if ( minq10[id] < MIN_RATIO_10 ) continue;
               if ( UNIQUE_ONLY && !aligns_index[id].solo( ) ) continue;
               int min_errors = 1000000000;
               for ( int j = 0; j < aligns_index[id].isize( ); j++ ) {
		   const look_align& la = aligns[ aligns_index[id][j] ];
		   if ( align_direction == 1 && la.Rc1() ) continue;
		   if ( align_direction == 2 && ! la.Rc1() ) continue;
		   min_errors = Min( min_errors, la.Errors( ) );
	       }
               if ( min_errors == 1000000000 ) continue;
               if ( MAX_ERRORS > 0 && min_errors > MAX_ERRORS ) continue;
               if ( (PAIRED == True) && (pairalignOK[ih][id] == 0) ) continue ;
               static vec<int> at_min;
               at_min.clear( );
               for ( int j = 0; j < aligns_index[id].isize( ); j++ ) {
		   const look_align& la = aligns[ aligns_index[id][j] ];
		   if ( align_direction == 1 && la.Rc1() ) continue;
		   if ( align_direction == 2 && ! la.Rc1() ) continue;
		   if ( la.Errors( ) == min_errors ) at_min.push_back(j);
	       }
               int u = randomx( ) % at_min.isize( );
               int theAtMin = at_min[u];
               if (PAIRED == True) {
                 if (ih % 2 == 0)
                     theAtMin = pairWhichAlignA[ih][id];
                 else
                     theAtMin = pairWhichAlignB[ih][id];
               }
               const look_align& la = aligns[ aligns_index[id][ theAtMin ] ];

                // JRM
                if (SKIP_INDELS && (la.indels != 0)) { continue; }

	       if ( align_direction == 1 && la.Rc1() ) { cout << "ERROR! WRONG ALIGN DIRECTION!"; exit(1); }
	       if ( align_direction == 2 && !la.Rc1() ) { cout << "ERROR! WRONG ALIGN DIRECTION!"; exit(1); }
               int tig = la.target_id, p1 = la.pos1( ), p2 = la.pos2( );
               cov_this_tig[tig] += la.Pos2( ) - la.pos2( );
               FourBaseVec J;
               if (MIN_INTENSITY > 0) {
                  J = I[id];
                  if (la.rc1) J.ReverseMe( );
               }
               int n = bases[id].size( );
               qualvector q;
               if ( MIN_BASE_QUAL > 0.0 )
               {    q = quals[id];
                    q.resize(n);
                    if (la.rc1) q.ReverseMe( );    }
               int left_trim = ( !LEFT_TRIM ? 0 : read_length - n );

               for ( int j = 0; j < la.a.Nblocks( ); j++ )
               {    if ( la.a.Gaps(j) > 0 )
                    {    for ( int x = 0; x < la.a.Gaps(j); x++ )
                         {    int p1x = ( !la.rc1 ? p1 : n - p1 - 1 );
                              int p1y = p1x + left_trim;
                              if ( MIN_INTENSITY > 0 &&
                                   (J[p1].CallQuality( ) < MIN_RATIO
                                    || J[p1].MaxInt( ) < MIN_INTENSITY) )
                                 continue;
                              if ( cycle_qual[p1y][0] >= MIN_CYCLE_QUAL
                                   && ( MIN_BASE_QUAL == 0.0
                                        || q[p1] >= MIN_BASE_QUAL ) )
                              {    ++calls[4][tig][p2+x];
                                   if ( trace_base == 4 && trace_tig == tig
                                        && trace_pos == p2+x )
                                   {    cout << "trace: " << id
                                             << endl;    }    }    }
                         p2 += la.a.Gaps(j);    }
                    if ( la.a.Gaps(j) < 0 )
                    {    for ( int x = 0; x < -la.a.Gaps(j); x++ )
                         {    int p1x = ( !la.rc1 ? p1+x : n - (p1+x) - 1 );
                              int p1y = p1x + left_trim;
                              if ( MIN_INTENSITY > 0 &&
                                   (J[p1+x].CallQuality( ) < MIN_RATIO
                                    || J[p1+x].MaxInt( ) < MIN_INTENSITY) )
                                 continue;
                              if ( cycle_qual[p1y][0] >= MIN_CYCLE_QUAL
                                   && ( MIN_BASE_QUAL == 0.0
                                        || q[p1] >= MIN_BASE_QUAL ) )
                              {    ++calls[5][tig][p2];
                                   if ( trace_base == 5 && trace_tig == tig
                                        && trace_pos == p2 )
                                   {    cout << "trace: " << id
                                             << endl;    }    }    }
                         p1 -= la.a.Gaps(j);    }
                    for ( int x = 0; x < la.a.Lengths(j); x++ )
                    {    int p1x = ( !la.rc1 ? p1 : n - p1 - 1 );
                         int p1y = p1x + left_trim;
                         if ( MIN_INTENSITY > 0 &&
                              (J[p1].CallQuality( ) < MIN_RATIO
                               || J[p1].MaxInt( ) < MIN_INTENSITY) )
                            continue;
                         if ( cycle_qual[p1y][0] >= MIN_CYCLE_QUAL
                              && ( MIN_BASE_QUAL == 0.0
                                   || q[p1] >= MIN_BASE_QUAL ) )
                         {    char base;
                              if ( !la.rc1 ) base = bases[id][p1];
                              else base = basesrc[id][p1];
                              if ( p2 >= ref[tig].isize( ) )
                              {    cout << "Alignment of read " << id
                                        << " appears to be illegal.\n";
                                   cout << "Could it be that the reference "
                                        << "has changed?\n";
                                   TracebackThisProcess( );    }
                              ++calls[base][tig][p2];


//                              printf("%s id: %9d pos: %2d base: %c rc: %s \n", bases[id].ToString().c_str(), id, p1x, as_base(base), (la.rc1?"Yes":"No"));
                              // mark the position as being used or not against the original forward read position.  E.g. the
                              // last based used in an rc read will be marked as base 0 in this vecbitvector
                              if (COV_BASES_USED != "") basesUsed[id].Set(p1x,True);


                              if ( trace_base == base && trace_tig == tig
                                   && trace_pos == p2 )
                              {    cout << "trace: " << id << endl;    }
                              if (SHOW_INTENSITY)
                              {    if ( base == ref[tig][p2] )
                                   {    Icalls[base][tig][p2] +=
                                             J[p1].base(base);    }    }    }
                         ++p1;
                         ++p2;
                }
            }
        }


        // write out vecbitvector of basesUsed if requested
        if (COV_BASES_USED != "") basesUsed.WriteAll(COV_BASES_USED);
    }


     // Assess and print coverage.

     int tig_digits = 0;
     for ( size_t tig = 0; tig < ref.size( ); tig++ )
          tig_digits = Max( tig_digits, ToString( ref[tig].size( ) - 1 ).isize( ) );
     vec< vec<ho_interval> > holes( ref.size( ) );
     longlong total = 0, wrong = 0;
     longlong wrong_AC = 0, wrong_AG = 0, wrong_AT = 0;
     longlong wrong_CA = 0, wrong_CG = 0, wrong_CT = 0;
     longlong wrong_GA = 0, wrong_GC = 0, wrong_GT = 0;
     longlong wrong_TA = 0, wrong_TC = 0, wrong_TG = 0;
     for ( size_t tig = 0; tig < ref.size( ); tig++ )
     {    vec<int> observed, expected, covs;
          for ( int j = 0; j < ref[tig].isize( ); j++ )
          {    int count = 0;
               for ( int x = 0; x < 4; x++ )
                    count += calls[x][tig][j];
               covs.push_back(count);    }
          Sort(covs);
          double Q3Q1 = double( covs[ ( 3 * covs.size( ) ) / 4 ] )
               / double( covs[ covs.size( ) / 4 ] );
          if (!TSV_ONLY) {
             out << "\nReference contig "  << tig << " has size " << ref[tig].size( )
                  << "; high-quality coverage = " << ToString( Mean(covs), 2 )
                  << ", Q3/Q1 = " << ToString( Q3Q1, 2 ) << "\n";
             out << "Legend:\n";
             int field = 1;
             if ( GC_WINDOW > 0 )
             {    out << "Field " << field++ << ". GC% in " << GC_WINDOW
                       << "-base window centered at base\n";    }
             if (SHOW_INTENSITY)
             {    out << "Field " << field++ << ". Mean intensity of bases called "
                       << "as reference base\n";    }
             out << "Field " << field++ << ". Reference contig\n";
             out << "Field " << field++ << ". Position on reference\n";
             out << "Field " << field++ << ". Reference base[instances]\n";
             out << "Field " << field++ << ". Other read bases\n";
             out << "Field " << field++ << ". Flag if different from reference\n";
             out << "\n";
          }
          int gc = 0;
          if ( GC_WINDOW > 0 )
          {    ForceAssert( ref[tig].isize( ) >= GC_WINDOW );
               for ( int i = 0; i < GC_WINDOW; i++ )
                    if ( IsGC( ref[tig][i] ) ) ++gc;    }
          for ( int j = 0; j < ref[tig].isize( ); j++ )
          {    Bool empty = True;
               for ( int x = 0; x < 6; x++ )
                    if ( calls[x][tig][j] > 0 ) empty = False;
               if ( !empty ) continue;
               int k;
               for ( k = j + 1; k < ref[tig].isize( ); k++ )
               {    Bool empty = True;
                    for ( int x = 0; x < 6; x++ )
                         if ( calls[x][tig][k] > 0 ) empty = False;
                    if ( !empty ) break;    }
               if ( k - j >= 100 ) holes[tig].push_back( ho_interval( j, k ) );
               j = k - 1;    }
          for ( int j = 0; j < ref[tig].isize( ); j++ )
          {    if ( GC_WINDOW > 0 )
               {    if ( j >= GC_WINDOW/2 && j < ref[tig].isize( ) - GC_WINDOW/2 )
                    {    if ( IsGC( ref[tig][j+GC_WINDOW/2] ) ) ++gc;
                         if ( IsGC( ref[tig][j-GC_WINDOW/2] ) ) --gc;    }
                    if (!TSV_ONLY) {
                       RightPrecisionOut( out,
                            100.0 * double(gc) / double(GC_WINDOW), 1 );
                       out << "% ";
                    }
               }
               if (SHOW_INTENSITY)
               {    int count = 0;
                    float Isum = 0;
                    for ( int x = 0; x < 4; x++ )
                    {    count += calls[x][tig][j];
                         Isum += Icalls[x][tig][j];    }
                    if ( count == 0 ) out << "   N/A";
                    else if (!TSV_ONLY)
                    {    out << setiosflags(ios::fixed) << setprecision(1)
                              << setw(6) << ( Isum / float(count) ) / 1000.0
                              << resetiosflags(ios::fixed);    }    }
               int blanks = tig_digits + 1 - ToString(tig).isize( );
               if (!TSV_ONLY) out <<      String( blanks, ' ' ) << tig;
               else out << tig;
               blanks = tig_digits + 1 - ToString(j).size( );
               if (!TSV_ONLY) out <<      String( blanks, ' ' ) << j << "   ";
               else out << " \t " << j;
               int A = calls[0][tig][j];
               int C = calls[1][tig][j];
               int G = calls[2][tig][j];
               int T = calls[3][tig][j];
               int D = calls[4][tig][j];
               int I = calls[5][tig][j];
               dumbcalls[tig][j] = dumbcall( A, C, G, T, D, I );
               int c = 0, Ref = 0;
               char refbase = ref[tig][j];
               for ( int x = 0; x < 4; x++ )
               {    c += calls[x][tig][j];
                    if ( x == refbase ) Ref += calls[x][tig][j];    }
               observed.push_back(c);
               for ( int x = 4; x < 6; x++ )
                    c += calls[x][tig][j];
               total += c;
               wrong += c - Ref;
               if ( as_base(refbase) == 'A' )
               {    wrong_AC += C;
                    wrong_AG += G;
                    wrong_AT += T;    }
               if ( as_base(refbase) == 'C' )
               {    wrong_CA += A;
                    wrong_CG += G;
                    wrong_CT += T;    }
               if ( as_base(refbase) == 'G' )
               {    wrong_GA += A;
                    wrong_GC += C;
                    wrong_GT += T;    }
               if ( as_base(refbase) == 'T' )
               {    wrong_TA += A;
                    wrong_TC += C;
                    wrong_TG += G;    }
               if (!TSV_ONLY) {
                  out << as_base(refbase) << "[" << Ref << "]"
                      <<      String( 6 - ToString(Ref).isize( ), ' ' );
               } else {
                  out << " \t " << as_base(refbase) << " \t " << Ref;
               }
               if ( as_base(refbase) == 'A' ) A = 0;
               if ( as_base(refbase) == 'C' ) C = 0;
               if ( as_base(refbase) == 'G' ) G = 0;
               if ( as_base(refbase) == 'T' ) T = 0;
               if (!TSV_ONLY) {
                  if ( A > 0 ) out << "A[" << A << "]";
                  if ( C > 0 ) out << "C[" << C << "]";
                  if ( G > 0 ) out << "G[" << G << "]";
                  if ( T > 0 ) out << "T[" << T << "]";
                  if ( D > 0 ) out << "D[" << D << "]";
                  if ( I > 0 ) out << "I[" << I << "]";
               } else {
                  out << " \t " << A << " \t " << C << " \t " << G << " \t " << T;
               }

               ++cov[Ref];
               double xr = DIFF_MIN_FRAC * double(Ref);
               if ( ( double(A) > xr && A >= DIFF_MIN_COVER )
                    || ( double(C) > xr && C >= DIFF_MIN_COVER )
                    || ( double(G) > xr && G >= DIFF_MIN_COVER )
                    || ( double(T) > xr && T >= DIFF_MIN_COVER )
                    || ( double(D) > xr && D >= DIFF_MIN_COVER )
                    || ( double(I) > xr && I >= DIFF_MIN_COVER )  )
               { if (!TSV_ONLY)
                   out << " DIFFERENT";
                 else
                   out << " \t DIFFERENT";
               } else if (TSV_ONLY)
                  out << " \t ";
               out << "\n";    }
          if ( ref.size( ) != 1 ) continue; // REST NOT IMPLEMENTED FOR > 1 CONTIG.

          // Considering 21-base windows centered at a given base, how does coverage
          // vary as a function of GC content of the window?

          vec<double> GC( ref[0].size( ), 0.0 );
          for ( int j = 10; j < ref[0].isize( ) - 10; j++ )
          {    int gc = 0;
               for ( int k = j - 10; k <= j+10; k++ )
                    if ( IsGC(ref[0][k]) )
                         ++gc;
               GC[j] = 100.0 * double(gc) / 21.0;    }
          for ( int pass = 0; pass < 20; pass++ )
          {    double low = pass * 5.0;
               double high = low + 5.0;
               int x1 = 0, y1 = 0;
               for ( int j = 10; j < ref[0].isize( ) - 10; j++ )
               {    if ( low <= GC[j] && GC[j] < high )
                    {    for ( int u = 0; u < 6; u++ )
                              x1 += calls[u][tig][j];
                         y1++;    }    }
               if ( y1 > 0 && !TSV_ONLY )
               {    out << "coverage for windows having GC of " << low << " - "
                    << high << "%: " << setprecision(3) << double(x1)/double(y1)
                    << " (N = " << y1 << " windows)\n";    }    }    }

     // Generate coverage histogram.

     if (COVERAGE_HISTOGRAM)
     {    int maxcov;
          for ( maxcov = cov.isize( ) - 1; maxcov >= 0; maxcov-- )
               if ( cov[maxcov] != 0 ) break;
          if (!TSV_ONLY) {
             out << "\nCoverage histogram:\n";
             for ( int i = 0; i <= maxcov; i++ )
                  out << i << " " << cov[i] << "\n";    }    }

     // Generate coverability analysis.

     if (COVERABILITY_ANALYSIS && !TSV_ONLY)
     {    vecKmerPath paths, pathsrc;
          vec<big_tagged_rpint> pathsdb;
          ReadsToPathsCoreY( ref, 100, paths, pathsrc, pathsdb );
          int nref = ref.size( );
          size_t total = ref.sumSizes();
          vec< vec<ho_interval> > duplicated(nref);
          for ( int tig = 0; tig < nref; tig++ )
          {    int pos = 0;
               const KmerPath& g = paths[tig];
               vec<ho_interval> dups;
               for ( int i = 0; i < g.NSegments( ); i++ )
               {    for ( longlong j = g.Start(i); j <= g.Stop(i); j++ )
                    {    static vec<longlong> con;
                         Contains( pathsdb, j, con );
                         if ( con.size( ) > 1 )
                              dups.push_back( ho_interval( pos, pos + 100 ) );
                         ++pos;    }    }
               for ( int i = 0; i < dups.isize( ); i++ )
               {    int j;
                    for ( j = i + 1; j < dups.isize( ); j++ )
                         if ( dups[j].Start( ) > dups[j-1].Stop( ) + 1 ) break;
                    int start = dups[i].Start( ), stop = dups[j-1].Stop( );
                    duplicated[tig].push( start, stop );
                    i = j - 1;    }    }
          size_t dupcount = 0, holecount = 0, dupholecount = 0;
          for ( int tig = 0; tig < nref; tig++ )
          {    dupcount += TotalCovered( duplicated[tig] );
               holecount += TotalCovered( holes[tig] );
               vec<ho_interval> d;
               ExtractGivenCoverage( ref[tig].size( ), 1, duplicated[tig], d );
               for ( int i = 0; i < holes[tig].isize( ); i++ )
                    dupholecount += Overlap( holes[tig][i], d );    }
          out << "\nintervals of genome in duplicated 100-mers:\n";
          for ( int tig = 0; tig < nref; tig++ )
          {    for ( int j = 0; j < duplicated[tig].isize( ); j++ )
                    out << tig << "." << duplicated[tig][j] << " ("
                         << duplicated[tig][j].Length( ) << ")\n";    }
          out << "\nintervals of genome in 100-base uncovered stretches:\n";
          for ( int tig = 0; tig < nref; tig++ )
          {    for ( int j = 0; j < holes[tig].isize( ); j++ )
                    out << tig << "." << holes[tig][j] << " ("
                         << holes[tig][j].Length( ) << ")\n";    }
          out << "\nfraction of genome in duplicated 100-mers: "
               << PERCENT_RATIO( 3, dupcount, total ) << "\n";
          out << "fraction of genome in 100-base uncovered stretches: "
               << PERCENT_RATIO( 3, holecount, total ) << "\n";
          out << "fraction of genome shared between these two: "
               << PERCENT_RATIO( 3, dupholecount, total ) << "\n";    }

     double Q = -10.0 * log10( double(wrong)/double(total) );
     if (!TSV_ONLY) {
        out << "\nTotal trusted bases = " << ToString( double(total)/1000000.0, 1 )
            << " Mb\n";
        out << "\nInferred base quality: " << ToString( Q, 1 ) << endl;
        longlong wrong_AN = wrong_AC + wrong_AG + wrong_AT;
        longlong wrong_CN = wrong_CA + wrong_CG + wrong_CT;
        longlong wrong_GN = wrong_GA + wrong_GC + wrong_GT;
        longlong wrong_TN = wrong_TA + wrong_TC + wrong_TG;
        PRINT4_TO( out, wrong_AC, wrong_AG, wrong_AT, wrong_AN );
        PRINT4_TO( out, wrong_CA, wrong_CG, wrong_CT, wrong_CN );
        PRINT4_TO( out, wrong_GA, wrong_GC, wrong_GT, wrong_GN );
        PRINT4_TO( out, wrong_TA, wrong_TC, wrong_TG, wrong_TN );
     }
     if (!NO_DUMB_CALLS)
         dumbcalls.WriteAll( (outfile + ".binary").c_str() );
}
