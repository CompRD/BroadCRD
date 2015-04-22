///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "CoreTools.h"
#include "PairsManager.h"
#include "ParseSet.h"
#include "Qualvector.h"
#include "lookup/LookAlign.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/fosmid/Ebv.h"

// The following code might have been under svn control.  Found in file having
// time stamp July 6, 2004.
//
// Determine if a given read is likely composed entirely of human alpha satellite
// 171-bp repeats, 
//
// AGCATTCTCAGAAACTTCTTTGTGATGTGTGTATTCAACTCACAGAGTTGAACATTTCTTTTGATAGAGCAGTTTGG
// AAACACTCTTTTTGTAGAATCTGCAAGTGGATATTTGGAGCGCTTTGAGGATTATGGTGGAAAAGGGAATATCTTCA
// TATAAAAACTAGACAGA
// (consensus from Waye and Willard, Nucleic Acids Research 14 (1986), 6915-6927)
//
// References:
// 1. J S Waye and H F Willard
// Nucleic Acids Res.1986 September 11;14 (17): 6915-6927
// "Molecular analysis of a deletion polymorphism in alpha satellite of human
// chromosome 17: evidence for homologous unequal crossing-over and
// subsequent fixation."
// 2.  Dirk Schindelhauer and Tobias Schwarz
// Genome Research, Vol. 12, Issue 12, 1815-1826, December 2002
// "Evidence for a Fast, Intrachromosomal Conversion Mechanism From Mapping of
// Nucleotide Variants Within a Homogeneous -Satellite DNA Array".
// (Asserts that "Human alpha-satellite DNA comprises 2%-5% of the genome",
// but I'm not sure what the basis for this claim is.)
// 3. Choo, K.H., Vissel, B., Nagy, A., Earle, E., and Kalitsis, P. 1991. A survey 
// of the genomic distribution of satellite DNA on all the human chromosomes, ad
// derivation of a new consensus sequence. Nucleic Acids Res. 19: 1179-1182.
//
// Note: some human alpha sequences in Genbank e.g. L08556 seem to be substantially
// different from the consensus.  I'm not sure what this means.

#include "Alignment.h"
#include "math/Functions.h"
#include "PrintAlignment.h"
#include "pairwise_aligners/SmithWatFree.h"

Bool IsAlphaHumanRepeat( const basevector& b, const double max_score, double& score )
{
     // Set up rotated versions of human repeat and its reverse complement.

     basevector rep(684 * 2), repr( 684 * 2 );
     Bool first_call(True);
     if (first_call)
     {    basevector r;
          r.SetFromString(
"AGCATTCTCAGAAACTTCTTTGTGATGTGTGTATTCAACTCACAGAGTTGAACATTTCTTTTGATAGAGCAGTTTGG"
"AAACACTCTTTTTGTAGAATCTGCAAGTGGATATTTGGAGCGCTTTGAGGATTATGGTGGAAAAGGGAATATCTTCA"
"TATAAAAACTAGACAGA"
"AGCATTCTCAGAAACTTCTTTGTGATGTGTGTATTCAACTCACAGAGTTGAACATTTCTTTTGATAGAGCAGTTTGG"
"AAACACTCTTTTTGTAGAATCTGCAAGTGGATATTTGGAGCGCTTTGAGGATTATGGTGGAAAAGGGAATATCTTCA"
"TATAAAAACTAGACAGA"
"AGCATTCTCAGAAACTTCTTTGTGATGTGTGTATTCAACTCACAGAGTTGAACATTTCTTTTGATAGAGCAGTTTGG"
"AAACACTCTTTTTGTAGAATCTGCAAGTGGATATTTGGAGCGCTTTGAGGATTATGGTGGAAAAGGGAATATCTTCA"
"TATAAAAACTAGACAGA"
"AGCATTCTCAGAAACTTCTTTGTGATGTGTGTATTCAACTCACAGAGTTGAACATTTCTTTTGATAGAGCAGTTTGG"
"AAACACTCTTTTTGTAGAATCTGCAAGTGGATATTTGGAGCGCTTTGAGGATTATGGTGGAAAAGGGAATATCTTCA"
"TATAAAAACTAGACAGA" );
          for ( int i = 0; i < 2; i++ )
               for ( int j = 0; j < 679; j++ )
                    rep.Set( 679 * i + j, r[j] );
          repr.ReverseComplement(rep);
          first_call = False;    }

     // Smith-Waterman.

     Bool use_fw = True, use_rc = True;
     ForceAssertLt( b.isize( ), 1358 );
     alignment a, ar;
     int best_loc;
     vector<int> mgg, mggr;
     int m = 1000000, mr = 1000000, i = 1000000, ir = 1000000;

     if (use_fw)
     {    SmithWatFree( b, rep, best_loc, a );
          mgg = a.MutationsGap1Gap2( b, rep );
          m = mgg[0];
          i = mgg[1] + mgg[2];    }

     if (use_rc)
     {    SmithWatFree( b, repr, best_loc, ar );
          mggr = ar.MutationsGap1Gap2( b, repr );
          mr = mggr[0];
          ir = mggr[1] + mggr[2];    }

     Bool reversed = False;
     if ( mr < m )
     {    reversed = True;
          swap( m, mr );
          swap( i, ir );    }

     int n = b.size( );
     float mp = 100.0 * float(m) / float(n);
     float mrp = 100.0 * float(mr) / float(n);
     float ip = 100.0 * float(i) / float(n);
     float irp = 100.0 * float(ir) / float(n);
     score = mp + (3 * ip);

     return score <= max_score;     }

     /*
     PRINT(score);
     cout << setprecision(3) << mp << " " << ip << " " << mrp << " " << irp 
          << " " << b.Length( ) << "\n";
     align al, alr;
     al.UnpackFrom(a);
     alr.UnpackFrom(ar);
     if ( !reversed ) PrintVisualAlignment( True, cout, b, rep, al );
     else PrintVisualAlignment( True, cout, b, repr, alr );
     return False;    }
     */

void CleanSatelliteReads( const String& TMP, String targets,
     const double max_alpha_score, Bool verbose )
{    

     vec<String> tar;
     ParseStringSet( targets, tar );

     const String sat_A( "CCATTCCATTCCATTCCATTCCATTCCATTCCATTCCATTCCATTCCATT" );
     const int max_sat_A_mis = 8;

     const basevector GenBank_JX174276(
     "GTCGACCAAATATATATTATATACTGTGCATAAAATATGAAGGTACATCAAATATATACTATATGCTGCA"
     "CATAAAATATCAAAGTACCCAAAGTATGTATTATATACTGTACATAAAATATCAAAGTACCCAAAGTATA"
     "TATTATATACTGTACATAAAATATAAAATACCCCAAATATATATTATATGCTGTACATAAAATATCAAAG"
     "TACACAAAATATTTATTATATACTGTACGTAAAATACCAAAGGACCCAAAATATATATTATATACTGTAG"
     "ATAAAATATCAAAATTCCCAAAATATATATTATATACTGTACATAAAATACCAAAGTACCCAAAGTATAT"
     "ATTATATACCGTACAAAAAATATCAAAGTAGCCACAATATGTATTATACACTGTACATAATATATGAAAG"
     "TACCCAAACATTTATAATAAACTGTACATAAAATATCAATGTACTCCAACTATATATTATGTGCGGTACA"
     "TAAGATATCAAAGTACCCATACTATATAGTATATACTGTACTTAAAATATCAAAATACCCAAAATATATA"
     "TTATATACGTTACATAAACTATCAAATACACAAAGTATGCATTATATACTGTACATAAAATAGCAAAGTT"
     "CCCAAAATATGTATTATATACTGTACATAAAATATGAAAGTACCCAACCATATTTAATAAACTCTACATA"
     "AAATACCAAAGTACCCAAACTGTGTATTATATACTGTCCATAAAATATCAAAGTACCCAAAATATATATT"
     "GTAAACTCTACATAAAATATCAAAGTAGCCAAAATACATATTGTATACTGTACATAAAATATAAAAGTAC"
     "CCCAAATATATATTATATACTGTACATAAAATATCAAAGGACCCAAAATAAACATTATATACTGTACATA"
     "AAATATCCAAATACTCAAAATATATATTACATACTGTACATAAAATATCAAAGTACCCAAATTATGTATT"
     "ATATACTGTACGTAAAATATCAGAGTACCCAAAATATGCATTATATACTGTACATAAAATATGAAAGTAA"
     "CAAAACATTTATAATAAACTGTACATAAAATATCAAAGTACTCAAACTATATATCATATACTGTACATAA"
     "AATATCAAAGTACCCAAACTGTCTATTATTATATACTGTATATAAAATATCAAAGTACCTGAAATATATA"
     "CTGTATACTGTACATAAAATATCAAAGTACCCAAAATATATACTATATACTGTACATAAAATATGAAGGT"
     "ACATCAAATATATATTATATCCTGTACATAAAATATCAAAGTACACCTATTATATATTGTGTACTGTACA"
     "TAATATTTTAAATCGCCCAAATATATATTGTATACTGTATATAAAATATGAAGGTACTTCAAATATATAT"
     "TATATGCTGCACATAAAATATCAAAGTACCCAAAGTATGCATTACATACTGTACATAAAATATCAAAGTA"
     "CACAAATTATATATTATATACTGTACATAAAATATGTAGGTTCATCAAATATATTTTATATTCTGTACAT"
     "AATATATCAAAGTACATCAAATATATATTATATACGGTACATAAAATATCAAAGTACCCAAAATATGTAT"
     "TATATACTGTACATAAAACATGAAAGTAACCAAACATTTATAATGAACTGTACATAAAATATCAAAGTAC"
     "TCAATCTACATATTATATACACTATATAAAATATCAAAGTACACAAACTGTATATTATATACTGTACATA"
     "ATATATCAAAGTACCCAGAATATATAGCATACTGTACATAAAATATCAAAGTACCCAAAATATATATTGT"
     "ATACTGTACATAAAATATGAAGGTACATCAAAGATATATTATATTCTGTACATAAAATATCAAAGTACAC"
     "CAAATATATATTATATGCTGTACATAATATATTAAAGTTGCCCATATATACATTCTATACTGTACATAAA"
     "ATATGAAGGTACATCAAATATATAGTATATGCTGCACATAAAATATCAAAGTACCCAAAGTATGTATTAT"
     "ACACTGTACATAAAATATCAAAGTACCAAAAGTATATATTACATACTGTAAATAAAATATGAAGGTACCT"
     "CAAATATATATTATGTTCTGTACATAAAATATCAAAGTACACCAAGTATATACTATATACTCTGCATAAA"
     "ATATTAAATTACCCAAAATATGTATTGTATTGTGTACATAAAATATGAAAGTAACCAAACATTTATAATA"
     "AACTGTACATAAGATATCAAAGTACTCAAACTGTATATTATATACTGTACATAAAGTATCAAAGTACCCA"
     "AACTGTGTAGTATATACTGTACATAAAATATCAAGGTACCCAAAATATATATTGTAAACTGTAGATAAAA"
     "TATCAAAGTACCCAAAATATATATTAAATACTGTACATAAAATATGAAGGCACATCCAATATATATTATA"
     "TTCTGTACATAAAATATCAAACTACACCAAATGTATATTATATACTGTACATAAAATATCAAAGTACCCA"
     "AATGTATATTATATACTGTACATAAAATATCAAAGTCCCCAAAGTGTATATTATATACTGTACATAAAAT"
     "ATGAAGGTTCATCAAATATATTTTATATTCTGCACATAATATATCAAAGTACAGCAAATATATATTATAT"
     "ACTGCACATAAAATATCAAAGTACCCAAAATATGTATTATATACTGTACATAAAATATCCAAATACTCAA"
     "AATATATATTATATACTGTACATAAAATATCAAAGTACCCAAATTATGTATTATATACTGTACGTAAAAT"
     "ATCAGAGTACCCAAAATATGCATTATATACTGTACATAAAATATGAAAGTAACAAAACATTTATAATACA"
     "CTGTACATAAAATATCAAAGTACTCAAACTATATATCATATACTGTACATAAAATATCAAAGTACCCAAA"
     "CTGTGTATTATTACATACTGTATATAAAATATCAAAGTACCTGAAATATATATTGTATACTGTACATAAA"
     "ATATCAAAGTACCCAAAATATATACTATATACTGTACATAAAATATGAAGGTACATCCAATATATATTAT"
     "ATCCTGTACATGAAATATCAAAGTACACCTATTATATATTGTGTACTGTACATAATATATTAAATCGCCC"
     "AAATATATATTGTATACTGTATATAAAATATGAAGGTACTTCAAATTCATATTACATGCTGCACATAAAA"
     "ATCAAAGTACCCAAAGTATGCATTACATACTGTACATAAAATATCAAAGTACACAAATTATATATTATAT"
     "ACTGTACATAAAATATGAAGGTTCATCAAATATATTTTATATTCTGTACATAATATATCAAAGTACATCA"
     "AATATATATTATATACTGTACATAAAATATCAAAGTACCAAAAATATGTGTTATATACTGTATATAAAAT"
     "ATGAAAGTAACCAAACATTTATAATAAATTGTACATAAAATATGAAAGTACTCAATCTACATATTATATA"
     "CAGTATATAAAATATCAAAGTACCCAAACTGTATATTATATACTGTACATAATATATCAAAGTACCCAGA"
     "ATATATAGCATACTGTACATACAATATCAAAGTACCCAGAATGTATATTGTATACTGTACATAAAATATG"
     "AAGGAACATCAAAGATATATTATATTCTTTACATAAAATATCAAAGTACACCAAATATATATTATATGCT"
     "GTACATAATATATTAAACTTGCCCATGTATAGATTTTATACTGTACATAAAATATGAAGGTACATCAAAT"
     "ATATAGTATATGCTGCACATAAAATATCAAAGTACCCAAAGTATGTATTATATACTGTACATAAAATATC"
     "AAAGTACCAAAATTATATATTACCTACTGTAAATAAAATATGAAGATACCTCAAATATATATTATATTCT"
     "GTACATAAAATGTCAAAGTACACCAAGTATATACTATATACTCTACATAAAATATTAAATTACCCAAAAT"
     "ATGTATTGTATTGTGTACATAAAATATGAAAGTAACCAAACATTTATAATAAACTGTACATAAAATGTCA"
     "AAGTACTCAAACTGTATATTATATACTGTACATAAAGTATCAAAGTACCCAAACTGTGTAGTATATACTG"
     "TACATAAAATATCAAGGTACCCAAAATATATATTGTAAACTGTAGATAAAATATCAAAGTACCCAAAATA"
     "TATATTAAATACTGTACATAAAATATGAAGGCACATCCAATATATATTATATTCTGTACATAAAATATCA"
     "AAGTACACCAAATGTATATTATATACTGTACATAAAATATCAAAGTCCCCAAAGTGTATATTATATACTG"
     "TACATAAAATATGAAGGTTCATCAAATATATTTTATATTCTGCACATAATATATCAAAGTACAGCAAATA"
     "TATATTATATACTGTACATAAAATATCAAAGTACCCAAAATATGTATTATATACTGTACATAAAATGGGA"
     "AAGTAACCAAACATTTATAATAAACTGTACATAAAATATCAAAGTACTCAAACTACATATTATATACAGT"
     "ATATAAAATGTCAAAGTACCCAAACTGTATATTATATACTGTACATAATAAATCAAAGTACCCGAATATA"
     "TAGCATACTGTACATAAAATATCAAAGTACCCAAAATATATATTGTATACTGTACATAAAATATGAAGGT"
     "ACATCAAAGATATATTATATTCTGTACATAAAATATCAAAGTACACCAAATATATATTATATGCTGTACA"
     "GAATATATTAAAGTCGAC" );

     // Homo sapiens clone hsat2cln40 satellite 2 sequence.
     const basevector GenBank_AF361742(
     "ATGGAAATGAAAGGGGTCATCATCTAATGGAATCGCATGGAATCATCATCAAATGGAATCGAATGGAATC"
     "ATCATCAAATGGAATCTAATGGAATCATTGAACAGAATTGAATGGAATGGTCATCGTATGAATTGAATGC"
     "AATCATCGAATGGTCTCGAATGGAATCATCTTCTAATGGAAAGGAATGGAATCATCACATAGAATCGAAT"
     "GGAATTATCATCGAATGGAATCGAATGGTATCAACACCAAACGGAAAAAAAAGGAATTATCGAATGGAAT"
     "CGAAGAGAATCTTCGAACGGACCCGAATGGAATCATCTAATGGAATGGAATGGAATAATCCATAGACTCG"
     "AATGCAATCATCATTGTATAGAATCGAATGGAATCATCGAATGGACTCGAAT" );

     const basevector GenBank_JN211119_part(
     "GAATTCATGCCTAGGCTTTGCCTACAGGGGACATTGTGACATATCCCTGCACTGATCACCCAGGTGATGC"
     "AACGCTTCTCTATGCTCTGCCTACAGCGGACATTGTGACATATCTCTGCACTGATCACCCAGTTGATGCA"
     "ACTCTTCTCTATGCTCTGCCTACAGGGGGCATTGTGACATATCTCTGCACTGATCAACGAGGTGATGTAA"
     "CTCTTTTCTAGTCTCTGCCGACAGAGGGCGTTGTGACATCACTCTGCACGGATCACACGGTTTATGTAAC"
     "TCTTGACTAGGCTCTGCCTACGGGGGCATTTTCACATATCACTGCACTGATCACCGAGATGATGTAACTC"
     "TTGTATAGGCTTCGCCGACAGGGGGCATTGAGACATATCTCTTCACTGATCACCGAGGTGATGCAACTGT"
     "TGTCTGGGATCTGCTTACAGGCGGCATTGTGACATATCTCTGCCCTGATCACCCAGGTGATGTAACTCTT"
     "GTCTAGGCTCTGCCTACTGGAGACATTGTGACATATCTCCGCACTGATCACCCAGTTGATGTCACTATTG"
     "TCAAGGATATGGCTATAGAGATATTGTGACATGTCACTGCACTGATCACAGAGCTGATATAACTCTTGTC"
     "TAGGCTCTGGCAACAGGGGGCTAGTGACACATCTCTGCACTGATCACACAGGTGATGTAACTCTGTTCTA"
     "AGCTCTGCCTAAAGGGGCATTGTGACAGATCTCTGCACTGAGCACTCAGGTGGTGTAACTATTGTCTAGG"
     "CTCTGCTTCAAGGGGCCTTGTCACATATCTCTGCACTGATCACCCAGGTGATGTAACTCTTGTCTCGGCT"
     "CTGCTTACAGGGGGTATTGTGACATAACCCTGCACTGATCACCTAAGTGACGTAACACTTATGTAGTCTC"
     "TGCCTACAGTGGCATTTTGACATATCTCTGCACTGTTAACCGAGGTGATGAAACTCGTGTCTAGTCTGTG"
     "CCCACAGGGGGATTGAGACATATCTCTGCACTGATCCCGAGGTGATCCAACTCTTGCCCGGTCTCTGCCT"
     "ACTGGGGACATTGTGACATACCTCTGCTCTGATCACCCAGGTGCTGTAACTTTAGTGTAGGCTCTGGCTA"
     "CACGGCATTGTGACATATCACTGCACTGATTACCCAGGTGATATAACTCTTGTCTAGGCTCTGCCTACAG"
     "GGGGCTTGTGACATATCTCTGCACTGATCACCCAGGTGATATAACTCTTCTCTAGGATCTGCCTACAGGG" );

     String head = TMP + "/frag_reads_orig";
     vecbasevector bases( head + ".fastb" );
     vecqualvector quals( head + ".qualb" );
     PairsManager pairs;
     pairs.Read( head + ".pairs" );

     vec<Bool> pairs_to_delete( pairs.nPairs( ), False );
     vec<Bool> reads_to_delete( bases.size( ), False );

     // Check for EBV.

     if ( Member( tar, String("ebv") ) )
     {    vecbasevector all;
          all.push_back( EBV( ) );
          all.Append(bases);
          const int K = 60;
          vec< triple<kmer<K>,int,int> > kmers_plus;
          MakeKmerLookup2( all, kmers_plus );
          for ( int64_t i = 0; i < (int64_t) kmers_plus.size( ); i++ )
          {    int64_t j;
               for ( j = i + 1; j < (int64_t) kmers_plus.size( ); j++ )
                    if ( kmers_plus[j].first != kmers_plus[i].first ) break;
               Bool valid = False;
               for ( int64_t k = i; k < j; k++ )
                    if ( kmers_plus[k].second == 0 ) valid = True;
               if (valid)
               for ( int64_t k = i; k < j; k++ )
               {    if ( kmers_plus[k].second > 0 ) 
                         reads_to_delete[ kmers_plus[k].second - 1 ] = True;    }
               i = j - 1;    }
          for ( int pid = 0; pid < (int) pairs.nPairs( ); pid++ )
          {    vec<int> id(2);
               id[0] = pairs.ID1(pid), id[1] = pairs.ID2(pid);
               if ( reads_to_delete[ id[0] ] || reads_to_delete[ id[1] ] )
               {    pairs_to_delete[pid] = True;
                    reads_to_delete[ id[0] ] = True;
                    reads_to_delete[ id[1] ] = True;    }    }    }
          
     // Main loop.

     #pragma omp parallel for
     for ( int pid = 0; pid < (int) pairs.nPairs( ); pid++ )
     {    vec<int> id(2);
          id[0] = pairs.ID1(pid), id[1] = pairs.ID2(pid);

          // Check for total junk, reads having little except Q <= 2 bases.

          const int min_qual = 1000;
          const int min_qual_count = 3;
          int qual_sum = 0;
          for ( int j = 0; j < 2; j++ )
          {    int max_qual = 0;
               const basevector& b = bases[ id[j] ];
               const qualvector& q = quals[ id[j] ];
               for ( int l = 0; l < b.isize( ); l++ )
                    if ( q[l] >= min_qual_count ) qual_sum += q[l];    }
          if ( qual_sum < min_qual )
          {    if (verbose) cout << id[0] << ": low quality" << endl;
               if (verbose) cout << id[1] << ": low quality" << endl;
               pairs_to_delete[pid] = True;
               reads_to_delete[ id[0] ] = reads_to_delete[ id[1] ] = True;
               continue;    }

          vec<Bool> junky(2, False);
          for ( int j = 0; j < 2; j++ )
          {    if ( j == 1 && junky[0] ) break;
               basevector b = bases[ id[j] ];

               // Check for alpha satellite.

               if ( Member( tar, String("alpha") ) || Member( tar, String("all") ) )
               {    double score;
                    if ( IsAlphaHumanRepeat( b, max_alpha_score, score ) ) 
                    {    junky[j] = True;
                         if (verbose) 
                              cout << id[j] << ": alpha, score = " << score << endl;
                         break;    }    }

               if ( !Member( tar, String("all") ) && !Member( tar, String("two") ) ) 
                    continue;

               for ( int pass = 1; pass <= 2; pass++ )
               {    if ( pass == 2 ) b.ReverseComplement( );

                    // Check for "satellite 1" and other stuff.

                    const int maxsat1 = 20;
                    int best_loc;
                    alignment a;
                    basevector b0;
                    if ( pass == 1 ) b0.SetToSubOf( b, 0, 100 );
                    else b0.SetToSubOf( b, b.isize( ) - 100, 100 );
                    int errs;

                    if ( Member( tar, String("all") ) )
                    {    errs = SmithWatFree( b0, GenBank_JX174276, best_loc, a,
                              false, false, 1, 1 );
                         if ( errs <= maxsat1 ) 
                         {    junky[j] = True;
                              if (verbose) 
                              {    cout << id[j] << ": JX174276, errs = " 
                                        << errs << endl;    }
                              break;    }    }

                    if ( Member( tar, String("all") ) 
                              || Member( tar, String("two") ) )
                    {    errs = SmithWatFree( b0, GenBank_AF361742, best_loc, a,
                              false, false, 1, 1 );
                         if ( errs <= maxsat1 ) 
                         {    junky[j] = True;
                              if (verbose) 
                              {    cout << id[j] << ": satellite two, errs = " 
                                        << errs << endl;    }
                              break;    }    }

                    if ( Member( tar, String("all") ) )
                    {    errs = SmithWatFree( b0, GenBank_JN211119_part, best_loc, a,
                              false, false, 1, 1 );
                         if ( errs <= maxsat1 ) 
                         {    junky[j] = True;
                              if (verbose) 
                              {    cout << id[j] << ": JN211119, errs = " 
                                        << errs << endl;    }
                              break;    }    }
                    
                    // Test for (CCATT)^10.

                    if ( Member( tar, String("all") ) )
                    {    String s = b.ToString( );
                         for ( int i = 0; i <= s.isize( ) - sat_A.isize( ); i++ )
                         {    int mis = 0;
                              for ( int l = 0; l < sat_A.isize( ); l++ )
                              {    if ( sat_A[l] != s[i+l] )
                                   {    mis++;
                                        if ( mis > max_sat_A_mis ) break;    }    }
                              if ( mis <= max_sat_A_mis )
                              {    junky[j] = True;
                                   if (verbose) cout << id[j] << ": (CCATT)" << endl;
                                   break;    }    }    }
                    if ( junky[j] ) break;    }    }

          if ( junky[0] || junky[1] ) 
          {    pairs_to_delete[pid] = True;
               reads_to_delete[ id[0] ] = reads_to_delete[ id[1] ] = True;    }    }

     cout << Date( ) << ": identified satellite DNA in "
          << Sum(pairs_to_delete) << " pairs" << endl;

     vec<int> to_new( bases.size( ) );
     int new_id = 0;
     for ( int id = 0; id < (int) bases.size( ); id++ )
     {    to_new[id] = new_id;
          if ( !reads_to_delete[id] ) new_id++;    }
     pairs.removePairs(pairs_to_delete);
     for ( int pid = 0; pid < (int) pairs.nPairs( ); pid++ )
          pairs.SetIDs( pid, to_new[ pairs.ID1(pid) ], to_new[ pairs.ID2(pid) ] );
     bases.EraseIf(reads_to_delete);
     quals.EraseIf(reads_to_delete);
     bases.WriteAll( head + ".fastb" );
     quals.WriteAll( head + ".qualb" );
     pairs.Write( head + ".pairs" );
     if ( IsRegularFile( head + ".qltout" ) )
     {    vec<look_align> aligns; 
          LoadLookAligns( head + ".qltout", aligns );
          vec<Bool> to_delete( aligns.size( ), False );
          for ( int id = 0; id < aligns.isize( ); id++ )
               if ( reads_to_delete[ aligns[id].query_id ] ) to_delete[id] = True;
          EraseIf( aligns, to_delete );
          WriteLookAligns( head + ".qltout", aligns );    }    }
