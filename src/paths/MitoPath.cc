///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MitoPath.  A simple mitochondrial genome assembler.  It takes as input a BAM
// file containing reads, selects those reads that are aligned to the mitochondrial
// genome, and assembles them.  They should provide high coverage.
//
// Assumptions:
// - The program samtools is in your path.
// - All the programs listed below as dependencies are in your path.
// - All reads have the same length.
//
// Known defects:
// - Pairing information is not used.
// - Output format may not be what is needed.

// MakeDepend: dependency CleanCorrectedReads
// MakeDepend: dependency Col
// MakeDepend: dependency FindErrors
// MakeDepend: dependency RandomLines 
// MakeDepend: dependency SAM2CRDDump 

#include "Alignment.h"
#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "PackAlign.h"
#include "graph/Digraph.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "paths/HyperKmerPath.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/Unipath.h"

void SmithWatFreeSym( const basevector& b1, const basevector& b2, align& a )
{    alignment al;
     int best_loc;
     if ( b1.size( ) <= b2.size( ) )
     {    SmithWatFree( b1, b2, best_loc, al, false, false, 1, 1 );
          a.UnpackFrom(al);    }
     else
     {    SmithWatFree( b2, b1, best_loc, al, false, false, 1, 1 );
          a.UnpackFrom(al);
          a.Flip( );    }    }

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(BAM, "name of bam file");
     CommandArgument_String_OrDefault_Doc(MITO_NAME, "MT", 
          "name of mitochondrial genome in bam file");
     CommandArgument_String_Doc(HEAD, "head for output; generates HEAD.edges.fasta, "
          "HEAD.dot and many intermediate files");
     CommandArgument_Int_OrDefault_Doc(PHASE, 1,
          "If PHASE=1, begin with read extraction, which is slow.\n"
          "If PHASE=2, begin with error correction.\n"
          "If PHASE=3, begin with graph generation.");
     CommandArgument_Double_OrDefault_Doc(TARGET_COVERAGE, 1000,
          "in phase 1, aim to extract roughly this level of coverage");
     CommandArgument_Int_OrDefault_Doc(MIN_COMPONENT, 500,
          "in phase 3, graph components smaller than this size are discarded");
     CommandArgument_Bool_OrDefault_Doc(REMOVE_HANGING_ENDS, True,
          "in phase 3, clean the graph");
     CommandArgument_Double_OrDefault_Doc(MIN_FRAC, 0.1,
          "in phase 3, at branches, delete those edges whose coverage is below "
          "this fraction of the total");
     CommandArgument_Int_OrDefault_Doc(K, 40, 
          "in phase 3, minimum overlap in graph");
     CommandArgument_Int_OrDefault_Doc(DELTA_KEEP, 0,
          "given two parallel edges, delete the least covered one if its number "
          "of differences with the more covered one is at most this");
     EndCommandArguments;

     // Begin phase 1.

     if ( PHASE == 1 )
     {
          // Get stats for bam file.  Duplicates are excluding using the samtools
          // view flag '-F 1028'.

          if ( !IsRegularFile(BAM) )
          {    cout << "I can't find your BAM file." << endl;
               cout << "Abort." << endl;
               exit(1);    }
          fast_pipe_ifstream bin( "samtools idxstats " + BAM );
          String line, recname;
          vec<int64_t> gsizes;
          while(1)
          {    getline( bin, line );
               if ( bin.fail( ) ) break;
               int64_t gsize;
               istrstream iline( line.c_str( ) );
               iline >> recname >> gsize;
               if ( recname == MITO_NAME ) 
               {    if ( gsize == 0 )
                    {    cout << "The mitochondrial genome in your bam file is "
                              << "listed as having a genome size of zero.  This "
                              << "doesn't make sense.  I'm not sure what to do and "
                              << "am giving up." << endl;
                         exit(1);    }
                    gsizes.push_back(gsize);    }    }
          if ( gsizes.empty( ) )
          {    cout << "I couldn't find the mitochondrial genome in your bam file.\n"
                    << "Perhaps you need to specify a different name using "
                    << "MITO_NAME.  Meanwhile, I am giving up." << endl;
               exit(1);    }
          if ( !gsizes.solo( ) )
          {    cout << "I found " << gsizes.size( ) << " mitochondrial reference "
                    << "sequences in your bam file, rather than one, as required.\n"
                    << "Either you need to fix your bam file, or this program "
                    << "needs to be made smarter." << endl;
               cout << "Abort." << endl;
               exit(1);    }
          int64_t gsize = gsizes[0];
          cout << Date( ) << ": counting non-duplicate reads" << endl;
          int readlength = StringOfOutput( "samtools view " + BAM 
               + " | head -1 | Col 10 | wc -c" ).Int( );
          int64_t nreads = LineOfOutput( "samtools view -F 1028 " + BAM + " " 
               + MITO_NAME + " | wc --lines" ).Int( );

          // Extract the reads.

          double bam_coverage = double(readlength) * double(nreads) / double(gsize);
          double frac_to_use = Min( 1.0, TARGET_COVERAGE / bam_coverage );
          cout << Date( ) << ": estimated coverage by non-duplicates in BAM file = " 
               << bam_coverage << endl;
          cout << Date( ) << ": using " << setprecision(3)
               << 100.0 * frac_to_use << "% of reads" << endl;
          cout << Date( ) << ": extracting reads" << endl;
          SystemSucceed( "samtools view -F 1028 " + BAM + " " + MITO_NAME 
               + " | RandomLines FRAC=" + ToString(frac_to_use) 
               + " NH=True | SAM2CRDDump OUT_HEAD=" + HEAD 
               + " PRESERVE_ORIENTATION=True LOG_TO_CERR=False " 
               + "> " + HEAD + ".extraction.log" );    }

     // Begin phase 2, error correction.

     if ( PHASE <= 2 )
     {    cout << Date( ) << ": correcting errors" << endl;
          SystemSucceed( "FindErrors HEAD_IN=" + HEAD + " HEAD_OUT="
               + HEAD + "_edit > " + HEAD + ".FindErrors.log" );
          Cp( HEAD + ".pairs", HEAD + "_edit.pairs" );
          SystemSucceed( "CleanCorrectedReads DELETE=True HEAD_IN=" + HEAD + "_edit"
                  + " HEAD_OUT=" + HEAD + "_corr" );    }

     // Begin phase 3, graph generation.

     cout << Date( ) << ": building graph" << endl;
     vecbasevector genome( HEAD + "_corr.fastb" );
     cout << Date( ) << ": have " << genome.size( ) << " error-corrected reads"
          << endl;
     if ( genome.size( ) == 0 )
     {    cout << Date( ) << ": no reads survived error correction, giving up." 
               << endl;
          exit(1);    }
     int min_length = 1000000000, max_length = 0;
     for ( size_t i = 0; i < genome.size( ); i++ )
     {    min_length = Min( min_length, genome[i].isize( ) );
          max_length = Max( max_length, genome[i].isize( ) );    }
     cout << Date( ) << ": reads have length between " << min_length
          << " and " << max_length << endl;
     if ( max_length < K )
     {    cout << Date( ) << ": your reads are all shorter than K" << endl;
          cout << Date( ) << ": please try using a smaller K, perhaps 20" << endl;
          cout << Date( ) << ": giving up" << endl;
          return 1;    }
     vecKmerPath paths, paths_rc, unipaths;
     vec<big_tagged_rpint> pathsdb, unipathsdb;
     ReadsToPathsCoreY( genome, K, paths );
     CreateDatabase( paths, paths_rc, pathsdb );
     Unipath( paths, paths_rc, pathsdb, unipaths, unipathsdb );

     // Generate initial graph and prune it.

     Ofstream( out, HEAD + ".dot" );
     digraph A;
     BuildUnipathAdjacencyGraph( paths, paths_rc, pathsdb, unipaths, unipathsdb, A );
     HyperKmerPath h;
     BuildUnipathAdjacencyHyperKmerPath( K, A, unipaths, h );
     if (REMOVE_HANGING_ENDS) RemoveHangingEnds( h, &KmerPath::KmerCount, 250, 5.0 );
     if ( MIN_COMPONENT > 0 ) h.RemoveSmallComponents(MIN_COMPONENT);
     h.RemoveDeadEdgeObjects( );
     h.RemoveEdgelessVertices( );
     h.RemoveUnneededVertices( );
     KmerBaseBrokerBig kbb( K, paths, paths_rc, pathsdb, genome );

     // Remove edges consisting entirely of Cs, which are almost certainly
     // artifacts.

     for ( int e = 0; e < h.EdgeObjectCount( ); e++ )
     {    Bool all_C = True;
          basevector b = kbb.Seq( h.EdgeObject(e) );
          for ( int i = 0; i < b.isize( ); i++ )
               if ( as_base( b[i] ) != 'C' ) all_C = False;
          if (all_C)
          {    vec<int> d;
               d.push_back(e);
               h.DeleteEdges(d);
               break;    }    }

     // Prune branches using MIN_FRAC.

     vec<double> cov;
     Coverage( h, pathsdb, cov );
     vec<int> to_delete;
     for ( int v = 0; v < h.N( ); v++ )
     {    for ( int j1 = 0; j1 < h.From(v).isize( ); j1++ )
          {    for ( int j2 = 0; j2 < h.From(v).isize( ); j2++ )
               {    int e1 = h.EdgeObjectIndexByIndexFrom( v, j1 );
                    int e2 = h.EdgeObjectIndexByIndexFrom( v, j2 );
                    double c1 = cov[e1], c2 = cov[e2];
                    if ( c1 < int( floor( MIN_FRAC * (c1+c2) ) ) )
                         to_delete.push_back(e1);    }    }
          for ( int j1 = 0; j1 < h.To(v).isize( ); j1++ )
          {    for ( int j2 = 0; j2 < h.To(v).isize( ); j2++ )
               {    int e1 = h.EdgeObjectIndexByIndexTo( v, j1 );
                    int e2 = h.EdgeObjectIndexByIndexTo( v, j2 );
                    double c1 = cov[e1], c2 = cov[e2];
                    if ( c1 < int( floor( MIN_FRAC * (c1+c2) ) ) )
                         to_delete.push_back(e1);    }    }    }
     UniqueSort(to_delete);
     h.DeleteEdges(to_delete);
     if (REMOVE_HANGING_ENDS) RemoveHangingEnds( h, &KmerPath::KmerCount, 250, 5.0 );
     if ( MIN_COMPONENT > 0 ) h.RemoveSmallComponents(MIN_COMPONENT);
     h.RemoveDeadEdgeObjects( );
     h.RemoveEdgelessVertices( );
     h.RemoveUnneededVertices( );
     h.RemoveDeadEdgeObjects( );

     // Remove edges that are highly similar to other edges.

     if ( DELTA_KEEP > 0 )
     {    while(1)
          {    vec<double> cov;
               Coverage( h, pathsdb, cov );
               vec<int> to_delete;
               for ( int v = 0; v < h.N( ); v++ )
               {    for ( int j1 = 0; j1 < h.From(v).isize( ); j1++ )
                    {    for ( int j2 = 0; j2 < h.From(v).isize( ); j2++ )
                         {    if ( j1 == j2 ) continue;
                              int e1 = h.EdgeObjectIndexByIndexFrom( v, j1 );
                              int e2 = h.EdgeObjectIndexByIndexFrom( v, j2 );
                              if ( h.From(v)[j1] != h.From(v)[j2] ) continue;
                              double c1 = cov[e1], c2 = cov[e2];
                              if ( c1 > c2 ) continue;
                              basevector b1 = kbb.Seq( h.EdgeObject(e1) );
                              basevector b2 = kbb.Seq( h.EdgeObject(e2) );
                              align a;
                              SmithWatFreeSym( b1, b2, a );
                              int errs = ActualErrors( b1, b2, a, 1, 1 );
                              if ( errs <= DELTA_KEEP )
                                   to_delete.push_back(e1);    }    }    }
               UniqueSort(to_delete);
               h.DeleteEdges(to_delete);
               if (REMOVE_HANGING_ENDS) 
                    RemoveHangingEnds( h, &KmerPath::KmerCount, 250, 5.0 );
               if ( MIN_COMPONENT > 0 ) h.RemoveSmallComponents(MIN_COMPONENT);
               h.RemoveDeadEdgeObjects( );
               h.RemoveEdgelessVertices( );
               h.RemoveUnneededVertices( );
               h.RemoveDeadEdgeObjects( );
               if ( to_delete.empty( ) ) break;    }    }

     // Generate output.

     vec<String> edge_labels_extra( h.EdgeObjectCount( ) );
     Coverage( h, pathsdb, cov );
     for ( int e = 0; e < h.EdgeObjectCount( ); e++ )
     {    const KmerPath& p = h.EdgeObject(e);
          ostringstream out;
          if ( p.KmerCount( ) < 1000 ) out << "(" << p.KmerCount( ) << " b)";
          out << "[" << setiosflags(ios::fixed) << setprecision(1) << cov[e] << "]";
          edge_labels_extra[e] = out.str( );    }
     Ofstream( fout, HEAD + ".edges.fasta" )
     for ( int i = 0; i < h.EdgeObjectCount( ); i++ )
          kbb.Seq( h.EdgeObject(i) ).Print( fout, "edge_" + ToString(i) );
     h.PrintSummaryDOT0w( out, False, False,
          True, NULL, False, &edge_labels_extra );    }
