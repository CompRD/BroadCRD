///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Basevector.h"
#include "Histogram.h"
#include "PairsManager.h"
#include "lookup/LookAlign.h"
#include "paths/reporting/ReftigUtils.h"
#include <omp.h>
// MakeDepend: library OMP

/**
 * AlignPairedReads
 *
 * Align a set of paired reads onto a given lookup reference, and
 * generate the histogram of the observed insert lengths (with
 * reasonable stretch). There will be separate output files (eps,
 * data) for each library.
 *
 * READS: base name for reads (loads <READS>.fastb and <READS>.pairto)
 * LOOKUP: full path name for lookup reference
 * OUTDIR: full path name for output directory
 * K: kmer size, used by GetAlignsFast
 * NUM_THREADS: use all threads available if 0
 * MAX_STRETCH: used to define valid inserts
 * FORCE: do not use cached aligns
 * DIGITS_TO_RIGHT: how many decimal places to show in the tics of the eps
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( READS );
  CommandArgument_String( LOOKUP );
  CommandArgument_String_OrDefault( OUTDIR, READS + ".AlignPairedReads" );
  CommandArgument_Int_OrDefault( K, 12 );
  CommandArgument_Int_OrDefault( NUM_THREADS, 0 );
  CommandArgument_Double_OrDefault( MAX_STRETCH, 10. );
  CommandArgument_Bool_OrDefault( FORCE, False );
  CommandArgument_Int_OrDefault( DIGITS_TO_RIGHT, 0 );
  EndCommandArguments;

 // Thread control.
  NUM_THREADS = configNumThreads( NUM_THREADS );
  omp_set_num_threads( NUM_THREADS );
  
 // Dir and file names.
  String fastb_file = READS + ".fastb";
  String pairs_file = READS + ".pairs";

  String aligns_file = OUTDIR + "/all_aligns.qlt";
    
  Mkpath( OUTDIR );

  // Load.
  size_t n_reads = MastervecFileObjectCount( fastb_file );

  cout << Date( ) << ": loading pairs" << endl;
  PairsManager pairs( pairs_file );

  vec<String> lib_names = pairs.getLibraryNames( );
  vec<String> data_files;
  vec<String> eps_files;
  for (int ii=0; ii<lib_names.isize( ); ii++) {
    lib_names[ii] = OUTDIR + "/" + lib_names[ii];
    data_files.push_back( lib_names[ii] + ".data" );
    eps_files.push_back( lib_names[ii] + ".eps" );
  }
  
  // Generate or load aligns.
  vec<look_align> hits;
  GetAlignsFast( K, fastb_file, LOOKUP, aligns_file, hits, !FORCE, OUTDIR );
  
  cout << Date( ) << ": indexing aligns" << endl;
  vec<int> idx( n_reads, -1 );
  for (size_t ii=0; ii<hits.size( ); ii++) {
    int rid = hits[ii].query_id;
    if ( idx[rid] == -1 ) idx[rid] = ii;
    else idx[rid] = -2;
  }

  // One histogram per library.
  vec< histogram<int> > histos( lib_names.size( ) );
  for (int lib_id=0; lib_id<lib_names.isize( ); lib_id++) {
    int sep = pairs.getLibrarySep( lib_id );
    int dev = pairs.getLibrarySD( lib_id );
    
    int low_bin = 0;
    int n_bins = 60;
    int bin_size = 0;
    float max_sep = sep + ( MAX_STRETCH * dev );
    for (int ii=1; ii<=20; ii++) {
      if ( max_sep <= ii * 3000 ) {
	bin_size = ii * 50;
	break;
      }
    }
    ForceAssertGt( bin_size, 0 );
    histos[lib_id].AddLinearProgressionOfBins( low_bin, bin_size, n_bins );
  }

  // Loop over all pairs.
  cout << Date( ) << ": generating histograms" << endl;

  vec<size_t> lib_valid( lib_names.size( ), 0 );
  for (size_t pid=0; pid<pairs.nPairs( ); pid++) {
    int id1 = pairs.ID1( pid );
    int id2 = pairs.ID2( pid );
    
    // Not uniquely aligned.
    int idx1 = idx[id1];
    int idx2 = idx[id2];
    if ( idx1 < 0 || idx2 < 0 ) continue;

    // Not a logical pair.
    if ( ! IsLogicalPair( hits[idx1], hits[idx2] ) ) continue;

    // Pair is too stretched.
    float given_sep = pairs.sep( pid );
    float given_dev = pairs.sd( pid );
    float obs_sep = ObservedSeparation( hits[idx1], hits[idx2] );
    float stretch = ( obs_sep - given_sep ) / given_dev;
    if ( obs_sep < -MAX_STRETCH || stretch > MAX_STRETCH ) continue;
    
    // Insert size is negative.
    int len1 = hits[idx1].query_length;
    int len2 = hits[idx2].query_length;
    int insert_size = obs_sep + len1 + len2;
    if ( insert_size < 0 ) continue;

    // Ok, add insert size.
    int lib_id = pairs.libraryID( pid );
    histos[lib_id].AddDatum( insert_size );
    lib_valid[lib_id]++;
  }
  
  // Generating output files.
  cout << Date( ) << ": generating output files" << endl;

  vec<size_t> lib_sizes = pairs.getLibrarySizes( );
  for (int lib_id=0; lib_id<lib_names.isize( ); lib_id++) {
    histogram<int> &histo = histos[lib_id];
    
    // Data output.
    ofstream out( data_files[lib_id].c_str( ) );
    histo.PrintAsColumns( out, true, true );
    out.close( );

    // Eps output.
    using namespace ns_psplot;

    vec<freetext> labels;
    String title1 = "Histogram of insert sizes for";
    String title2 = lib_names[lib_id] ;
    String title3 = ToString( 2 * lib_sizes[lib_id]) + " total reads in input";
    String title4 = ToString( 2 * lib_valid[lib_id] ) + " reads in valid pairs";
    if ( lib_sizes[lib_id] > 0 ) {
      float ratio = float( lib_valid[lib_id] ) / float( lib_sizes[lib_id]);
      title4 += " (" + ToString( 100.0 * ratio, 1 ) + "% of the total)";
    }

    labels.push_back( freetext( title1, black, 12 ) );
    labels.push_back( freetext( title2, black, 12 ) );
    labels.push_back( freetext( title3, black, 10 ) );
    labels.push_back( freetext( title4, black, 10 ) );
    
    ofstream eps_out( eps_files[lib_id].c_str( ) );
    histo.PrintAsEps( eps_out, labels, DIGITS_TO_RIGHT );
    eps_out.close( );
  }
  
  // Done.
  cout << Date( ) << ": AlignPairedReads done" << endl;

}

