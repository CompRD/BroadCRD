// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology
//

#include "Basevector.h"
#include "MainTools.h"
#include "LocsHandler.h"
#include "SupersHandler.h"
#include "graphics/BarGraph.h"
#include "graphics/BuildHistogram.h"
#include "reporting/SnapToGrid.h"

/**
 * CoverageDistribution
 *
 * Find the distribution of coverages for contigs or supers, using a naive
 * algorithm: for every contig (super) determine the assembly coverage, and
 * assign that coverage to all the bases in the contig (super). Plot the
 * histogram.
 *
 * SUPERS: find distribution for supers (or contigs, if False)
 * INCLUDE_GAPS: if SUPERS=True, use gapped super length as the denominator
 * MAX_COV: max coverage to plot
 * INCREMENT: for histogram bins, as a fraction of coverage
 * BASE_OUT: base for output (BASE_OUT.txt, BASE_OUT.eps in SUBDIR)
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String_OrDefault( SUBDIR, "" );
  CommandArgument_Bool_OrDefault( SUPERS, True );
  CommandArgument_Bool_OrDefault( INCLUDE_GAPS, True );
  CommandArgument_Double_OrDefault( MAX_COV, 12 );
  CommandArgument_Double_OrDefault( INCREMENT, 0.05 );
  CommandArgument_String_OrDefault( BASE_OUT, "CoverageDistribution" );
  EndCommandArguments;
  
  // Dir names.
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/" + SUBDIR;

  // File names.
  String reads_bases_file = run_dir + "/reads.fastb";

  String contigs_bases_file = sub_dir + "/mergedcontigs.fastb";
  String supers_file = sub_dir + "/mergedcontigs.superb";
  String locs_file = sub_dir + "/mergedcontigs_orig.locs";

  String out_txt_file = sub_dir + "/" + BASE_OUT + ".txt";
  String out_eps_file = sub_dir + "/" + BASE_OUT + ".eps";

  // Load.
  int n_reads = MastervecFileObjectCount( reads_bases_file );
  int n_contigs = MastervecFileObjectCount( contigs_bases_file );

  cout << Date( ) << ": loading locs" << endl;
  lhandler locs( n_reads, n_contigs, locs_file );

  cout << Date( ) << ": loading supers" << endl;
  shandler supers( n_contigs, supers_file );

  // Histogram bins.
  vec<float> x_bars( 1, 0 );
  while ( x_bars.back( ) < MAX_COV )
    x_bars.push_back( (float)x_bars.size( ) * (float)INCREMENT );
  vec<float> y_bars( x_bars.size( ), 0 );
  ForceAssert( x_bars.size( ) > 0 );

  // Calculate coverages.
  cout << Date( ) << ": generating coverages" << endl;

  longlong tot_bases = 0;
  if ( SUPERS ) {
    for (int super_id=0; super_id<supers.Size( ); super_id++) {
      const superb &sup = supers[super_id];
      int super_len = INCLUDE_GAPS ? sup.TrueLength( ) : sup.ReducedLength( );
      if ( super_len < 1 )
	continue;
      tot_bases += super_len;
      int n_bases = 0;
      for (int cgpos=0; cgpos<sup.Ntigs( ); cgpos++) {
	int contig_id = sup.Tig( cgpos );
	int floc = locs.FirstLoc( contig_id );
	for (int locpos=floc; locpos<locs.Size( ); locpos++) {
	  if ( locs[locpos].Contig( ) != contig_id )
	    break;
	  n_bases += locs[locpos].LengthOfRead( );
	}
      }
      float cov = SafeQuotient( n_bases, super_len );
      int pos = PlaceInBin( cov, x_bars );
      y_bars[pos] += (float)n_bases;
    }
  }
  else {
    for (int contig_id=0; contig_id<n_contigs; contig_id++) {
      int floc = locs.FirstLoc( contig_id );
      int contig_len = locs[floc].LengthOfContig( );
      if ( contig_len < 1 )
	continue;
      tot_bases += contig_len;
      int n_bases = 0;
      for (int locpos=floc; locpos<locs.Size( ); locpos++) {
	if ( locs[locpos].Contig( ) != contig_id )
	  break;
	n_bases += locs[locpos].LengthOfRead( );
      }
      float cov = SafeQuotient( n_bases, contig_len );
      int pos = PlaceInBin( cov, x_bars );
      y_bars[pos] += (float)n_bases;
    }
  }
  
  // Normalize.
  for (int ii=0; ii<(int)y_bars.size( ); ii++)
    y_bars[ii] = 100.0 * SafeQuotient( (longlong)y_bars[ii], tot_bases );
  
  // Save txt.
  String title = ( SUPERS ? "Super " : "Contig " );
  title +=  "Coverage Histogram";

  ofstream out_txt( out_txt_file.c_str( ) );
  out_txt << title << "\n" << sub_dir << "\n\n";
  for (int ii=0; ii<(int)x_bars.size( ); ii++)
    out_txt << x_bars[ii] << "\t" << y_bars[ii] << "\n";
  out_txt.close( );

  // Generate and save eps.
  using namespace ns_psplot;

  float max_bar_height = 150;
  float bar_width = 6;
  float bar_sep = 1.5;
  color col_label = black;
  color col_bars = gray;
  short font_size = 7;

  vec<freetext> labels( 2 );
  labels[0] = freetext( title, col_label, 12 );
  labels[1] = freetext( sub_dir, col_label, 10 );

  float low = *min_element( x_bars.begin( ), x_bars.end( ) );
  float high = *max_element( x_bars.begin( ), x_bars.end( ) );
  vec<float> bot_tics;

  bargraph graph = BarPlot( y_bars, low, high, bar_width, bar_sep,
			    col_bars, font_size, max_bar_height, labels,
			    bot_tics, true);
  
  ofstream out_eps( out_eps_file.c_str() );
  out_eps << graph << "\n";
  out_eps.close( );
  
  // Done.
  cout << Date( ) << ": done" << endl;
  
}
