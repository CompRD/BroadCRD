// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology
//

#include "Loader.h"
#include "MainTools.h"
#include "ParseSet.h"
#include "ReadLocation.h"
#include "String.h"
#include "graphics/BarGraph.h"

/*
 * CoverageHistogram
 *
 * Generate a coverage histogram, where bin[ii]=jj means that there
 * are jj bases having coverage exactly ii. It only loads a locs
 * file (so the coverage is computed contig by contig).
 *
 * OUTDIR: full path name for output dir
 * CONTIG_IDS: if non empty, use only specified contigs
 * NBINS: number of bins to show
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( LOCS );
  CommandArgument_String( OUTDIR );
  CommandArgument_String_OrDefault( CONTIG_IDS, "" );
  CommandArgument_UnsignedInt_OrDefault( NBINS, 32 );
  EndCommandArguments;

  // File names.
  String raw_file = OUTDIR + "/raw.data";
  String histo_eps_file = OUTDIR + "/histo.eps";
  String histo_txt_file = OUTDIR + "/histo.txt";
  Mkdir777( OUTDIR );
  
  // Load.
  cout << Date( ) << ": loading" << endl;
  READ( LOCS, vec<read_location>, locs );
  sort( locs.begin( ), locs.end( ) );

  vec<int> flocs;
  LoadFirstLocs( locs, flocs );

  // Select contigs.
  int n_contigs = 1 + locs.back( ).Contig( );
  vec<Bool> select;
  if ( CONTIG_IDS == "" )
    select.resize( n_contigs, True );    
  else {
    select.resize( n_contigs, False );
    vec<int> ids;
    ParseIntSet( CONTIG_IDS, ids );
    for (int ii=0; ii<(int)ids.size( ); ii++)
      select[ids[ii]] = True;
  }
  
  // Status bar.
  int dotter = 1000;
  cout << n_contigs << " contigs found (.=" << dotter << " contigs):" << endl;

  // Generate bins.
  vec<longlong> rawbins;

  for (int contig_id=0; contig_id<n_contigs; contig_id++) {
    if ( 0 == contig_id % dotter )
      Dot( cout, contig_id / dotter );
    if ( !select[contig_id] )
      continue;
    int floc = flocs[contig_id];
    if ( floc < 0 )
      continue;
    
    int cg_len = locs[floc].LengthOfContig( );
    vec<longlong> loc_cov( cg_len, 0 );
    for (int ii=floc; ii<(int)locs.size( ); ii++) {
      if ( locs[ii].Contig( ) != contig_id )
	break;
      int beg = Max( 0, locs[ii].StartOnContig( ) );
      int end = Min( cg_len, 1 + locs[ii].StopOnContig( ) );
      for (int jj=beg; jj<end; jj++)
	loc_cov[jj] += 1;
    }
    
    for (int ii=0; ii<(int)loc_cov.size( ); ii++) {
      int cov = loc_cov[ii];
      if ( (int)rawbins.size( ) <= cov )
	rawbins.resize( cov + 1, 0 );
      rawbins[cov] += 1;
    }
  }
  cout << "\n";
  
  // Labels.
  String str_histo = "Histogram of coverage for";
  
  // Save raw data to file.
  longlong sum = 0;
  for (int ii=0; ii<(int)rawbins.size( ); ii++)
    sum += rawbins[ii];
  
  ofstream raw_out( raw_file.c_str( ) );
  raw_out << "Raw data for" << LOCS << "\n\n"
	  << "Total bases: " << sum << "\n\n";
  for (int ii=0; ii<(int)rawbins.size( ); ii++)
    raw_out << ii << "\t" << rawbins[ii] << "\n";
  raw_out.close( );
  
  // Normalize.
  cout << Date( ) << ": generating histogram" << endl;
  vec<float> bins( rawbins.size( ), 0 );
  for (int ii=0; ii<(int)bins.size( ); ii++)
    bins[ii] = 100.0 * SafeQuotient( rawbins[ii], sum );
  vec<float> cutbins = bins;
  cutbins.resize( NBINS );

  // Generate figure.
  using namespace ns_psplot;

  float low = 0;
  float high = Max( bins );
  vec<float> bot_tics( (int)NBINS, 0 );
  for (int ii=0; ii<(int)bot_tics.size( ); ii++)
    bot_tics[ii] = (float)ii;
   
  float max_bar_height = 150;
  float bar_width = 6;
  float bar_sep = 1.5;
  short font_size = 7;
  color col = gray;

  vec<ns_psplot::freetext> h_labels( 2 );
  h_labels[0] = freetext( str_histo, black, 12 );
  h_labels[1] = freetext( LOCS, black, 9 );

  bargraph h_graph = BarPlot( cutbins, low, high, bar_width, bar_sep,
			      col, font_size, max_bar_height, h_labels,
			      bot_tics, true );
  
  // Save histogram.
  cout << Date( ) << ": saving histogram" << endl;
  
  ofstream histo_eps_out( histo_eps_file.c_str( ) );
  histo_eps_out << h_graph << "\n";
  histo_eps_out.close( );
  
  ofstream histo_txt_out( histo_txt_file.c_str( ) );
  histo_txt_out << str_histo << "\n" << LOCS << "\n\n";
  for (int ii=0; ii<(int)cutbins.size( ); ii++)
    histo_txt_out << ii << "\t" << ToString( cutbins[ii], 2 ) << "%\n";
  histo_txt_out.close( );
  
  // Done.
  cout << Date( ) << ": done" << endl;
  
}
