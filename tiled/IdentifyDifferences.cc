// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology
//

#include "tiled/CKDiff.h"
#include "tiled/IdentifyDifferences.h"
#include "tiled/Tiling.h"
#include "tiled/TilingPrinter.h"

typedef pair<int,int> win;



bool Confirms( const vec<char> &master,
	       const vec<char> &read_b,
	       const vec<int> &read_q )
{
  int min_acceptable_qual_mean = 15;

  bool they_match = true;
  for (int ii=0; ii<(int)master.size( ); ii++) {
    if ( master[ii] != read_b[ii] ) {
      they_match = false;
      break;
    }
  }
  
  if ( they_match ) {
    float q_mean = Mean( read_q );
    if ( q_mean < min_acceptable_qual_mean )
      they_match = false;
  }
  
  return they_match;
}



ck_diff IdentifyDifferences( int radius,
			     int tile_id,
			     tiling_analyzer &analyzer,
			     ostream &out )
{
  const bool verbose = true;

  const tiling &tile_set = analyzer.Tiling( );
  const padded_seq &cpads = tile_set.Contig( ).PaddedSeq( );
  const padded_seq &kpads = tile_set.Known( ).PaddedSeq( );
  int w_begin = cpads.Begin( ) + cpads.AlBegin( );
  int w_end = w_begin + cpads.PaddedLength( );
  int w_len = w_end - w_begin;
  int n_reads = tile_set.ReadsCount( );
  int contig_id = tile_set.ContigId( );
  int known_id = tile_set.KnownId( );
  
  ForceAssert( contig_id > -1 && known_id > -1 );
  int width = 150;
  
  vec<int> posts;
  int pos_on_master = w_begin;
  while( pos_on_master < w_end ) {
    posts.push_back( pos_on_master );
    pos_on_master += width;
  }

  // Basic statistics. Notice in general diff.n_total <= w_len.
  ck_diff diff;
  diff.n_total = 0;

  // Bases for which consensus != known are called bad.
  vec<int> bad;          // list of bad bases
  vec<win> bad_win;      // regions containing bad bases.
  vec<Bool> is_bad;      // has size w_len (boolean if base is bad or not)
  
  // Loop over all chunks (detect differences).
  bad.reserve( w_len );
  for (int chunk_id=0; chunk_id<(int)posts.size( ); chunk_id++) {
    int post = posts[chunk_id];
    
    vec<int> read_pos;
    read_pos.push_back( -2 );
    read_pos.push_back( -1 );

    vec< vec<char> > rbases;
    vec< vec<int> > rquals;
    
    analyzer.TamedRectangles( post, width, read_pos, rbases, rquals );
    ForceAssert( rbases.size( ) > 0 );

    for (int jj=0; jj<(int)rbases[0].size( ); jj++) {
      if ( post + jj > w_end - 1 )
	break;
      if ( rbases.size( ) < 2 )
	continue;
      if ( rbases[0][jj] == empty_base || rbases[0][jj] == empty_base )
	continue;
      diff.n_total += 1;
      if ( rbases[0][jj] != rbases[1][jj] ) {
	bad.push_back( post + jj - w_begin );
	if ( rbases[0][jj] == gap_base )
	  diff.n_pads_k++;
	if ( rbases[1][jj] == gap_base )
	  diff.n_pads_c++;
      }
    }
  }

  if ( bad.size( ) < 1 )
    return diff;
  
  // Build is_bad.
  is_bad.resize( w_len, False );
  for (int ii=0; ii<(int)bad.size( ); ii++)
    is_bad[bad[ii]] = True;

  // Build bad_win.
  int first_pos = 0;
  int pos = 0;
  while ( pos < (int)bad.size( ) ) {
    while ( pos + 1 < (int)bad.size( ) && bad[pos+1] - bad[pos] <= radius )
      pos++;
    int last_pos = pos;
    int wbegin = Max( 0, bad[first_pos] - radius );
    int wend = Min( w_len, bad[last_pos] + radius + 1 );
    bad_win.push_back( make_pair( wbegin, wend ) );
    
    pos++;
    first_pos = pos;
  }

  // Analyze bad intervals.
  for (int ii=0; ii<(int)bad_win.size( ); ii++) {
    int print_beg = bad_win[ii].first;
    int print_end = bad_win[ii].second;
    int post = print_beg;
    int width = print_end - print_beg;
    int n_bad = 0;
    for (int jj=print_beg; jj<print_end; jj++)
      if ( is_bad[jj] )
	n_bad++;

    out << "tile_" << tile_id
	<< "  known_" << known_id
	<< "  contig_" << contig_id
	<< "  [" << print_beg
	<< ", " << print_end
	<< ")  n_diff=" << n_bad;
    
    vec<int> read_pos;
    vec< vec<char> > rbases;
    vec< vec<int> > rquals;
    analyzer.Rectangles( post, width, read_pos, rbases, rquals );
    
    // Check if reads confirm known and/or contig.
    bool confirms_known = false;
    bool confirms_contig = false;
    for (int jj=2; jj<(int)rbases.size( ); jj++) {
      const vec<char> &read_b = rbases[jj];
      const vec<int> &read_q = rquals[jj];

      if ( !confirms_known ) {
	const vec<char> &known_b = rbases[0];
	if ( Confirms( known_b, read_b, read_q ) )
	  confirms_known = true;
      }

      if ( !confirms_contig ) {
	const vec<char> &contig_b = rbases[1];
	if ( Confirms( contig_b, read_b, read_q ) )
	  confirms_contig = true;
      }

      if ( confirms_known && confirms_contig )
	break;
    }

    if ( confirms_known && confirms_contig ) {
      diff.n_confirm_both += n_bad;
      out << "  confirm_c, and confirm_k";
    }
    else if ( confirms_known ) {
      out << "  confirm_k";
      diff.n_confirm_k += n_bad;
    }
    else if ( confirms_contig ) {
      out << "  confirm_c";
      diff.n_confirm_c += n_bad;
    }
    else {
      out << "  neither c nor k confirmed";
      diff.n_confirm_none += n_bad;
    }
    out << "\n";

    if ( verbose && !confirms_contig ) {
      tiling_printer printer ( &analyzer );
      printer.SetPrintWings( false );
      printer.SetWindowToPrint( print_beg, print_end );
      printer.ToStream( out );
    }

  }

  // One line summary for this tile.
  out << "Summary for tile_" << tile_id << ":  ";
  diff.PrintCompact( out );
  out << "\n\n";

#ifdef BAD_WINS_HISTO
  vec<int> lens;
  for (int ii=0; ii<(int)bad_win.size( ); ii++) {
    int raw_len = bad_win[ii].second - bad_win[ii].first;
    if ( raw_len >= (int)lens.size( ) )
      lens.resize( raw_len + 1, 0 );
    lens[raw_len] += 1;
  }

  for (int ii=0; ii<(int)lens.size( ); ii++)
    cout << ii << "\t" << lens[ii] << "\n";
#endif
  
  return diff;
  
}



