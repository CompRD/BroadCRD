/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "Qualvector.h"
#include "String.h"
#include "tiled/CharBaseUtil.h"
#include "tiled/PaddedSeq.h"
#include "tiled/Tiling.h"
#include "tiled/TilingAnalyzer.h"
#include "tiled/TilingPrinter.h"

/**
 * tiling_printer
 * Constructor
 */
tiling_printer::tiling_printer( tiling_analyzer *analyzer ) :
  analyzer_ ( analyzer )
{
  this->ResetOptions( );
}

/**
 * tiling_printer
 * ResetOptions
 */
void tiling_printer::ResetOptions( )
{
  width_ = 150;
  to_print_begin_ = -1;
  to_print_end_ = -1;
  print_wings_ = true;
  print_quals_ = true;
  print_discrep_ = true;
}

/**
 * tiling_printer
 * SetWidth
 */
void tiling_printer::SetWidth ( int width )
{
  width_ = width;
}

/**
 * tiling_printer
 * SetWindowToPrint
 */
void tiling_printer::SetWindowToPrint( int begin, int end )
{
  to_print_begin_ = begin;
  to_print_end_ = end;
}

/**
 * tiling_printer
 * SetPrintWings
 */
void tiling_printer::SetPrintWings( bool print_wings )
{
  print_wings_ = print_wings;
}

/**
 * tiling_printer
 * SetPrintQuals
 */
void tiling_printer::SetPrintQuals( bool print_quals )
{
  print_quals_ = print_quals;
}

/**
 * tiling_printer
 * SetPrintDiscrep
 */
void tiling_printer::SetPrintDiscrep( bool print_discrep )
{
  print_discrep_ = print_discrep;
}

/**
 * tiling_printer
 * ToStream
 *
 * If sel is given, print only the specified sequences (see method
 * TamedRectangle in TilingAnalyzer.cc).
 */
void tiling_printer::ToStream( ostream &out, vec<int> *sel ) const
{
  // Fetch actual window to be printed.
  pair<int, int> actual_win;
  this->ActualWindow( actual_win );

  // Define posts (partition actual window in width_ sized chunks).
  vec<int> posts;
  
  int pos_on_master = actual_win.first;
  while ( pos_on_master < actual_win.second ) {
    posts.push_back( pos_on_master );
    pos_on_master += width_;
  }

  // Print each chunk.
  const tiling &the_tiling = analyzer_->Tiling( );

  for (int chunk_id=0; chunk_id<(int)posts.size( ); chunk_id++) {
    int post = posts[chunk_id];
    int actual_width = Min( width_, actual_win.second - post );
    
    // The rectangles of bases and quals
    vec< vec<char> > rbases;
    vec< vec<int> > rquals;
    
    vec<int> rpos;
    if ( sel )
      analyzer_->TamedRectangles( post, actual_width, *sel, rbases, rquals );
    else
      analyzer_->Rectangles( post, actual_width, rpos, rbases, rquals );
    const vec<int> &read_pos = sel ? *sel : rpos;

    // Line with discrepancies tags.
    bool tagKC =( read_pos.size( ) > 1 && read_pos[0] < 0 && read_pos[1] < 0 );
    vec<char> discrepancies = this->DiscrepTag( rbases, tagKC );
    
    // Vector of names.
    vec<String> names;
    names.reserve( read_pos.size( ) );
    for (int ii=0; ii<(int)read_pos.size( ); ii++) {
      const tile *the_tile;
      if ( read_pos[ii] == -2 )
	the_tile = &( analyzer_->Tiling( ).Known( ) );
      else if ( read_pos[ii] == -1 )
	the_tile = &( analyzer_->Tiling( ).Contig( ) );
      else
	the_tile = &( analyzer_->Tiling( ).Read( read_pos[ii] ) );

      names.push_back( the_tile->PrettyName( ) );
    }

    // Coordinates to be printed.
    out << "  " << chunk_id << "." << posts.size( ) - 1 << " ";
    this->PrintCoordinates( out, post, post + actual_width );
    
    // Print rectangle.
    if ( print_discrep_ && names.size( ) > 0 ) {
      for (int jj=0; jj<(int)names[0].size( ); jj++)
	out << " ";
      for (int jj=0; jj<(int)discrepancies.size( ); jj++)
	out << discrepancies[jj];
      out << "\n";
    }

    for (int ii=0; ii<(int)names.size( ); ii++) {
      out << names[ii];
      for (int jj=0; jj<(int)rbases[ii].size( ); jj++)
	out << rbases[ii][jj];
      out << "\n";
    }

    if ( print_quals_ ) {
      for (int ii=0; ii<(int)names.size( ); ii++) {
	out << names[ii];
	for (int jj=0; jj<(int)rquals[ii].size( ); jj++)
	  out << analyzer_->QualToChar( rquals[ii][jj] );
	out << "\n";
      }
    }

    out << "\n";
  } // next chunk.

  // Flush stream.
  out << flush;

}

/**
 * tiling_printer
 * ToStream454
 *
 * Compactify the view of a deep coverage contig.
 */
void tiling_printer::ToStream454( ostream &out ) const
{
  // Fetch actual window to be printed.
  pair<int, int> actual_win;
  this->ActualWindow( actual_win );

  // Define posts (partition actual window in width_ sized chunks).
  vec<int> posts;
  
  int pos_on_master = actual_win.first;
  while ( pos_on_master < actual_win.second ) {
    posts.push_back( pos_on_master );
    pos_on_master += width_;
  }

  // Print each chunk.
  const tiling &the_tiling = analyzer_->Tiling( );

  for (int chunk_id=0; chunk_id<(int)posts.size( ); chunk_id++) {
    int post = posts[chunk_id];
    int actual_width = Min( width_, actual_win.second - post );
    
    // The rectangles of bases and quals
    vec< vec<char> > rbases;
    vec< vec<int> > rquals;
    
    vec<int> read_pos;    
    analyzer_->Rectangles( post, actual_width, read_pos, rbases, rquals );

    for (int ii=0; ii<(int)rbases.size( ); ii++)
      for (int jj=0; jj<(int)rbases[ii].size( ); jj++)
	ToUpperCase( rbases[ii][jj] );
    
    // Select representatives.
    vec< pair<int,int> > sel2count;
    for (int ii=0; ii<(int)read_pos.size( ); ii++) {
      // Known and/or contig always selected.
      if ( read_pos[ii] < 0 ) {
	sel2count.push_back( make_pair( ii, 1 ) );
	continue;
      }
      
      // Match read with current placed reads.
      int bag_id = -1;
      const vec<char> &curr = rbases[ii];
      for (int jj=0; jj<(int)sel2count.size( ); jj++) {
	if ( read_pos[sel2count[jj].first] < 0 ) continue;
	const vec<char> &plac = rbases[sel2count[jj].first];
	bool match = true;
	for (int kk=0; kk<(int)plac.size( ); kk++) {
	  if ( plac[kk] != curr[kk] ) {
	    match = false;
	    break;
	  }
	}
	if ( match ) {
	  bag_id = jj;
	  break;
	} 
      }
      if ( bag_id < 0 ) {
	sel2count.push_back( make_pair( ii, 0 ) );
	bag_id = (int)sel2count.size( ) - 1;
      }
      sel2count[bag_id].second++;
    }
   
    // Run TamedRectanlges.
    vec< pair<int,int> > tempsel;
    for (int ii=0; ii<(int)sel2count.size( ); ii++) {
      int pos = read_pos[sel2count[ii].first];
      int count = sel2count[ii].second;
      tempsel.push_back( make_pair( pos, count ) );
    }
    swap( sel2count, tempsel );
    sort( sel2count.begin( ), sel2count.end( ) );

    read_pos.clear( );
    for (int ii=0; ii<(int)sel2count.size( ); ii++)
      read_pos.push_back( sel2count[ii].first );

    analyzer_->TamedRectangles( post, actual_width, read_pos, rbases, rquals );
    
    // All upper case.
    for (int ii=0; ii<(int)rbases.size( ); ii++)
      for (int jj=0; jj<(int)rbases[ii].size( ); jj++)
	if ( ! IsEmpty (rbases[ii][jj] ) )
	  ToUpperCase( rbases[ii][jj] );
    
    // Line with discrepancies tags.
    bool tagKC =( read_pos.size( ) > 1 && read_pos[0] < 0 && read_pos[1] < 0 );
    vec<char> discrepancies = this->DiscrepTag( rbases, tagKC );
    
    // Vector of names.
    vec<String> names;
    names.reserve( read_pos.size( ) );
    for (int ii=0; ii<(int)read_pos.size( ); ii++) {
      const tile *the_tile;
      if ( read_pos[ii] == -2 )
	the_tile = &( analyzer_->Tiling( ).Known( ) );
      else if ( read_pos[ii] == -1 )
	the_tile = &( analyzer_->Tiling( ).Contig( ) );
      else
	the_tile = &( analyzer_->Tiling( ).Read( read_pos[ii] ) );

      names.push_back( the_tile->PrettyName( ) );
    }

    // Adjust names (multiplicity).
    int old_size = names[0].size( );
    int new_size = old_size + 12;
    for (int ii=0; ii<(int)names.size( ); ii++) {
      int pos = sel2count[ii].first;
      vec<int>::iterator it = find( read_pos.begin( ), read_pos.end( ), pos );
      int sel_id = distance( read_pos.begin( ), it );
      
      if ( pos == -2 ) names[ii] += " -K-";
      else if ( pos == -1 ) names[ii] += " -C-";
      else names[ii] += " +" + ToString( sel2count[sel_id].second - 1 );

      while( (int)names[ii].size( ) < new_size ) names[ii] += " ";
    }
    
    // Coordinates to be printed.
    out << "  " << chunk_id << "." << posts.size( ) - 1 << " ";
    this->PrintCoordinates( out, post, post + actual_width );

    // Print rectangle.
    if ( print_discrep_ && names.size( ) > 0 ) {
      for (int jj=0; jj<(int)names[0].size( ); jj++)
	out << " ";
      for (int jj=0; jj<(int)discrepancies.size( ); jj++)
	out << discrepancies[jj];
      out << "\n";
    }

    for (int ii=0; ii<(int)names.size( ); ii++) {
      out << names[ii];
      for (int jj=0; jj<(int)rbases[ii].size( ); jj++)
	out << rbases[ii][jj];
      out << "\n";
    }

    if ( print_quals_ ) {
      for (int ii=0; ii<(int)names.size( ); ii++) {
	out << names[ii];
	for (int jj=0; jj<(int)rquals[ii].size( ); jj++)
	  out << analyzer_->QualToChar( rquals[ii][jj] );
	out << "\n";
      }
    }

    out << "\n";
  } // next chunk.

  // Flush stream.
  out << flush;

}

/**
 * tiling_printer
 * ActualWindow
 *
 * The actual window that will be printed. It may differ from
 * [to_print_begin_, to_print_end_] by virtue of print_wings_.
 */
void tiling_printer::ActualWindow( pair<int, int> &win ) const
{
  const tiling &the_tiling =analyzer_->Tiling( );

  int known_id = the_tiling.KnownId( );
  int contig_id = the_tiling.ContigId( );
  
  int master_len = 0;
  if ( known_id < 0 )
    master_len = analyzer_->PaddedLengthContig( );
  else
    master_len = analyzer_->PaddedLengthKnown( );

  win.first  = ( to_print_begin_ < 0 ) ? 0 : to_print_begin_;
  win.second = ( to_print_end_ < 0 ) ? master_len : to_print_end_;
  
  // Take into account (roughly) only regions for which there are reads.
  if ( the_tiling.ReadsCount( ) < 1 )
    return;
  
  pair<int, int> reads_win = analyzer_->RWin( 0 );
  for (int ii=0; ii<the_tiling.ReadsCount( ); ii++) {
    pair<int, int> r_win = analyzer_->RWin( ii );
    int aa = win.first;
    int bb = win.second;
    int cc = r_win.first;
    int dd = r_win.second;
    if ( IntervalOverlap( aa, bb, cc, dd ) ) {
      reads_win.first = Min( reads_win.first, r_win.first );
      reads_win.second = Max( reads_win.second, r_win.second );
    }
  }
  
  // Check print_wings_ only if a precise window is not specified.
  if ( to_print_begin_ < 0 && to_print_end_ < 0 ) {
    if ( print_wings_ ) {
      win.first = Min( reads_win.first, win.first );
      win.second = Max( reads_win.second, win.second );
    }
  }
  
}

/**
 * tiling_printer
 * PrintCoordinates
 */
void tiling_printer::PrintCoordinates( ostream &out, int begin, int end ) const
{
  const tiling &the_tiling = analyzer_->Tiling( );
  int known_id = the_tiling.KnownId( );
  int contig_id = the_tiling.ContigId( );

  const padded_seq *kpads = 0;
  const padded_seq *cpads = 0;
  if ( known_id > -1 ) kpads = &( the_tiling.Known( ).PaddedSeq( ) );
  if ( contig_id > -1 ) cpads = &( the_tiling.Contig( ).PaddedSeq( ) );
  
  const padded_seq *pads = kpads ? kpads : cpads;
  int unpadded_begin = pads->ToUnpaddedPos( begin );
  int unpadded_end = pads->ToUnpaddedPos( end );
  
  out << "  padded window: [" << begin << "," << end
      << ")  unpadded window: [" << unpadded_begin << "," << unpadded_end
      << ")\n";

}

/**
 * tiling_printer
 * DiscrepTag
 *
 * If tagKC is set to true, tag only events where contig and known do not
 * agree.
 */
vec<char> tiling_printer::DiscrepTag( const vec< vec<char> > &rbases,
				      bool tagKC ) const
{
  vec<char> discrepancies;

  // Tag only differences between contig and known.
  if ( tagKC ) {
    for (int ii=0; ii<(int)rbases[0].size( ); ii++) {
      bool match = ( rbases[0][ii] == rbases[1][ii] );
      char column_tag = match ? ' ' : 'X';
      discrepancies.push_back( column_tag );
    }
  }
  
  // Tag any column difference.
  if ( ! tagKC ) {
    int n_rows = rbases.size( );
    if ( print_discrep_ && n_rows > 0 ) {
      int n_cols = rbases[0].size( );
      discrepancies.reserve( n_cols );
      discrepancies.clear( );
      for (int jj=0; jj<n_cols; jj++) {
	char column_tag = ' ';
	char lead_symbol = ' ';
	for (int ii=0; ii<n_rows; ii++) {
	  if ( rbases[ii][jj] != empty_base ) {
	    lead_symbol = rbases[ii][jj];
	    break;
	  }
	}
	if ( lead_symbol == ' ' )
	  break;
	for (int ii=0; ii<n_rows; ii++) {
	  char base = rbases[ii][jj];
	  if ( base == empty_base )
	    continue;
	  if ( IsLowerCase( base ) )
	    continue;
	  if ( base != lead_symbol ) {
	    column_tag = '|';
	  }
	}
	discrepancies.push_back( column_tag );
      }
    }
  }

  // Return.
  return discrepancies;
}

