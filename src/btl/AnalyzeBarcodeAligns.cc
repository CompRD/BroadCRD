///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "lookup/LookAlign.h"
#include "util/RunCommand.h"

/**
 * AddCounter
 */
void AddCounter( int ii, vec<int> &counters )
{
  if ( counters.isize( ) <= ii ) counters.resize( ii+1, 0 );
  counters[ii]++;
}

/**
 * TestPair
 */
bool TestPair( const int MIN_SEP,
	       const int MAX_SEP,
	       const String &barcode,
	       const vec< pair<look_align,int> > &h2c,
	       const int ii,
	       const int jj,
	       int &sep,
	       ostream &out )
{
  sep = -666;

  const look_align &al1 = h2c[ii].first;
  const look_align &al2 = h2c[jj].first;
  
  if ( al1.rc1 == al2.rc1 ) return false;
  
  const look_align &left_al = al1.rc1 ? al1 : al2;
  const look_align &right_al = al1.rc1 ? al2 : al1;
  const int left_m = al1.rc1 ? h2c[ii].second : h2c[jj].second;
  const int right_m = al1.rc1 ? h2c[jj].second : h2c[ii].second;
  
  sep = right_al.pos2( ) - left_al.Pos2( );
  if ( sep < MIN_SEP || sep > MAX_SEP ) return false;
  
  out << barcode << "\n";
  for (int ii=0; ii<2; ii++) {
    const int mult = ( ii == 0 ) ? left_m : right_m;
    const look_align &hit = ( ii == 0 ) ? left_al : right_al;

    out << " m" << mult << "\t"
	<< ( hit.rc1 ? "rc on " : "" )
	<< hit.target_id << "\t"
	<< "[" << hit.a.pos2( ) << ", " << hit.a.Pos2( ) << ")\n";
  }
  out << "\n";

  return true;  
  
}

/**
 * AnalyzeBarcodeAligns
 *
 * Load output from AlignRank4Barcodes and generate simple stats.
 *
 * Duplication rate: number of reads that are identical to each other.
 *
 * Multiplicity: after removing duplicated reads, count number of
 *               distinct reads in each barcode.
 */
int main( int argc, char* argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( BARCODES_FILE );
  CommandArgument_Int_OrDefault( SLACK, 10 );

  CommandArgument_Int_OrDefault( MIN_SEP, 1000 );
  CommandArgument_Int_OrDefault( MAX_SEP, 100000 );
  
  EndCommandArguments;

  String dup_file = BARCODES_FILE.Before( ".qlt" ) + ".dup";
  String mult_file = BARCODES_FILE.Before( ".qlt" ) + ".mult";
  String cand_file = BARCODES_FILE.Before( ".qlt" ) + ".candidates";
  String seps_file = BARCODES_FILE.Before( ".qlt" ) + ".seps";

  vec<String> needed;
  needed.push_back( BARCODES_FILE);
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;
  
  // Candidate pairs stream.
  ofstream cand_out( cand_file.c_str( ) );

  // dup: duplication rate; mult: distinct reads count per barcode.
  vec<int> dup;
  vec<int> mult;
  
  // Vector of candidate seps found.
  vec<int> seps;

  // Parse barcodes file.
  String line;
  String barcode;
  ifstream in( BARCODES_FILE );
  while ( in ) {
    
    // Goto start of block.
    bool start_found = false;
    while ( in ) {
      getline( in, barcode );
      if ( ! in ) break;
      if ( ! barcode.Contains( "r", 0 ) ) continue;
      start_found = true;
      break;
    }
    if ( ! start_found ) break;
    
    // Load aligns.
    vec< pair<look_align,int> > h2c;
    while ( in ) {
      getline( in, line );
      if ( ! in ) break;
      if ( line == "" ) break;
      ForceAssert( line.Contains( "QUERY", 0 ) );

      look_align hit;
      hit.ReadParseable( line );

      bool is_new = true;
      for (int ii=0; ii<h2c.isize( ); ii++) {
	const look_align &prev = h2c[ii].first;
	if ( ( prev.target_id == hit.target_id )
	     && ( ( Abs( prev.pos2( ) - hit.pos2( ) ) <= SLACK ) 
		  || ( Abs( prev.Pos2( ) - hit.Pos2( ) ) <= SLACK ) ) ) {
	  is_new = false;
	  break;
	}
      }

      if ( is_new ) h2c.push_back( make_pair( hit, 1 ) );
      else h2c[ h2c.size( )-1 ].second += 1;
    }
    
    // Update counters.
    AddCounter( h2c.isize( ), mult );
    for (int ii=0; ii<h2c.isize( ); ii++)
      AddCounter( h2c[ii].second, dup );
    
    // Print candidates.
    if ( h2c.size( ) > 1 ) {
      for (int ii=0; ii<h2c.isize( ); ii++) {
	for (int jj=ii+1; jj<h2c.isize( ); jj++) {
	  int sep = -666;
	  if ( TestPair(MIN_SEP,MAX_SEP,barcode,h2c,ii,jj,sep,cand_out) )
	    seps.push_back( sep );
	}
      }
    }      
    
  }

  in.close( );
  cand_out.close( );

  // Save seps.
  sort( seps.begin( ), seps.end( ) );
  WRITE( seps_file, seps );

  // Report.
  ofstream dup_out( dup_file.c_str( ) );
  dup_out << "DUPLICATION RATE TABLE\n\n";
  for (int ii=0; ii<dup.isize( ); ii++)
    dup_out << ii << "\t" << dup[ii] << "\n";
  dup_out << "\n";
  dup_out.close( );

  ofstream mult_out( mult_file.c_str( ) );
  mult_out << "MULTIPLICITY FOUND\n\n";
  for (int ii=0; ii<mult.isize( ); ii++)
    mult_out << ii << "\t" << mult[ii] << "\n";
  mult_out << endl;
  mult_out.close( );

  // Done.
  cout << Date( ) << ": AnalyzeBarcodeAligns done" << endl;

  
}
