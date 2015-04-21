///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Basevector.h"
#include "Qualvector.h"

/**
 * FishMonomersFromReads
 *
 * Find putative monomers by looking at reads sequenced from
 * constructs, where a construct is a sequence of primers and
 * monomers, as in:
 *
 *    -->[============]
 *
 * the --> part is the primer (usually 4 bases long), and the
 * [============] part is the actual monomer (MONOMER_SIZE bases
 * long).
 */
int main( int argc, char *argv[] )
{
  RunTime( );
 
  BeginCommandArguments;
  CommandArgument_String( READS );
  CommandArgument_String( OUT_DIR );
  CommandArgument_String_OrDefault( BAIT, "CCAGAGCAGGTCGTGGCA" );
  CommandArgument_Int_OrDefault( MONOMER_SIZE, 98 );
  CommandArgument_Int_OrDefault( BAIT_START, 3 );
  CommandArgument_Int_OrDefault( MIN_TO_REPORT, 4 );
  EndCommandArguments;

  // Dir and file names.
  String reads_file = READS + ".fastb";

  Mkpath( OUT_DIR );

  // Clean from old runs.
  vec<String> all_files = AllFiles( OUT_DIR );
  for (int ii=0; ii<all_files.isize( ); ii++) {
    String file = OUT_DIR + "/" + all_files[ii];
    if ( IsRegularFile( file ) ) Remove( file );
  }

  // Load.
  bvec bait( BAIT );
  const int K = bait.size( );
  ForceAssertLt( K, MONOMER_SIZE );

  vecbvec reads( reads_file.c_str( ) );

  // Find matches in orig reads (and sort them).
  vec<String> matches;
  for (size_t ii=0; ii<reads.size( ); ii++) {
    bvec rc_read;
    if ( ii%2 == 1 ) {
      rc_read = reads[ii];
      rc_read.ReverseComplement( );
    }
    const bvec &read = ( ii%2==1 ) ? rc_read : reads[ii];
    const int read_len = read.size( );
    for (int pos=BAIT_START; pos<read_len-MONOMER_SIZE+BAIT_START; pos++) {
      bvec chunk( read, pos, K );
      if ( chunk == bait ) {
	int start = pos - BAIT_START;
	bvec candidate( read, start, MONOMER_SIZE );
	matches.push_back( candidate.ToString( ) );
      }
    }
  }
  sort( matches.begin( ), matches.end( ) );

  // Each pairs contains multiplicity, sequence of monomer (core part only).
  vec< pair<int,String> > n2mono;
  for (int ii=0; ii<matches.isize( ); ii++) {
    if ( ii == 0 || matches[ii] != n2mono.back( ).second ) 
      n2mono.push_back( make_pair( 1, matches[ii] ) );
    else
      n2mono[n2mono.size( )-1].first += 1;
  }
  sort( n2mono.rbegin( ), n2mono.rend( ) );
  
  // Remove putatives with low frequency.
  {
    vec< pair<int,String> > temp;
    for (int ii=0; ii<n2mono.isize( ); ii++)
      if ( n2mono[ii].first >= MIN_TO_REPORT )
	temp.push_back( n2mono[ii] );
    swap( temp, n2mono );
  }
  if ( n2mono.size( ) < 1 ) {
    cout << "No monomer candidate detected.\n" << endl;
    return 0;
  }
  
  // Generate line to visualize differences.
  String diff;
  for (int jj=0; jj<n2mono[0].second.isize( ); jj++) {
    vec<char> col;
    for (int ii=0; ii<n2mono.isize( ); ii++)
      col.push_back( n2mono[ii].second[jj] );
    sort( col.begin( ), col.end( ) );
    col.erase( unique( col.begin( ), col.end( ) ), col.end( ) );

    diff += col.size( ) == 1 ? ' ' : '*';
  }
  cout << diff << "\n";

  // Print monomers, and save as fastas.
  int ntot = 0;
  for (int ii=0; ii<n2mono.isize( ); ii++) {
    ntot += n2mono[ii].first;
    cout << n2mono[ii].second << "\t" << n2mono[ii].first << "\n";

    String name = "monomer_" + ToString( ii );
    String out_file = OUT_DIR + "/" + name + ".fasta";
    ofstream out( out_file.c_str( ) );
    bvec bmono( n2mono[ii].second );
    bmono.Print( out, name );
    out.close( );
  }
  cout << endl;
  
}

