// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
//

#include <strstream>

#include "Basevector.h"
#include "FastIfstream.h"
#include "Qualvector.h"
#include "String.h"
#include "Vec.h"
#include "tiled/PaddedSeq.h"
#include "tiled/ReadsTiling.h"
#include "tiled/TilesContigParser.h"



/*
 * tc_parser
 * Reserve
 */
void tc_parser::Reserve( int n_reads, longlong n_bases )
{
  bases_.clear( );
  quals_.clear( );
  pads_.clear( );

  bases_.Reserve( n_bases / 16 + n_reads, n_reads );
  quals_.Reserve( n_bases, n_reads );
  pads_.reserve( n_reads );
}



/*
 * tc_parser
 * Append
 */
void tc_parser::Append( const String &in_fasta, const String &in_qual )
{
  // A non empty fasta file contains at least two lines.
  if ( LineCount( in_fasta ) < 2 ) {
    ForceAssert( LineCount( in_qual ) < 2 );
    return;
  }

  // Fetch bases as vec of char and quals as vec of int.
  fast_ifstream in_bases( in_fasta.c_str( ) );
  fast_ifstream in_quals( in_qual.c_str( ) );

  char a_char;
  int a_qual;
  basevector a_fastb;
  qualvector a_qualb;
  String read_name;
  String a_line;

  in_bases.peek( a_char );
  ForceAssert( a_char == '>' );

  in_quals.peek( a_char );
  ForceAssert( a_char == '>' );

  vec<char> read_bases;
  vec<int> read_quals;

  while ( 1 ) {

    // End of files.
    if ( in_bases.fail( ) ) {
      ForceAssert( in_quals.fail( ) );
      break;
    }

    getline( in_bases, read_name );
    if ( in_bases.fail( ) ) {
      in_quals.get( a_char );
      ForceAssert( in_quals.fail( ) );
      break;
    }

    ForceAssert( read_name.Contains( ">", 0 ) );

    // Load bases first.
    read_bases.resize( 0 );

    while( 1 ) {
      in_bases.peek( a_char );
      if ( in_bases.fail( ) || a_char == '>' )
	break;
      getline( in_bases, a_line );
      istrstream inbline( a_line.c_str( ) );
      while ( 1 ) {
	inbline >> a_char;
	if ( !inbline )
	  break;
	read_bases.push_back( a_char );
      }
    }

    // Then load quals.
    read_quals.resize( 0 );
    read_quals.reserve( read_bases.size( ) );

    getline( in_quals, a_line );
    ForceAssert( a_line.Contains( ">", 0 ) );

    while ( 1 ) {
      in_quals.peek( a_char );
      if ( in_quals.fail( ) || a_char == '>' )
	break;
      getline( in_quals, a_line );
      istrstream inqline( a_line.c_str( ) );
      while ( 1 ) {
	inqline >> a_qual;
	if ( !inqline )
	  break;
	read_quals.push_back( a_qual );
      }
    }

    // Check, match and append to vecbasevector and vecqualvector.
    int n_valid_bases = 0;
    ForceAssert( read_bases.size( ) == read_quals.size( ) );
    for (int ii=0; ii<(int)read_bases.size( ); ii++) {
      if ( read_bases[ii] != gap_base ) {
	ForceAssert( read_quals[ii] != gap_qual );
	n_valid_bases++;
      }
    }

    a_fastb.Setsize( n_valid_bases );
    a_qualb.resize( n_valid_bases );

    int pos = 0;
    for (int ii=0; ii<(int)read_bases.size( ); ii++) {
        char base = read_bases[ii];
        if ( Base::isCanonicalBase(base) )
        {
            a_fastb.Set(pos,Base::char2Val(base));
            a_qualb[pos++] = read_quals[ii] > -1 ? read_quals[ii] : 0;
        }
    }

    bases_.push_back( a_fastb );
    quals_.push_back( a_qualb );

    // Create a padded_seq with the contig's pads.
    padded_seq read_pads( 0, 0, a_fastb.size( ) );
    for (int ii=0; ii<(int)read_bases.size( ); ii++)
      if ( read_bases[ii] == gap_base )
	read_pads.AddPad( ii );

    pads_.push_back( read_pads );
  }

}



/*
 * tc_parser
 * Write
 */
void tc_parser::Write( const String &out_fastb,
		       const String &out_qualb,
		       const String &out_pads ) const
{
  bases_.WriteAll( out_fastb );

  quals_.WriteAll( out_qualb );

  ofstream pout( out_pads.c_str( ) );
  pout << pads_.size( ) << "\n";
  for (int ii=0; ii<(int)pads_.size( ); ii++)
    pout << pads_[ii] << "\n";
}



