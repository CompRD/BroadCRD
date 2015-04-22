///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "Fastavector.h"
#include "Superb.h"
#include "Vec.h"
#include "math/Functions.h"
#include "paths/AssemblyIntegrity.h"

/**
 * AssemblyIntegrity
 */
int AssemblyIntegrity( ostream *log,
		       const vec<superb> *superbs,
		       const vec<fastavector> *fastas,
		       const vecbasevector *fastbs,
		       const vecbasevector *efastas )
{
  ofstream devnull ( "/dev/null" );
  ostream &out = log ? *log : devnull;
  
  // Total errors count.
  int n_errors = 0;

  // Check fasta and fastb sizes match.
  if ( fastas && fastbs ) {
    out << "Testing fasta sizes vs fastb sizes" << endl;
    if ( fastas->size( ) != fastbs->size( ) ) {
      out << " ERROR: sizes of fasta (" << fastas->size( )
	  << " entries), and fastb ("  << fastbs->size( )
	  << " entries) do not match\n";
      n_errors++;
    }
    else out << fastas->size( ) << " entries found\n";

    vec< vec<String> > table;
    vec<String> line;
    {
      line.push_back( "  contig_id" );
      line.push_back( "fasta_length" );
      line.push_back( "fastb_length" );
      table.push_back( line );
    }
    for (size_t ii=0; ii<Min( fastas->size( ), fastbs->size( ) ); ii++) {
      if ( (*fastas)[ii].size( ) == (*fastbs)[ii].size( ) ) continue;
      line.clear( );
      line.push_back( ToString( ii ) );
      line.push_back( ToString( (*fastas)[ii].size( ) ) );
      line.push_back( ToString( (*fastbs)[ii].size( ) ) );
      table.push_back( line );
    }
    if ( table.size( ) > 1 ) {
      out << " ERROR: some fasta/fastb entries have different length:\n";
      PrintTabular( out, table, 3, "rrr" );
      n_errors += table.isize( ) - 1;
    }

    out << endl;

  } // check fasta and fastb sizes
  
  // Check superb lengths macth fasta (or fastb).
  if ( superbs && ( fastas || fastbs ) ) {
    out << "Testing contig lengths match lengths in the superb\n"
	<< superbs->size( ) << " supers found" << endl;

    vec<int> clens;
    if ( fastbs ) {
      clens.reserve( fastbs->size( ) );
      for (size_t ii=0; ii<fastbs->size( ); ii++)
	clens.push_back( int( (*fastbs)[ii].size( ) ) );
    }
    else {
      clens.reserve( fastas->size( ) );
      for (int ii=0; ii<fastas->isize( ); ii++)
	clens.push_back( int( (*fastas)[ii].size( ) ) );
    }
    
    vec<String> missing;
    for (int ii=0; ii<superbs->isize( ); ii++)
      for (int jj=0; jj<(*superbs)[ii].Ntigs( ); jj++)
	if ( (*superbs)[ii].Tig( jj ) >= clens.isize( ) )
	  missing.push_back( "  " + ToString( (*superbs)[ii].Tig( jj ) )
			     + " = s" + ToString( ii )
			     + "." + ToString( jj )
			     + "/" + ToString( (*superbs)[ii].Ntigs( ) - 1 )
			     + " (" + ToString( (*superbs)[ii].Len( jj ) )
			     + " bp)" );
    if ( missing.size( ) > 0 ) {
      out << " ERROR: some contig ids in the superb are not found:\n";
      for (int ii=0; ii<missing.isize( ); ii++)
	out << missing[ii] << "\n";
      n_errors += missing.isize( );
    }

    vec< vec<String> > table;
    vec<String> line;
    {
      line.push_back( "  contig_id" );
      line.push_back( "super_pos" );
      line.push_back( fastbs ? "fastb_length" : "fasta_length" );
      line.push_back( "super_length" );
      table.push_back( line );
    }
    for (int ii=0; ii<superbs->isize( ); ii++) {
      String nsup = ToString( (*superbs)[ii].Ntigs( ) );
      for (int jj=0; jj<(*superbs)[ii].Ntigs( ); jj++) {
	int c_id = (*superbs)[ii].Tig( jj );
	int b_len = (*superbs)[ii].Len( jj );
	if ( c_id < clens.isize( ) && b_len == clens[c_id] ) continue;
	line.clear( );
	line.push_back( ToString( c_id ) );
	line.push_back( ToString( ii ) + "." + ToString( jj ) + "/" + nsup );
	line.push_back( c_id < clens.isize() ? ToString( clens[c_id] ) : "na" );
	line.push_back( ToString( b_len ) );
	table.push_back( line );
      }
    }
    if ( table.size( ) > 1 ) {
      out << "ERROR: mismatches found between fast* and super lengths:\n";
      PrintTabular( out, table, 3, "rrrr" );
      n_errors += table.isize( ) - 1;
    }
    
    out << endl;

  }

  // Check fastbs and efastas match base by base.
  if ( fastbs && efastas ) {
    out << "Base by base testing fastbs vs efastas" << endl;
    if ( efastas->size( ) != fastbs->size( ) ) {
      out << " ERROR: sizes of efasta (" << efastas->size( )
	  << " entries), and fastb ("  << fastbs->size( )
	  << " entries) do not match\n";
      n_errors++;
    }
    else out << efastas->size( ) << " entries found\n";
    
    vec< vec<String> > table;
    vec<String> line;
    {
      line.push_back( "  contig_id" );
      line.push_back( "efasta_length" );
      line.push_back( "fastb_length" );
      table.push_back( line );
    }
    for (size_t ii=0; ii<Min( efastas->size( ), fastbs->size( ) ); ii++) {
      if ( (*efastas)[ii] == (*fastbs)[ii] ) continue;
      line.clear( );
      line.push_back( ToString( ii ) );
      line.push_back( ToString( (*efastas)[ii].size( ) ) );
      line.push_back( ToString( (*fastbs)[ii].size( ) ) );
      table.push_back( line );
    }
    if ( table.size( ) > 1 ) {
      out << " ERROR: some efasta/fastb entries are different:\n";
      PrintTabular( out, table, 3, "rrr" );
      n_errors += table.isize( ) - 1;
    }

    out << endl;
  }

  // Done.
  out << n_errors << " total error(s) found\n" << endl;
  return n_errors;

}

