///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "String.h"
#include "Vec.h"
#include "btl/MonomerUtilsTALENs.h"
#include "lookup/LookAlign.h"

/**
 * ShortName
 */
String ShortName( const String &name )
{
  if ( name.Contains( "monomer" ) ) return name.After( "monomer" );
  if ( name.Contains( "termin" ) ) return name.After( "termin" ) + ".t";
  return name;
}

/**
 * MonomersOverlap
 */
bool MonomersOverlap( const look_align_plus &al1,
		      const look_align_plus &al2 )
{
  if ( al1.target_id != al2.target_id ) return false;

  int begin = Max( al1.pos2( ), al2.pos2( ) );
  int end = Min( al1.Pos2( ), al2.Pos2( ) );
  return ( end - begin > 0 );
}

/**
 * AlignsChains
 */
void AlignsToChains( const vec<look_align_plus> &all_hits,
		     const vec<bool> &is_termin,
		     vec<look_align_plus> &chains )
{
  chains.clear( );

  // Minimal filtering.
  vec<look_align_plus> hits;
  hits.reserve( all_hits.size( ) );
  for (size_t ii=0; ii<all_hits.size( ); ii++) {
    const look_align_plus &al = all_hits[ii];
    const int tid = al.target_id;
    
    // Monomer/terminator not fully embedded.
    if ( al.pos1( ) != 0 || al.Pos1( ) != (int)al.query_length ) continue;
    
    // Only accept fw hits on read 0 and 2, rc hits on read 1.
    if ( ( al.Fw1( ) && tid != 1 ) || ( al.Rc1( ) && tid == 1 ) )
      hits.push_back( al );
  }
  
  // Sort aligns by error count.
  vec< pair<int,int> > err2ids;
  err2ids.reserve( hits.size( ) );
  for (int ii=0; ii<hits.isize( ); ii++)
    err2ids.push_back( make_pair( hits[ii].Errors( ), ii ) );
  sort( err2ids.begin( ), err2ids.end( ) );
  
  // Generate a list of non-overlapping monomers (first) / terminators (later).
  chains.reserve( hits.size( ) );
  for (int pass=0; pass<2; pass++) {
    for (int ii=0; ii<err2ids.isize( ); ii++) {
      int hit_id = err2ids[ii].second;
      const look_align_plus &hit = hits[hit_id];
      const int mono_id = hit.query_id;
      if ( pass == 0 && is_termin[mono_id] ) continue;
      if ( pass == 1 && ! is_termin[mono_id] ) continue;

      bool overlaps = false;
      for (int jj=0; jj<chains.isize( ); jj++) {
	if ( MonomersOverlap( chains[jj], hit ) ) {
	  overlaps = true;
	  break;
	}
      }
      if ( ! overlaps ) chains.push_back( hit );
    }
  }

  // Sort chains.
  vec< pair<int,int> > starts;
  starts.reserve( chains.size( ) );
  for (int ii=0; ii<chains.isize( ); ii++)
    starts.push_back( make_pair( chains[ii].target_id, chains[ii].pos2( ) ) );
  SortSync( starts, chains );
  
}

/**
 * PrintChains
 */
void PrintChains( const vec<look_align_plus> &chains,
		  const vec<bool> &is_termin,
		  const vecString &parts_ids,
		  ostream &log )
{
  // HEURISTICS.
  int max_errors = 2;

  // Build printable tables (fw reads first, then rc).
  vec< vec<String> > table;

  table.push_back( MkVec( String( "  read" ),
			  String( "range" ),
			  String( "name" ),
			  String( "errors" ) ) );

  String ws = "";
  vec<String> empty_line = MkVec( ws, ws, ws, ws );
  table.push_back( empty_line );

  bool found2 = false;
  for (int ii=0; ii<chains.isize( ); ii++) {
    const look_align_plus &hit = chains[ii];
    if ( hit.target_id == 2 ) found2 = true;
    if ( hit.target_id != 0  ) continue;

    String str_p2 = ToString( hit.pos2( ) );
    String str_P2 = ToString( hit.Pos2( ) );
    table.push_back( MkVec( String( "fw1" ),
			    String( "[" + str_p2 + ", " + str_P2 + ")" ),
			    ShortName( parts_ids[hit.query_id] ),
			    ToString( hit.Errors( ) ) ) );
  }
  table.push_back( empty_line );
  
  if ( found2 ) {
    for (int ii=0; ii<chains.isize( ); ii++) {
      const look_align_plus &hit = chains[ii];
      if ( hit.target_id != 2  ) continue;
      
      String str_p2 = ToString( hit.pos2( ) );
      String str_P2 = ToString( hit.Pos2( ) );
      table.push_back( MkVec( String( "fw2" ),
			      String( "[" + str_p2 + ", " + str_P2 + ")" ),
			      ShortName( parts_ids[hit.query_id] ),
			      ToString( hit.Errors( ) ) ) );
    }
    table.push_back( empty_line );
  }
  
  for (int ii=chains.isize( )-1; ii>=0; ii--) {
    const look_align_plus &hit = chains[ii];
    if ( hit.target_id != 1 ) continue;
    
    String str_p2 = ToString( hit.pos2( ) );
    String str_P2 = ToString( hit.Pos2( ) );
    table.push_back( MkVec( String( "rc" ),
			    String( "[" + str_P2 + ", " + str_p2 + ")" ),
			    ShortName( parts_ids[hit.query_id] ),
			    ToString( hit.Errors( ) ) ) );
  }

  // Print table.
  PrintTabular( log, table, 3, "rrlr" );
  
}

/**
 * PrintError
 */
void PrintError( const vec< pair<int,int> > &to_putative,
		 const vecString &parts_ids,
		 const String &descrip,
		 const int qq,
		 const int pos,
		 ostream &log )
{
  bool is_primer = to_putative[pos].first < 0;
  int id = is_primer ? - 1 - to_putative[pos].first : to_putative[pos].first;
  int mono_pos = to_putative[pos].second;
  
  String name = "";
  if ( is_primer ) 
    name = char( 97 + (id%3) );
  else {
    name = parts_ids[id];
    if ( name.Contains( "monomer" ) ) name = name.After( "monomer" );
    if ( name.Contains( "termin" ) ) name = name.After( "termin" );
  }
    
  log << pos << "\t"
      << name << "."
      << mono_pos << "\t"
      << descrip << "(q=" 
      << qq << ")\n";
}

