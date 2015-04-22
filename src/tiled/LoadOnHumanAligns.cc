// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#include "system/System.h"
#include "Vec.h"
#include "lookup/LookAlign.h"
#include "tiled/LoadOnHumanAligns.h"



/*
 * LoadOnHumanAligns
 */
int LoadOnHumanAligns( vec<look_align_plus> &hits,
		       const String &hits_file,
		       const vec<int> *individual,
		       ostream *input_log,
		       ostream *excluded_log )
{
  ofstream devnull( "/dev/null" );
  ostream &inlog = input_log ? *input_log : devnull; 
  ostream &exlog = excluded_log ? *excluded_log : devnull;
  
  // Reserve memory.
  int n_queries = LineCount( hits_file );
  hits.clear( );
  hits.reserve( n_queries );
  
  // Load.
  ifstream in( hits_file.c_str( ) );
  String a_line;
  static look_align_plus la;
  while( 1 ) {
    getline( in, a_line );
    if ( in.fail( ) )
      break;

    if ( a_line.Contains( "QUERY", 0 ) ) {
      la.ReadParseable( a_line );
      int read_id = la.query_id;

      inlog << la.target_id << "\t" << la.query_id << "\n";

      // Discard reads which do not belong to individual 0.
      if ( individual && 0 != (*individual)[read_id] ) {
	exlog << la.target_id << "\t"
	      << la.query_id << "\tindividual_"
	      << (*individual)[read_id] << "\n";
	continue;
      }
      
      // Add hit.
      hits.push_back( la );
    }
  }
  in.close( );

  // Sort.  
  order_lookalign_TargetLookAlignOffset order;
  sort( hits.begin( ), hits.end( ), order );
  
  // Return.
  return n_queries;
}
