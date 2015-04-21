// Copyright (c) 2003 Broad Institute/Massachusetts Institute of Technology

#include <time.h>

#include "math/Functions.h"
#include "Qualvector.h"
#include "system/System.h"
#include "Vec.h"
#include "select/GenerateQ20Pluses.h"



/*
 * GenerateQ20Pluses
 */
void GenerateQ20Pluses( vec<int> &q20pluses,
			const String &quals_file,
			ostream *log)
{
  ofstream devnull( "/dev/null" );
  ostream &out = ( log ) ? *log : devnull;

  // The q20pluses file name.
  String plain_qfile = quals_file;
  while( plain_qfile.Contains( "/" ) )
    plain_qfile = plain_qfile.After( "/" );
  String q20_file = quals_file.Before( plain_qfile ) + "reads.q20pluses";

  // The q20pluses file exists, and it is younger than the quals file.
  if ( IsRegularFile( q20_file ) ) {
    struct stat file_stat;

    ForceAssert( stat( q20_file.c_str(), &file_stat ) == 0 );
    time_t q20_time = file_stat.st_mtime;

    ForceAssert( stat( quals_file.c_str(), &file_stat ) == 0 );
    time_t quals_time = file_stat.st_mtime;

    if ( q20_time > quals_time ) {
      READX( q20_file, q20pluses );
      return;
    }
  }
  
  // Generate q20pluses. Page in chunk_size quals at a time.
  int n_quals = MastervecFileObjectCount( quals_file );

  q20pluses.clear( );
  q20pluses.resize( n_quals, 0 );
  
  int chunk_begin = 0;
  int chunk_size = 1000000;
  while ( chunk_begin < n_quals ) {
    int chunk_end = Min( chunk_begin + chunk_size, n_quals );
    out << Date( ) << ": generating q20pluses for quals "
	<< chunk_begin << " to "
	<< chunk_end - 1
	<< endl;

    vecqualvector quals;
    vec<int> ids( chunk_end - chunk_begin, -1 );
    for (int ii=0; ii<chunk_end-chunk_begin; ii++)
      ids[ii] = chunk_begin + ii;
    quals.SparseRead( quals_file, ids, 0 );

    for (int ii=0; ii<(int)ids.size( ); ii++) {
      int id = ids[ii];
      const qualvector &the_qual = quals[id];
      for (int jj=0; jj<(int)the_qual.size( ); jj++)
	if( int( the_qual[jj] ) >= 20 )
	  q20pluses[id] += 1;
    }
    
    chunk_begin += chunk_size;
  }
  
  // Save.
  out << Date( ) << ": saving q20pluses onto " << q20_file << endl;
  WRITE( q20_file, q20pluses );

}
