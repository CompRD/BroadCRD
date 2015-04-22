// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#include "DeletedReads.h"
#include "String.h"
#include "system/System.h"



/*
 * deleted_reads
 * Constructor
 */
deleted_reads::deleted_reads(const String &deleted_reads_file,
			     const String &unknown_reads_file) :
  deleted_reads_file_ ( deleted_reads_file ),
  unknown_reads_file_ ( unknown_reads_file )
{
  del_stream_.open( deleted_reads_file_.c_str( ), ios::app );
  unk_stream_.open( unknown_reads_file_.c_str( ), ios::app );
  
  Assert ( !del_stream_.bad() || !unk_stream_.bad() );
}



/*
 * deleted_reads
 * CreateFiles
 */
void deleted_reads::CreateFiles()
{
  del_stream_.close();
  unk_stream_.close();

  del_stream_.open( deleted_reads_file_.c_str( ) );
  unk_stream_.open( unknown_reads_file_.c_str( ) );
}



/*
 * deleted_reads
 * GetDeletedReads
 */
vec< del_datum > deleted_reads::GetDeletedReads()
{
  // I will read the file twice: the first time to get an idea of
  // its size, and the second to gather the data.
  vec< del_datum > del_reads;

  int n_lines = LineCount( deleted_reads_file_ );
  del_reads.reserve( n_lines );

  // Get data.
  ifstream in( deleted_reads_file_.c_str( ), ios::in );
  while ( 1 ) {
    String read_name;
    int read_id;
    int n_del_code;

    in >> n_del_code >> read_name >> read_id;
    if ( !in )
      break;
    
    del_datum deleted(read_name, read_id, (deletion_code)n_del_code);
    del_reads.push_back( deleted );
  }
  
  return del_reads;
}



/*
 * deleted_reads
 * GetUnknownReads
 */
vec<String> deleted_reads::GetUnknownReads()
{
  vec<String> unknown_reads;
  
  // Determine size of file.
  ifstream in( unknown_reads_file_.c_str( ) );
  if ( in.bad() )
    return unknown_reads;

  int n_lines = LineCount( unknown_reads_file_ );
  unknown_reads.reserve( n_lines );

  // Load data.
  while ( 1 ) {
    String read;

    in >> read;
    if ( !in )
      break;
    
    unknown_reads.push_back(read);
  }

  return unknown_reads;
}

/**
 * deleted_reads
 * NoIdsCount
 */
int deleted_reads::NoIdsCount( )
{
  int tot = 0;
  
  ifstream in( deleted_reads_file_.c_str( ), ios::in );
  while ( 1 ) {
    String read_name;
    int read_id;
    int n_del_code;
    in >> n_del_code >> read_name >> read_id;
    if ( !in ) break;
    if ( read_id < 0 ) tot++;
  }
  in.close( );
  
  return tot;
}

/**
 * deleted_reads
 * NoIdsNames
 */
void deleted_reads::NoIdsNames( vec<String> &names )
{
  names.clear( );

  vec<del_datum> all_deleted = this->GetDeletedReads( );
  names.reserve( all_deleted .size( ) );
  for (int ii=0; ii<(int)all_deleted.size( ); ii++)
    if ( all_deleted[ii].GetReadId( ) < 0 )
      names.push_back( all_deleted[ii].GetName( ) );

  sort( names.begin( ), names.end( ) );
}

