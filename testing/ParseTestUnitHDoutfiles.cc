/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "CoreTools.h"
#include "PrettyPrintTable.h"
#include "TokenizeString.h"
#include "math/Functions.h"

void CollectData( const String &in_file, vec< vec<String> > &table );

void DigestData( const vec< vec<String> > &data, vec< vec<String> > &digest );

/**
 * ParseTestUnitHDoutfiles
 *
 * Recursively parse the given BASE_DIR and digest all the log files from
 * TestUnitHD. Output sent to files in BASE_DIR.
 *
 * BASE_DIR: full path name of root dir where TestUnitHD was run
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( BASE_DIR );
  EndCommandArguments;

  vec< vec<String> > table;
  CollectData( BASE_DIR, table );
  BeautifyTable( table );

  vec< vec<String> > digest;
  DigestData( table, digest );
  BeautifyTable( digest);

  for (int ii=0; ii<2; ii++) {
    vec< vec<String> > &toprint = ( ii == 0 ) ? table : digest;
    String outfile = BASE_DIR + ( ii == 0 ? "/full.out" : "/digest.out" );
    ofstream out( outfile.c_str( ) );

    for (int ii=0; ii<(int)toprint.size( ); ii++) {
      for (int jj=0; jj<(int)toprint[ii].size( ); jj++)
	out << toprint[ii][jj] << " ";
      out << "\n";
    }

    out.close( );
  }

}

/**
 * CollectData
 *
 * The recursive core function to harvest data. Find all outputs from
 * TestUnitHD runs and save the content as vec<String>.
 */
void CollectData( const String &in_file, vec< vec<String> > &table )
{
  if ( table.size( ) < 1 ) {
    vec<String> row;
    row.push_back( "wga," );
    row.push_back( "scr," );
    row.push_back( "n_files," );
    row.push_back( "tot_size_Gb," );
    row.push_back( "garbage_size_Gb," );
    row.push_back( "buffer_size_b," );
    row.push_back( "writing," );
    row.push_back( "garbage_in," );
    row.push_back( "garbage_out," );
    row.push_back( "read_1," );
    row.push_back( "read_2," );
    row.push_back( "delete_files" );
    
    table.push_back( row );
  }
  
  vec<String> all_files = AllFiles( in_file );
  for (int ii=0; ii<(int)all_files.size( ); ii++) {
    const String full_file = in_file + "/" + all_files[ii];
    if ( IsDirectory( full_file ) ) CollectData( full_file, table );
    if ( all_files[ii] != "TestUnitHD.log" ) continue;
    
    vec<String> row( 12, "" );

    String aline;
    vec<String> tokens;
    ifstream in( full_file.c_str( ) );
    while ( in ) {
      getline( in, aline ); 
      if ( !in ) break;
      Tokenize( aline, tokens );
      
      if ( aline.Contains( "Testing", 0 ) ) {
	row[0] = tokens[3] +  ",";
	row[1] = tokens[1] + ",";
	continue;
      }
      if ( aline.Contains( "number of files", 0 ) ) {
	row[2] = tokens[3] + ",";
	continue;
      }
      if ( aline.Contains( "total size", 0 ) ) {
	row[3] = tokens[2] + ",";
	continue;
      }
      if ( aline.Contains( "size of garbage", 0 ) ) {
	row[4] = tokens[4] + ",";
	continue;
      }
      if ( aline.Contains( "buffer size", 0 ) ) {
	row[5] = tokens[2] + ",";
	continue;
      }
      if ( tokens.size( ) > 5 && tokens[5] == "writing..." ) {
	row[6] = tokens[6] + ",";
	continue;
      }
      if ( tokens.size( ) > 6 && tokens[6] == "garbage" ) {
	if ( tokens[5] == "write" ) row[7] = tokens[10] + ",";
	else row[8] = tokens[10] + ",";
	continue;
      }
      if ( tokens.size( ) > 5 && tokens[5] == "reading" ) {
	if ( tokens[6] == "(first" ) row[9] = tokens[8] + ",";
	else row[10] = tokens[8] + ",";
	continue;
      }
      if ( tokens.size( ) > 5 && tokens[5] == "deleting" ) {
	row[11] = tokens[7];
	continue;
      }
    }
    in.close( );

    table.push_back( row );
  }

}

/**
 * DigestData
 *
 * Digest collected data (runs with the same number of files will be
 * pooled together and averaged).
 */
void DigestData( const vec< vec<String> > &data, vec< vec<String> > &digest )
{
  vec<int> nfiles;                       // n_files
  vec< vec< vec<double> > > coredata;    // stats for all runs with n_files

  // Find stats.
  for (int ii=0; ii<(int)data.size( ); ii++) {
    if ( data[ii][0].Contains( "wga," ) )
      continue;

    int nf = data[ii][2].Before( "," ).Int( );
    double twrite = data[ii][6].Before( "," ).Double( );
    double tread1 = data[ii][9].Before( "," ).Double( );
    double tread2 = data[ii][10].Before( "," ).Double( );
    double tdel = data[ii][11].Double( );

    vec<double> newset;
    newset.push_back( twrite );
    newset.push_back( tread1 );
    newset.push_back( tread2 );
    newset.push_back( tdel );
    
    int pos = -1;
    vec<int>::iterator it = find( nfiles.begin( ), nfiles.end( ), nf );
    if ( it == nfiles.end( ) ) {
      pos = nfiles.size( );
      nfiles.push_back( nf );
      coredata.resize( coredata.size( ) + 1 );
    }
    else
      pos = distance( nfiles.begin( ), it );

    coredata[pos].push_back( newset );
  }

  // Average on all runs.
  digest.clear( );

  vec<String> summary;
  summary.push_back( "n_files," );
  summary.push_back( "write," );
  summary.push_back( "read1," );
  summary.push_back( "read2," );
  summary.push_back( "delete," );
  digest.push_back( summary );

  for (int ii=0; ii<(int)coredata.size( ); ii++) {
    summary.clear( );

    summary.push_back( ToString( nfiles[ii] ) + "," );

    for (int kk=0; kk<4; kk++) {
      vec<double> tempdata;
      for (int jj=0; jj<(int)coredata[ii].size( ); jj++)
	tempdata.push_back( coredata[ii][jj][kk] );
      sort( tempdata.begin( ), tempdata.end( ) );
      double tmed = Median( tempdata, 0 );
      summary.push_back( ToString( tmed ) + "," );
    }
    digest.push_back( summary );
    
  }
  
}

