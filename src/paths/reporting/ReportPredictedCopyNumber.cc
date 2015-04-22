/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "VecUtilities.h"
#include "math/Functions.h"
#include "paths/KmerPath.h"
#include "paths/PdfEntry.h"

/**
 * ReportPredictedCopyNumber
 *
 * Load <HEAD>.predicted_count.k<K> and print vartious statistics.
 */
int main( int argc, char *argv[] )
{
  BeginCommandArguments;
  CommandArgument_Int( K );
  CommandArgument_String( HEAD );
  CommandArgument_Int_OrDefault( MAX_CN, 12 );
  EndCommandArguments;

  String strK = "k" + ToString ( K );
  String copynum_file = HEAD + ".predicted_count." + strK;
  String unipaths_file = HEAD + "." + strK;

  cout << Date( ) << ": loading copy numbers" << endl;
  VecPdfEntryVec copynum( copynum_file );
  
  cout << Date( ) << ": loading unipaths" << endl;
  vecKmerPath unipaths( unipaths_file );
  
  // The vector of lengths for each CN.
  cout << Date( ) << ": parsing unipaths" << endl;
  vec< vec<int> > lens( 1 + MAX_CN );
  for (int ii=1; ii<lens.isize( ); ii++)
    lens[ii].reserve( unipaths.size( ) );

  longlong tot_unipaths_len = 0;
  for (size_t ii=0; ii<copynum.size( ); ii++) {
    tot_unipaths_len += (longlong)unipaths[ii].TotalLength( );
    int cn = copynum[ii][0].first;
    double prob = copynum[ii][0].second;
    for (size_t jj=1; jj<copynum[ii].size( ); jj++) {
      if ( copynum[ii][jj].second > prob ) {
	cn = copynum[ii][jj].first;
	prob = copynum[ii][jj].second;	
      }
    }
    cn = Min( MAX_CN, cn );
    lens[cn].push_back( unipaths[ii].TotalLength( ) );
  }
  
  vec<longlong> tot_lens( 1 + MAX_CN );
  for (int ii=0; ii<lens.isize( ); ii++)
    tot_lens[ii] = BigSum( lens[ii] );

  // The table.
  vec< vec<String> > table;
  vec<String> aline;
  {
    aline.push_back( "CN" );
    aline.push_back( "n_unipaths" );
    aline.push_back( "total_length" );
    aline.push_back( "N50" );
    aline.push_back( "%_length" );
    aline.push_back( "%_cumul" );
    table.push_back( aline );
  }

  for (int ii=1; ii<=MAX_CN; ii++) {
    aline.clear( );
    sort( lens[ii].begin( ), lens[ii].end( ) );
    longlong cumul_len = 0;
    for (int jj=ii; jj<=MAX_CN; jj++)
      cumul_len += tot_lens[jj];
    int n_unipaths = lens[ii].size( );
    int n50_len = lens[ii].size( ) < 1 ? 0 : N50( lens[ii] );
    double ratioL = 100.0 * SafeQuotient( tot_lens[ii], tot_unipaths_len );
    double ratioC = 100.0 * SafeQuotient( cumul_len, tot_unipaths_len );

    aline.push_back( ToString( ii ) );
    aline.push_back( ToString( n_unipaths ) );
    aline.push_back( ToString( tot_lens[ii] ) );
    aline.push_back( ToString( n50_len ) );
    aline.push_back( ToString( ratioL, 2 ) );
    aline.push_back( ToString( ratioC, 2 ) );
    table.push_back( aline );
  }

 
  // Print.
  cout << "\n";
  PrintTabular( cout, table, 3, "rrrrrr" );
  cout << "\n";


  // Report copy number for 20 largest unipaths
  vec<pair<uint64_t, uint64_t> > CN_summary(unipaths.size());
  for (size_t i = 0; i < unipaths.size(); i++) {
    int cn = copynum[i][0].first;
    double prob = copynum[i][0].second;
    for (size_t jj=1; jj<copynum[i].size( ); jj++) {
      if ( copynum[i][jj].second > prob ) {
	cn = copynum[i][jj].first;
	prob = copynum[i][jj].second;	
      }
    }
    CN_summary[i] = make_pair(unipaths[i].TotalLength(),cn);
  }
  
  sort( CN_summary.begin(), CN_summary.end(), greater< pair< uint64_t,uint64_t > > () );
  cout << "Copy numbers for 20 largest unipaths:" << endl;
  cout << "Length, CN" << endl;
  for (size_t i = 0; i < 20; i++) {
    cout << CN_summary[i].first << ", " << CN_summary[i].second << endl;
  }    

  // Done.
  cout << Date( ) << ": done" << endl;
  
}
