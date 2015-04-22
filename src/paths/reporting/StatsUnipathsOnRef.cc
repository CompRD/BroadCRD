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
#include "Basevector.h"
#include "Qualvector.h"
#include "SupersHandler.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "paths/UnibaseUtils.h"

/**
 * StatsUnipathsOnRef
 *
 * Simple stats on coverage of reference (Sanger) assembly in terms of
 * given look_align_pluses (these are minimally filtered). Output is
 * sent to cout.
 *
 * HITS: full path name of look_aligns (aligns on contigs)
 * MAX_ERROR_RATE: do not accept aligns with high error rate.
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String( SUBDIR );
  CommandArgument_String( HITS );
  CommandArgument_Double_OrDefault( MAX_ERROR_RATE, 0.1 );
  EndCommandArguments;
  
  // Dir and file names.
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/" + SUBDIR;

  String cgbases_file = sub_dir + "/mergedcontigs.fastb";
  String cgquals_file = sub_dir + "/mergedcontigs.qualb";
  String supers_file = sub_dir + "/mergedcontigs.superb";

  // Load.
  int n_contigs = MastervecFileObjectCount( cgbases_file );
  
  cout << Date( ) << ": loading supers" << endl;
  shandler supers( n_contigs, supers_file );
  vec<int> clens( n_contigs, 0 );
  for (int ii=0; ii<supers.Size( ); ii++)
    for (int jj=0; jj<supers[ii].Ntigs( ); jj++)
      clens[ supers[ii].Tig( jj ) ] = supers[ii].Len( jj );

  cout << Date( ) << ": loading quals" << endl;
  vecqualvector quals;
  quals.ReadAll( cgquals_file );

  cout << Date( ) << ": loading hits (fw only)" << endl;
  vec<look_align_plus> hits;
  vec<int> fhits( n_contigs, -1 );
  {
    vec<look_align_plus> full_hits;
    LoadLookAlignPlus( HITS, full_hits );
    hits.reserve( full_hits.size( ) / 2 );
    for (int ii=0; ii<full_hits.isize( ); ii++)
      if ( ! full_hits[ii].rc1 )
	hits.push_back( full_hits[ii] );

    order_lookalign_TargetBegin sorterT;
    if ( !is_sorted( hits.begin( ), hits.end( ), sorterT ) ) {
      cout << Date( ) << ": sorting hits" << endl;
      sort( hits.begin( ), hits.end( ), sorterT );
    }

    cout << Date( ) << ": generating fhits" << endl;
    for (int ii=hits.isize( )-1; ii>=0; ii--)
      fhits[ hits[ii].target_id ] = ii;
  }
  
  // Counters.
  vec<int> cglen( n_contigs, 0 );     // length of each contig
  vec<int> covlen( n_contigs, 0 );    // covered portion of contig
  vec<int> covlenA( n_contigs, 0 );   // with unipaths with <= 1 error
  vec<int> cglen50( n_contigs, 0 );   // length of q50 portion of contig
  vec<int> covlen50( n_contigs, 0 );  // covered portion of q50 part of contig
  vec<int> covlenA50( n_contigs, 0 ); // with unipaths with <= 1 error

  // Loop over all contigs.
  int dotter = 1000;

  cout << Date( ) << ": parsing "
       << n_contigs << " contigs (.="
       << dotter << " contigs):\n"
       << endl;

  for (int contig_id=0; contig_id<n_contigs; contig_id++) {
    if ( contig_id % dotter == 0 ) Dot( cout, contig_id / dotter );
    int fhit = fhits[contig_id];
    if ( fhit < 0 ) continue;
      
    const qualvector &cgqual = quals[contig_id];
    vec<int> cov( cgqual.size( ), 0 );
    vec<int> covA( cgqual.size( ), 0 );

    for (int hit_id=fhit; hit_id<hits.isize( ); hit_id++) {
      if ( hits[hit_id].target_id != contig_id ) break;
      const look_align_plus &lap = hits[hit_id];

      // Not a proper align.
      int qlen = lap.query_length;
      if  ( lap.a.pos1( ) != 0 && lap.a.Pos1( ) != qlen ) continue;
      
      // Too many errors.
      if ( lap.ErrorRate( ) > MAX_ERROR_RATE ) continue;
      
      // Ok
      int beg = lap.a.pos2( );
      int end = lap.a.Pos2( );
      for (int ii=beg; ii<end; ii++)
	cov[ii] += 1;
      if ( hits[hit_id].mutations < 2 && hits[hit_id].indels < 1 ) 
	for (int ii=beg; ii<end; ii++)
	  covA[ii] += 1;
    }
    
    cglen[contig_id] = cgqual.size( );
    for (int ii=0; ii<cov.isize( ); ii++) {
      if ( cov[ii] > 0 ) covlen[contig_id] += 1;
      if ( covA[ii] > 0 ) covlenA[contig_id] += 1;
      if ( (int)cgqual[ii] >= 50 ) {
	cglen50[contig_id] += 1;
	if ( cov[ii] > 0 ) covlen50[contig_id] += 1;
	if ( covA[ii] > 0 ) covlenA50[contig_id] += 1;
      }
    }
  }
  cout << endl;

  // Print results.
  longlong tot_cglen = BigSum( cglen );
  longlong tot_covlen = BigSum( covlen );
  longlong tot_covlenA = BigSum( covlenA );
  longlong tot_cglen50 = BigSum( cglen50 );
  longlong tot_covlen50 = BigSum( covlen50 );
  longlong tot_covlenA50 = BigSum( covlenA50 );

  double ratio = 100.0 * SafeQuotient( tot_covlen, tot_cglen );
  double ratioA = 100.0 * SafeQuotient( tot_covlenA, tot_cglen );
  double ratio50 = 100.0 * SafeQuotient( tot_covlen50, tot_cglen50 );
  double ratioA50 = 100.0 * SafeQuotient( tot_covlenA50, tot_cglen50 );

  cout << "\n"
       << tot_cglen << "\ttotal contig length\n"
       << tot_covlen << "\ttotal covered length ("
       << ToString( ratio, 2 ) << "%)\n"
       << tot_covlenA << "\ttotal covered length, 1 mismatch or less ("
       << ToString( ratioA, 2 ) << "%)\n"
       << "\n"
       << tot_cglen50 << "\ttotal q50 contig length\n"
       << tot_covlen50 << "\ttotal q50 covered length ("
       << ToString( ratio50, 2 ) << "%)\n"
       << tot_covlenA50 << "\ttotal q50 covered length, 1 mismatch or less ("
       << ToString( ratioA50, 2 ) << "%)\n"
       << endl;

  // Done.
  cout << Date( ) << ": done" << endl;


}
