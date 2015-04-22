/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"
#include "CoverageAnalyzer.h"
#include "Loader.h"
#include "ParseSet.h"
#include "graphics/CContigPlotter.h"
#include "lookup/LookAlign.h"
#include "util/PickTic.h"

/**
 * SimpleACGTPlot
 *
 * Plot ACGT content on sliding windows (the plot represents the
 * percent of combined As and Ts).
 *
 * FASTB: full path name of target fastb
 * OUTDIR: full path name for output dir
 * TARGET_IDS: all if empty (caution! all means a lot of output)
 * WIN_SIZE: sliding windows size
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( FASTB );
  CommandArgument_String( OUTDIR );
  CommandArgument_String_OrDefault( TARGET_IDS, "" );
  CommandArgument_UnsignedInt_OrDefault( WIN_SIZE, 1000 );
  EndCommandArguments;

  Mkdir777( OUTDIR );
  
  // Load.
  cout << Date( ) << ": loading target fastb" << endl;
  vecbasevector targets;
  targets.ReadAll( FASTB );

  vec<int> tlens;
  LoadContigLengths( targets, tlens );
  
  vec<int> select;
  if ( TARGET_IDS == "" ) {
    select.reserve( tlens.size( ) );
    for (int ii=0; ii<(int)tlens.size( ); ii++)
      select.push_back( ii );
  }
  else
    ParseIntSet( TARGET_IDS, select );

  // Generate data.
  cout << Date( ) << ": looking at " << select.size( ) << " contigs" << endl;
  for (int sel_id=0; sel_id<(int)select.size( ); sel_id++) {
    Dot( cout, sel_id );

    // Files and txt stream.
    int target_id = select[sel_id];
    int target_len = targets[target_id].size( );
    int target_nwin = Max( 0, target_len - (int)WIN_SIZE );
    
    String base_out = OUTDIR + "/contig_" + ToString( target_id );
    String txt_out = base_out + ".txt";
    String eps_out = base_out + ".eps";

    ofstream out( txt_out.c_str( ) );

    if ( target_nwin < 1 ) {
      out << "contig length (" << target_len
	  << " bp) is smaller than WIN_SIZE (" << WIN_SIZE
	  << " bp). No data generated\n";
      out.close( );
      continue;
    }

    // AT coverage.
    vec<float> cov_sol;
    
    const basevector &target = targets[target_id];
    CDataPlotter data( red, 0.4 );
    data.Reserve( target_nwin );
    
    vec<int> ATvec;
    for (int ii=0; ii<(int)target.size( ); ii++) {
      char aschar = as_base( target[ii] );
      ATvec.push_back( ( aschar == 'A' || aschar == 'T' ) ? 1 : 0 );
    }

    longlong cov_amount = 0;
    for (int pos=0; pos<target_nwin; pos++) {
      if ( pos == 0 ) {
	for (int ii=0; ii<(int)WIN_SIZE; ii++)
	  cov_amount += ATvec[ii];
      }
      else
	cov_amount += ATvec[-1+pos+(int)WIN_SIZE] - ATvec[pos-1];
      float xpt = pos;
      float ypt = SafeQuotient( cov_amount, WIN_SIZE );
      data.AddPoint( xpt, ypt );
      cov_sol.push_back( ypt );
    }
    
    // Tics.
    vec<CTic> xtics;
    StockXTics( xtics, target_nwin );
  
    vec<CTic> ytics;
    StockYTics( ytics, 1 );

    // Generate eps.
    CContigPlotter plotter( target_nwin );
    plotter.AddData( data );
    plotter.AddHorizontalTics( xtics );
    plotter.AddVerticalTics( ytics );
  
    plotter.GoFigure( eps_out );
    
    // Save as text.
    out << "pos\tATcoverage\n";
    for (int ii=0; ii<(int)cov_sol.size( ); ii++)
      out << ii << "\t" << cov_sol[ii] << "\n";
    out.close( );

  }
  cout << endl;

  // Done.
  cout << Date( ) << ": done" << endl;

}
