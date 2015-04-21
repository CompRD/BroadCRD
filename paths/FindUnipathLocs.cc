///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/**
   Program: FindUnipathLocs

   Finds locations of unipaths on the reference.
*/


#include "Alignment.h"
#include "MainTools.h"
#include "lookup/LookAlign.h"
#include "lookup/QueryLookupTableCore.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "pairwise_aligners/PerfectAlignerLG.h"
#include "paths/ImproperMerge.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"
#include "paths/Unipath.h"
#include "paths/simulation/Placement.h"

// MakeDepend: dependency QueryLookupTable

int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_Int(K);
  CommandArgument_Int_OrDefault_Doc( NUM_THREADS, 8, "Number of threads for PerfectAlignerLG" );
  CommandArgument_Bool_OrDefault(SHOW_PLACEMENTS, False);
  CommandArgument_String_OrDefault(READS, "reads");
  EndCommandArguments;
  
  String data_dir = PRE + "/" + DATA;
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String KS = ToString(K);
  
  
  cout << Date() << ": Loading input files" << endl;
  String seqs_file = run_dir + "/" + READS + ".unibases.k" + KS;
  vecbasevector seqs(seqs_file);
  int num_unibases = seqs.size();
  seqs.ReadAll( data_dir + "/genome.fastb", True );
  
  // Find locations of unipaths on reference, using PerfectAlignerLG.
  cout << Date() << ": First attempt at aligning (using PerfectAlignerLG)" << endl;
  vec<alignment_plus> aligns;
  PerfectAlignerLG( K, PerfectAlignerLG::findProperOnly )
    .Align( seqs, aligns, NUM_THREADS, num_unibases );
  
  // Try to align harder, using QueryLookupTable.
  cout << Date() << ": Second attempt at aligning (using QueryLookupTable)" << endl;
  vec< vec<int> > aligns_index(num_unibases);
  for ( int i = 0; i < aligns.isize( ); i++ ) {
    const alignment_plus& ap = aligns[i];
    aligns_index[ ap.Id1( ) ].push_back(i);
  }
  Mkdir777( run_dir + "/tmp" );
  temp_file qlt( run_dir + "/tmp/QueryLookupTable.out.XXXXXXX" );
  QueryLookupTableCore( "K=12 L=" + data_dir + "/genome.lookup " + " MM="
			+ ToString(K-1) + " SEQS=" + seqs_file + " SMITH_WAT=True"
			+ " SEQS_IS_FASTB=True PARSEABLE=True NH=True OUTFILE=" + qlt );
  vec<look_align> aligns2;
  vec< vec<int> > aligns_index2;
  LoadLookAligns( qlt, aligns2, aligns_index2, num_unibases );
  VecPlacementVec locs0( num_unibases );
  for ( int i = 0; i < aligns.isize( ); i++ ) {
    const alignment_plus& ap = aligns[i];
    if ( !ap.Rc2( ) )
      locs0[ ap.Id1( ) ].push_back( placement( ap.Id2( ), 
					       ap.a.pos2( ), ap.a.Pos2( ), ap.Rc2( ) ) );
    else {
      int g = seqs[ num_unibases + ap.Id2( ) ].size( );
      locs0[ ap.Id1( ) ].push_back( placement( ap.Id2( ), 
					       g - ap.a.Pos2( ), g - ap.a.pos2( ), ap.Rc2( ) ) );
    }
  }
  
  for ( int u = 0; u < num_unibases; u++ ) {
    if ( aligns_index[u].nonempty( ) ) continue;
    for ( int i = 0; i < aligns_index2[u].isize( ); i++ ) {
      const look_align& la = aligns2[ aligns_index2[u][i] ];
      locs0[la.query_id].push_back( placement( la.target_id,
					       la.pos2( ), la.Pos2( ), la.Rc1( ) ) );
    }
  }
  
  if (SHOW_PLACEMENTS) cout << "\n";
  
  // Sort the set of placements of each unipath.
  for ( VecPlacementVec::size_type i = 0; i < locs0.size( ); i++ ) {
    PlacementVec& locs0i = locs0[i];
    using std::sort;
    sort( locs0i.begin(), locs0i.end() );
    if ( SHOW_PLACEMENTS ) {
      for ( PlacementVec::size_type j = 0; j < locs0i.size( ); j++ ) {
	const placement& p = locs0i[j];
	cout << "unipath " << i << " (l="
	     << ( seqs[i].size( ) - K + 1 ) << ") placed at "
	     << p.GenomeId( ) << "." << p.pos( ) << "-"
	     << p.Pos( ) << " (" << ( p.Fw( ) ? "fw" : "rc" )
	     << ")\n";
      }
    }
  }
  
  locs0.WriteAll( run_dir + "/"+READS+".unipaths.k" + KS + ".locs" );
  
  cout << Date( ) << ": Done with FindUnipathLocs!" << endl;
}
