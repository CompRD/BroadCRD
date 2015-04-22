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
#include "SupersHandler.h"
#include "lookup/LookAlign.h"
#include "paths/UnibaseUtils.h"

/**
 * PrintUnipathsOnRef
 *
 * Load a set of look_align_plus and sort them by start on target
 * (where targets are the contigs of the given Sanger assembly). Is
 * SUPERS is given, print output one super at a time. Output is sent
 * to cout.
 *
 * HITS: full path name of objects aligning the contigs
 * REFERENCE: full path name of rereference (fastb)
 * SUPERS: if given, load superb structure (where REFERENCE are the contigs)
 * FW_ONLY: do not print rc aligns
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( HITS );
  CommandArgument_String( REFERENCE );
  CommandArgument_String_OrDefault( SUPERS, "" );
  CommandArgument_Bool_OrDefault( FW_ONLY, True );
  EndCommandArguments;
  
  // Load.
  int n_contigs = MastervecFileObjectCount( REFERENCE );

  shandler supers( n_contigs );
  if ( SUPERS != "" ) {
    cout << Date( ) << ": loading supers" << endl;
    supers.LoadFromFile( SUPERS );
  }

  cout << Date( ) << ": loading hits" << endl;
  vec<look_align_plus> hits;
  {
    vec<look_align_plus> full_hits;
    LoadLookAlignPlus( HITS, full_hits );
    if ( FW_ONLY ) hits.reserve( full_hits.size( ) / 2 );
    else hits.reserve( full_hits.size( ) );
    for (int ii=0; ii<full_hits.isize( ); ii++) {
      if ( FW_ONLY && full_hits[ii].rc1 ) continue;
      hits.push_back( full_hits[ii] );
    }
  }
  
  // Sort hits by target_id, start_on_target.
  cout << Date( ) << ": sorting hits\n" << endl;
  order_lookalign_TargetBegin sorterT;
  sort( hits.begin( ), hits.end( ), sorterT );

  // fhits map.
  vec<int> fhits( n_contigs, -1 );
  for (int ii=hits.isize( )-1; ii>=0; ii--)
    fhits[ hits[ii].target_id ] = ii;

  // Print info (with or without supers info).
  if ( SUPERS != "" ) {  
    for (int super_id=0; super_id<supers.Size( ); super_id++) {
      const superb &sup = supers[super_id];

      // Loop over contigs in super.
      for (int cgpos=0; cgpos<sup.Ntigs( ); cgpos++) {
	int contig_id = sup.Tig( cgpos );
	cout << "s" << super_id
	     << "_" << cgpos
	     << "." << sup.Ntigs( )-1
	     << " = c" << contig_id
	     << "\n";

	int fhit = fhits[contig_id];
	if ( fhit < 0 ) {
	  cout << "  no hits on this contig\n";
	  continue;
	}

	// Loop over hits on contig.
	for (int hitid=fhit; hitid<hits.isize( ); hitid++) {
	  if ( hits[hitid].TargetId( ) != contig_id ) break;
	  cout << "  ";
	  hits[hitid].PrintParseable( cout );
	}
      }
    }
  }
  else {
    for (int contig_id=0; contig_id<n_contigs; contig_id++) {
      cout << "c" << contig_id << "\n";
	   
      int fhit = fhits[contig_id];
      if ( fhit < 0 ) {
	cout << "  no hits on this contig\n";
	continue;
      }
      
      // Loop over hits on contig.
      for (int hitid=fhit; hitid<hits.isize( ); hitid++) {
	if ( hits[hitid].TargetId( ) != contig_id ) break;
	cout << "  ";
	hits[hitid].PrintParseable( cout );
      }
    }  
  }
  
  // Done.
  cout << "\n" << Date( ) << ": done" << endl;

}
