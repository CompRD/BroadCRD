///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Intvector.h"
#include "paths/ReadLoc.h"
#include "util/RunCommand.h"

/**
 * ReadlocsCoverage
 *
 * Generate a quala-like object with base by base coverage (stored as
 * int), as implied from a set of readlocs.
 *
 * INPUT:
 *    <ASSEMBLY>.{contigs.fastb,readlocs}
 *
 * OUTPUT:
 *    <ASSEMBLY>.readlocs.cov (as a quala file)
 *
 * SCORES_PER_LINE: how many scores (coverages) printed per line
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( ASSEMBLY );
  CommandArgument_Int_OrDefault( SCORES_PER_LINE, 20 );
  EndCommandArguments;

  // Dir and file names.
  String run_dir = Dirname( ASSEMBLY );
  String contigs_file = ASSEMBLY + ".contigs.fastb";
  String readlocs_file = ASSEMBLY + ".readlocs";
  String out_file = ASSEMBLY + ".readlocs.cov";

  vec<String> needed;
  needed.push_back( contigs_file );
  needed.push_back( readlocs_file );
  if ( ! CheckFilesExist( needed, &cout )) return 1;

  // Load.
  read_locs_on_disk locs_file( ASSEMBLY, run_dir );

  cout << Date( ) << ": loading contigs" << endl;
  vecbvec contigs( contigs_file );
  size_t n_contigs = contigs.size( );

  // The output.
  VecIntVec coverages( n_contigs );
  for (size_t ii=0; ii<n_contigs; ii++)
    coverages[ii].resize( contigs[ii].size( ), 0 );
  
  // Loop over all contigs.
  size_t dotter = 1000;
  cout << Date( ) << ": parsing " << n_contigs
       << " contigs (.=" << dotter
       << " contigs)" << endl;
  for (size_t cid=0; cid<n_contigs; cid++) {
    if ( cid % dotter == 0 ) Dot( cout, cid / dotter );
    vec<read_loc> locs;
    locs_file.LoadContig( cid, locs );

    // Loop over all reads in contig.
    for (size_t lid=0; lid<locs.size( ); lid++) {
      int begin = Max( locs[lid].Start( ), 0 );
      int end = Min( locs[lid].Stop( ), (int)contigs[cid].size( ) );
      for (int ii=begin; ii<end; ii++)
	coverages[cid][ii]++;
    }
  }
  cout << endl;

  // Save.
  cout << Date( ) << ": saving output" << endl;
  ofstream out( out_file.c_str( ) );
  for (size_t cid=0; cid<n_contigs; cid++) {
    out << ">contig_" << cid;
    for (size_t ii=0; ii<contigs[cid].size( ); ii++) {
      out << ( ii % SCORES_PER_LINE == 0 ? "\n" : " " )
	  << coverages[cid][ii];
    }
    out << '\n';
  }
  out.close( );

  // Done.
  cout << Date( ) << ": ReadlocsCoverage done" << endl;
  
}
