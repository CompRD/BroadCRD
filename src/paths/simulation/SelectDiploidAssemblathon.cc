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
#include "pairwise_aligners/AlignTwoBasevectors.h"

/**
 * SelectDiploidAssemblathon
 *
 * Select a diploid interval from a range (eg "0:500000-600000"). The
 * code assumes that the intput fastb contains n "chromosomes", where
 * chromosomes i and n+i are matched.
 */ 
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( FASTB_IN );
  CommandArgument_String( FASTB_OUT );
  CommandArgument_String( RANGE );
  EndCommandArguments;
  
  int id = RANGE.Before( ":" ).Int( );
  int begin = RANGE.After( ":" ).Before( "-" ).Int( );
  int end = RANGE.After( "-" ).Int( );
  int length = end - begin;

  vecbvec bases( FASTB_IN );
  int n_copies = (int)bases.size( ) / 2;

  // Haplotype 1: matches RANGE exactly.
  vecbvec select;
  select.reserve( 2 );
  select.push_back( bvec( bases[id], begin, length ) );

  // Haplotype 2: look for region (try to match edges).
  int stitch_size = 3000;
  int slop = 100000;

  bvec anchorL( bvec( select[0], 0, stitch_size ) );
  bvec anchorR( bvec( select[0], length - stitch_size, stitch_size ) );

  const bvec &full_target= bases[id+n_copies];
  const int fulltsize = (int)full_target.size( );
  int targetL_beg = Max( 0, begin - slop );
  int targetL_end = Min( fulltsize, begin + slop + stitch_size );
  int targetL_len = targetL_end - targetL_beg;
  int targetR_beg = Max( 0, end - stitch_size - slop );
  int targetR_end = Min( fulltsize, end + slop );
  int targetR_len = targetR_end - targetR_beg;
  bvec targetL( full_target, targetL_beg, targetL_len );
  bvec targetR( full_target, targetR_beg, targetR_len );
  
  ofstream devnull( "/dev/null" );
  int minOv = 0;
  int maxOv = 2 * stitch_size;
  float maxErrRate = 0.35;
  ostream &a2bLog = devnull;
  int a2bRC = false;
  int mode = 5;
  int K = 24;
  int stretch = 12;
  int nstretch = 6;
  const qvec q1( 0 );
  const qvec q2( 0 );
  float maxScore = 1000000.0;
  float maxErr = 5000;
  ostream &badlog = devnull;
  bool avoidProm = False;
  int maxCliq8 = 1000;
  int maxAligns8 = 10000;
  int maxErr8 = 1000;
  int locMaxErr = 150;
  bool altMethod = false;
  int band = 0;
  int minMutmer = 0;
  float ambigThreshold = 3.0;
  double minPolyscoreRatio = 100.0;
  bool swMethod = True;

  align alignL;
  AlignTwoBasevectors( anchorL, targetL, alignL, minOv, maxOv, maxErrRate,
		       &a2bLog, a2bRC, mode, K, stretch, nstretch, q1, q2,
		       maxScore, maxErr, badlog, avoidProm, maxCliq8,
		       maxAligns8, maxErr8, locMaxErr, altMethod, band,
		       minMutmer, ambigThreshold, minPolyscoreRatio,
		       swMethod );
  
  align alignR;
  AlignTwoBasevectors( anchorR, targetR, alignR, minOv, maxOv, maxErrRate,
		       &a2bLog, a2bRC, mode, K, stretch, nstretch, q1, q2,
		       maxScore, maxErr, badlog, avoidProm, maxCliq8,
		       maxAligns8, maxErr8, locMaxErr, altMethod, band,
		       minMutmer, ambigThreshold, minPolyscoreRatio,
		       swMethod );

  bool is_bad = false;
  if ( alignL.pos1( ) != 0 || alignL.Pos1( ) != stitch_size )
    is_bad = true;
  else if ( alignR.pos1( ) != 0 || alignR.Pos1( ) != stitch_size )
    is_bad = true;

  int beg2 = Max( 0, alignL.pos2( ) - slop + begin );
  int end2 = Min( fulltsize, alignR.Pos2( ) - slop - stitch_size + end );
  if ( is_bad ) {
    beg2 = begin;
    end2 = end;
    cout << "WARNING! No proper alignment found. Stitches align as:\n"
	 << "  left: " << alignL.pos1( ) << " to " << alignL.Pos1( )
	 << "\t" << alignL.pos2( ) << " to " << alignL.Pos2( ) << "\n"
    	 << "  right: " << alignR.pos1( ) << " to " << alignR.Pos1( )
	 << "\t" << alignR.pos2( ) << " to " << alignR.Pos2( ) << "\n"
	 << "Using fallback ranges:\n";
  }
  else
    cout << "Alignment found:\n";    
  int len2 = end2 - beg2;
  
  cout << "h1: [" << begin << ", " << end << ")_"
       << bases[id].size( ) << " on id_" << id << "\n"
       << "h2: [" << beg2 << ", " << end2 << ")_"
       << bases[id+n_copies].size( ) << " on id_" << id+n_copies << "\n"
       << endl;

  select.push_back( bvec( bases[id+n_copies], beg2, len2 ) );

  // Save.
  select.WriteAll( FASTB_OUT );

}

