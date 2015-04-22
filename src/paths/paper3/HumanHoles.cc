///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// HumanHoles.  Quantify coverage holes in the human data set.  Note that the
// jump BAM files include alignments for a low fraction of the reads.

#include "Bitvector.h"
#include "FeudalMimic.h"
#include "MainTools.h"
#include "Charvector.h"
#include "math/Functions.h"
#include "lookup/LookAlign.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_Int(MIN_COV);
     EndCommandArguments;

     // Define directory.

     String dir = "/wga/scr1/ALLPATHS/H.sapiens.NA12878";

     // Load coverage map.

     VecUCharVec cov1( dir + "/bams/human.frag.cov" );
     VecUCharVec cov2( dir + "/bams/human.jump.cov" );
     int ntigs = 24;

     // Load genome Ns.

     vecbitvector genome_amb( dir + "/other_ref/build18.fastamb" );

     // Find bases having zero coverage.

     vecbitvector cov;
     Mimic( genome_amb, cov );
     int64_t total = 0, uncovered = 0;
     for ( int i = 0; i < ntigs; i++ )
     {    for ( size_t j = 0; j < cov1[i].size( ); j++ )
          {    if ( genome_amb[i][j] ) continue;
               total++;
               if ( cov1[i][j] + cov2[i][j] < MIN_COV )
               {    // cout << i << "." << j << endl;    
                    uncovered++;    }
               else cov[i].Set( j, True );    }    }
     PRINT2( uncovered, total );
     vec<int> covints;
     for ( int i = 0; i < ntigs; i++ )
     {    for ( size_t j = 0; j < cov[i].size( ); j++ )
          {    if ( !cov[i][j] ) continue;
               size_t start = j;
               while( cov[i][j] && j < cov[i].size( ) ) j++;
               covints.push_back( j - start );    }    }
     Sort(covints);
     cout << "N50 covered stretch = " << N50(covints) << endl;    }
