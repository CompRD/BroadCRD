///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Simulation to test ability to phase SNPs with long reads.

#include "Basevector.h"
#include "MainTools.h"
#include "random/Random.h"

int main( )
{    
     // Genome size.  More or less irrelevant.

     const int gsize = 10000;

     // SNP positions.  Only distance between them relevant.  And distance from
     // ends needs to exceed read length.

     const int pos1 = 4000;
     const int pos2 = 5000;

     // Read properties.

     const int rl = 3000;
     const double cov = 20.0;
     const double err_rate = 0.08;

     // Other.

     const int min_ratio = 5;
     const int sample = 1000000;

     int nreads = int( round( (gsize*cov) / rl ) );

     int good = 0, bad = 0;
     for ( int pass = 1; pass <= sample; pass++ )
     {    int AA = 0, TT = 0, AT = 0, TA = 0;
          for ( int i = 0; i < nreads; i++ )
          {    int start = randomx( ) % ( gsize - rl );
               if ( start > pos1 ) continue;
               int stop = start + rl;
               if ( stop <= pos2 ) continue;
               int allele = randomx( ) % 2;
               int b1, b2;
               if ( allele == 0 )
               {    b1 = 0, b2 = 0;    }
               else
               {    b1 = 3, b2 = 3;    }
               if ( randomx( ) % 1000000 < err_rate * 1000000 )
                    b1 = ( b1 + 1 + ( randomx( ) % 3 ) ) % 4;
               if ( randomx( ) % 1000000 < err_rate * 1000000 )
                    b2 = ( b2 + 1 + ( randomx( ) % 3 ) ) % 4;
               if ( b1 == 0 && b2 == 0 ) AA++;
               if ( b1 == 3 && b2 == 3 ) TT++;
               if ( b1 == 0 && b2 == 3 ) AT++;
               if ( b1 == 3 && b2 == 0 ) TA++;    }
          if ( AA >= 2 && TT >= 2 && AA + TT >= min_ratio * (AT+TA) ) good++;
          else if ( AT >= 2 && TA >= 2 && AT + TA >= min_ratio * (AA+TT) ) 
          {    bad++;
               PRINT4( AA, TT, AT, TA );    }    }

     PRINT3( good, bad, sample );    }
