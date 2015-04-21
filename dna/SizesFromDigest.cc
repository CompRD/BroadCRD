/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/// SizesFromDigest.  Given sequence reads from a size-selected restriction 
/// digest, assay the true sizes in the mix.  Output the inferred distribution.
///
/// REF: the reference file to which the reads were aligned.  It is assumed that
/// this was generated in a certain way, via Restrict.
///
/// ALIGNS: the alignments.  They should all be forward and all start at zero on
/// a reference contig.
///
/// READS: fastb file for the reads.  We only use the number of reads.

#include "Basevector.h"
#include "FetchReads.h"
#include "MainTools.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"

int main( int argc, char *argv[] ) 
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(ALIGNS);
     CommandArgument_String(REF);
     CommandArgument_String(READS);
     CommandArgument_Int_OrDefault(MAX_ERROR_DIFF, 2);
     CommandArgument_Int_OrDefault(MAX_DIFF_IN_ALTS, 2);
     CommandArgument_Int_OrDefault(MAX_MIN_ERRORS, 1);
     CommandArgument_Int_OrDefault(SMOOTH_RADIUS, 10);
     EndCommandArguments;

     vecbasevector reads(READS);

     vecbasevector ref;
     vecString refnames;
     FetchReads( ref, refnames, REF );

     vec<look_align> aligns;
     vec< vec<int> > aligns_index;
     LoadLookAligns( ALIGNS, aligns, aligns_index, reads.size( ) );

     vec<int> potential_sizes;
     for ( size_t i = 0; i < refnames.size( ); i++ )
     {    int start = refnames[i].Between( "bases ", "-" ).Int( );
          int stop = refnames[i].After( "bases " ).Between( "-", "," ).Int( );

          // Ignore fragments that cross circular juncture.

          if ( stop <= start ) continue;

          potential_sizes.push_back( stop - start );    }

     vec<int> observed_sizes;
     for ( size_t id = 0; id < reads.size( ); id++ )
     {    int min_errors = 1000000000;
          for ( int j = 0; j < aligns_index[id].isize( ); j++ )
          {    const look_align& la = aligns[ aligns_index[id][j] ];
               min_errors = Min( min_errors, la.Errors( ) );    }
          if ( min_errors > MAX_MIN_ERRORS ) continue;
          vec<int> lengths;
          for ( int j = 0; j < aligns_index[id].isize( ); j++ )
          {    const look_align& la = aligns[ aligns_index[id][j] ];
               if ( la.Errors( ) > min_errors + MAX_ERROR_DIFF ) continue;
               int id = la.query_id;
               int start = refnames[la.target_id].Between( "bases ", "-" ).Int( );
               int stop = refnames[la.target_id].
                    After( "bases " ).Between( "-", "," ).Int( );

               // Ignore fragments that cross circular juncture.

               if ( stop <= start ) continue;

               lengths.push_back( stop - start );    }
          if ( lengths.empty( ) ) continue;
          Sort(lengths);
          int n = lengths.size( );
          if ( n >= 2 && lengths[n-1] - lengths[0] > MAX_DIFF_IN_ALTS ) continue;
          observed_sizes.push_back( lengths[n/2] );    }

     int m = Min(observed_sizes), M = Max(observed_sizes);
     vec<int> observed( M-m+1, 0 ), potential( M-m+1, 0 );
     for ( int i = 0; i < potential_sizes.isize( ); i++ )
     {    if ( potential_sizes[i] < m || potential_sizes[i] > M ) continue;
          ++potential[ potential_sizes[i] - m ];    }
     for ( int i = 0; i < observed_sizes.isize( ); i++ )
          ++observed[ observed_sizes[i] - m ];
     vec<double> rate(M-m+1), srate(M-m+1);
     for ( int j = 0; j <= M-m; j++ )
     {    if ( potential[j] == 0 ) rate[j] = -1;
          else rate[j] = double(observed[j])/double(potential[j]);    }
     for ( int j = m; j <= M; j++ )
     {    vec<double> x;
          for ( int k = -SMOOTH_RADIUS; k <= SMOOTH_RADIUS; k++ )
          {    if ( j+k < m || j+k > M ) continue;
               if ( rate[j+k-m] < 0 ) continue;
               x.push_back( rate[j+k-m] );    }
          if ( x.empty( ) ) srate[j-m] = -1;
          else srate[j-m] = Mean(x);    }
     double srate_sum = 0;
     for ( int j = 0; j <= M-m; j++ )
          if ( srate[j] >= 0 ) srate_sum += srate[j];
     for ( int j = 0; j <= M-m; j++ )
          srate[j] /= srate_sum;
     for ( int j = 0; j <= M-m; j++ )
     {    if ( srate[j] >= 0 ) 
               cout << m+j << " " << setprecision(3) << srate[j] << "\n";    }    }
