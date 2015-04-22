/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// MutationSensitivity.  Predict sensitivity of somatic mutation calling, either
// at the coverage provided by a given data set, or at higher coverage.

/*

Example of typical usage:

MakeRandomRanges REF_AMB=/wga/scr1/LOOKUP/human36/hg18.fastamb N=10K S=1 NH=True \
     SEED=543219876 > random_ranges_10K

ExtractRangeFromCLFs \
     CLFS_DIR=/home/radon01/kiran/tmp/tcga/tcga-tumor-freeze1-081121/clfs/wg.e4/ \
     REFERENCE=/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta.lookuptable.fastb \
     RANGES=@random_ranges_10K SILENT=True BIN_QUALS=random_ranges_10K_quals_tumor

ExtractRangeFromCLFs \
     CLFS_DIR=/home/radon01/kiran/tmp/tcga/tcga-normal-freeze1-081121/clfs/wg.e4/ \
     REFERENCE=/seq/references/Homo_sapiens_assembly18/v0/Homo_sapiens_assembly18.fasta.lookuptable.fastb \
     RANGES=@random_ranges_10K SILENT=True BIN_QUALS=random_ranges_10K_quals_normal

MutationSensitivity TUMOR_QUALS=random_ranges_10K_quals_tumor \
     NORMAL_QUALS=random_ranges_10K_quals_normal

*/

#include "MainTools.h"
#include "Qualvector.h"
#include "random/Random.h"
#include "util/AlleleLikelihood.h"

int main(int argc, char **argv)
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(TUMOR_QUALS);
     CommandArgument_String(NORMAL_QUALS);
     CommandArgument_Double_OrDefault(TUMOR_COV, -1.0);
     CommandArgument_Double_OrDefault(NORMAL_COV, -1.0);
     CommandArgument_Int_OrDefault(TUMOR_COV_FLOOR, 0);
     CommandArgument_Int_OrDefault(NORMAL_COV_FLOOR, 0);
     EndCommandArguments;

     // Load the data.

     vecqualvector TQ(TUMOR_QUALS), NQ(NORMAL_QUALS);
     ForceAssertEq( TQ.size()/4, NQ.size()/4 );

     // Compute coverage.

     double tumor_coverage_mean = 0, normal_coverage_mean = 0;
     vecqvec::size_type n = TQ.size( ) / 4;
     for ( vecqvec::size_type i = 0; i < 4*n; i++ )
     {    tumor_coverage_mean += TQ[i].size( );
          normal_coverage_mean += NQ[i].size( );    }
     tumor_coverage_mean /= double(n), normal_coverage_mean /= double(n);
     PRINT2( tumor_coverage_mean, normal_coverage_mean );

     // Simulate calling of somatic mutations.

     vec<String> alleles;
     vec<double> likelihoods;
     int have_power = 0;
     char refbase = 'A', altbase = 'C';
     int count = 0;
     for ( vecqvec::size_type i = 0; i < n; i++ )
     {    genotype_likelihoods TG, NG;
          int tumor_coverage = 0, normal_coverage = 0;
          for ( vecqvec::size_type j = 0; j < 4; j++ )
          {    for ( qvec::size_type k = 0; k < TQ[ 4*i + j ].size( ); k++ )
               {    char call = as_base( randomx( ) % 2 );
                    genotype_likelihoods TGk( 'A', call, TQ[ 4*i + j ][k] );
                    TG += TGk;
                    ++tumor_coverage;    }    }
          if ( tumor_coverage < TUMOR_COV_FLOOR ) continue;

          // Add coverage to tumor.

          if ( TUMOR_COV >= 0 )
          {    if ( TUMOR_COV > tumor_coverage_mean )
               {    int tc = tumor_coverage;
                    double extras = double(tc)
                         * (TUMOR_COV - tumor_coverage_mean) / tumor_coverage_mean;
                    for ( int u = 0; u < ceil(extras); u++ )
                    {    if ( u == ceil(extras) - 1 )
                         {    double r = double( randomx( ) % 1000000 ) / 1000000.0;
                              if ( extras - floor(extras) < r ) break;    }
                         qvec::size_type k = randomx( ) % tc;
                         for ( int j = 0; j < 4; j++ )
                         {    if ( k < TQ[ 4*i + j ].size( ) )
                              {    char call = as_base( randomx( ) % 2 );
                                   genotype_likelihoods
                                        TGk( 'A', call, TQ[ 4*i + j ][k] );
                                   TG += TGk;
                                   ++tumor_coverage;
                                   break;    }
                              k -= TQ[ 4*i + j ].size( );    }    }    }
               else
               {    cout << "Reduction in tumor coverage not implemented.\n";
                    exit(1);    }    }

          for ( int j = 0; j < 4; j++ )
          {    for ( qvec::size_type k = 0; k < NQ[ 4*i + j ].size( ); k++ )
               {    genotype_likelihoods NGk( 'A', 'A', NQ[ 4*i + j ][k] );
                    NG += NGk;
                    ++normal_coverage;    }    }
          if ( normal_coverage < NORMAL_COV_FLOOR ) continue;
          ++count;

          // Add coverage to normal.

          if ( NORMAL_COV >= 0 )
          {    if ( NORMAL_COV > normal_coverage_mean )
               {    int nc = normal_coverage;
                    double extras = double(nc)
                         * (NORMAL_COV - normal_coverage_mean)/normal_coverage_mean;
                    for ( int u = 0; u < ceil(extras); u++ )
                    {    if ( u == ceil(extras) - 1 )
                         {    double r = double( randomx( ) % 1000000 ) / 1000000.0;
                              if ( extras - floor(extras) < r ) break;    }
                         qvec::size_type k = randomx( ) % nc;
                         for ( int j = 0; j < 4; j++ )
                         {    if ( k < NQ[ 4*i + j ].size( ) )
                              {    char call = as_base( randomx( ) % 2 );
                                   genotype_likelihoods
                                        NGk( 'A', 'A', NQ[ 4*i + j ][k] );
                                   NG += NGk;
                                   ++normal_coverage;
                                   break;    }
                              k -= NQ[ 4*i + j ].size( );    }    }    }
               else
               {    cout << "Reduction in normal coverage not implemented.\n";
                    exit(1);    }    }

          TG.GetSorted( alleles, likelihoods );
          int AA = Position( alleles, String("AA") );
          int AC = Position( alleles, String("AC") );
          int CC = Position( alleles, String("CC") );
          double AC_CC =
               log10( pow( 10, likelihoods[AC] ) + pow( 10, likelihoods[CC] ) );
          double tumor_score = AC_CC - likelihoods[AA];

          NG.GetSorted( alleles, likelihoods );
          AA = Position( alleles, String("AA") );
          AC = Position( alleles, String("AC") );
          CC = Position( alleles, String("CC") );
          AC_CC = log10( pow( 10, likelihoods[AC] ) + pow( 10, likelihoods[CC] ) );
          double normal_score = likelihoods[AA] - AC_CC;

          cout << "\n[" << i << "], " << "tumor cov = " << tumor_coverage
               << ", tumor score = " << tumor_score
               << ", normal cov = " << normal_coverage
               << ", normal score = " << normal_score << "\n";

          if ( tumor_score >= 6.3 && normal_score >= 2.3 ) ++have_power;    }

     cout << "\n";
     PRINT2( count, have_power );
     cout << PERCENT_RATIO( 3, have_power, count ) << endl;    }
