/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// CancerCopy.  Experimental code to look at cancer copy number.

#include "Basevector.h"
#include "Bitvector.h"
#include "FeudalMimic.h"
#include "Intvector.h"
#include "MainTools.h"
#include "PackAlign.h"
#include "feudal/BinaryStream.h"
#include "lookup/LookupTable.h"
#include "math/Functions.h"
#include "random/Poisson.h"
#include "solexa/SolexaPipeline.h"

class cov_record {

     public:

     String chr;
     double start;
     int N1, N2;
     double ratio_low, ratio, ratio_high;
     double gcp;

     cov_record( const String& chr, const double start, const int N1, const int N2,
          const double ratio_low, const double ratio, const double ratio_high, 
          const double gcp ) : chr(chr), start(start), N1(N1), N2(N2), 
          ratio_low(ratio_low), ratio(ratio), ratio_high(ratio_high), gcp(gcp) { }

     void Print( ostream& out, const int prec )
     {    out << chr << ":" 
               << setiosflags(ios::fixed) << setprecision(prec) << setw(6)
               << double(start)/1000000 << resetiosflags(ios::fixed) << "M,  "
               << "N = " << setiosflags(ios::fixed) << setw(5)
               << N1 << resetiosflags(ios::fixed)
               << " / " << setiosflags(ios::fixed) << setw(5)
               << N2 << resetiosflags(ios::fixed) << ",  "
               << "ratio = " << setiosflags(ios::fixed) << setprecision(2)
               << setw(6) << ratio_low << resetiosflags(ios::fixed) << " -"
               << setiosflags(ios::fixed) << setprecision(2)
               << setw(5) << ratio << resetiosflags(ios::fixed) << " -"
               << setiosflags(ios::fixed) << setprecision(2)
               << setw(5) << ratio_high << resetiosflags(ios::fixed) << ",  "
               << "GC = " << setiosflags(ios::fixed) << setprecision(1)
               << setw(4) << gcp << resetiosflags(ios::fixed) 
               << "%";
          if ( ratio_high < 1.0 ) out << " LOW";
          if ( ratio_low > 1.0 ) out << " HIGH";
          out << "\n";    }

};

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_Int(WINDOW);
     CommandArgument_Double_OrDefault(GENOME_SIZE_MULTIPLIER, -1.0);
     CommandArgument_Int_OrDefault(OFFSET,0);
     EndCommandArguments;

     // #1 = tumor
     // #2 = paired normal

     Bool female = True;

     const int readlength = 76;

     String rdir = SOLEXAPIPELINE;
     String gdir = "/wga/scr1/LOOKUP/human36";
     
     vecbasevector genome( gdir + "/hg18.fastb" );
     vecbitvector ambig( gdir + "/hg18.fastamb" );
     lookup_table look( gdir + "/hg18.lookup" );
     VecIntVec cov1, cov2;
     Mimic( genome, cov1 );
     Mimic( genome, cov2 );

     // Load placement marks.

     vec<placement_mark> M1, M2;
     String mdir = "/wga/scr1/cancer_Illumina_WGS/marks_tmp";
     vec<String> all = AllFiles(mdir);
     for ( int i = 0; i < all.isize( ); i++ )
     {    PRINT2( i, all[i] );
          vec<placement_mark> X;
          BinaryReader::readFile( mdir + "/" + all[i], &X );
          if ( all[i].Contains( "tumor", 0 ) ) M1.append(X);
          if ( all[i].Contains( "normal", 0 ) ) M2.append(X);    }
     cout << M1.size( ) << " placed reads from tumor sample" << endl;
     cout << M2.size( ) << " placed reads from paired normal\n" << endl;

     // Generate table of uniformly distributed random numbers.

     vec<double> unif;
     for ( int i = 0; i < 1000000; i++ )
          unif.push_back( drand48( ) );
     int uptr = 0;

     // Compute coverage.

     for ( int i = 0; i < M1.isize( ); i++ )
     {    int tig = M1[i].Tig( ), pos = M1[i].Pos( );
          for ( int j = 0; j < readlength; j++ )
               ++cov1[tig][pos+j];    }
     for ( int i = 0; i < M2.isize( ); i++ )
     {    int tig = M2[i].Tig( ), pos = M2[i].Pos( );
          for ( int j = 0; j < readlength; j++ )
               ++cov2[tig][pos+j];    }

     // Generate window records.

     vec<cov_record> C;
     for ( size_t tig = 0; tig < genome.size( ); tig++ )
     {    const basevector& g = genome[tig];
          const bitvector& a = ambig[tig];
          String chr = look.ContigName(tig);
          if ( female && chr == "chrY" ) continue;
          if ( chr.Contains( "random" ) ) continue;
          for ( int start = 0 + OFFSET; start <= g.isize( ) - WINDOW; start += WINDOW )
          {    int count1 = 0, count2 = 0, gc = 0, amb = 0;
               for ( int j = 0; j < WINDOW; j++ )
               {    count1 += cov1[tig][start+j];
                    count2 += cov2[tig][start+j];
                    if ( g[start+j] == 1 || g[start+j] == 2 ) ++gc;
                    if ( a[start+j] ) ++amb;    }
               count1 /= readlength;
               count2 /= readlength;

               // Ignore windows containing 50% or more ambiguous bases.

               if ( double(amb)/double(WINDOW) > 0.5 ) continue;

               // Compute ratio.

               double ratio = ( double(count1)/double(count2) )
                    * ( double( M2.size( ) ) / double( M1.size( ) ) );

               // Generate confidence interval for ratio.

               int sample = 10000;
               double ratio_low, ratio_high;
               PoissonRatioConfidenceInterval( count1, count2, 0.99, sample, 
                    unif, uptr, ratio_low, ratio_high );
               ratio_low *= double( M2.size( ) ) / double( M1.size( ) );
               ratio_high *= double( M2.size( ) ) / double( M1.size( ) );

               // Make record.

               double gcp = 100.0 * double(gc) / double(WINDOW);
               C.push( chr, double(start), count1, count2, 
                    ratio_low, ratio, ratio_high, gcp );    }    }    

     // Determine genome size multiplier.

     if ( GENOME_SIZE_MULTIPLIER < 0 )
     {    vec< pair<int,double> > goods_mult;
          for ( double mult = 1.0; mult <= 1.15; mult += 0.005 )
          {    int goods = 0;
               for ( int i = 0; i < C.isize( ); i++ )
               {    if ( mult * C[i].ratio_low <= 1.0 
                         && 1.0 <= mult * C[i].ratio_high ) 
                    {    ++goods;    }    }
               goods_mult.push( goods, mult );    }
          ReverseSort(goods_mult);
          GENOME_SIZE_MULTIPLIER = goods_mult[0].second;
          cout << "genome size multiplier = " << setprecision(4) 
               << GENOME_SIZE_MULTIPLIER << "\n\n";    }
     for ( int i = 0; i < C.isize( ); i++ )
     {    C[i].ratio_low *= GENOME_SIZE_MULTIPLIER;
          C[i].ratio *= GENOME_SIZE_MULTIPLIER;
          C[i].ratio_high *= GENOME_SIZE_MULTIPLIER;    }

     // Print windows.

     int prec;
     if ( WINDOW >= 1000000 ) prec = 0;
     else if ( 100000 <= WINDOW && WINDOW <= 500000 ) prec = 1;
     else prec = 2;
     for ( int i = 0; i < C.isize( ); i++ )
          C[i].Print( cout, prec );    }
