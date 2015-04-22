///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// XFP: estimate variant calling false positives using male vs female X.
// Currently hardwired to use GATK calls from 100 or 250 base reads, and
// Cortex calls.

#include "FastIfstream.h"
#include "MainTools.h"
#include "TokenizeString.h"
#include "math/Functions.h"
#include "random/Random.h"

int main(int argc, char** )
{    RunTime( );

     Bool dump = False;
     if ( argc > 1 ) dump = True;

     const int MIN_GQ = 0;
     const double MIN_QUAL_SNP = 180;
     const double MIN_QUAL_INDEL = 100;
     String FILTER_PLUS = "PASS";
     // String FILTER_PLUS = "{PASS,99.00to,99.90to}";
     // String FILTER_PLUS = "{PASS,99.90to}";
     // String FILTER_MINUS = "Low";
     String FILTER_MINUS = "";
     vec<String> filter_plus, filter_minus;
     ParseStringSet( FILTER_PLUS, filter_plus );
     ParseStringSet( FILTER_MINUS, filter_minus );

     String female = "NA12878";
     String male = "NA12891"; // NA12878's dad
     // String male = "NA12877"; // NA12878's husband
     // String male = "NA12889"; // NA12878's husband's dad
     // String male = "NA12882"; // a son
     // String male = "NA12883"; // a son
     // String male = "NA12884"; // a son
     // String male = "NA12886"; // a son
     // String male = "NA12888"; // a son
     // String male = "NA12893"; // a son

     String calls = "/wga/scr4/NA12878_calls";

     for ( int c = 1; c <= 4; c++ )
     {    String caller;
          int rl;
          if ( c == 1 ) 
          {    caller = "GATK";
               rl = 100;    }
          if ( c == 2 ) 
          {    caller = "GATK";
               rl = 250;    }
          if ( c == 3 ) 
          {    caller = "CORTEX";
               rl = 250;    }
          if ( c == 4 ) 
          {    caller = "DISCOVAR";
               rl = 250;    }
          cout << "\n";
          PRINT2( caller, rl );

          String female_250, male_250;

          if ( caller == "GATK" )
          {    female_250 = calls + "/mem4/haplotype-caller-pcr-none/"
                    "reverted.12.aligned.wholegenome.sorted.indel_cleaned_local."
                    "recal.unfiltered.recal_snp_recal_indel.vcf";
               male_250 = calls + "/H03N7ADXX+H05F1ADXX/gatk/3.hc-no.pcr/"
                    "Solexa-135851.aligned.sorted.merged.indel_cleaned_local."
                    "recal.unfiltered.recal_snp_recal_indel.vcf";    }

          if ( caller == "CORTEX" )
          {    female_250 = calls + "/cortex/"
                    "cortex-NA12878-again/vcf/NA12878.decomp.processed.vcf";
               male_250 = calls + "/cortex/"
                    "cortex-NA12891-again/vcf/NA12891.decomp.processed.vcf";    }

          if ( caller == "DISCOVAR" )
          {    female_250 = "/wga/scr4/human_assemblies/1/v9/v9_combined.filtered.vcf";
               male_250 = 
                    "/wga/scr4/human_assemblies/5/v2/v2_chrX_combined.filtered.vcf";    }

          const double div = 1.8;

          const int start = 10000000;
          const int stop = 110000000;
          vec<vec<double>> hets( 2, vec<double>(100,0) );
          for ( int pass = 0; pass < 2; pass++ )
          {    String sample = ( pass == 0 ? male : female );
               String vcf = "/wga/scr4/human_data/CEPH/" + sample + "_S1.genome.vcf";
               if ( rl == 250 )
               {    if ( pass == 0 ) vcf = male_250;
                    else vcf = female_250;    }
               fast_ifstream in(vcf);
               String line;
               vec<String> fields;
               while(1)
               {    getline( in, line );
                    if ( in.fail( ) ) break;
                    if ( rl == 100 )
                    {    if ( !line.Contains( "chrX" , 0 ) ) continue;    }
                    else
                    {    if ( !line.Contains( "X" , 0 ) ) continue;    }
                    Bool OK = False;
     
                    for ( int i = 0; i < filter_plus.isize( ); i++ )
                         if ( line.Contains( filter_plus[i] ) ) OK = True;
     
                    for ( int i = 0; i < filter_minus.isize( ); i++ )
                         if ( line.Contains( filter_minus[i] ) ) OK = False;

                    if ( !OK ) continue;
                    line.GlobalReplaceBy( "|", "/" );
                    if ( !line.Contains( "/" ) ) continue;

                    String p1 = line.Before( "/" );
                    if ( p1.Contains( "\t" ) ) p1 = p1.RevAfter( "\t" );
                    // if ( p1.Contains( " " ) ) p1 = p1.RevAfter( " " );
                    String p2 = line.After( "/" ).Before( ":" );
                    if ( p1 == p2 ) continue;

                    Tokenize( line, {'\t'}, fields );
     
                    if ( caller == "GATK" && rl == 250 )
                    {    vec<String> gfields;
                         Tokenize( fields[9], {':'}, gfields );
                         int GQ = gfields[3].Int( );
                         if ( GQ < MIN_GQ ) continue;
                         Bool SNP 
                              = ( fields[3].size( ) == 1 && fields[4].size( ) == 1 );
                         if ( SNP && fields[5].Double( ) < MIN_QUAL_SNP ) continue;
                         if ( !SNP && fields[5].Double( ) < MIN_QUAL_INDEL ) 
                              continue;    }
     
                    int pos = fields[1].Int( );
                    if ( pos < start || pos >= stop ) continue;

                    if ( dump ) cout << ((pass==0)?"XY ":"XX ") << caller << "-" << rl << " " << pos << endl;

                    hets[pass][ (pos-start) / 1000000 ]++;    }    }
          const int bootstrap = 1000;
          vec<double> FP;
          for ( int pass = 0; pass < bootstrap; pass++ )
          {    int m = 0, f = 0; 
               for ( int j = 0; j < 100; j++ )
               {    int s = randomx( ) % 100;
                    m += hets[0][s]; 
                    f += hets[1][s];    }
               FP.push_back( 100.0 * double(m) / ( div * double(f) ) );    }
          cout << "male hets = " << int(round(Sum(hets[0]))) << endl;
          cout << "female hets = " << int(round(Sum(hets[1]))) << endl;
          cout << "FP rate = " << PERCENT_RATIO( 3, Sum(hets[0]), div * Sum(hets[1]) )
               << " +/- " << StdDev( FP, Mean(FP) ) << endl;    }    }
