///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// LookMisses.  Look at some HapMap3 variants present in GATK-250 but not found
// by DISCOVAR.  Assemble the 20 kb region centered at the site and see if it it
// appears to contain the variant.

#include "FastIfstream.h"
#include "MainTools.h"
#include "VecUtilities.h"

int main( )
{
     RunTime( );

     const int pcount = 20;
     const int sep = 2500;

     fast_ifstream in1( "/wga/scr4/NA12878_calls/genotype_concordance/"
          "our-concordance-calculations/dataset1.v5/dataset1.v5.misses" );
     String line, chr, ref, alt, junk;
     int pos;
     int count = 0;
     vec< pair< pair<String,int>, pair<String,String> > > misses;
     vec<String> mlines;
     vec<int> mcounts;
     while(1)
     {    getline( in1, line );
          if ( in1.fail( ) ) break;
          if ( ++count % sep != 0 ) continue;
          istringstream iline( line.c_str( ) );
          iline >> chr >> pos >> junk >> ref >> alt;
          if ( alt.Contains( "," ) )
          {    count--;
               continue;    }
          misses.push( make_pair( chr, pos ), make_pair( ref, alt ) );
          mlines.push_back(line), mcounts.push_back(count);
          if ( misses.isize( ) == pcount ) break;    }
     SortSync( misses, mlines, mcounts );

     vec<Bool> gatk( pcount, False );
     String gatk_250 = "/wga/scr4/NA12878_calls/mem4/haplotype-caller-pcr-none/"
          "reverted.12.aligned.wholegenome.sorted.indel_cleaned_local.recal."
          "unfiltered.recal_snp_recal_indel.vcf";

     fast_ifstream in2(gatk_250);
     while(1)
     {    getline( in2, line );
          if ( in2.fail( ) ) break;
          istringstream iline( line.c_str( ) );
          iline >> chr >> pos >> junk >> ref >> alt;
          pair< pair<String,int>, pair<String,String> > x;
          if ( alt.Contains( "," ) ) continue;
          x = make_pair( make_pair(chr,pos), make_pair(ref,alt) );
          int p = BinPosition( misses, x );
          if ( p >= 0 ) gatk[p] = True;    }

     for ( int i = 0; i < pcount; i++ )
     {    if ( !gatk[i] ) continue;
          cout << "\nmisses2." << sep * (i+1) << endl;
          cout << mlines[i] << endl;
          String chr = misses[i].first.first;
          int pos = misses[i].first.second;

          const int flank = 10000;
          SystemSucceed( "LongProto SAMPLE=human READS=#picard TMP=tmp.xxx "
               "OUT_INT_HEAD=/wga/dev/jaffe/BroadCRD/aaa "
               "LOGGING=REFTRACE_VARIANTS=True "
               "X=" + chr + ":" + ToString(pos-flank) + "-" + ToString(pos+flank)
               + " > xxx.woof" );
          System( "cat /wga/dev/jaffe/BroadCRD/aaa.final.variant | grep "
               + ToString(pos) );    }    }
               
