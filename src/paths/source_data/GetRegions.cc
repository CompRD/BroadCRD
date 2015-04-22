///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// This code makes trunk/references/Homo_sapiens/WIBR_Fosmid_Pool.fastb.  The ends
// of the reference are trimmed back to avoid including vector in read sets.  The
// amount is specified in LongProtoTools.h.
//
// You need to first run the code shown below.  This creates an edited sequence
// for the entire genome, which is stored in a temporary location.  It takes a 
// couple of hours.


/*

crd8% cd /wga/dev/jaffe
crd8% mkdir GATK2
crd8% cd GATK2
crd8% (got GenomeAnalysisTK-2.1-12.tar.bz2 from the GATK site on the web)
crd8% bunzip2 GenomeAnalysisTK-2.1-12.tar.bz2
crd8% cat GenomeAnalysisTK-2.1-12.tar | tar xf -

crd8% cd /local/scratch/jaffe/BroadCRD/fixed

crd8% java -jar /wga/dev/jaffe/GATK2/GenomeAnalysisTK-2.1-12-ga99c19d/GenomeAnalysisTK.jar -R /wga/scr2/bigrefs/human19/genome.fasta -T UnifiedGenotyper -I frag.list -o raw2.vcf -baq CALCULATE_AS_NECESSARY -nt 48

crd8% java -jar /wga/dev/jaffe/GATK2/GenomeAnalysisTK-2.1-12-ga99c19d/GenomeAnalysisTK.jar -R /wga/scr2/bigrefs/human19/genome.fasta -T SelectVariants --variant raw2.vcf -o select.vcf -restrictAllelesTo BIALLELIC

crd8% java -jar /wga/dev/jaffe/GATK2/GenomeAnalysisTK-2.1-12-ga99c19d/GenomeAnalysisTK.jar -R /wga/scr2/bigrefs/human19/genome.fasta -T VariantFiltration -V select.vcf -o var.vcf --filterExpression "QD<5.0||AC<2||DP<6" --filterName junk

crd8% cat var.vcf | grep -v junk > var.clean.vcf

crd8% r43253:EditRefUsingVcf SAMPLE="High GC Fosmid Pool" VCFS=var.clean.vcf     \
      REF_IN=/wga/scr2/bigrefs/human19/genome.fasta REF_OUT=genome_fixed.fasta

crd8% Fastb PRE=. FILE=genome_fixed.fasta

*/


#include "Basevector.h"
#include "FastIfstream.h"
#include "FetchReads.h"
#include "MainTools.h"
#include "paths/long/LoadAndCorrect.h"

int main( )
{    
     RunTime( );

     // String pool = "WIBR";
     String pool = "NA12878";

     vec<String> regions;
     fast_ifstream in( "/wga/dev/references/Homo_sapiens/" + pool 
          + "_Fosmid_Pool.regions" );
     String line;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          regions.push_back(line);    }
     // vecbasevector genome( "/wga/scr2/bigrefs/human19/genome.lookup.fastb" );
     vecbasevector genome( 
          "/local/scratch/jaffe/BroadCRD/fixed/genome_fixed.fastb" );
     vecbasevector X;
     for ( int i = 0; i < regions.isize( ); i++ )
     {    String chr = regions[i].Before( ":" );
          int g = ( chr == "X" ? 22 : chr.Int( ) - 1 );
          int start = regions[i].Between( ":", "-" ).Int( ) + Fosmid_trim_back;
          int stop = regions[i].After( "-" ).Int( ) - Fosmid_trim_back;
          basevector x( genome[g], start, stop - start );
          X.push_back_reserve(x);    }

     // Load validated assemblies.

     if ( IsRegularFile( "/wga/dev/references/Homo_sapiens/" + pool 
          + "_Fosmid_Pool.validated.fasta" ) )
     {    vecbasevector validated;
          vecString vnames;
          FetchReads( validated, vnames, "/wga/dev/references/Homo_sapiens/" + pool 
               + "_Fosmid_Pool.validated.fasta" );
          for ( int i = 0; i < (int) validated.size( ); i++ )
          {    int id = vnames[i].Int( );
               X[id] = validated[i];    }    }

     // Write assemblies.

     X.WriteAll( "/wga/dev/references/Homo_sapiens/" + pool 
          + "_Fosmid_Pool.fastb" );    }
