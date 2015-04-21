///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Horrible one-off code to generate some plots.  They show coverage at the first
// exon of the RB1 gene in two human NA12878 data sets, one from the 1000 Genomes 
// Project and one consisting of fragment reads in our human assembly in the PNAS
// paper.  GC content is also shown.

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "dna/Bases.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"

// % cd /wga/scr1/ALLPATHS/H.sapiens.NA12878/bams
// % build_mini_bam chr13:47200000-48300000 wiggle

int main( int argc, char *argv[] )
{
     RunTime( );

     vecbasevector ref;
     ref.ReadOne( 
          "/wga/scr1/ALLPATHS/H.sapiens.NA12878/other_ref/build18.fastb", 13 );

     String TMP = "tmp";
     Mkdir777(TMP);

     Ofstream( out, "hhh" );
     for ( int pass = 1; pass <= 2; pass++ )
     {    String head;

          if ( pass == 1 )
          {    head = "/humgen/1kg/DCC/ftp/data/NA12878/alignment/NA12878.chrom13.ILLUMINA.bwa.CEU.high_coverage.20100311";    }
          if ( pass == 2 )
          {    head = "/wga/scr1/ALLPATHS/H.sapiens.NA12878/bams/wiggle";    }

          int start = 0, stop = 0;
          if ( pass == 1 ) // hg19
          {    start = 48877882; 
               stop = 48878185;    }
          if ( pass == 2 ) // hg18;
          {    start = 47775883;
               stop = 47776186;    }
          start -= 1000;
          stop += 1000;
          int center = ( start + stop ) / 2;
          int radius = ( stop - start ) / 2;
          String tig;
          if ( pass == 1 ) tig = "13";
          else tig = "chr13";

          // Get mean coverage for 1 Mb window.

          double meancov = 0.0;
          int big_radius = 500000;
          SystemSucceed( "/seq/dirseq/samtools/current/samtools view " + head 
               + ".bam \"" + tig + "\":" + ToString(center-big_radius) + "-" 
               + ToString(center+big_radius) + " | SAM2CRDDump OUT_HEAD=" + TMP 
               + "/regionx SEP=-15 DEV=10 NH=True" );
          vec<look_align> big_aligns;
          LoadLookAligns( "tmp/regionx.qltout", big_aligns );
          int64_t total_cov = 0;
          for ( int i = 0; i < big_aligns.isize( ); i++ )
          {    const look_align& la = big_aligns[i];
               int start = la.pos2( ) - (center-big_radius);
               int stop = la.Pos2( ) - (center-big_radius);
               for ( int j = Max(0, start); j < Min(big_radius*2, stop); j++ )
                    total_cov++;    }
          double mean_cov = double(total_cov) / double(2*big_radius);
          PRINT(mean_cov);

          SystemSucceed( "/seq/dirseq/samtools/current/samtools view " + head 
               + ".bam \"" + tig + "\":" + ToString(center-radius) + "-" 
               + ToString(center+radius) + " | SAM2CRDDump OUT_HEAD=" + TMP 
               + "/regionx SEP=-15 DEV=10 NH=True" );

          vec<look_align> aligns;
          LoadLookAligns( "tmp/regionx.qltout", aligns );
          vec<int> cov( radius*2, 0 );
          for ( int i = 0; i < aligns.isize( ); i++ )
          {    const look_align& la = aligns[i];
               int start = la.pos2( ) - (center-radius);
               int stop = la.Pos2( ) - (center-radius);
               for ( int j = Max(0, start); j < Min(radius*2, stop); j++ )
                    ++cov[j];    }
          for ( int i = 0; i < radius*2; i++ )
          {    out << i << " " << double(cov[i]) / mean_cov;
               if ( pass == 1 ) out << " 1 0 0\n";
               else out << " 0 0 1\n";    }    

          // Get GC composition.

          if ( pass == 2 )
          {    int window = 100;
               for ( int i = 0; i < radius*2; i++ )
               {    int gc = 0;
                    for ( int j = 0; j < window; j++ )
                    {    char c = as_base( ref[0][ 
                              center + i - radius + j - window/2 ] );
                         if ( c == 'G' || c == 'C' ) gc++;    }
                    out << i << " " << double(gc)/double(window) << " 0 1 0\n";    }
                    }    }

     // Generate plots.

     flush(out);
     SystemSucceed( "PlotPoints IN=hhh OUT=hhh.png CONNECT=True POINTSIZE=0.1 " 
          "X_AXIS_EXTEND=0 Y_AXIS_EXTEND=0 X_AXIS_OFFSET=0 " 
          "COLOR_BY_POINT=True MIN_Y=0 MAX_Y=1.7 POSTSCALE=0.2" );
     System( "cat hhh | grep ' 1 0 0' > hhh1" );
     SystemSucceed( "PlotPoints IN=hhh1 OUT=hhh1.png CONNECT=True POINTSIZE=0.1 " 
          "X_AXIS_EXTEND=0 Y_AXIS_EXTEND=0 X_AXIS_OFFSET=0 " 
          "COLOR_BY_POINT=True MIN_Y=0 MAX_Y=1.7 POSTSCALE=0.2" );
     System( "cat hhh | egrep ' 1 0 0| 0 0 1' > hhh2" );
     SystemSucceed( "PlotPoints IN=hhh2 OUT=hhh2.png CONNECT=True POINTSIZE=0.1 " 
          "X_AXIS_EXTEND=0 Y_AXIS_EXTEND=0 X_AXIS_OFFSET=0 " 
          "COLOR_BY_POINT=True MIN_Y=0 MAX_Y=1.7 POSTSCALE=0.2" );
     System( "cat hhh | grep ' 0 1 0' > hhh3" );
     SystemSucceed( "PlotPoints IN=hhh3 OUT=hhh3.png CONNECT=True POINTSIZE=0.1 " 
          "X_AXIS_EXTEND=0 Y_AXIS_EXTEND=0 X_AXIS_OFFSET=0 " 
          "COLOR_BY_POINT=True MIN_Y=0 POSTSCALE=0.2" );    }
