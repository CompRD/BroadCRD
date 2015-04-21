///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Compute human assembly N50 as a function of %PF.
// Should add our favorite flowcell and the one Illumina sent.
//
// - N50PF > ddd
//
// - cat ddd | Col 12 9 | tr -d ',' > data
//
// - PlotPoints IN=data OUT=data.png PR=1 TITLE="N50 contig size in kb as a function of %PF" POINTSIZE=2 XMAX=0.9 MAX_Y=110

#include "FastIfstream.h"
#include "MainTools.h"
#include "TokenizeString.h"

int main( )
{
     vec<String> sx = { "HG00096:H7AGFADXX", "HG00268:H7A1KADXX",
          "HG00419:H7A26ADXX", "HG00759:H7A29ADXX", "HG01051:H7A22ADXX",
          "HG01112:H7A2BADXX", "HG01500:H7A83ADXX", "HG01565:H7A63ADXX",
          "HG01583:H7A2NADXX", "HG01595:H7A7BADXX", "HG01879:H7A2JADXX",
          "HG02568:H7A2DADXX", "HG02922:H7A2VADXX", "HG03006:H7CLPADXX",
          "HG03052:H7B6EADXX", "HG03642:H7A2EADXX", "HG03742:H7A81ADXX",
          "NA18525:H7A33ADXX", "NA18939:H781JADXX", "NA19017:H7A7YADXX",
          "NA19625:H7A1NADXX", "NA19648:H7AGDADXX", "NA20502:H7A1WADXX",
          "NA20845:H7A30ADXX", "NA12878:H01UJADXX" };

     for ( int i = 0; i < sx.isize( ); i++ )
     {    String sample = sx[i].Before( ":" ), fc = sx[i].After( ":" );
          String id = sample;
          if ( sample == "NA12878" ) id = "50642";
          if ( !IsRegularFile( "GapToy." + id ) ) continue;
          fast_ifstream in1( "GapToy." + id );
          String line;
          double N50 = -1;
          while(1)
          {    getline( in1, line );
               if ( in1.fail( ) ) break;
               if ( line.Contains( "contig line N50" ) )
                    N50 = line.After( "contig line N50 = " ).Double( )/1000;    }
          if ( N50 < 0 ) continue;
          fast_pipe_ifstream in2( "find -L /seq/picard/" + fc 
               + " -name \"*.alignment_summary_metrics\" -print" );
          vec<String> lines;
          while(1)
          {    getline( in2, line );
               if ( in2.fail( ) ) break;
               lines.push_back(line);    }
          ForceAssertEq( lines.isize( ), 2 );
          vec<double> pfs;
          for ( int j = 0; j < 2; j++ )
          {    fast_ifstream in( lines[j] );
               while(1)
               {    getline( in, line );
                    if ( in.fail( ) ) break;
                    if ( !line.Contains( "FIRST_OF_PAIR", 0 ) ) continue;
                    vec<String> x;
                    Tokenize( line, '\t', x );
                    pfs.push_back( x[3].Double( ) );    }    }
          double mpf = Mean(pfs);
          PRINT4( sample, fc, N50, mpf );    }    }
