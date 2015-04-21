///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// For some positions that appear to be homozygous in Fosmid 40 and around it,
// compute the multiplicity of the corresponding 21-mer in two NA12878 data sets.
// Possibly useful for determing if Fosmid 40 is correctly placed.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "FastIfstream.h"
#include "FetchReads.h"
#include "MainTools.h"
#include "TokenizeString.h"

// Assemble Fosmid region 40 plus 20 kb to left.
//
// LongProto SAMPLE=human READS=#picard TMP=tmp.xxx X=16:15012318-15067606
//     OUT_INT_HEAD= ~/crd/ccc OUT_GENOME=genome.fastb

int main( )
{    RunTime( );
     vecbasevector ref( "genome.fastb" );
     vecbasevector A;
     FetchReads( A, 0, "/wga/dev/jaffe/BroadCRD/ccc.final.fasta" );
     vecbasevector bases( "tmp.xxx/frag_reads_orig.fastb" );

     vec<int> core;
     core.push_back( 100, 97, 52, 21, 143, 141, 4, 119, 116, 34 );
     core.push_back( 64, 63, 105, 104, 101, 59, 55, 26, 123, 74 );
     core.push_back( 37, 14, 7, 40, 44, 11, 128, 127, 18, 15 );
     core.push_back( 145, 144, 113, 48, 138, 76, 39, 90, 88 );
     core.push_back( 85, 82 );

     const int K = 21;
     vec<triple<basevector,int,int>> stuff;
     for ( int t = 0; t < core.isize( ); t++ )
     {    const basevector& x = A[ core[t] ];
          for ( int j = 0; j <= x.isize( ) - K; j++ )
          {    basevector b( x, j, K );
               stuff.push( b, t, j );    }    }
     Sort(stuff);
     vec<basevector> kmers;
     for ( int i = 0; i < stuff.isize( ); i++ )
          kmers.push_back( stuff[i].first );

     vec<vec<int>> count( 2, vec<int>( stuff.size( ), 0 ) );

     for ( int spass = 0; spass < 2; spass++ )
     {    vec<String> flowcell, picard_run, lanes, libs;
          if ( spass == 0 )
          {    flowcell.push_back( "H01UJADXX", "H01UJADXX" );
               picard_run.push_back( 
                    "C1-508_2012-11-01_2012-11-04", "C1-508_2012-11-01_2012-11-04" );
               lanes.push_back( "1", "2" );
               libs.push_back( "Solexa-125532", "Solexa-125532" );    }
          if ( spass == 1 )
          {    flowcell.push_back( "H06HDADXX", "H06HDADXX", "H06JUADXX" );
               picard_run.push_back( "C1-508_2013-01-10_2013-01-13",
                    "C1-508_2013-01-10_2013-01-13", "C1-508_2013-01-10_2013-01-13" );
               lanes.push_back( "1", "2", "1" );
               libs.push_back( "Solexa-135852", "Solexa-135852", "Solexa-135852" );
                    }
          // String R = "16"; // ***************************************************
          String R = "";
          String shead = "samtools view -X /seq/picard/";
          String dbam = ".aligned.duplicates_marked.bam";
          String look = " " + R;
          for ( int i = 0; i < flowcell.isize( ); i++ )
          {    fast_pipe_ifstream in( shead + flowcell[i] + "/" + picard_run[i] + "/"
                    + lanes[i] + "/" + libs[i] + "/" + flowcell[i]
                    + "." + lanes[i] + dbam + look );
               String line;
               vec<String> fields;
               vec<basevector> bases;
               scan:
               while(1)
               {    getline( in, line );
                    if ( in.fail( ) ) break;
                    Tokenize( line, {'\t'}, fields );
                    String s = fields[9];
                    s.GlobalReplaceBy( "N", "A" );
                    bases.push_back( basevector(s) );
                    if ( bases.size( ) == 10000 ) break;    }
               #pragma omp parallel for
               for ( int id = 0; id < bases.isize( ); id++ )
               {    basevector r = bases[id];
                    for ( int pass = 1; pass <= 2; pass++ )
                    {    if ( pass == 2 ) r.ReverseComplement( );
                         for ( int j = 0; j <= r.isize( ) - K; j++ )
                         {    basevector b( r, j, K );
                              int low = LowerBound( kmers, b );
                              int high = UpperBound( kmers, b );
                              if ( low < high )
                              {
                                   #pragma omp critical
                                   {    for ( int m = low; m < high; m++ )
                                             count[spass][m]++;    
                                                  }    }    }    }    }
               bases.clear( );
               if ( !in.fail( ) ) goto scan;    }    }

     vec<vec<vec<int>>> xcount( 2, vec<vec<int>>( core.isize( ) ) );
     for ( int spass = 0; spass < 2; spass++ )
     {    for ( int t = 0; t < core.isize( ); t++ )
               xcount[spass][t].resize( A[ core[t] ].isize( ) - K + 1, 0 );
          for ( int i = 0; i < count[0].isize( ); i++ )
          {    xcount[spass][ stuff[i].second ][ stuff[i].third ] 
                    = count[spass][i];    }    }
     for ( int i = 0; i < core.isize( ); i++ )
     {    cout << "\nCOVERAGE OF EDGE " << core[i] << endl;
          for ( int j = 0; j < xcount[0][i].isize( ); j++ )
          {    if ( xcount[0][i][j] < 60 || xcount[1][i][j] < 60 )
               {    cout << j << "  " << xcount[0][i][j] 
                         << "," << xcount[1][i][j] << "  ";
                    const int flank = 10;
                    for ( int l = 0; l < flank; l++ )
                         cout << as_base( A[ core[i] ][j+l] );
                    cout << "/" << as_base( A[ core[i] ][j+flank] ) << "/";
                    for ( int l = 0; l < flank; l++ )
                         cout << as_base( A[ core[i] ][j+l+flank+1] );
                    cout << "\n";    }    }    }    }
