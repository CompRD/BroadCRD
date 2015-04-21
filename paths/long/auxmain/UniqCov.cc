///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// For a sample of 10000 21-mers x that occur exactly once in the human reference 
// sequence, determine the number of times that x occurs in two NA12878 datasets.
// Possibly needed as part of analyses for DISCOVAR manuscript.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "Bitvector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "math/Functions.h"
#include "random/Random.h"

int main( )
{    RunTime( );
     vecbasevector ref( "/wga/scr4/bigrefs/human19/genome.fastb" );
     vecbitvector amb( "/wga/scr4/bigrefs/human19/genome.lookup.fastamb" );
     vec<basevector> x21;
     const int K = 21;
     for ( int i = 0; i < (int) ref.size( ); i++ )
     for ( int j = 0; j <= ref[i].isize( ) - K; j++ )
     {    Bool un = False;
          for ( int l = 0; l < K; l++ )
               if ( amb[i][j+l] ) un = True;
          if (un) continue;
          basevector b( ref[i], j, K );
          x21.push_back(b);    }
     ParallelSort(x21);
     const int sample = 10000;
     vec<basevector> samples;
     while( samples.isize( ) < sample )
     {    int64_t i = big_random( ) % x21.size( );
          if ( i > 0 && x21[i] == x21[i-1] ) continue;
          if ( i < x21.isize( ) - 1 && x21[i] == x21[i+1] ) continue;
          if ( Member( samples, x21[i] ) ) continue;
          basevector b = x21[i];
          b.ReverseComplement( );
          if ( LowerBound( x21, b ) < UpperBound( x21, b ) ) continue;
          samples.push_back( x21[i] );    }
     Sort(samples);

     vec<vec<int>> count( 2, vec<int>( samples.size( ), 0 ) );

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
                              int low = LowerBound( samples, b );
                              int high = UpperBound( samples, b );
                              if ( low < high )
                              {
                                   #pragma omp critical
                                   {    for ( int m = low; m < high; m++ )
                                             count[spass][m]++;    
                                                  }    }    }    }    }
               bases.clear( );
               if ( !in.fail( ) ) goto scan;    }    }

     const int flank = 10;
     for ( int i = 0; i < samples.isize( ); i++ )
     {    for ( int l = 0; l < flank; l++ )
               cout << as_base( samples[i][l] );
          cout << "/" << as_base( samples[i][flank] ) << "/";
          for ( int l = 0; l < flank; l++ )
               cout << as_base( samples[i][l+flank+1] );
          cout << "  " << count[0][i] << "," << count[1][i] << "\n";    }

     cout << "\n";
     vec<double> m(2);
     for ( int i = 0; i < 2; i++ )
     {    m[i] = Mean( count[i] );
          cout << "mean coverage for dataset " << i+1 << " = " << m[i] << endl;    }

     for ( int spass = 0; spass < 3; spass++ )
     {    vec<int> depth(1000, 0);

          if ( spass < 2 )
          {    for ( int j = 0; j < samples.isize( ); j++ )
               {    double d = count[spass][j]/m[spass];
                    depth[ int(round(10*d)) ]++;    }    }
          else
          {    for ( int j = 0; j < samples.isize( ); j++ )
               {    double d = ( count[0][j] + count[1][j] ) / ( m[0] + m[1] );
                    depth[ int(round(10*d)) ]++;    }    }

          cout << "\n";
          double n = Sum(depth);
          int top = 20;
          for ( int i = 0; i < top; i++ )
          {    if ( depth[i] > 0 )
               {    cout << i/10.0 << "-" << (i+1)/10.0 << "  " << depth[i] 
                         << "  " << 100 * depth[i]/n << "%" << endl;    }    }
          double tail = 0;
          for ( int i = top; i < depth.isize( ); i++ )
               tail += depth[i];
          cout << top/10.0 << "-infinity  " << tail << "  " << 100 * tail/n 
               << "%" << endl;    }    }
