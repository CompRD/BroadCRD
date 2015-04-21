///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Code used only a couple times, to generate the coordinates of NA12878 Fosmid
// regions.

#include "FastIfstream.h"
#include "MainTools.h"

int main( )
{
     // String lib = "Solexa-127359"; // hpool2
     String lib = "Solexa-127365"; // hpool3

     fast_pipe_ifstream in( "samtools view "
          "/wga/scr4/picard/A2925/C1-508_2012-11-12_2012-11-14/1/"
          + lib + "/A2925.1.aligned.duplicates_marked.bam | Col 3 4" );

     String line;
     vec< pair<String,int> > hits;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          istringstream iline( line.c_str( ) );
          String x;
          int y;
          iline >> x >> y;
          hits.push( x, y );    }

     const int max_dist = 300;
     const int min_group = 15000;

     for ( int i = 0; i < hits.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < hits.isize( ); j++ )
          {    if ( hits[j].first != hits[i].first ) break;
               if ( hits[j].second - hits[j-1].second > max_dist ) break;    }
          int count = j - i;
          if ( count >= min_group )
          {    int start = hits[i].second, stop = hits[j-1].second;
               double cov = double(count) / double(stop-start);
               cout << hits[i].first << "." << start << "-" << stop
                    << " [len=" << stop - start << ",reads=" << count
                    << ",cov=" << cov << "]" << endl;    }
          i = j - 1;    }    }
