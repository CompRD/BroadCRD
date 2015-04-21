///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "paths/HyperBasevector.h"

void MakeLexLookup( const HyperBasevector& hb, vec< pair<int,int> >& locs )
{    vec<int64_t> starts;
     starts.push_back(0);
     int K = hb.K( );
     for ( int e = 0; e < hb.E( ); e++ )
     {    const basevector& E = hb.EdgeObject(e);
          starts.push_back( starts.back( ) + Max( 0, E.isize( ) - K + 1 ) );    }
     locs.resize( starts.back( ) );
     #pragma omp parallel for
     for ( int e = 0; e < hb.E( ); e++ )
     {    const basevector& u = hb.EdgeObject(e);
          for ( int j = 0; j <= u.isize( ) - K; j++ )
          {    int64_t r = starts[e] + j;
               locs[r].first = e, locs[r].second = j;    }    }
     double clock = WallClockTime( ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     ParallelSort( locs, [K, &hb](const pair<int,int>& x1, const pair<int,int>& x2)
     {    auto itr1 = hb.EdgeObject(x1.first).begin(x1.second);
          auto itr2 = itr1+hb.K();
          auto itr3 = hb.EdgeObject(x2.first).begin(x2.second);
          auto itr4 = itr3+hb.K();
          return std::lexicographical_compare(itr1,itr2,itr3,itr4);    }
               );    
     cout << TimeSince(clock) << " used in sort" << endl; // XXXXXXXXXXXXXXXXXXXXXXX

}

int main( )
{

     String INSTANCE = "50718";

     String SAMPLE = "NA12878";
     String work_dir = "/wga/scr4/jaffe/GapToy/" + INSTANCE;

     HyperBasevector hb;
     BinaryReader::readFile( work_dir + "/a.final/a.hbv", &hb );

     vec< pair<int,int> > locs;
     MakeLexLookup( hb, locs );    }

