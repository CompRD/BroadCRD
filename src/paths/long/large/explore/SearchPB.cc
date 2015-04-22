///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Find PacBio reads that bridge a gap, defined by two edges.  Look for a read
// that has a k-mer match on both side, where k = 28 and the kmer is unique in
// the genome.
//
// Hardcoded for the moment.

#include "MainTools.h"
#include "paths/RemodelGapTools.h"

int main( )
{    RunTime( );
     String dir = "/wga/scr4/jaffe/GapToy/51881.bigger/a.final";
     vecbasevector tigs( dir + "/a.fastb" );
     vec<int> inv;
     BinaryReader::readFile( dir + "/a.inv", &inv );
     int id1 = 6306734;
     int id2 = 6338016;
     for ( int p = 1; p <= 3; p++ )
     {    String fn = "/wga/scr4/human_data/NA12878/Schadt/chemistry_"
               + ToString(p) + ".fastb";
          int N = MastervecFileObjectCount(fn);
          for ( int z = 0; z < 10; z++ )
          {    PRINT2( p, z );
               int start = (z * N) / 10, stop = ((z+1) * N) / 10;
               vecbasevector pb;
               pb.ReadRange( "/wga/scr4/human_data/NA12878/Schadt/chemistry_"
                    + ToString(p) + ".fastb", start, stop );
               vecbasevector all(tigs);
               all.Append(pb);
               cout << Date( ) << ": making lookup" << endl;
               const int K = 28;
               vec< triple<kmer<K>,int,int> > kmers_plus;
               MakeKmerLookup0( all, kmers_plus );
               cout << Date( ) << ": searching" << endl;
               vec<Bool> have1( pb.size( ), False );
               vec<Bool> have2( pb.size( ), False );
               for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
               {    int64_t j, m;
                    for ( j = i + 1; j < kmers_plus.jsize( ); j++ )
                         if ( kmers_plus[j].first != kmers_plus[i].first ) break;
                    for ( m = i; m < j; m++ )
                         if ( kmers_plus[m].second >= (int) tigs.size( ) ) break;
                    if ( m - i == 1 )
                    {    int tig = kmers_plus[i].second;
                         if ( tig == id1 || inv[tig] == id1 )
                         {    for ( int64_t k = m; k < j; k++ )
                              {    int rid = kmers_plus[k].second 
                                        - (int) tigs.size( );
                                   have1[rid] = True;    }    }
                         if ( tig == id2 || inv[tig] == id2 )
                         {    for ( int64_t k = m; k < j; k++ )
                              {    int rid = kmers_plus[k].second 
                                        - (int) tigs.size( );
                                   have2[rid] = True;    }    }    }
                    i = j - 1;   }
               for ( int i = 0; i < (int) pb.size( ); i++ )
               {    if ( have1[i] && have2[i] ) 
                         cout << p << "." << start+i << endl;    }    }    }    }
