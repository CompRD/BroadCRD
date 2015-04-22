///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Count number of edges that each read appears in.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "kmers/KmerRecord.h"
#include "paths/long/MakeKmerStuff.h"
#include "paths/long/large/tenx/TenxDirs.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_Int_Doc(N, "1 or 2 or 3 or 4 or 5 or 6");
     EndCommandArguments;

     // Hardcoded directories.

     String dir, odir, tdir;
     SetTenxDirs( N, dir, odir, tdir );

     // Load assembly.

     cout << Date( ) << ": loading assembly" << endl;
     vecbasevector tigs( dir + "/a.fastb" );
     vec<int> inv;
     BinaryReader::readFile( dir + "/a.inv", &inv );
     for ( int e = 0; e < (int) tigs.size( ); e++ )
          if ( inv[e] < e ) tigs[e].resize(0);

     // Load 10X data.

     cout << Date( ) << ": loading reads" << endl;
     cout << Date( ) << ": loading 10X data" << endl;
     vecbasevector bases( tdir + "/10X.fastb" );

     // Make lookup table.

     cout << Date( ) << ": making lookup table" << endl;
     const int K = 88;
     vecbasevector all(tigs);
     all.Append(bases);
     vec< triple<kmer<K>,int,int> > kmers_plus;
     MakeKmerLookup2( all, kmers_plus );

     // Count.

     cout << Date( ) << ": counting" << endl;
     vec<int> count( bases.size( ), 0 );
     for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
     {    int64_t j, m;
          for ( j = i + 1; j < kmers_plus.jsize( ); j++ )
               if ( kmers_plus[j].first != kmers_plus[i].first ) break;
          for ( m = i + 1; m < j; m++ )
               if ( kmers_plus[m].second >= (int) tigs.size( ) ) break;
          if ( m > i )
          {    for ( int64_t l = m; l < j; l++ )
               {    int rid = kmers_plus[l].second - (int) tigs.size( );
                    count[rid] += m - i;    }    }
          i = j - 1;    }

     // Write.

     cout << Date( ) << ": writing" << endl;
     BinaryWriter::writeFile( odir + "/10X.count", count );
     cout << Date( ) << ": done" << endl;
     cout << "peak memory usage = " << ToStringAddCommas( PeakMemUsageBytes( ) )
          << endl;
     Scram(0);    }
