///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "Qualvector.h"
#include "kmers/KmerRecord.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(INSTANCE);
     EndCommandArguments;

     cout << Date( ) << ": importing reads" << endl;
     SystemSucceed( "LongProto SAMPLE=unknown "
          "READS=/wga/scr4/jaffe/fos_filter/sub." + INSTANCE + ".bam "
          "TMP=/wga/scr4/jaffe/Glumper/tmp.xxx." + INSTANCE 
          + " EXIT=LOAD > lp.outx." + INSTANCE );

     cout << Date( ) << ": calling ReadQGrapher" << endl;
     SystemSucceed( "ReadQGrapher IN_HEAD=/wga/scr4/jaffe/Glumper/tmp.xxx." + INSTANCE
          + "/frag_reads_orig "
          "OUT_HEAD=fido." + INSTANCE + " MIN_QUAL=20 MIN_FREQ=5 FILL_GAPS=False "
          "JOIN_OVERLAPS=False > rqg.outx." + INSTANCE );
     cout << Date( ) << ": loading RQG edges" << endl;
     vecbasevector edges( "fido." + INSTANCE + ".60.hbv.fastb" );
     cout << Date( ) << ": loading Fosmid reads" << endl;
     vecbasevector fbases(
          "/wga/scr4/jaffe/fos_filter/tmp.fos/frag_reads_orig.fastb" );

     cout << Date( ) << ": building RQG kmer list" << endl;
     const int K = 60;
     vec< kmer<K> > strong_wgs;
     for ( int64_t e = 0; e < (int64_t) edges.size( ); e++ )
     {    kmer<K> x;
          for ( int j = 0; j <= edges[e].isize( ) - K; j++ )
          {    x.SetToSubOf( edges[e], j );
               strong_wgs.push_back(x);    }    }
     cout << Date( ) << ": sorting kmers" << endl;
     Sort(strong_wgs);
     PRINT( strong_wgs.size( ) );

     cout << Date( ) << ": traversing Fosmid reads" << endl;
     vec<Bool> infos( strong_wgs.size( ), False );
     #pragma omp parallel for
     for ( int i = 0; i < (int) fbases.size( ); i++ )
     {    kmer<K> x;
          basevector b = fbases[i];
          for ( int pass = 1; pass <= 2; pass++ )
          {    if ( pass == 2 ) b.ReverseComplement( );
               for ( int j = 0; j < (int) b.isize( ) - K; j++ )
               {    x.SetToSubOf( b, j );
                    int p = BinPosition( strong_wgs, x );
                    if ( p >= 0 && !infos[p] ) 
                    {
                         #pragma omp critical
                         {    infos[p] = True;    }    }    }    }    }
     PRINT( Sum(infos) );

     cout << Date( ) << ": loading wgs reads" << endl;
     vecbasevector bases( "/wga/scr4/jaffe/Glumper/tmp.xxx." + INSTANCE 
          + "/frag_reads_orig.fastb" );
     vecqualvector quals( "/wga/scr4/jaffe/Glumper/tmp.xxx." + INSTANCE 
          + "/frag_reads_orig.qualb" );
     PairsManager pairs;
     pairs.Read( "/wga/scr4/jaffe/Glumper/tmp.xxx." + INSTANCE + "/frag_reads_orig.pairs" );
     PRINT( bases.size( ) );

     vec<Bool> pairs_to_delete( pairs.nPairs( ), False );
     vec<Bool> reads_to_delete( bases.size( ), False );
     cout << Date( ) << ": looking for wgs reads with illegal kmers" << endl;
     #pragma omp parallel for
     for ( int64_t pid = 0; pid < (int64_t) pairs.nPairs( ); pid++ )
     {    int64_t id1 = pairs.ID1(pid), id2 = pairs.ID2(pid);
          kmer<K> x;
          for ( int pass = 1; pass <= 2; pass++ )
          {    int id = ( pass == 1 ? id1 : id2 );
               for ( int j = 0; j <= bases[id].isize( ) - K; j++ )
               {    x.SetToSubOf( bases[id], j );
                    int p = BinPosition( strong_wgs, x );
                    if ( p >= 0 && !infos[p] )
                    {    reads_to_delete[id1] = reads_to_delete[id2] = True;
                         pairs_to_delete[pid] = True;
                         break;    }    }    }    }
     PRINT2( pairs.nPairs( ), Sum(pairs_to_delete) );

     cout << Date( ) << ": writing modified data" << endl;
     vec<int> to_new( bases.size( ) );
     int new_id = 0;
     for ( int id = 0; id < (int) bases.size( ); id++ )
     {    to_new[id] = new_id;
          if ( !reads_to_delete[id] ) new_id++;    }
     pairs.removePairs(pairs_to_delete);
     for ( int pid = 0; pid < (int) pairs.nPairs( ); pid++ )
          pairs.SetIDs( pid, to_new[ pairs.ID1(pid) ], to_new[ pairs.ID2(pid) ] );
     bases.EraseIf(reads_to_delete);
     quals.EraseIf(reads_to_delete);
     String head = "/wga/scr4/jaffe/Glumper/tmp.xxx." + INSTANCE + "/frag_reads_orig";
     bases.WriteAll( head + ".fastb" );
     quals.WriteAll( head + ".qualb" );
     pairs.Write( head + ".pairs" );
     cout << Date( ) << ": done" << endl;    }
