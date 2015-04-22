///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// KmerRepeatContent.  Compute the k-mer repeat content of a given genome, defined
// to be the fraction of positions on the genome at which the k-mer starting there 
// is an exact duplicate (or the reverse complement) of a k-mer starting at 
// another position (ignoring kmers containing ambiguous bases).

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "FetchReads.h"
#include "FetchReadsAmb.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "kmers/KmerRecord.h"

template<int K> void KmerRepeatContent( const vecbasevector& genome, 
     const vecbitvector& genome_amb )
{    vec<int64_t> starts;
     starts.push_back(0);
     for ( size_t i = 0; i < genome.size( ); i++ )
     {    const basevector& u = genome[i];
          starts.push_back( starts.back( ) + u.isize( ) - K + 1 );    }
     vec< kmer<K> > kmers( starts.back( ) );
     vec<Bool> to_remove( starts.back( ), False );
     #pragma omp parallel for schedule( dynamic, 1 )
     for ( size_t i = 0; i < genome.size( ); i++ )
     {    const basevector& u = genome[i];
          kmer<K> x, xrc;
          for ( int j = 0; j <= u.isize( ) - K; j++ )
          {    int64_t r = starts[i] + j;
               x.SetToSubOf( u, j );
               xrc = x;
               xrc.ReverseComplement( );
               kmers[r] = ( x < xrc ? x : xrc );
               Bool amb = False;
               for ( int m = 0; m < K; m++ )
               {    if ( genome_amb[i][j+m] ) 
                    {    amb = True;
                         break;    }    }
               if (amb) to_remove[r] = True;    }    }
     EraseIf( kmers, to_remove );
     ParallelSort(kmers);
     int64_t total = 0, repetitive = 0;
     for ( size_t i = 0; i < kmers.size( ); i++ )
     {    size_t j = kmers.NextDiff(i);
          total += j - i;
          if ( j - i > 1 ) repetitive += j - i;
          i = j - 1;    }
     cout << PERCENT_RATIO( 4, repetitive, total ) << endl;    }


int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(F, "fasta file for genome");
     CommandArgument_Int(K);
     EndCommandArguments;

     vecbasevector genome;
     vecbitvector genome_amb;
     FetchReads( genome, 0, F );
     FetchReadsAmb( genome_amb, F );

     if ( K == 4 ) KmerRepeatContent<4>( genome, genome_amb );
     else if ( K == 100 ) KmerRepeatContent<100>( genome, genome_amb );
     else
     {    cout << "K value unsupported" << endl;
          exit(1);    }    }

