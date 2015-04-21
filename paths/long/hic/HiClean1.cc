///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// HiClean1.  Load read pairs from parallel files.  Truncate each read at the first 
// base having quality < MINQ = 20.  Kill pairs if either read has length < K = 60.  
// Kill pairs that have the same starting kmers, saving the pair that has the most 
// kmers.  Write modified pairs.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "Qualvector.h"
#include "kmers/KmerRecord.h"

template<int K> void MarkDups( const vec<vecbasevector>& bases, vec<Bool>& to_delete )
{    cout << Date( ) << ": finding kmers" << endl;
     vec< quad< kmer<K>, kmer<K>, int, int > > kmers_plus( bases[0].size( ) );
     #pragma omp parallel for
     for ( size_t i = 0; i < bases[0].size( ); i++ )
     {    if ( to_delete[i] ) continue;
          kmers_plus[i].first.SetToSubOf( bases[0][i], 0 );
          kmers_plus[i].second.SetToSubOf( bases[1][i], 0 );
          kmers_plus[i].third = - bases[0][i].isize( ) - bases[1][i].isize( );
          kmers_plus[i].fourth = i;    }
     ParallelSort(kmers_plus);
     cout << Date( ) << ": marking duplicates" << endl;
     for ( int i = 0; i < kmers_plus.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < kmers_plus.isize( ); j++ )
          {    if ( kmers_plus[j].first != kmers_plus[i].first ) break;
               if ( kmers_plus[j].second != kmers_plus[i].second ) break;    }
          for ( int k = i + 1; k < j; k++ ) 
               to_delete[ kmers_plus[k].fourth ] = True;
          i = j - 1;    }    }

int main(int argc, char *argv[])
{    RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(HEAD, "inputs HEAD{_R1,_R2}.{fastb,qualb}, "
          "outputs HEAD{_R1,_R2}.clean1.fastb");
     CommandArgument_Int_OrDefault_Doc(MINQ, 20, "minimum quality score");
     CommandArgument_Int_OrDefault_Doc(K, 60, "K value");
     EndCommandArguments;

     // Trim reads.

     cout << Date( ) << ": loading" << endl;
     vec<vecbasevector> bases(2);
     bases[0].ReadAll( HEAD + "_R1.fastb" ), bases[1].ReadAll( HEAD + "_R2.fastb" );
     vec<vecqualvector> quals(2);
     quals[0].ReadAll( HEAD + "_R1.qualb" ), quals[1].ReadAll( HEAD + "_R2.qualb" );
     vec<Bool> to_delete( bases[0].size( ), False );
     cout << Date( ) << ": trimming" << endl;
     for ( int j = 0; j < 2; j++ )
     {    
          #pragma omp parallel for
          for ( int id = 0; id < (int) bases[j].size( ); id++ )
          {    for ( int p = 0; p < bases[j][id].isize( ); p++ )
               {    if ( quals[j][id][p] < MINQ )
                    {    bases[j][id].resize(p);
                         break;    }    }    }    }
     for ( int id = 0; id < (int) bases[0].size( ); id++ )
     {    if ( bases[0][id].isize( ) < K || bases[1][id].isize( ) < K )
               to_delete[id] = True;    }
     
     // Find kmers and mark duplicates.

     if ( K == 60 ) MarkDups<60>( bases, to_delete );
     else 
     {    cout << "Not implemented for that K value." << endl;
          Scram(1);    }

     // Remove killed reads.

     cout << Date( ) << ": removing killed reads" << endl;
     for ( int p = 0; p < 2; p++ )
          bases[p].EraseIf(to_delete);

     // Write files.
          
     cout << Date( ) << ": writing files" << endl;
     bases[0].WriteAll( HEAD + "_R1.clean1.fastb" );
     bases[1].WriteAll( HEAD + "_R2.clean1.fastb" );
     cout << Date( ) << ": done" << endl;    }
