/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// TryLock.  For a given DNA sequence, and for each n <= 10, compute the maximum 
// melting temperature that can be obtained by locking n nucleotides, and output the 
// sequence of the modified DNA.
//
// The implementation sucks.

#include "math/Functions.h"
#include "MainTools.h"
#include "dna/DNAHybridization.h"

void LockPrint( const String& S, const vec<Bool>& locked )
{    for ( int i = 0; i < S.isize( ); i++ )
     {    if ( locked[i] ) cout << "+";
          cout << S[i];    }    }
     
int main( int argc, char *argv[] )
{
     RunTime( );
     
     BeginCommandArguments;
     CommandArgument_String(S);
     EndCommandArguments;

     double T = Tm_NearestNeighbor( S, 0.00000025, 0.05 );
     cout << "+0 locked: " << T << "\n";

     vec<Bool> best_locked;

     {
     double max_T = 0.0;
     for ( int i1 = 2; i1 < S.isize( ) - 2; i1++ )
     {    vec<Bool> locked( S.size( ), False );
          locked[i1] = True;
          double T = Tm_NearestNeighbor( S, 0.00000025, 0.05, locked );
          if ( T > max_T ) best_locked = locked;
          max_T = Max( max_T, T );    }
     cout << "+1 locked: " << max_T << " ";
     LockPrint( S, best_locked );
     cout << "\n";
     }

     {
     double max_T = 0.0;
     for ( int i1 = 2; i1 < S.isize( ) - 2; i1++ )
     for ( int i2 = i1 + 2; i2 < S.isize( ) - 2; i2++ )
     {    vec<Bool> locked( S.size( ), False );
          locked[i1] = True;
          locked[i2] = True;
          double T = Tm_NearestNeighbor( S, 0.00000025, 0.05, locked );
          if ( T > max_T ) best_locked = locked;
          max_T = Max( max_T, T );    }
     cout << "+2 locked: " << max_T << " ";
     LockPrint( S, best_locked );
     cout << "\n";
     }

     {
     double max_T = 0.0;
     for ( int i1 = 2; i1 < S.isize( ) - 2; i1++ )
     for ( int i2 = i1 + 2; i2 < S.isize( ) - 2; i2++ )
     for ( int i3 = i2 + 2; i3 < S.isize( ) - 2; i3++ )
     {    vec<Bool> locked( S.size( ), False );
          locked[i1] = True;
          locked[i2] = True;
          locked[i3] = True;
          double T = Tm_NearestNeighbor( S, 0.00000025, 0.05, locked );
          if ( T > max_T ) best_locked = locked;
          max_T = Max( max_T, T );    }
     cout << "+3 locked: " << max_T << " ";
     LockPrint( S, best_locked );
     cout << "\n";
     }

     {
     double max_T = 0.0;
     for ( int i1 = 2; i1 < S.isize( ) - 2; i1++ )
     for ( int i2 = i1 + 2; i2 < S.isize( ) - 2; i2++ )
     for ( int i3 = i2 + 2; i3 < S.isize( ) - 2; i3++ )
     for ( int i4 = i3 + 2; i4 < S.isize( ) - 2; i4++ )
     {    vec<Bool> locked( S.size( ), False );
          locked[i1] = True;
          locked[i2] = True;
          locked[i3] = True;
          locked[i4] = True;
          double T = Tm_NearestNeighbor( S, 0.00000025, 0.05, locked );
          if ( T > max_T ) best_locked = locked;
          max_T = Max( max_T, T );    }
     cout << "+4 locked: " << max_T << " ";
     LockPrint( S, best_locked );
     cout << "\n";
     }

     {
     double max_T = 0.0;
     for ( int i1 = 2; i1 < S.isize( ) - 2; i1++ )
     for ( int i2 = i1 + 2; i2 < S.isize( ) - 2; i2++ )
     for ( int i3 = i2 + 2; i3 < S.isize( ) - 2; i3++ )
     for ( int i4 = i3 + 2; i4 < S.isize( ) - 2; i4++ )
     for ( int i5 = i4 + 2; i5 < S.isize( ) - 2; i5++ )
     {    vec<Bool> locked( S.size( ), False );
          locked[i1] = True;
          locked[i2] = True;
          locked[i3] = True;
          locked[i4] = True;
          locked[i5] = True;
          double T = Tm_NearestNeighbor( S, 0.00000025, 0.05, locked );
          if ( T > max_T ) best_locked = locked;
          max_T = Max( max_T, T );    }
     cout << "+5 locked: " << max_T << " ";
     LockPrint( S, best_locked );
     cout << "\n";
     }

     {
     double max_T = 0.0;
     for ( int i1 = 2; i1 < S.isize( ) - 2; i1++ )
     for ( int i2 = i1 + 2; i2 < S.isize( ) - 2; i2++ )
     for ( int i3 = i2 + 2; i3 < S.isize( ) - 2; i3++ )
     for ( int i4 = i3 + 2; i4 < S.isize( ) - 2; i4++ )
     for ( int i5 = i4 + 2; i5 < S.isize( ) - 2; i5++ )
     for ( int i6 = i5 + 2; i6 < S.isize( ) - 2; i6++ )
     {    vec<Bool> locked( S.size( ), False );
          locked[i1] = True;
          locked[i2] = True;
          locked[i3] = True;
          locked[i4] = True;
          locked[i5] = True;
          locked[i6] = True;
          double T = Tm_NearestNeighbor( S, 0.00000025, 0.05, locked );
          if ( T > max_T ) best_locked = locked;
          max_T = Max( max_T, T );    }
     cout << "+6 locked: " << max_T << " ";
     LockPrint( S, best_locked );
     cout << "\n";
     }

     {
     double max_T = 0.0;
     for ( int i1 = 2; i1 < S.isize( ) - 2; i1++ )
     for ( int i2 = i1 + 2; i2 < S.isize( ) - 2; i2++ )
     for ( int i3 = i2 + 2; i3 < S.isize( ) - 2; i3++ )
     for ( int i4 = i3 + 2; i4 < S.isize( ) - 2; i4++ )
     for ( int i5 = i4 + 2; i5 < S.isize( ) - 2; i5++ )
     for ( int i6 = i5 + 2; i6 < S.isize( ) - 2; i6++ )
     for ( int i7 = i6 + 2; i7 < S.isize( ) - 2; i7++ )
     {    vec<Bool> locked( S.size( ), False );
          locked[i1] = True;
          locked[i2] = True;
          locked[i3] = True;
          locked[i4] = True;
          locked[i5] = True;
          locked[i6] = True;
          locked[i7] = True;
          double T = Tm_NearestNeighbor( S, 0.00000025, 0.05, locked );
          if ( T > max_T ) best_locked = locked;
          max_T = Max( max_T, T );    }
     cout << "+7 locked: " << max_T << " ";
     LockPrint( S, best_locked );
     cout << "\n";
     }

     {
     double max_T = 0.0;
     for ( int i1 = 2; i1 < S.isize( ) - 2; i1++ )
     for ( int i2 = i1 + 2; i2 < S.isize( ) - 2; i2++ )
     for ( int i3 = i2 + 2; i3 < S.isize( ) - 2; i3++ )
     for ( int i4 = i3 + 2; i4 < S.isize( ) - 2; i4++ )
     for ( int i5 = i4 + 2; i5 < S.isize( ) - 2; i5++ )
     for ( int i6 = i5 + 2; i6 < S.isize( ) - 2; i6++ )
     for ( int i7 = i6 + 2; i7 < S.isize( ) - 2; i7++ )
     for ( int i8 = i7 + 2; i8 < S.isize( ) - 2; i8++ )
     {    vec<Bool> locked( S.size( ), False );
          locked[i1] = True;
          locked[i2] = True;
          locked[i3] = True;
          locked[i4] = True;
          locked[i5] = True;
          locked[i6] = True;
          locked[i7] = True;
          locked[i8] = True;
          double T = Tm_NearestNeighbor( S, 0.00000025, 0.05, locked );
          if ( T > max_T ) best_locked = locked;
          max_T = Max( max_T, T );    }
     cout << "+8 locked: " << max_T << " ";
     LockPrint( S, best_locked );
     cout << "\n";
     }

     {
     double max_T = 0.0;
     for ( int i1 = 2; i1 < S.isize( ) - 2; i1++ )
     for ( int i2 = i1 + 2; i2 < S.isize( ) - 2; i2++ )
     for ( int i3 = i2 + 2; i3 < S.isize( ) - 2; i3++ )
     for ( int i4 = i3 + 2; i4 < S.isize( ) - 2; i4++ )
     for ( int i5 = i4 + 2; i5 < S.isize( ) - 2; i5++ )
     for ( int i6 = i5 + 2; i6 < S.isize( ) - 2; i6++ )
     for ( int i7 = i6 + 2; i7 < S.isize( ) - 2; i7++ )
     for ( int i8 = i7 + 2; i8 < S.isize( ) - 2; i8++ )
     for ( int i9 = i8 + 2; i9 < S.isize( ) - 2; i9++ )
     {    vec<Bool> locked( S.size( ), False );
          locked[i1] = True;
          locked[i2] = True;
          locked[i3] = True;
          locked[i4] = True;
          locked[i5] = True;
          locked[i6] = True;
          locked[i7] = True;
          locked[i8] = True;
          locked[i9] = True;
          double T = Tm_NearestNeighbor( S, 0.00000025, 0.05, locked );
          if ( T > max_T ) best_locked = locked;
          max_T = Max( max_T, T );    }
     cout << "+9 locked: " << max_T << " ";
     LockPrint( S, best_locked );
     cout << "\n";
     }

     {
     double max_T = 0.0;
     for ( int i1 = 2; i1 < S.isize( ) - 2; i1++ )
     for ( int i2 = i1 + 2; i2 < S.isize( ) - 2; i2++ )
     for ( int i3 = i2 + 2; i3 < S.isize( ) - 2; i3++ )
     for ( int i4 = i3 + 2; i4 < S.isize( ) - 2; i4++ )
     for ( int i5 = i4 + 2; i5 < S.isize( ) - 2; i5++ )
     for ( int i6 = i5 + 2; i6 < S.isize( ) - 2; i6++ )
     for ( int i7 = i6 + 2; i7 < S.isize( ) - 2; i7++ )
     for ( int i8 = i7 + 2; i8 < S.isize( ) - 2; i8++ )
     for ( int i9 = i8 + 2; i9 < S.isize( ) - 2; i9++ )
     for ( int i10 = i9 + 2; i10 < S.isize( ) - 2; i10++ )
     {    vec<Bool> locked( S.size( ), False );
          locked[i1] = True;
          locked[i2] = True;
          locked[i3] = True;
          locked[i4] = True;
          locked[i5] = True;
          locked[i6] = True;
          locked[i7] = True;
          locked[i8] = True;
          locked[i9] = True;
          locked[i10] = True;
          double T = Tm_NearestNeighbor( S, 0.00000025, 0.05, locked );
          if ( T > max_T ) best_locked = locked;
          max_T = Max( max_T, T );    }
     cout << "+10 locked: " << max_T << " ";
     LockPrint( S, best_locked );
     cout << "\n";
     }

          }
