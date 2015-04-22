#include "MainTools.h"
#include "kmers/SortKmers.h"
#include "TaskTimer.h"

#if __GNUC__ < 3
#include <numeric>
#else
#include <ext/numeric>
#endif

template <int K, int I>
void SortTheseKmers( unsigned int S, unsigned int passes, vecbasevector &bases, vec<int> &ids )
{
  dummy<1> d1;

  vec< kmer_record<K,I> > R(S);

  TaskTimer timer;
  for ( unsigned int pass = 0; pass < passes; pass++ )
  {
    TaskTimer onepass_timer;
    timer.Start();
    onepass_timer.Start();
    SortKmers( d1, bases, ids, pass, R, S );
    onepass_timer.Stop();
    timer.Stop();
    PRINT2( S, onepass_timer );
    PrintMemUsage();
  }
  PRINT( timer );
}

int main( int argc, char** argv )
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String(FASTB);
  CommandArgument_UnsignedInt_OrDefault( ACTUAL_PASSES, 1 );
  CommandArgument_UnsignedInt_OrDefault(PERCENT_TO_LOAD, 100);
  EndCommandArguments;

  ForceAssertGt( PERCENT_TO_LOAD,   0u );
  ForceAssertLe( PERCENT_TO_LOAD, 100u );

  vecbasevector bases;
  longlong numSeqs = MastervecFileObjectCount( FASTB );

  numSeqs *= PERCENT_TO_LOAD;
  numSeqs /= 100;

  PRINT( numSeqs );
  bases.ReadRange( FASTB, 0, numSeqs );

  vec<int> ids( bases.size() );
  iota( ids.begin(), ids.end(), 0 );

  const int K = 24;
  int I = 1;

  longlong S_guess = 0;   
  for ( size_t seq = 0; seq < bases.size( ); seq++ )
  {
    const int length = bases[seq].size( );
    if ( length >= K ) 
      S_guess += length - K + 1;
    if ( length >= 1024 )
      I = 2;
  }
  S_guess += S_guess/4;

  PRINT( S_guess );
  unsigned int S = S_guess;

  if ( K == 24 )
    if ( I == 1 )
      SortTheseKmers<24,1>( S, ACTUAL_PASSES, bases, ids );
    else
      SortTheseKmers<24,2>( S, ACTUAL_PASSES, bases, ids );

  PrintMemUsage();
}


