
#include "MainTools.h"
#include "IndexedAlignmentPlusVector.h"
#include "system/RunTime.h"
#include "system/System.h"

int main( int argc, char** argv )
{
  RunTime( );
  BeginCommandArguments;
  CommandArgument_UnsignedInt_OrDefault( START, 0 );
  CommandArgument_UnsignedInt_OrDefault( STOP, 0 );
  CommandArgument_Bool_OrDefault( PRINT, False );
  EndCommandArguments;

  String aligns_file = "/wga/scratch26/jbutler/Chimp.newaligns";
  String ids_file = "/wga/scratch23/Chimp/run/work/reads.ids";

  int n_reads = MastervecFileObjectCount( ids_file );
  if ( STOP == 0 )
    STOP = n_reads;

  VecAlignmentPlusReader aligns( aligns_file );

  const int reads_per_dot = 100000;

  cout << ( n_reads / reads_per_dot ) + 1 << " passes." << endl;

  vec<alignment_plus> subset_aligns;
  for (unsigned int id1=START; id1<STOP; id1++) 
  {
    if ( id1 % reads_per_dot == 0 )
      Dot( cout, id1 / reads_per_dot );
    aligns.ReadSubset( subset_aligns, id1 );
    if ( PRINT )
      for ( unsigned int align_idx = 0; align_idx < subset_aligns.size(); ++align_idx )
      {
        PRINT3( id1, align_idx, align(subset_aligns[align_idx].a).Nblocks() );
      }
  }

  cout << endl;
}
