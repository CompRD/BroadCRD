/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "math/Functions.h"
#include "SeqInterval.h"
#include "feudal/BinaryStream.h"
//#include "CoverageAnalyzer.h"

/// This loads a vec<seq_interval> pseudo-assembly and gives its stats,
/// and optionally will print out all contig sizes or coordinates.

int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String(PSEUDO);
  CommandArgument_Int_OrDefault(LONGEST_PERCENT, 0);
  CommandArgument_Bool_OrDefault(CONTIG_COORDS, False);
  EndCommandArguments;

  vec<seq_interval> pseudo;
  BinaryReader::readFile( PSEUDO, &pseudo );

  cout << "Loaded " << pseudo.size() << " contigs" << endl;

  vec<int> lengths(pseudo.size());
  for(int i=0; i<pseudo.isize(); i++)
    lengths[i] = pseudo[i].Length();

  sort(lengths.rbegin(), lengths.rend());

  vec<int> Nstats = NStatistics( lengths );
  for(int i=0; i<9; i++)
    cout << "N" << 10*(i+1) << " = " << ToStringAbbrev(Nstats[i]) << endl;

  if( LONGEST_PERCENT > 0 ) {
    int num_to_show = (lengths.isize() * LONGEST_PERCENT)/100;
    if( num_to_show == lengths.isize() )
      cout << "\nContig sizes:" << endl;
    else
      cout << "\nLargest " << num_to_show << " contig sizes:" << endl;

    for(int i=0; i<num_to_show; i++)
      cout << lengths[i] << endl;
  }

  if( CONTIG_COORDS ) {
    cout << "\nContig coordinates, in kmers:\n"
	 << endl;
    for( vec<seq_interval>::iterator contig = pseudo.begin();
	 contig != pseudo.end(); contig++ )
      cout << contig->SeqId() << ", kmers "
	   << contig->Begin() << " - "
	   << contig->End() << " ("
	   << contig->Length() << " kmers), in super "
	   << contig->IntervalId() << endl;
//       cout << contig->SeqId() << '\t'
// 	   << contig->Begin() << '\t' 
// 	   << contig->End() << '\t'
// 	   << contig->Length()
// 	   << endl;
  }

}
