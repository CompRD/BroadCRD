///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "random/Shuffle.h"
#include "PairsManager.h"
#include "feudal/IncrementalWriter.h"
#include <map>


bool qual_ok(const qvec& qv, const u_int min_qual)
{
  u_int n = qv.size();
  u_int q = 0;
  for (size_t i = 0; i < n; ++i)
    q += qv[i];
  if (q/n < min_qual) return false;
  return true;
}



/**
 * SelectMinQualPairs
 *
 * Randomly pairs of minimum quality.
 *
 * READS_IN: it loads <READS_IN>.{fastb,qualb,pairs}
 * READS_OUT: it saves <READS_OUT>.{fastb,qualb,pairs,select}
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( READS_IN );
  CommandArgument_String( READS_OUT );
  CommandArgument_UnsignedInt( MIN_QUAL );
  EndCommandArguments;

  // File names.
  String in_bases_file = READS_IN + ".fastb";
  String in_quals_file = READS_IN + ".qualb";
  String in_pairs_file = READS_IN + ".pairs";

  String out_bases_file = READS_OUT + ".fastb";
  String out_quals_file = READS_OUT + ".qualb";
  String out_select_file = READS_OUT + ".select";
  
  // Output dir.
  String out_dir = out_bases_file;
  if ( ! out_dir.Contains( "/" ) ) out_dir = ".";
  else out_dir = out_dir.RevBefore( "/" );
  if ( out_dir != "." ) Mkpath( out_dir );

  cout << Date( ) << ": loading pairing info" << endl;
  PairsManager pairs( in_pairs_file );
  size_t n_pairs = pairs.nPairs();
    

  // Selected pairs.
  PairsManager sel_pairs(pairs.nReads());
  vec<size_t> select;
  
  { // scope to minimize qual footprint
    cout << Date( ) << ": loading quals" << endl;
    vecqvec quals(in_quals_file);

    size_t discarded = 0;
    size_t kept = 0;

    for (size_t i = 0; i < pairs.nPairs(); i++) {
	size_t id1 = pairs.ID1(i);
	size_t id2 = pairs.ID2(i);
	if (qual_ok(quals[id1], MIN_QUAL) && qual_ok(quals[id2], MIN_QUAL)) {
	  int new_id1 = select.size();
	  select.push_back(id1);
	  int new_id2 = select.size();
	  select.push_back(id2);
	  sel_pairs.addPair(new_id1, new_id2, pairs.sep(i), pairs.sd(i), pairs.libraryName(i));
	  ++kept;
	} else ++discarded;
    }

    cout << Date() << ": kept " << kept
	 << " pairs (" << round(100.0*kept/(kept+discarded)) << "%)" << endl;


    cout << Date( ) << ": saving pairs" << endl;
    sel_pairs.Write( READS_OUT + ".pairs" );
    cout << Date( ) << ": saving quals" << endl;
    IncrementalWriter<qvec> sel_quals(out_quals_file.c_str());
    for (size_t ii=0; ii < select.size(); ii++)
      sel_quals.add( quals[select[ii]] );
    sel_quals.close();
  }

  {
    //  bases - scoped for memory.
    // Load bases and pairs.
    cout << Date( ) << ": loading bases" << endl;
    vecbvec bases( in_bases_file );
    
    cout << Date( ) << ": saving bases" << endl;
    IncrementalWriter<bvec> sel_bases(out_bases_file.c_str());
    for (size_t ii=0; ii < select.size(); ii++)
      sel_bases.add( bases[select[ii]] );
    sel_bases.close();
  }
  // Done.
  cout << Date( ) << ": done" << endl;
  
}
