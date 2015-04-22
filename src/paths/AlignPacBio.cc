/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"
#include "pairwise_aligners/LongReadsAligner.h"
#include "random/Shuffle.h"
#include "util/RunCommand.h"
#include <omp.h>
// MakeDepend: library OMP



void distrib_to_text_file(const map<size_t, size_t> & distrib, 
                          const String fn)
{
  ofstream os(fn.c_str());

  os << "# size      counts" << endl;
  for (map<size_t, size_t>::const_iterator it = distrib.begin();
       it != distrib.end(); it++) {
    os << setw(12) << it->first << " " 
       << setw(12) << it->second << endl;
  }
  os.close();

}


void distrib_from_look_aligns(const vec<look_align> & look_aligns,
                              map<size_t, size_t> * distrib_p)
{
  // The windows (aligned portions).
  map<size_t, int> al_lens;
  const size_t na = look_aligns.size();
  for (size_t ii = 0; ii < na; ii++) {
    const look_align & al = look_aligns[ii];
    size_t qid = al.query_id;
    int win = al.a.Pos2() - al.a.pos2();
    if (win > al_lens[qid])
      al_lens[qid] = win;
  }


  for (map<size_t, int>::const_iterator it = al_lens.begin();
       it != al_lens.end(); it++)
    (*distrib_p)[it->second]++;

}  






/**
 * AlignPacBio
 *
 * Align long reads to a reference by first seeding on shared
 * k-mers. See LongReadsAligner for details.
 *
 * READS: full path name of reads (fastb)
 * REFERENCE: full path name of reference (fastb)
 * QLTOUT: output file (look_aligns)
 * K: used to seed aligns (must be <= 12)
 * NUM_THREADS: use all threads if 0
 * QUERY_IDS: optionally align these ids only (overrides SAMPLE_SIZE)
 * SAMPLE_SIZE: randomly select these many reads to align (do all if -1)
 * SEED: used with SAMPLE_SIZE to randomize reads
 * VERBOSE: verbose log flag (sent to cout), see LongReadsAligner
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( READS );
  CommandArgument_String( REFERENCE );
  CommandArgument_String( QLTOUT );
  CommandArgument_Int_OrDefault( K, 11 );
  CommandArgument_Int_OrDefault( NUM_THREADS, 0 );
  CommandArgument_String_OrDefault( QUERY_IDS, "" );
  CommandArgument_Int_OrDefault( SAMPLE_SIZE, -1 );
  CommandArgument_Int_OrDefault( SEED, 666 );
  CommandArgument_Int_OrDefault( VERBOSE, 0 );
  EndCommandArguments;
  
  Mkpath( Basename( QLTOUT ) );

  // Tread control.
  NUM_THREADS = configNumThreads(NUM_THREADS);
  omp_set_num_threads( NUM_THREADS );

  // Needed files.
  vec<String> needed;
  needed.push_back( READS );
  needed.push_back( REFERENCE );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;

  // Load
  cout << Date( ) << ": loading reference" << endl;
  vecbvec reference( REFERENCE );

  cout << Date( ) << ": loading reads" << endl;
  vecbvec reads( READS );

  // Select reads to align.
  vec<int> qselect;
  if ( QUERY_IDS != "" )
    ParseIntSet( QUERY_IDS, qselect );
  else if ( SAMPLE_SIZE > -1 ) {
    Shuffle( (int)reads.size( ), qselect, SEED );
    qselect.resize( SAMPLE_SIZE );
  }
  else
    qselect.resize( 0 );

  const vec<int> *psel = qselect.size( ) > 0 ? &qselect : 0;
  const size_t nsel = (qselect.size( ) > 0 ? qselect.size() : reads.size()); 
  // Align queries (logging within).
  vec<look_align> look_aligns;
  LongReadsAligner(K, reference, reads, &look_aligns, 
                   psel, VERBOSE, &cout);
  

  // Save aligns (and clean up).
  cout << Date() << ": saving " << look_aligns.size() << " aligns\n" << endl;
  WriteLookAligns(QLTOUT, look_aligns);

  // save distributions
  map<size_t, size_t> size_distrib_aligned;
  distrib_from_look_aligns(look_aligns, &size_distrib_aligned);
  distrib_to_text_file(size_distrib_aligned, QLTOUT + ".aligned.distrib");

  // Generate coverages.
  cout << "n_sel= " << nsel << endl;
  LongReadsCoverages(reference, reads, look_aligns, nsel, cout);

  
  // Done.
  cout << Date( ) << ": AlignPacBio done" << endl;

}


