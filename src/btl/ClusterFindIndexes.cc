///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"

#include "Basevector.h"
#include "Fastavector.h"
#include "Histogram.h"
#include "SeqInterval.h"
#include "VecUtilities.h"
#include "btl/CBarcodes.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "paths/reporting/ReftigUtils.h"
#include "util/RunCommand.h"
#include <omp.h>
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

typedef pair<seq_interval,int> pseudo_t;

/**
 * PrintPseudoRead
 *
 * A one-liner pseudo_read printer.
 */
void PrintPseudoRead( const seq_interval &si, const int &mult, ostream &out )
{
  int tid = si.SeqId( ) < 0 ? -1 - si.SeqId( ) : si.SeqId( );
  char orient = si.SeqId( ) < 0? '-' : '+';
  out << tid << orient << "\t"
      << si.Begin( ) << ", " << si.End( ) << "\t"
      << si.Length( ) << "\t"
      << mult << "\n";
}

/**
 * FindValidSets
 *
 * Look for all possible combinations of valid pairs.
 */
void FindValidSets( const int &MIN_SEP,
		    const int &MAX_SEP,
		    const bool &CIRCULAR,
		    const vec<pseudo_t> &pseudos,
		    const vec<int> &tlens,
		    vec<int> &vlens,
		    vec<int> &pcounts,
		    vec<int> &plens,
		    ostream *log = 0 )
{
  vlens.clear( );
  pcounts.clear( );
  plens.clear( );
  
  for (int ii=0; ii<pseudos.isize( ); ii++) {
    const pseudo_t &left = pseudos[ii];
    const seq_interval &lsi = left.first;
    if ( lsi.SeqId( ) > -1 ) continue;
    
    for (int jj=ii+1; jj<pseudos.isize( ); jj++) {
      const pseudo_t &right = pseudos[jj];
      const seq_interval &rsi = right.first;
      if ( rsi.SeqId( ) < 0 ) continue;

      if ( - 1 - lsi.SeqId( ) != rsi.SeqId( ) ) continue;

      int vlen = rsi.End( ) - lsi.Begin( );
      if(CIRCULAR){
          vlen = (vlen+tlens[0])%tlens[0];
          ForceAssert(vlen>=0);
      }
      if ( vlen < MIN_SEP || vlen > MAX_SEP ) continue;
      
      vlens.push_back( vlen );
      
      if ( log ) {
	pcounts.push_back( left.second );
	pcounts.push_back( right.second );
	plens.push_back( lsi.Length( ) );
	plens.push_back( rsi.Length( ) );
	
	*log << "insert size: " << vlen << "\n";
	PrintPseudoRead( lsi, left.second, *log );
	PrintPseudoRead( rsi, right.second, *log );
	*log << "\n";
      }
    }
  }

}

/**
 * ClusterIndexedReads
 *
 * Create pseudo-reads from fragment with the same indexes. This code
 * must run after FindIndexes, and it's supposed to just generate a
 * bunch of stats on how paired pseudo-reads look like. A pseudo-read
 * is what you get when you merge reads (with the same orientation),
 * on the reference.
 */
int main( int argc, char* argv[] )
{
  RunTime( );
 
  BeginCommandArguments;
  CommandArgument_String( IN_DIR );
  CommandArgument_String( REF_HEAD );
  CommandArgument_String_OrDefault( HEAD, "indexed_reads" );

  // If reference is circular.
  CommandArgument_Bool_OrDefault( CIRCULAR, True );
  
  // Max gap size.
  CommandArgument_Int_OrDefault( MAX_GAP, 500 );
  
  // Min and max separation size, and histogram related args
  CommandArgument_Int_OrDefault( MIN_SEP, 0 );
  CommandArgument_Int_OrDefault( MAX_SEP, 75000 );
  CommandArgument_Int_OrDefault( BIN_SIZE, 500 );

  // Do not use cached data, if FORCE=True (and other aligners args).
  CommandArgument_Int_OrDefault( NUM_THREADS, 0 );
  CommandArgument_Bool_OrDefault( FORCE, True );

  EndCommandArguments;
  
  // Dir and file names.
  String lookup = REF_HEAD + ".lookup";
  String ref_fastb = REF_HEAD + ".fastb";

  String full_head = IN_DIR + "/" + HEAD;
  String rfasta_file = full_head + ".fasta";
  String rfastb_file = full_head + ".fastb";
  String names_file = full_head + ".ids";
  String aligns_file = full_head + ".qlt";
  String sis_file = full_head + ".sis";
  String barcounts_file = full_head + ".barcounts";

  String histo_eps_file = full_head + ".ClusterFindIndexes.eps";
  String histo_txt_file = full_head + ".ClusterFindIndexes.txt";
  String valid_sets_file =  full_head + ".ClusterFindIndexes.ValidSets.txt";
    
  vec<String> needed;
  needed.push_back( lookup );
  needed.push_back( ref_fastb );
  needed.push_back( rfasta_file );
  needed.push_back( barcounts_file );
  
  // Thread control (needed by FastAlignShortReads).
  NUM_THREADS = configNumThreads( NUM_THREADS );
  omp_set_num_threads( NUM_THREADS );

  // Load references lengths.
  cout << Date( ) << ": loading reference" << endl;
  vec<int> tlens;
  {
    bool bMaxLen=MAX_SEP<0;
    vecbvec ref( ref_fastb );
    tlens.reserve( ref.size( ) );
    for (size_t ii=0; ii<ref.size( ); ii++){
      if(bMaxLen){MAX_SEP=max(MAX_SEP,ref[ii].isize());}
      tlens.push_back( ref[ii].size( ) );
    }
  }

  // Load barcounts file.
  cout << Date( ) << ": loading barcounts" << endl;
  vec<CBarcodes> indexes;
  {
    String name = "";
    int cnt = 0;
    ifstream in( barcounts_file.c_str( ) );
    while( in ) {
      in >> name >> cnt;
      if ( !in ) break;
      CBarcodes idx;
      idx.SetIndexes( name );
      indexes.push_back( idx );
    }
    in.close( );
    sort( indexes.begin( ), indexes.end( ) );
  }

  // Convert reads to fastb.
  if ( FORCE || ! IsRegularFile( rfastb_file ) ) {
    cout << Date( ) << ": converting genomic chunks to fastb" << endl;
    
    vec<fastavector> fastas;
    vec<String> names;
    LoadFromFastaFile( rfasta_file, fastas, names );
    
    vecbvec bases;
    bases.reserve( fastas.size( ) );
    for (size_t ii=0; ii<fastas.size( ); ii++)
      bases.push_back( fastas[ii].ToBasevector( ) );
    
    bases.WriteAll( rfastb_file );
    WRITE( names_file, names );
  }

  // Align reads (or load cached aligns ), and generate seq_intervals.
  if ( FORCE || ! IsRegularFile( sis_file ) ) {
    const int K = 26;
    vec<look_align> hits;
    GetAlignsFast( K, rfastb_file, lookup, aligns_file, hits, !FORCE, IN_DIR );
    
    // Load read names (needed to map reads to barcodes).
    cout << Date( ) << ": loading read names" << endl;
    READ( names_file, vec<String>, names );
    
    // Map read ids to align ids.
    cout << Date( ) << ": mapping reads to aligns" << endl;
    vec<int> to_hit( names.size( ), -1 );
    for (int ii=0; ii<hits.isize( ); ii++) {
      int read_id = hits[ii].query_id;
      to_hit[read_id] = ( to_hit[read_id]  == -1 ) ? ii : -2;
    }

    // Convert aligns to seq_intervals.
    cout << Date( ) << ": converting to seq_intervals" << endl;
    vec<seq_interval> sis;
    sis.reserve( names.size( ) );
    vec<CBarcodes>::iterator it;
    for (int rid=0; rid<names.isize( ); rid++) {
      if ( to_hit[rid] < 0 ) continue;
      const look_align &hit = hits[ to_hit[rid] ];
      bool rc = hit.Rc1( );
      
      // seq_id (= target_id, signed for orientation), begin, and end.
      int seq_id = rc ? - hit.target_id - 1 : hit.target_id;
      int beg = hit.StartOnTarget( );
      int end = hit.EndOnTarget( );

      // interval_id (= barcode_id).
      CBarcodes barc;
      barc.SetIndexes( names[rid].After( "_" ).Before( "_" ) );
      it = lower_bound( indexes.begin( ), indexes.end( ), barc );
      int int_id = distance( indexes.begin( ), it );
      ForceAssert( it != indexes.end( ) );
      ForceAssert( barc == indexes[int_id] );

      // Add seq_interval to list.
      sis.push_back( seq_interval( int_id, seq_id, beg, end ) );
    }

    // Sort and save (note we sort by interval_id first).
    cout << Date( ) << ": sorting and saving seq_intervals" << endl;
    seq_interval::OrderByIntervalId sorter;
    sort( sis.begin( ), sis.end( ), sorter );
    ofstream out( sis_file.c_str( ) );
    for (int ii=0; ii<sis.isize( ); ii++)
      out << sis[ii] << "\n";
    out.close( );
    
  }

  // Loading (or reloading) seq_interval.
  cout << Date( ) << ": reloading seq_intervals" << endl;
  vec<seq_interval> sis;
  vec<int> first_si;
  {
    int nreads = MastervecFileObjectCount( rfastb_file );
    sis.reserve( nreads );
    ifstream in( sis_file.c_str( ) );
    while( in ) {
      seq_interval si;
      in >> si;
      if ( !in ) break;
      sis.push_back( si );
    }
    in.close( );

    first_si.resize( indexes.size( ), -1 );
    for (int ii=sis.isize( )-1; ii>=0; ii--)
      first_si[ sis[ii].IntervalId( ) ] = ii;
  }
  
  // The list of (valid) insert lengths (with log stream).
  ofstream vout( valid_sets_file.c_str( ) );
  vout << "OVERALL STATS AT END OF FILE\n\n";

  vec<int> ilens;
  vec<int> pcounts;
  vec<int> plens;

  // Loop over all barcodes.
  for (int bar_id=0; bar_id<indexes.isize( ); bar_id++) {
    int range_begin = first_si[bar_id];
    if ( range_begin < 0 ) continue;

    int range_end = range_begin;
    for (int ii=range_begin; ii<=sis.isize( ); ii++) {
      range_end = ii;
      if ( ii == sis.isize( ) || sis[ii].IntervalId( ) != bar_id ) break;
    }

    // The pseudo reads (these allow for gaps), with count.
    vec<pseudo_t> pseudo;

    // Loop over all intervals for this barcode.
    for (int ii=range_begin; ii<range_end; ii++) {
      const seq_interval &si = sis[ii];
      const seq_interval &plast = pseudo.back( ).first;

      bool append_new
	= ( pseudo.size( ) < 1 )
	|| ( plast.SeqId( ) != si.SeqId( ) )
	|| ( si.Begin( ) - plast.End( ) > MAX_GAP );
      
      if ( append_new ) {
	pseudo.push_back( make_pair( si, 1 ) );
	continue;
      }

      int int_id = si.IntervalId( );
      int seq_id = si.SeqId( );
      int begin = plast.Begin( );
      int end = si.End( );

      int new_count = 1 + pseudo.back( ).second;
      seq_interval new_si( int_id, seq_id, begin, end );
      
      pseudo[ pseudo.size( )-1 ] = make_pair( new_si, new_count );
    }

    // List all possible valid sets.
    vec<int> vl;
    vec<int> rc;
    vec<int> rl;
    
    FindValidSets(MIN_SEP, MAX_SEP, CIRCULAR, pseudo, tlens, vl, rc, rl, &vout);
    copy( vl.begin( ), vl.end( ), back_inserter( ilens ) );
    copy( rc.begin( ), rc.end( ), back_inserter( pcounts ) );
    copy( rl.begin( ), rl.end( ), back_inserter( plens ) );
  }

  // Overall stats.
  cout << Date( ) << ": sorting valid insert lenghts" << endl;

  sort( ilens.begin( ), ilens.end( ) );
  sort( pcounts.begin( ), pcounts.end( ) );
  sort( plens.begin( ), plens.end( ) );

  vout << "OVERALL STATS\n\n"
       << "Number of genomic fragments in pseudo-reads:   "
       << ToStringAddCommas( BigSum(pcounts) ) << "\n"
       << "N50 of genomic fragments in each pseudo-read:  "
       << ToStringAddCommas( N50(pcounts) ) << "\n"
       << "N50 of pseudo-reads:                           "
       << ToStringAddCommas( N50(plens) ) << "\n"
       << "\n";
  
  vout.close( );
  
  // Build histogram of separations.
  histogram<int> histo;

  vec<int> bins;
  int step = 0;
  while( 1 ) {
    int bin = MIN_SEP + step * BIN_SIZE;
    if ( bin > MAX_SEP ) break;
    bins.push_back( bin );
    step++;
  }
  histo.AddBins( bins );
  
  for (int ii=0; ii<ilens.isize( ); ii++)
    histo.AddDatum( ilens[ii] );
  
  String str_nval = ToStringAddCommas( ilens.size( ) );
  String str_min = ToStringAddCommas( MIN_SEP );
  String str_max = ToStringAddCommas( MAX_SEP );
  String m1 = "Histogram of " + str_nval + " valid insert lengths";
  String m2 = "MIN_SEP = " + str_min;
  String m3 = "MAX_SEP = " + str_max;
  
  ofstream txt_out( histo_txt_file.c_str( ) );
  txt_out << m1 << "\n"
	  << m2 << "\n"
	  << m3 << "\n"
	  << "\n";
  histo.PrintAsColumns( txt_out, true, true, false );
  txt_out << endl;
  txt_out.close( ) ;

  {
    using namespace ns_psplot;
    
    vec<ns_psplot::freetext> labels( 4 );
    labels[0] = freetext( m1, black, 16 );
    labels[1] = freetext( m2, black, 12 );
    labels[2] = freetext( m3, black, 12 );

    ofstream eps_out( histo_eps_file.c_str( ) );
    histo.PrintAsEps( eps_out, labels, 0 );
    eps_out.close( );
  }
  
  // Done.
  cout << Date( ) << ": ClusterFindIndexes done" << endl;

}
