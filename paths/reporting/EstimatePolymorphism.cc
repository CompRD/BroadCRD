///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"
#include "Intvector.h"
#include "ParseSet.h"
#include "VecUtilities.h"
#include "paths/ReadLoc.h"

/**
 * IsHighQuality
 *
 * Decide if this base consensus is high quality, and return base(s)
 * call.
 */
bool IsHighQuality( const dumbcall &raw,
		    String &hq_call,
		    const int MIN_COV,
		    const double MIN_RATIO,
		    const double NEXT_RATIO,
		    const double MAX_DISCREP,
		    const bool VERBOSE = False )
{
  // Defaults.
  hq_call = "";
  
  // Find best and second best bases.
  vec< pair<int,char> > n2base;
  n2base.reserve( 6 );
  n2base.push_back( make_pair( raw.A( ), 'A' ) );
  n2base.push_back( make_pair( raw.C( ), 'C' ) );
  n2base.push_back( make_pair( raw.G( ), 'G' ) );
  n2base.push_back( make_pair( raw.T( ), 'T' ) );
  n2base.push_back( make_pair( raw.D( ), 'D' ) );
  n2base.push_back( make_pair( raw.I( ), 'I' ) );
  sort( n2base.rbegin( ), n2base.rend( ) );
  
  int n0 = n2base[0].first;
  int n1 = n2base[1].first;
  int n2 = n2base[2].first;
  int ntotal =0;
  for (int ii=0; ii<6; ii++)
    ntotal += n2base[ii].first;
  
  char base0 = n2base[0].second;
  char base1 = n2base[1].second;

  // Log verbose info (one line).
  if ( VERBOSE ) {
    cout << n0 << " " << base0 << "s "
      	 << n1 << " " << base1 << "s (total: " << ntotal << ")";
    if ( n0 > 0 ) {
      cout << " ratio0: " << ToString( SafeQuotient( n0, ntotal ), 2 )
	   << " ratio1: " << ToString( SafeQuotient( n0 + n1, ntotal ), 2 )
	   << " discrep: " << ToString( SafeQuotient( n1, n0 ), 2 )
	   << "   ";
    }
  }
  
  // High quality call (one haplotype).
  if ( n0 >= MIN_COV &&
       SafeQuotient( n0, ntotal ) >= MIN_RATIO &&
       ( n1 == 0 || SafeQuotient( n0, n1 ) >= NEXT_RATIO ) &&
       base0 != 'D' ) {
    hq_call = String( n2base[0].second );
    if ( VERBOSE ) cout << " hq, haploid\n";
    return true;
  }
  
  // High quality call (two haplotypes).
  if ( n0 + n1 >= MIN_COV &&
       SafeQuotient( n0 + n1, ntotal ) >= MIN_RATIO &&
       ( n2 == 0 || SafeQuotient( n1, n2 ) >= NEXT_RATIO ) && 
       SafeQuotient( n1, n0 ) >= MAX_DISCREP && 
       ( base0 != 'D' && base1 != 'D' ) ) {
    hq_call = String( n2base[0].second ) + String( n2base[1].second );
    if ( VERBOSE ) cout << " hq, diploid\n";
    return true;
  }

  // Everything else is not trusted.
  if ( VERBOSE ) cout << " BAD\n";
  return false;
}

/**
 * EstimatePolymorphism
 *
 * Estimate polymorphism of a genome by looking at pile ups of aligned
 * reads.
 *
 * RADIUS: defines a NQS-like neighborhood
 * MIN_COV, MIN_RATIO, NEXT_RATIO, MAX_DISCREP: args of IsHighQuality
 * USE_JUMPS: use also jump (and long jump) reads in the Pileups
 * ANNOTATIONS: optionally, save an annotation file (as a VecIntVec)
 * VERBOSE: toggle verbose mode
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_String_OrDefault(ASSEMBLY, "linear_scaffolds0.clean.patched");
  CommandArgument_String_OrDefault( TIGS, "" );
  CommandArgument_Int_OrDefault( RADIUS, 5 );
  CommandArgument_Int_OrDefault( MIN_COV, 20 );
  CommandArgument_Double_OrDefault( MIN_RATIO, 0.85 );
  CommandArgument_Double_OrDefault( NEXT_RATIO, 2.0 );
  CommandArgument_Double_OrDefault( MAX_DISCREP, 0.65 );
  CommandArgument_Bool_OrDefault( USE_JUMPS, True );
  CommandArgument_String_OrDefault( ANNOTATIONS, "" );
  CommandArgument_Bool_OrDefault( VERBOSE, False );
  EndCommandArguments;
  
  // Dir and file names.
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR; 

  String frag_reads_head = run_dir + "/frag_reads_filt_cpd";
  String jump_reads_head = run_dir + "/jump_reads_filt_cpd";
  String long_jump_reads_head = run_dir + "/long_jump_reads_filt";

  String assembly_head = sub_dir + "/" + ASSEMBLY;
  String contigs_file = assembly_head + ".contigs.fastb";
  
  // Load.
  vecbvec contigs( contigs_file );
  size_t n_tigs = contigs.size ( );

  read_locs_on_disk locs_file( assembly_head, run_dir );

  bool get_frags
    = IsRegularFile( frag_reads_head + ".fastb" );
  bool get_jumps
    = USE_JUMPS ? IsRegularFile( jump_reads_head + ".fastb" ) : false;
  bool get_long_jumps
    = USE_JUMPS ? IsRegularFile( long_jump_reads_head + ".fastb" ) : false;

  cout << Date( ) << ": using "
       << ( get_frags ? "frags" : "" )
       << ( get_jumps ? ", jumps" : "" )
       << ( get_long_jumps ? ", long_jumps" : "" )
       << "\n" << endl;

  // Selected contigs.
  vec<bool> select;
  if ( TIGS == "" ) select.resize( n_tigs, true );
  else {
    vec<int> ids;
    ParseIntSet( TIGS, ids );
    select.resize( n_tigs, false );
    for (int ii=0; ii<ids.isize( ); ii++)
      select[ ids[ii] ] = true;
  }

  // Core counters.
  vec<longlong> n_bases( n_tigs, 0 );   // bases in contig
  vec<longlong> n_hqs( n_tigs, 0 );     // high quality bases in contig
  vec<longlong> n_snps( n_tigs, 0 );    // snps in contig
  
  // Annotations.
  //   0: untrusted (default)
  //   1: trusted, but not a snp
  //   2: trusted, and a snp
  VecIntVec annotations;
  if ( ANNOTATIONS != "" ) {
    annotations.resize( contigs.size( ) );
    for (size_t ii=0; ii<contigs.size( ); ii++)
      annotations[ii].resize( contigs[ii].size( ), 0 );
  }

  // Loop over all contigs.
  for (size_t contig_id=0; contig_id<n_tigs; contig_id++) {
    if ( ! select[contig_id] ) continue;
    cout << contig_id << "\t";
    if ( VERBOSE ) cout << endl;
    else cout << flush;

    // Get locs and dumbcalls.
    const bvec &tig = contigs[contig_id];
    vec<read_loc> locs;
    vec<dumbcall> calls;
    locs_file.LoadContig( contig_id, locs );
    Pileup( tig, locs, run_dir, calls, get_frags, get_jumps, get_long_jumps );

    // Parse dumbcalls to collect high quality bases.
    vec<String> hq_bases( tig.size( ), "" );
    for (size_t pos=0; pos<calls.size( ); pos++) {
      if ( VERBOSE ) cout << pos << "\t";
      bool is_hq
	= IsHighQuality( calls[pos], hq_bases[pos],
			 MIN_COV, MIN_RATIO, NEXT_RATIO, MAX_DISCREP, VERBOSE );
      if ( is_hq && ANNOTATIONS != "" ) annotations[contig_id][pos] = 1;
	
    }
    
    // Fill counters.
    n_bases[contig_id] = tig.size( );
    for (int ii=RADIUS; ii<(int)tig.size( )-RADIUS-1; ii++) {
      bool hq_stretch = true;
      for (int jj=ii-RADIUS; jj<ii+RADIUS;jj++) {
	if ( hq_bases[jj].size( ) < 1 ) {
	  hq_stretch = false;
	  break;
	}
      }
      if ( ! hq_stretch ) continue;

      n_hqs[contig_id]++;
      if ( hq_bases[ii].size( ) > 1 ) {
	n_snps[contig_id]++;
	if ( ANNOTATIONS != "" ) annotations[contig_id][ii] = 2;
      }
    }
    
    // Log contig.
    if ( VERBOSE ) cout << contig_id << " done:\t";
    cout << n_bases[contig_id] << "\t"
	 << n_hqs[contig_id] << "\t"
	 << n_snps[contig_id] << endl;
  }

  // Save annotations.
  if ( ANNOTATIONS != "" ) {
    cout << "\n" << Date( ) << ": saving annotations" << endl;
    annotations.WriteAll( ANNOTATIONS );
  }
  
  // Overall stats.
  longlong tot_bases = BigSum( n_bases );
  longlong tot_hqs = BigSum( n_hqs );
  longlong tot_snps = BigSum( n_snps);
  double ratio1 = tot_bases<1 ? .0 : 100.0 * SafeQuotient( tot_hqs, tot_bases );
  double ratio2 = tot_hqs<1 ? .0 : 100.0 * SafeQuotient( tot_snps, tot_hqs );

  cout << "\n"
       << "OVERALL STATISTICS\n"
       << "\n"
       << "Total bases in input:  " << tot_bases << "\n"
       << "Total high qual bases: " << tot_hqs
       << " (" << ToString( ratio1, 2 )
       << "% of the total bases)\n"
       << "Total snp called:      " << tot_snps
       << " (" << ToString( ratio2, 2 )
       << "% of the total high qual bases)\n"
       << "\n";
  
  // Done.
  cout << Date( ) << ": done" << endl;

}

