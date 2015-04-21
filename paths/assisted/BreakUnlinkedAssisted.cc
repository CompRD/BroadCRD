///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Equiv.h"
#include "Fastavector.h"
#include "PairsManager.h"
#include "SupersHandler.h"
#include "Vec.h"
#include "VecUtilities.h"
#include "lookup/LookAlign.h"
#include "paths/Alignlet.h"
#include "paths/SaveScaffoldGraph.h"
#include "util/RunCommand.h"

/**
 * class CSARead (Scaffold Assisted Read)
 *
 * Helper class to deal with paired reads. Beware! Read ids are stored
 * as int.
 */
class CSARead {

public:

  CSARead( ) :
    cid_ ( -1 ), rid_ ( -1 ), beg_ ( 0 ), end_ ( 0 ) { }
  
  CSARead( int cid, int rid, int beg, int end ) : 
    cid_ ( cid ), rid_ ( rid ), beg_ ( beg ), end_ ( end ) { }
    
  CSARead( const alignlet &al, const int rid ) {
    cid_ = al.TargetId( );
    rid_ = rid;
    beg_ = al.pos2( );
    end_ = al.Pos2( );
    ForceAssertLt( beg_, end_ );
    if ( ! al.Fw1( ) ) swap( beg_, end_ );
  }

  int ContigId( ) const { return cid_; }
  int ReadId( ) const { return rid_; }
  bool Fw( ) const { return beg_ < end_; }
  bool Rc( ) const { return ! this->Fw( ); }
  int Begin( ) const { return ( this->Fw( ) ? beg_ : end_ ); }
  int End( ) const { return ( this->Fw( ) ? end_ : beg_ ); }
  
  friend bool operator< ( const CSARead &left, const CSARead &right ) {
    if ( left.ContigId( ) < right.ContigId( ) ) return true;
    if ( left.ContigId( ) > right.ContigId( ) ) return false;
    if ( left.Begin( ) < right.Begin( ) ) return true;
    if ( left.Begin( ) > right.Begin( ) ) return false;
    if ( left.End( ) < right.End( ) ) return true;
    if ( left.End( ) > right.End( ) ) return false;
    return ( left.ReadId( ) < right.ReadId( ) );
  }
  
  
private:

  int cid_;
  int rid_;
  int beg_;  // read is rc if end_ < beg_
  int end_;
  
};

/**
 * BreakUnlinkedAssisted
 *
 * Break scaffolds into connected components (note there is another
 * Arachne-specific module called BreakUnlinked). This module should
 * be run after MergeContigsOnReference, which builds scaffolds based
 * on aligned of contigs on a reference.
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_String_OrDefault( READS, "scaffold_reads" );
  CommandArgument_String_OrDefault( ALIGNS, "scaffold_reads_filtered" );
  CommandArgument_String_OrDefault( HEAD_IN, "initial_scaffolds0.ref.applied" );
  CommandArgument_String_OrDefault( HEAD_OUT, "assisted_scaffolds" );

  // A link is consistent iff the separation between the two end reads
  // implied by the given gap satisfies:
  //
  //    MIN_SEP <= sep <= LIB_MULTIPLIER * lib_size
  //
  // where lib_size is the library mean separation.
  CommandArgument_Int_OrDefault( MIN_SEP, 10 );
  CommandArgument_Double_OrDefault( LIB_MULTIPLIER, 5.0 );

  // Two contigs are linked iff they own >= MIN_LINKS consistent links.
  CommandArgument_Int_OrDefault( MIN_LINKS, 3 );
  
  // Gap dev args: dev = Max( DEV_RATIO * gap_size, MIN_DEV ).
  CommandArgument_Double_OrDefault( DEV_RATIO, 0.25 );
  CommandArgument_Int_OrDefault( MIN_DEV, 20 );

  // Toggle verbose log.
  CommandArgument_Bool_OrDefault( VERBOSE, False );

  EndCommandArguments;
  
  // Dir and file names.
  String data_dir = PRE + "/" + DATA;
  String run_dir = data_dir + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;

  String pairs_file = run_dir + "/" + READS + ".pairs";
  String aligns_file = sub_dir + "/" + ALIGNS + ".qltoutlet";
  String index_file =  sub_dir + "/" + ALIGNS + ".qltoutlet.index";

  String contigs_file = sub_dir + "/" + HEAD_IN + ".contigs.fasta";
  String supers_file = sub_dir + "/" + HEAD_IN + ".superb";

  String final_head =  sub_dir + "/" + HEAD_OUT;
  
  // Needed.
  vec<String> needed;
  needed.push_back( pairs_file );
  needed.push_back( aligns_file );
  needed.push_back( index_file );
  needed.push_back( contigs_file );
  needed.push_back( supers_file );
  if ( ! CheckFilesExist( needed, &cout ) ) return 1;
  
  // Load.
  cout << Date( ) << ": loading pairs" << endl;
  PairsManager pairs( pairs_file );
  size_t n_pairs = pairs.nPairs( );
  
  cout << Date( ) << ": loading aligns" << endl;
  vec<alignlet> aligns;
  BinaryReader::readFile( aligns_file, &aligns );
  
  cout << Date( ) << ": loading index" << endl;
  vec<int> index;
  BinaryReader::readFile( index_file, &index );
  const int n_reads = index.size( );
  
  cout << Date( ) << ": loading contigs" << endl;
  vec<fastavector> contigs;
  LoadFromFastaFile( contigs_file, contigs );
  const int n_contigs = contigs.size( );

  cout << Date( ) << ": loading supers" << endl;
  shandler supers( contigs.isize( ), supers_file );
  const int n_supers = supers.Size( );
    
  // Build maps.
  cout << Date( ) << ": digesting and sorting aligns" << endl;
  vec<CSARead> raligns;
  raligns.reserve( n_reads );
  for (int ii=0; ii<n_reads; ii++) {
    if ( index[ii] < 0 ) continue;
    raligns.push_back( CSARead( aligns[ index[ii] ], ii ) );
  }
  sort( raligns.begin( ), raligns.end( ) );

  vec<int> first_ralign( n_contigs, -1 );
  vec<int> n_aligns( n_contigs, 0 );
  vec<int> rindex( n_reads, -1 );
  for (int ii=raligns.isize( )-1; ii>=0; ii--) {
    int cid = raligns[ii].ContigId( );
    int rid = raligns[ii].ReadId( );
    n_aligns[cid] += 1;
    first_ralign[cid] = ii;
    rindex[rid] = ii;
  }

  vec< pair<int,int> > sep_ranges;
  vec<String> lib_names = pairs.getLibraryNames( );
  vec<PM_LibraryStats> lib_stats = pairs.getLibraryStats( );
  sep_ranges.reserve( lib_names.size( ) );
  for (int ii=0; ii<lib_names.isize( ); ii++) {
    int max_sep = LIB_MULTIPLIER * (double)lib_stats[ii].sep;
    sep_ranges.push_back( make_pair( MIN_SEP, max_sep ) );
  }

  vec< vec<String> > table;
  vec<String> line = MkVec( String( "lib_name" ),
			    String( "min_sep" ),
			    String( "max_sep" ) );
  table.push_back( line );
  for (int ii=0; ii<lib_names.isize( ); ii++) {
    line = MkVec( lib_names[ii],
		  ToString( sep_ranges[ii].first ),
		  ToString( sep_ranges[ii].second ) );
    table.push_back( line );
  }
  cout << "\n";
  PrintTabular( cout, table, 3, "rrr" );

  // Keep track of unplaced contigs.
  vec<bool> placed( n_contigs, false );

  // The final supers.
  vec<superb> out_supers;

  // Loop over all supers.
  cout << "\n" << Date( ) << ": loop over " << n_supers << " supers" << endl;
  for (int super_id=0; super_id<n_supers; super_id++) {
    
    // Estimate reads and pairs (needed to reserve links below).
    const int n_tigs = supers[super_id].Ntigs( );
    int n_reads_estim = 0;
    for (int tpos=0; tpos<n_tigs; tpos++)
      n_reads_estim += n_aligns[ supers[super_id].Tig( tpos ) ];
    int n_pairs_estim = ( ( 1 + n_reads_estim ) / 2 );
    
    // All consistent links.
    vec< pair<int,int> > links;
    links.reserve( n_pairs_estim );

    // Loop over all contigs in super.
    for (int tpos1=0; tpos1<n_tigs; tpos1++ ) {
      const int tig1 = supers[super_id].Tig( tpos1 );
      int first_ral = first_ralign[tig1];
      
      // Repeat contigs have no links, force-link them to adjacent contigs.
      if ( first_ral < 0 ) {
	if ( tpos1 > 0 )
	  for (int jj=0; jj<MIN_LINKS; jj++)
	    links.push_back( make_pair( tpos1-1, tpos1 ) );
	if ( tpos1 < n_tigs - 1 )
	  for (int jj=0; jj<MIN_LINKS; jj++)
	    links.push_back( make_pair( tpos1, tpos1+1 ) );
	continue;
      }
      
      // Loop over all reads in contig.
      for (int ral_id=first_ral; ral_id<raligns.isize( ); ral_id++) {
	const CSARead &ral1 = raligns[ral_id];
	if ( ral1.ContigId( ) != tig1 ) break;
      	if ( ral1.Rc( ) ) continue;

	int id1 = ral1.ReadId( );
	int id2 = pairs.getPartnerID( (longlong)id1 );
	int libid = pairs.libraryID( pairs.getPairID( (longlong)id1 ) );
	if ( rindex[id2] < 0 ) continue;

	const CSARead &ral2 = raligns[ rindex[id2] ];
	if ( ral2.Fw( ) ) continue;
	
	int tig2 = ral2.ContigId( );
	if ( tig2 == tig1 ) continue;
	if ( supers.ToSuper( tig2 ) != super_id ) continue;

	int tpos2 = supers.PosOnSuper( tig2 );
	int tig1len = supers[super_id].Len( tpos1 );
	int tig2len = supers[super_id].Len( tpos2 );
	int tig1start = supers.StartOnSuper( tig1 );
	int tig1end = tig1start + tig1len;
	int tig2start = supers.StartOnSuper( tig2 );
	int sep
	  =  tig1len - ral1.End( )
	  + tig2start - tig1end
	  + ral2.Begin( );

	bool out_of_range
	  = sep < sep_ranges[libid].first || sep > sep_ranges[libid].second;

	if ( VERBOSE )
	  cout << "s" << super_id << "\t"
	       << tpos1 << "-" << tpos2 << "\t"
	       << "sep: " << sep << "\t"
	       << ( out_of_range ? "out_of_range" : "ok" ) << "\n";
	
	if ( out_of_range ) continue;

	if ( tpos2 < tpos1 ) links.push_back( make_pair( tpos2, tpos1 ) );
	else links.push_back( make_pair( tpos1, tpos2 ) );
	
      } // loop over all reads in contig.

    } // loop over all contigs.
    
    sort( links.begin( ), links.end( ) );
    
    // Compactify links in bundles.
    vec< triple<int,int,int> > bundles;
    bundles.reserve( links.size( ) );
    for (int ii=0; ii<links.isize( ); ii++) {
      if ( bundles.size( ) < 1 ||
	   bundles.back( ).first != links[ii].first ||
	   bundles.back( ).second != links[ii].second ) {
	triple<int,int,int> newt( links[ii].first, links[ii].second, 0 );
	bundles.push_back( newt );
      }
      bundles[bundles.size( )-1].third += 1;
    }
    
    // Find connected components.
    equiv_rel er( n_tigs );
    for (int ii=0; ii<bundles.isize( ); ii++)
      if ( bundles[ii].third >= MIN_LINKS )
	er.Join( bundles[ii].first, bundles[ii].second );
    
    vec<int> reps;
    er.OrbitRepsAlt( reps );
    
    // One output super per orbit.
    for (int orbit_id=0; orbit_id<reps.isize( ); orbit_id++) {
      const superb &orig = supers[super_id];

      vec<int> select;
      er.Orbit( reps[orbit_id], select );
      sort( select.begin( ), select.end( ) );
      superb newsup;

      int tig = orig.Tig( select[0] );
      int len = orig.Len( select[0] );
      newsup.PlaceFirstTig( tig, len );
      placed[tig] = true;

      for (int ii=1; ii<select.isize( ); ii++) {
	int tig = orig.Tig( select[ii] );
	int len = orig.Len( select[ii] );
	int prev_tig = orig.Tig( select[ii-1] );
	int prev_len = orig.Len( select[ii-1] );
	int prev_end = prev_len + supers.StartOnSuper( prev_tig );
	int this_beg = supers.StartOnSuper( tig );
	int gap = this_beg - prev_end;
	int dev = Max( (int)( DEV_RATIO * (double)gap ), MIN_DEV );

	newsup.AppendTig( tig, len, gap, dev );
	placed[tig] = true;
      }

      out_supers.push_back( newsup );
    }

  } // loop over all supers.
  
  // Unplaced contigs.
  cout << Date( ) << ": adding unplaced contigs... " << flush;
  int n_unplaced = 0;
  for (int ii=0; ii<n_contigs; ii++) {
    if ( placed[ii] ) continue;
    n_unplaced++;
    superb newsup;
    int len = contigs[ii].size( );
    newsup.PlaceFirstTig( ii, len );
    out_supers.push_back( newsup );
  }
  cout << n_unplaced << " contigs added" << endl;
  
  // Save and leave.
  cout << Date( ) << ": saving" << endl;
  SaveScaffoldAssembly( final_head, supers.AllSupers( ), contigs );

  cout << Date( ) << ": BreakUnlinkedAssisted done" << endl;

}

