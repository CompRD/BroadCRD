///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"
#include "PairsManager.h"
#include "SupersHandler.h"
#include "VecUtilities.h"
#include "paths/ReadLoc.h"
#include "paths/assisted/SPlacem.h"
#include "random/RNGen.h"
#include "util/RunCommand.h"

/**
 * FindPlacements
 */
void FindPlacements( const int64_t read_id,
		     const shandler &supers,
		     const vec<int> &jlens,
		     const vec< vec<int> > &uid_to_cid,
		     const vec< triple<int64_t,int64_t,int> > &jaligns,
		     vec<SPlacem> &placements )
{
  placements.clear( );

  triple<int64_t,int64_t,int> bogus( read_id, -1 , 0 );
  vec< triple<int64_t,int64_t,int> >::const_iterator it
    = lower_bound( jaligns.begin( ), jaligns.end( ), bogus );

  // Unaligned read.
  if ( it == jaligns.end( ) || it->first != read_id  ) return;
  
  // This unibase is not a contig.
  int uid = it->second;
  int upos = ( it->third < 0 ) ? - it->third - 1 : it->third;
  bool urc = ( it->third < 0 );
  const vec<int> &cids = uid_to_cid[uid];
  const int n_placs = cids.size( );
  if ( n_placs < 1 ) return;

  // The placements.
  int read_len = jlens[read_id];
  placements.reserve( cids.size( ) );
  for (int ii=0; ii<cids.isize( ); ii++) {
    int32_t cid = ( cids[ii] < 0 ) ? - cids[ii] - 1 : cids[ii];
    bool crc =  ( cids[ii] < 0 );
    int sid = supers.ToSuper( cid );
    int beg = supers.StartOnSuper( cid ) + upos;
    int end = beg + read_len;
    bool rc = ( urc != crc );
    
    placements.push_back( SPlacem( read_id, cids[ii], sid, beg, end, rc ) );
  }
  
}

/**
 * FindLogicalMates
 *
 * It randomly select one logical pair from the set of all possible
 * logical pairs (or a pairs of null pointers, if none are found).
 */
void FindLogicalMates( const int MIN_SEP,
		       const int MAX_SEP,
		       const vec<SPlacem> &placs1,
		       const vec<SPlacem> &placs2,
		       RNGen &rgen,
		       pair<const SPlacem*, const SPlacem*> &mates )
{
  const SPlacem *answer1 = 0;
  const SPlacem *answer2 = 0;
  mates = make_pair( answer1, answer2 );

  vec< pair<const SPlacem*, const SPlacem*> > all_mates;
  for (int id1=0; id1<placs1.isize( ); id1++)
    for (int id2=0; id2<placs2.isize( ); id2++)
      if ( placs1[id1].IsLogicalPair( MIN_SEP, MAX_SEP, placs2[id2] ) )
	all_mates.push_back( make_pair( &( placs1[id1] ), &( placs2[id2] ) ) );
  
  if ( all_mates.size( ) < 1 ) return;
  
  int select = rgen.random( ) * all_mates.size( ) / RNGen::RNGEN_RAND_MAX;
  mates = all_mates[select];
}

/**
 * ReportMultiplicities
 *
 * Print a table with stats on contigs with multiplicities.
 */
void ReportMultiplicities( const vec< vec<int> > &uid_to_cid,
			   const vec<int> &clens,
			   ostream &out )
{
  int ntigs = clens.size( );
  int nunis = uid_to_cid.size( );
  longlong tot_clen = BigSum( clens );

  int max_mult = 0;
  for (int ii=0; ii<nunis; ii++)
    max_mult = Max( max_mult, uid_to_cid[ii].isize( ) );
  
  vec<int> mult_unicounts( max_mult + 1, 0 );
  vec<int> mult_counts( max_mult + 1, 0 );
  vec<longlong> mult_lens( max_mult + 1, 0 );
  for (int ii=0; ii<nunis; ii++) {
    int mult = uid_to_cid[ii].size( );
    mult_unicounts[mult] += 1;
    mult_counts[mult] += mult;
    if ( mult > 0 ) {
      int tig = uid_to_cid[ii][0];
      if ( tig < 0 ) tig = - tig - 1;
      mult_lens[mult] += ( mult * clens[tig] );
    }
  }
 
  vec< vec<String> > table;

  vec<String> line = MkVec( String( "mult" ),
			    String( "n_unis" ),
			    String( "n_tigs" ),
			    String( "tig_length" ),
			    String( "tig_length_%" ) );
  table.push_back( line );

  line = MkVec( String( "0" ),
		ToStringAddCommas( mult_lens[0] ),
		String( "na" ),
		String( "na" ),
		String( "na" ) );
  table.push_back( line );

  for (int mult=1; mult<mult_counts.isize( ); mult++) {
    double ratio = 100. * SafeQuotient( mult_lens[mult], tot_clen );
    line = MkVec( ToStringAddCommas( mult ),
		  ToStringAddCommas( mult_unicounts[mult] ),
		  ToStringAddCommas( mult_counts[mult] ),
		  ToStringAddCommas( mult_lens[mult] ),
		  ToString( ratio, 1 ) );
    table.push_back( line );
  }

  line = MkVec( String( "tot" ),
		ToStringAddCommas( nunis ),
		ToStringAddCommas( ntigs ),
		ToStringAddCommas( tot_clen ),
		String( "100.0" ) );
  table.push_back( line );

  out << "OVERALL CONTIGS MULTIPLICITY\n"
      << "NB: " << ToStringAddCommas( mult_unicounts[0] )
      << " unibases do not appear as contigs\n\n";
  PrintTabular( out, table, 3, "rrrrr" );
  out << endl;
}

/**
 * ReportPlacements
 *
 * Print a table with stats on placements of paired reads.
 */
void ReportPlacements( const PairsManager &pairs,
		       const vec<read_loc> &locs,
		       ostream &out )
{
  const int npairs = pairs.nPairs( );
  const int nlocs = locs.size( );

  vec<int> nplaced( npairs, 0 );
  for (int loc_id=0; loc_id<nlocs; loc_id++) {
    int id1 = locs[loc_id].ReadId( );
    int pid = pairs.getPairID( id1 );
    if ( ! locs[loc_id].PartnerPlaced( ) ) {
      nplaced[pid] = 1;
      continue;
    }
    nplaced[pid] = 2;
  }

  vec<int> ntots( 3, 0 );
  for (int ii=0; ii<npairs; ii++)
    ntots[ nplaced[ii] ]++;
  
  vec< vec<String> > table;

  vec<String> line = MkVec( String( "plac" ),
			    String( "n_pairs" ),
			    String( "%_pairs" ) );
  table.push_back( line );
  
  for (int ii=0; ii<ntots.isize( ); ii++) {
    double ratio = 100. * SafeQuotient( ntots[ii], npairs );
    line = MkVec( ToString( ii ),
		  ToStringAddCommas( ntots[ii] ),
		  ToString( ratio, 1 ) );
    table.push_back( line );
  }

  line = MkVec( String( "tot" ),
		ToStringAddCommas( npairs ),
		String( "100.0" ) );
  table.push_back( line );

  out << "OVERALL PLACEMENT RATE\n\n";
  PrintTabular( out, table, 3, "rrr" );
  out << endl;

}  

/**
 * BuildReadLocs
 *
 * Code specific to the assisting package (see "key assmptions" below):
 * build readlocs from a set of aligns (passed as triples, as from the
 * output of SearchFastb2).
 *
 * KEY ASSUMPTIONS:
 *
 *   1) each contig is a perfect copy of one unibase. Different
 *      contigs may be the copy of the same unibase (".orig.ids" maps
 *      each contig to its original unibase id). Not all unibases need
 *      to be "used" (ie, contigs are a subset of the unibases).
 *
 *   2) JUMPS are error corrected reads, and each read aligns on
 *      exactly one unibase.
 *
 * INPUT:
 *
 *   <RUN_DIR>/<JUMPS>.aligns
 *   <RUN_DIR>/<READS>.unibases.k<K>
 *   <RUN_DIR>/<ASSEMBLIES>/<SUB_DIR>/<ASSEMBLY>.contigs.fastb
 *   <RUN_DIR>/<ASSEMBLIES>/<SUB_DIR>/<ASSEMBLY>.orig.ids
 *   <RUN_DIR>/<ASSEMBLIES>/<SUB_DIR>/<ASSEMBLY>.superb
 *
 * OUTPUT:
 *
 *   <RUN_DIR>/<ASSEMBLIES>/<SUB_DIR>/<ASSEMBLY>.readlocs
 *
 */
int main( int argc, char *argv[] )
{
  RunTime( );
 
  BeginCommandArguments;
  CommandArgument_String( RUN_DIR );
  CommandArgument_String( SUB_DIR );
  CommandArgument_String( ASSEMBLY );
  CommandArgument_Int_OrDefault( K, 96 );
  CommandArgument_Int_OrDefault( MIN_SEP, 150 );
  CommandArgument_Int_OrDefault( MAX_SEP, 15000 );
  CommandArgument_String_OrDefault( READS, "all_reads" );
  CommandArgument_String_OrDefault( JUMPS, "jump_reads_ec" );
  CommandArgument_UnsignedInt_OrDefault( RANDOM_SEED, 666 );
  EndCommandArguments;

  String sub_dir = RUN_DIR + "/ASSEMBLIES/" + SUB_DIR;
  
  String strK = ToString( K );
  String jpairs_file = RUN_DIR + "/" + JUMPS + ".pairs";
  String jreads_file = RUN_DIR + "/" + JUMPS + ".fastb";
  String jaligns_file = RUN_DIR + "/" + JUMPS + ".aligns";
  String unibases_file = RUN_DIR + "/" + READS + ".unibases.k" + strK;
  
  String assembly_head = sub_dir + "/" + ASSEMBLY;
  String contigs_file = assembly_head + ".contigs.fastb";
  String orig_ids_file = assembly_head + ".orig.ids";
  String supers_file = assembly_head + ".superb";
  String locs_head = assembly_head + ".contigs";

  vec<String> needed;
  needed.push_back( jpairs_file );
  needed.push_back( jreads_file );
  needed.push_back( jaligns_file );
  needed.push_back( contigs_file );
  needed.push_back( orig_ids_file );
  needed.push_back( supers_file );
  
  // Load.
  size_t n_unibases = MastervecFileObjectCount( unibases_file );
  size_t n_contigs = MastervecFileObjectCount( contigs_file );
  size_t n_jreads = MastervecFileObjectCount( jreads_file );

  cout << Date( ) << ": loading aligns" << endl;
  vec< triple<int64_t,int64_t,int> > jaligns;
  BinaryReader::readFile( jaligns_file.c_str( ), &jaligns );
  const size_t njaligns_total = jaligns.size( );
  const String str_njaligns_total = ToStringAddCommas( njaligns_total );
  
  cout << Date( ) << ": loading reads lengths" << endl;
  vec<int> jlens;
  {
    jlens.reserve( n_jreads );
    vecbvec jreads( jreads_file );
    for (size_t ii=0; ii<jreads.size( ); ii++)
      jlens.push_back( jreads[ii].size( ) );
  }

  cout << Date( ) << ": loading orig.ids" << endl;
  READ( orig_ids_file, vec<int>, orig_ids );
  ForceAssertEq( orig_ids.size( ), n_contigs );

  cout << Date( ) << ": loading supers" << endl;
  shandler supers( (int)n_contigs, supers_file );
  const int nsupers = supers.Size( );
  const String str_nsupers = ToStringAddCommas( nsupers );
  
  cout << Date( ) << ": loading pairs" << endl;
  PairsManager jpairs( jpairs_file );
  const size_t n_pairs = jpairs.nPairs( );

  // Make sure jaligns is sorted.
  cout << Date( ) << ": sorting " << str_njaligns_total << " aligns" << endl;
  if ( ! is_sorted( jaligns.begin( ), jaligns.end( ) ) )
    sort( jaligns.begin( ), jaligns.end( ) );

  // Generate maps
  cout << Date( ) << ": building maps" << endl;
  vec< vec<int> > uid_to_cid( n_unibases );
  for (int cid=0; cid<(int)n_contigs; cid++) {
    int uid = ( orig_ids[cid] < 0 ) ? - orig_ids[cid] - 1 : orig_ids[cid];
    bool rc = ( orig_ids[cid] < 0 );
    uid_to_cid[uid].push_back( rc ? - cid - 1 : cid );
  }
  
  vec<int> clens( n_contigs, 0 );
  for (int ii=0; ii<supers.Size( ); ii++)
    for (int jj=0; jj<supers[ii].Ntigs( ); jj++)
      clens[supers[ii].Tig( jj )] = supers[ii].Len( jj );

  // Random generator.
  RNGen rgen( RANDOM_SEED );

  // The read locations.
  vec<read_loc> locs;
  locs.reserve( 2 * n_pairs );
  
  // Loop over all pairs.
  size_t dotter = 100000;
  cout << Date( ) << ": building read locs\n"
       << "\n"
       << "Parsing " << ToStringAddCommas( n_pairs )
       << " pairs (" << ToStringAddCommas( dotter )
       << " pairs per dot)" << endl;
  for (size_t pid=0; pid<n_pairs; pid++) {
    if ( pid % dotter == 0 ) Dot( cout, pid / dotter );

    // Read ids.
    const int64_t id1 = jpairs.ID1( pid );
    const int64_t id2 = jpairs.ID2( pid );

    // All placements.
    vec<SPlacem> plac1;
    FindPlacements( id1, supers, jlens, uid_to_cid, jaligns, plac1 );
    
    vec<SPlacem> plac2;
    FindPlacements( id2, supers, jlens, uid_to_cid, jaligns, plac2 );
    
    // No aligns found.
    if ( plac1.size( ) < 1 && plac2.size( ) < 1 ) continue;
    
    // One mate aligns, the other does not.
    if ( plac1.size( ) < 1 || plac2.size( ) < 1 ) {
      const vec<SPlacem> &placs = ( plac1.size( ) < 1 ) ? plac2 : plac1;      

      // Multiple placements, we do not know which to pick, leave.
      if ( placs.size( ) > 1 ) continue;

      // Add read_loc (single entry, partner unplaced), and continue.
      placs[0].AddReadLoc( jaligns, clens, jpairs, locs );
      continue;
    }

    // Both reads align, pick a random logical set.
    pair<const SPlacem*, const SPlacem*> mates;
    FindLogicalMates( MIN_SEP, MAX_SEP, plac1, plac2, rgen, mates );
    if ( mates.first != 0 && mates.second != 0 ) {
      const SPlacem &partner = *( mates.second );
      mates.first->AddReadLocs( jaligns, clens, jpairs, partner, locs );
    }
  }
  cout << "\n" << endl;

  // Sort and save read locs.
  cout << Date( ) << ": sorting and saving read locs\n" << endl;
  sort( locs.begin( ), locs.end( ) );
  WriteReadLocs( locs_head, n_contigs, locs );
  
  // Print some stats.
  ReportMultiplicities( uid_to_cid, clens, cout );
  ReportPlacements( jpairs, locs, cout );
  
  // Done.
  cout << Date( ) << ": BuildReadLocs done" << endl;
  
}

