///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"

#include "Fastavector.h"
#include "Histogram.h"
#include "PairsManager.h"
#include "VecUtilities.h"
#include "paths/Alignlet.h"
#include "paths/RemoveDuplicateAligns.h"
#include "paths/ScaffoldsUtils.h"
#include "feudal/BinaryStream.h"

/**
 * class CompactLibInfo
 */
class CompactLibInfo {
  
public:
  
  CompactLibInfo( ) : 
    n_pairs_ ( 0 ), n_aligned_ ( 0 ), n_innies_ ( 0 ) { }
  
  void SetNPairs( size_t n_pairs ) { n_pairs_ = n_pairs; }
  void SetNAligned( size_t n_aligned ) { n_aligned_ = n_aligned; }
  void SetNInnies( size_t n_innies ) { n_innies_ = n_innies; }
  
  size_t NPairs( ) const { return n_pairs_; }
  size_t NAligned( ) const { return n_aligned_; }
  size_t NInnies( ) const { return n_innies_; }
  size_t NValid( ) const { return n_aligned_ - n_innies_; }
  
  
private:
  
  size_t n_pairs_;
  size_t n_aligned_;
  size_t n_innies_;
  
};

/**
 * PlotLibsCore
 *
 * Generate eps and txt output files.
 */
void PlotLibsCore( const PairsManager &pairs,
		   const String &base_out,
		   const vec<int> &lens,
		   const int lib_id,
		   CompactLibInfo &lib_info )
{
  // Vector of lengths must be sorted.
  ForceAssert( is_sorted( lens.begin( ), lens.end( ) ) );

  // These are used to build histograms.
  const float MAX_MULT = 5.0;
  const vec<float> PROGRESSION = MkVec( float(25), float(50), float(100) );
  const int MIN_BINS = 60;
  const int MAX_BINS = 120;
  
  // Core data.
  String lib_name = pairs.getLibraryName( lib_id );
  int sep = pairs.getLibrarySep( lib_id );
  int dev = pairs.getLibrarySD( lib_id );

  size_t count = lens.size( );
  float radius = float( dev ) * MAX_MULT;
  float bin_begin = float( sep ) - radius;
  float bin_end = float( sep ) + radius;

  // File names.
  String raw_file = base_out + ".raw";
  String eps_file = base_out + ".eps";
  String histo_file = base_out + ".histo";

  // Raw data file.
  ofstream lout( raw_file.c_str( ) );
  for (size_t ii=0; ii<count; ii++)
    lout << lens[ii] << "\n";
  lout.close( );
  
  // Pick bin size.
  float bin_size = 0.0;
  float pow = 1.0;
  while ( 1 ) {
    for (int ii=0; ii<PROGRESSION.isize( ); ii++) {
      int nbins = int( ( bin_end - bin_begin ) / ( pow * PROGRESSION[ii] ) );
      if ( MIN_BINS <= nbins && nbins < MAX_BINS ) {
	bin_size = pow * PROGRESSION[ii];
	break;
      }
    }
    if ( bin_size > 0 ) break;
    pow *= 10;
  }
  
  // Generate x_bars.
  vec<float> x_bars;
  x_bars.reserve( radius / bin_size );
  int x_bar = 0;
  while ( x_bar < sep - radius )
    x_bar += bin_size;
  while ( x_bar <= sep + radius ) {
    x_bars.push_back( x_bar );
    x_bar += bin_size;
  }

  // AddDatum will silently discard innies, by storing them as Underflow.
  size_t n_innies = 0;
  histogram<float> histo;
  for (int ii=0; ii<x_bars.isize( ); ii++)
    histo.AddBin( x_bars[ii] );
  for (size_t  ii=0; ii<count; ii++) {
    if ( lens[ii] < 0 ) n_innies++;
    histo.AddDatum( lens[ii] );
  }
  
  // Overall stats.
  lib_info.SetNPairs( pairs.getLibrarySizes( )[lib_id] );
  lib_info.SetNAligned( lens.size( ) );
  lib_info.SetNInnies( n_innies );
    
  // Histogram's legend.
  color col = black;
  String str_mean = ToString( sep ) + " +/- " + ToString( dev );
  String str_sample = "(sample size: " + ToString( count ) + " inserts)";
  vec<ns_psplot::freetext> labels( 3 );
  labels[0] = ns_psplot::freetext( lib_name, col, 12 );
  labels[1] = ns_psplot::freetext( str_mean, col, 10 );
  labels[2] = ns_psplot::freetext( str_sample, col, 10 );

  // Generate eps file.
  ofstream eout( eps_file.c_str( ) );
  histo.PrintAsEps( eout, labels, 1, &sep, &dev );
  eout.close( );

  // Generate text histogram file.
  ofstream tout( histo_file.c_str( ) );
  histo.PrintAsColumns( tout, true, true );
  tout.close( );
  
}

/**
 * PlotLibsOnRef
 *
 * Realign reads onto a reference, and generate histograms of observed
 * insert lengths.  Alignments are cached, and duplicate molecules are
 * filtered out. Output is saved in <READS>.PlotLibsOnRef.
 * 
 * K: used by AlignReadsToContigs
 * READS: it loads <READS>.{fastb,pairs}
 * REF: it loads <REF>.fasta
 * REMOVE_DUPLICATES: run RemoveDuplicateAligns.h
 * FORCE: do not used cached alignments, regenerate them
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_Int( K );
  CommandArgument_String( READS );
  CommandArgument_String( REF );
  CommandArgument_Bool_OrDefault( REMOVE_DUPLICATES, True );
  CommandArgument_Bool_OrDefault( FORCE, False );
  EndCommandArguments;

  // Dir and file names.
  String out_dir = READS + ".PlotLibsOnRef";
  
  String reads_file = READS + ".fastb";
  String pairs_file = READS + ".pairs";
  String contigs_file = REF + ".fasta";

  String log_file = out_dir + "/main.log";
  String all_stats_file = out_dir + "/allstats.txt";
  String aligner_log_file = out_dir + "/AlignReadsToContigs.log";
  String aligns_file = out_dir + "/aligns.qltoutlet";
  String index_file = out_dir + "/aligns.index";
  
  // Output dir, log stream.
  Mkpath( out_dir );
  ofstream log( log_file.c_str( ) );
  PrintCommandPretty( log );

  cout << "Sending log to " << log_file << "\n" << endl;

  // Load.
  log << Date( ) << ": loading pairing info" << endl;
  PairsManager pairs( pairs_file );
  
  log << Date( ) << ": loading contigs fasta" << endl;
  vec<fastavector> contigs;
  LoadFromFastaFile( contigs_file, contigs );

  // Lib stats, as from the pairs file.
  {
    ofstream stats_out( all_stats_file.c_str( ) );
    stats_out << "Library statistics, as from the pairs file:\n\n";
    for (size_t ii=0; ii<pairs.nLibraries( ); ii++)
      stats_out << ii << "\t"
		<< pairs.getLibraryName( ii ) << "\t"
		<< pairs.getLibrarySep( ii ) << " +/- "
		<< pairs.getLibrarySD( ii ) << "\n";
    stats_out << endl;
    stats_out.close( );
  }

  // Align reads, or load aligns.
  vec<alignlet> aligns;
  vec<int> index;
  if ( FORCE || ! IsRegularFile( aligns_file ) ) {
    log << Date( ) << ": aligning read and saving aligns" << endl;
    ofstream alog( aligner_log_file.c_str( ) );

    alog << Date( ) << ": starting AlignReadsToContigs" << endl;
    AlignReadsToContigs( K, out_dir, reads_file, contigs, aligns, index, alog );

    if ( REMOVE_DUPLICATES )
      RemoveDuplicateAligns( pairs, aligns, index, alog );

    alog << Date( ) << ": saving aligns" << endl;
    BinaryWriter::writeFile( aligns_file, aligns );
    BinaryWriter::writeFile( index_file, index );

    alog << Date( ) << ": done" << endl;
    alog.close( );
  }
  else {
    log << Date( ) << ": loading aligns" << endl;
    BinaryReader::readFile( aligns_file, &aligns );
    BinaryReader::readFile( index_file, &index );
  }
  
  // The lengths of valid inserts.
  vec< vec<int> > lens( pairs.nLibraries( ) );

  // Reserve memory.
  vec<size_t> lib_sizes = pairs.getLibrarySizes( );
  for (size_t lib_id=0; lib_id<pairs.nLibraries( ); lib_id++)
    lens[lib_id].reserve( lib_sizes[lib_id] );

  // Loop over all pairs.
  log << Date( ) << ": parsing pairs" << endl;
  for (size_t pair_id=0; pair_id<pairs.nPairs( ); pair_id++) {
    int id1 = pairs.ID1( pair_id );
    int id2 = pairs.ID2( pair_id );
    if ( index[id1] < 0 || index[id2] < 0 ) continue;

    const alignlet &al1 = aligns[ index[id1] ];
    const alignlet &al2 = aligns[ index[id2] ];
    if ( al1.TargetId( ) != al2.TargetId( ) ) continue;
    if ( al1.Fw1( ) == al2.Fw1( ) ) continue;

    const alignlet &al_fw = al1.Fw1( ) ? al1 : al2;
    const alignlet &al_rc = al1.Fw1( ) ? al2 : al1;
    int lib_id = pairs.libraryID( pair_id );
    int sep = al_rc.pos2( ) - al_fw.Pos2( );

    lens[lib_id].push_back( sep );
  }
  
  // Sort.
  log << Date( ) << ": sorting" << endl;
  for (size_t lib_id=0; lib_id<pairs.nLibraries( ); lib_id++)
    sort( lens[lib_id].begin( ), lens[lib_id].end( ) );
  
  // Generating output files.
  vec<CompactLibInfo> libs_info( pairs.nLibraries( ) );
  for (size_t lib_id=0; lib_id<pairs.nLibraries( ); lib_id++) {
    String lib_name = pairs.getLibraryName( lib_id );
    String base_out = out_dir + "/" + lib_name;
    
    log << Date( ) << ": " << lib_name << endl;
    PlotLibsCore( pairs, base_out, lens[lib_id], lib_id, libs_info[lib_id] );
  }

  // Overall stats.
  vec< vec<String> > table;
  vec<String> line;
  
  vec<String> legend = MkVec( String( "lib_name" ),
			      String( "lib_stats" ),
			      String( "n_pairs" ),
			      String( "n_aligned" ),
			      String( "n_innies (%)" ),
			      String( "n_valid (%)" ) );
  table.push_back( legend );

  for (int ii=0; ii<libs_info.isize( ); ii++) {
    int sep = pairs.getLibrarySep( ii );
    int dev = pairs.getLibrarySD( ii );
    String sname = pairs.getLibraryName( ii );
    String sstats = ToString( sep ) + " +/- " + ToString( dev );
    size_t npairs = libs_info[ii].NPairs( );
    size_t nal = libs_info[ii].NAligned( );
    size_t ninn = libs_info[ii].NInnies( );
    size_t nval = libs_info[ii].NValid( );
    double rinn = 100.0 * double( ninn ) / double( nal );
    double rval = 100.0 * double( nval ) / double( nal );
    String snpairs = ToStringAddCommas( npairs );
    String snal = ToStringAddCommas( nal );
    String sinn = ToStringAddCommas( ninn ) + " (" + ToString( rinn, 1 ) + ")";
    String sval = ToStringAddCommas( nval ) + " (" + ToString( rval, 1 ) + ")";

    line = MkVec( sname, sstats, snpairs, snal, sinn, sval );
    table.push_back( line );
  }

  log << "\n";
  PrintTabular( log, table, 3, "rrrrrr" );
  log << endl;

  // Done.
  String str_done = Date( ) + ": done";
  cout << str_done << endl;
  log << str_done << endl;

}

