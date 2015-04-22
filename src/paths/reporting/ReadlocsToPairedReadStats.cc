///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "math/Functions.h"
#include "paths/ReadLoc.h"

/**
 * CSimpleHisto
 *
 * Small histogram class.
 */
class CSimpleHisto {
  
public:

  CSimpleHisto( const vec<int> &bins ) { this->Setup( bins ); }

  // Increment freq_[ii] by one, where ii is defined by:
  //   bins_[ii-1] <= datum < bins_[ii]
  // with bins_[-1] := -\infty, and bins_[bins_.size( )] := +\infty.
  void AddDatum( int datum ) {
    vec<int>::iterator it = upper_bound( bins_.begin( ), bins_.end( ), datum );
    freq_[ distance( bins_.begin( ), it ) ]++;
  }

  int NBins( ) const { return bins_.isize( ); }

  int TotalCount( ) const { return Sum( freq_ ); }
  
  // For 0 <= ii <= bins_.size( ).
  int Freq( int ii ) const { return freq_[ii]; }
  
  // For 0 <= ii <= bins_.size( ).
  double Ratio( int ii ) const {
    return SafeQuotient( freq_[ii], this->TotalCount( ) );
  }
  
private:
  
  void Setup( const vec<int> &bins ) {
    bins_ = bins;
    sort( bins_.begin( ), bins_.end( ) );
    freq_.resize( 1 + bins_.size( ), 0 );
  }
  
  
private:

  vec<int> bins_;
  vec<int> freq_;   // has size bins_.size( ) + 1
  
};

/**
 * ReadlocsToPairedReadStats
 *
 * Evaluate paired reads statistics from a set of read_locs.
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String( ASSEMBLY );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_String_OrDefault( READLOCS_PREFIX, "" );
  CommandArgument_Int_OrDefault( MIN_SEP, -1200 );
  CommandArgument_Int_OrDefault( MAX_SEP, 0 );
  CommandArgument_Int_OrDefault( BIN_SIZE, 25 );
  EndCommandArguments;

  // Dir and file names.
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
  
  String contigs_fastb_file = sub_dir + "/" + ASSEMBLY + ".contigs.fastb";

  String head = sub_dir + "/" + ASSEMBLY;
  if ( READLOCS_PREFIX != "" ) head += "." + READLOCS_PREFIX;
  
  // Create bins.
  vec<int> bins( 1, MIN_SEP );
  while ( bins.back( ) < MAX_SEP )
    bins.push_back( bins.back( ) + BIN_SIZE );
  
  // Create hystogram containers (one per class).
  vec< vec<CSimpleHisto> > histograms( 3 );
  
  // Parse contigs.
  read_locs_on_disk locs_file( head, run_dir );

  int ntigs = MastervecFileObjectCount( contigs_fastb_file );
  cout << Date( ) << ": " << ntigs << " contigs found" << endl;
  for (int contig_id=0; contig_id<ntigs; contig_id++) {
    Dot( cout, contig_id );
    vec<read_loc> locs;
    locs_file.LoadContig( contig_id, locs );

    for (int loc_id=0; loc_id<locs.isize( ); loc_id++) {
      const read_loc &rloc = locs[loc_id];
      if ( rloc.ReadClass( ) < 1 ) continue;   // skip fragments
      if ( ! rloc.PartnerPlaced( ) ) continue;
      if ( rloc.PartnerContigId( ) != rloc.ContigId( ) ) continue;
      if ( ! rloc.PartnerRc( ) ) continue;
      if ( ! rloc.Fw( ) ) continue;   // avoid duplication
      
      int class_id = rloc.ReadClass( );
      int lib_id = rloc.LibraryId( );
      if ( histograms[class_id].isize( ) <= lib_id )
	histograms[class_id].resize( lib_id + 1, CSimpleHisto( bins ) );
      
      int sep = rloc.PartnerStart( ) - rloc.Stop( );
      histograms[class_id][lib_id].AddDatum( sep );
    }
  }
  cout << "\n" << endl;
  
  // Generate printable table.
  vec< vec<String> > table;

  vec<String> legend( 1, "library" );
  for (int class_id=0; class_id<3; class_id++) {
    String str_lib;
    if ( class_id == 0 ) str_lib = "frag";
    else if ( class_id == 1 ) str_lib = "jump";
    else str_lib = "long_jump";
    for (int lib_id=0; lib_id<histograms[class_id].isize( ); lib_id++) {
      str_lib += ".[" + ToString( lib_id ) + "]";
      legend.push_back( str_lib );
    }
  }
  table.push_back( legend );

  vec<String> tot_count( 1, "# valid pairs" );
  for (int class_id=0; class_id<3; class_id++) {
    for (int lib_id=0; lib_id<histograms[class_id].isize( ); lib_id++) {
      int tot_valid = histograms[class_id][lib_id].TotalCount( );
      tot_count.push_back( ToStringAddCommas( tot_valid ) + " (100%)" );
    }
  }
  table.push_back( tot_count );
  
  for (int ii=0; ii<=bins.isize( ); ii++) {
    vec<String> line;
    String str_bin;
    str_bin += ( ii == 0 ) ? "(-inf" : "[" + ToString( bins[ii-1] );
    str_bin += ", ";
    str_bin += ( ii == bins.isize( ) ) ? "+inf" : ToString( bins[ii] );
    str_bin += ")";
    line.push_back( str_bin );
    for (int class_id=0; class_id<3; class_id++) {
      for (int lib_id=0; lib_id<histograms[class_id].isize( ); lib_id++) {
	int tot = histograms[class_id][lib_id].TotalCount( );
	if ( tot < 1 ) {
	  line.push_back( "na" );
	  continue;
	}
	int freq = histograms[class_id][lib_id].Freq( ii );
	double percent = 100.0 * histograms[class_id][lib_id].Ratio( ii );
	line.push_back( ToStringAddCommas( freq )
			+ " (" + ToString( percent, 1 ) + "%)" );
      }
    }
    table.push_back( line );
  }

  // Print table.
  String justif = "rr";
  for (int class_id=0; class_id<3; class_id++) {
    for (int lib_id=0; lib_id<histograms[class_id].isize( ); lib_id++) {
      justif += "r";
    }
  }

  PrintTabular( cout, table, 3, justif );
  cout << endl;
  
  // Done.
  cout << Date( ) << ": done" << endl;
  
}
