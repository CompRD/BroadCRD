/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Superb.h"
#include "paths/ReadLoc.h"

/**
 * Airiness
 *
 * Defined as: ( \sum Max( 0, gap_i ) ) / super_gapped_length.
 */
double Airiness( const superb &sup ) {
  double denumer = sup.TrueLength( );
  if ( denumer < 1.0 ) return 0.0;
  double numer = 0;
  for (int ii=0; ii<sup.Ntigs( )-1; ii++)
    numer += Max( 0, sup.Gap( ii ) );
  return numer / denumer;
}

/**
 * BreakAtGap
 *
 * Decide if the scaffold needs to be broken at gap <gap_id>.
 */
bool BreakAtGap( const int gap_id,
		 const superb &sup,
		 const vec<int> &to_super,
		 const vec<int> &to_superpos,
		 const vec<int> &start_on_super,
		 const vec< vec<read_loc> > &locs,
		 const int MIN_LINKS,
		 const double MAX_STRETCH,
		 ostream &log )
{
  // Number of valid links above the gap.
  int n_valid = 0;

  // Id of super.
  int super_id = -1;
  for (int ii=0; ii<locs.isize( ); ii++) {
    for (int jj=0; jj<locs[ii].isize( ); jj++) {
      super_id = to_super[ locs[ii][jj].ContigId( ) ];
      if ( super_id > -1 ) break;
    }
  }
  ForceAssertGe( super_id, 0 );

  // Loop over all contigs up to the gap.
  for (int cg_left=0; cg_left<=gap_id; cg_left++) {
    const vec<read_loc> &cg_locs = locs[cg_left];

    // Loop over all locs in contig.
    for (size_t loc_id=0; loc_id<cg_locs.size( ); loc_id++) {
      const read_loc &left_loc = cg_locs[loc_id];
      if ( left_loc.Rc( ) ) continue;
      if ( left_loc.PartnerFw( ) ) continue;
      int right_cg_id = left_loc.PartnerContigId( );
      if ( to_super[right_cg_id] != super_id ) continue;
      int cg_right = to_superpos[right_cg_id];
      if ( cg_right <= gap_id ) continue;

      int left_end = start_on_super[sup.Tig( cg_left )] + left_loc.Stop( );
      int right_beg = start_on_super[right_cg_id] + left_loc.PartnerStart( );
      int seen_sep = right_beg - left_end;
      double num_stretch = seen_sep - left_loc.Sep( );
      double denum_stretch = Max( 1, left_loc.Dev( ) );
      double stretch = Abs( num_stretch ) / denum_stretch;
      if ( stretch > MAX_STRETCH ) continue;
      
      ++n_valid;
      if ( n_valid >= MIN_LINKS ) return false;
    }
  }
  
  // No need to break.
  if ( n_valid >= MIN_LINKS ) return false;

  // Log breaking event.
  log << " s" << super_id
      << " (" << sup.ReducedLength( )
      << " bp, " << sup.Ntigs( )
      << " contigs) - break at tig " << gap_id
      << " (" << n_valid
      << " valid links)\n";
  
  // Done.
  return true;
  
}

/**
 * ReduceAiriness
 *
 * Reduce airiness in short scaffolds by breaking them at gaps spanned
 * by a low number of links.
 *
 * Only scaffolds with ungapped length < <MAX_UNGAPPED_LENGTH>, and
 * with airiness bigger than <MAX_AIRINESS> (see Airiness( ) above)
 * are considered for braking. A scaffold is broken at all gap which
 * are spanned by less than <MIN_LINKS> links (argument <MAX_STRETCH>
 * is used to define valid links).
 */ 
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_String_OrDefault( SCAFFOLDS_IN, "linear_scaffolds0.clean.patched" );
  CommandArgument_String_OrDefault( SCAFFOLDS_OUT, "linear_scaffolds0.clean.patched.reduced" );
  CommandArgument_Int_OrDefault( MAX_UNGAPPED_LENGTH, 200000 );
  CommandArgument_Double_OrDefault( MIN_AIRINESS, 0.15 );
  CommandArgument_Int_OrDefault( MIN_LINKS, 14 );
  CommandArgument_Double_OrDefault( MAX_STRETCH, 5.0 );
  EndCommandArguments;

  // Dir and file names.
  String full_run = PRE + "/" + DATA + "/" + RUN;
  String full_sub = full_run + "/ASSEMBLIES/" + SUBDIR;

  String head = full_sub + "/" + SCAFFOLDS_IN;
  String supers_file = head + ".superb";

  String supers_out_file = full_sub + "/" + SCAFFOLDS_OUT + ".superb";
  String summary_out_file = full_sub + "/" + SCAFFOLDS_OUT + ".summary";
  String log_file = full_sub + "/" + SCAFFOLDS_OUT + ".log";
  
  cout << " Logging to " << log_file << "\n" << endl;
  ofstream log( log_file.c_str( ) );
  PrintCommandPretty( log );

  // Load.
  read_locs_on_disk locs_file( head, full_run );

  vec<superb> supers;
  ReadSuperbs( supers_file, supers );

  int n_tigs = 0;
  for (int ii=0; ii<supers.isize( ); ii++)
    for (int jj=0; jj<supers[ii].Ntigs( ); jj++)
      n_tigs = Max( n_tigs, supers[ii].Tig( jj ) );
  ++n_tigs;
  
  vec<int> to_super( n_tigs, -1 );
  vec<int> to_superpos( n_tigs, -1 );
  vec<int> start_on_super( n_tigs, 0 );
  for (int ii=0; ii<supers.isize( ); ii++) {
    int start_pos = 0;
    for (int jj=0; jj<supers[ii].Ntigs( ); jj++) {
      to_super[ supers[ii].Tig( jj ) ] = ii;
      to_superpos[ supers[ii].Tig( jj ) ] = jj;
      start_on_super[ supers[ii].Tig( jj ) ] = start_pos;
      start_pos += supers[ii].Len( jj );
      if ( jj < supers[ii].Ntigs( )-1 )
	start_pos += supers[ii].Gap( jj );
    }
  }
  
  // Reduced supers (reserve memory, expecting a large number of breaks).
  vec<superb> reduced_supers;
  reduced_supers.reserve( 2 * supers.size( ) );

  // Parse selected supers.
  for (int super_id=0; super_id<supers.isize( ); super_id++) {
    const superb &sup = supers[super_id];
    double airiness = Airiness( sup );
    int ungapped_len = sup.ReducedLength( );
    if ( sup.Ntigs( ) < 2 
	 || airiness < MIN_AIRINESS
	 || ungapped_len > MAX_UNGAPPED_LENGTH ) {
      reduced_supers.push_back( sup );
      continue;
    }
    
    // Break points
    vec<int> breaks;

    // Locs for all contigs in super.
    vec< vec<read_loc> > locs( sup.Ntigs( ) );
    for (int ii=0; ii<sup.Ntigs( ); ii++)
      locs_file.LoadContig( sup.Tig( ii ), locs[ii] );
    
    // Parse each gap in super.
    for (int gap_id=0; gap_id<sup.Ntigs( )-1; gap_id++) {
      if ( BreakAtGap ( gap_id, sup, to_super, to_superpos, start_on_super,
			locs, MIN_LINKS, MAX_STRETCH, log ) ) {
	breaks.push_back( gap_id );
      }
    }

    // Create subsets.
    if ( breaks.size( ) < 1 ) {
      reduced_supers.push_back( sup );
      continue;
    }
    for (int ii=0; ii<=breaks.isize( ); ii++) {
      int begin = ( ii == 0 ) ? 0 : breaks[ii-1] + 1;
      int end = ( ii == breaks.isize( ) ) ? sup.Ntigs( ) : breaks[ii] + 1;
      vec<int> to_use;
      to_use.reserve( end - begin );
      for (int ii=begin; ii<end; ii++)
	to_use.push_back( ii );
      reduced_supers.push_back( sup.SubSuper( to_use ) );
    }
    
  }
  log << endl;

  // Save.
  cout << Date( ) << ": saving" << endl;
  WriteSuperbs( supers_out_file, reduced_supers );
  WriteSummary( summary_out_file, reduced_supers );
  
  // Done.
  cout << Date( ) << ": done" << endl;
  log << Date( ) << ": done" << endl;
  log.close( );
  
}
