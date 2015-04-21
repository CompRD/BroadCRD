///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "PairsManager.h"
#include "TokenizeString.h"
#include "paths/assisted/CShowLinks.h"

/**
 * struct SLibInfo
 * Set
 */
void SLibInfo::Set( const String &pairs_file )
{
  ReadPairsManagerLibInfo( pairs_file, nreads_, names_, seps_, devs_ );
}

/**
 * struct SLinkInfo
 * Constructor
 */
SLinkInfo::SLinkInfo( int rclass,
		      int lib,
		      bool rc1,
		      bool rc2,
		      int sep,
		      int dev ) :
  class_ ( rclass ),
  lib_ ( lib ),
  rc1_ ( rc1 ),
  rc2_ ( rc2 ),
  sep_ ( sep ),
  dev_ ( dev )
{
  ForceAssert( class_ == 1 || class_ == 2 );
}
  
/**
 * struct SLinkInfo
 * ToLine
 */
vec<String> SLinkInfo::ToLine( ) const {
  vec<String> line;
  
  String str_or1 = rc1_ ? "-" : "+";
  String str_or2 = rc2_ ? "-" : "+";
  line.push_back( "  " + str_or1 + str_or2 );
  
  String str_class = "[" + String( class_ == 1 ? "j" : "J" ) + "]";
  String str_lib =  "[" + ToString( lib_ ) + "]";
  line.push_back( str_class + str_lib );
  
  if ( dev_ > 0 ) line.push_back( ToString( sep_ ) );
  else line.push_back( "" );
  
  if ( dev_ > 0 ) line.push_back( ToString( dev_ ) );
  else line.push_back( "" );
  
  return line;
}

/**
 * struct SLinkInfo
 * operator<
 */
bool operator< ( const SLinkInfo &left, const SLinkInfo &right )
{
  if ( left.rc1_ < right.rc1_ ) return true;
  if ( left.rc1_ > right.rc1_ ) return false;
  
  if ( left.rc2_ < right.rc2_ ) return true;
  if ( left.rc2_ > right.rc2_ ) return false;
  
  if ( left.class_ < right.class_ ) return true;
  if ( left.class_ > right.class_ ) return false;
  
  if ( left.lib_ < right.lib_ ) return true;
  if ( left.lib_ > right.lib_ ) return false;
  
  return ( left.sep_ < right.sep_ );
}

/**
 * CShowLinks
 * Constructor
 */
CShowLinks::CShowLinks( const vecbvec &bases,
			const SLibInfo &ijumps,
			const SLibInfo &iJumps,
			const vec<look_align> &hits,
			read_locs_on_disk &locs_file ) :
  bases_ ( bases ),
  ijumps_ ( ijumps ),
  iJumps_ ( iJumps ),
  hits_ ( hits ),
  locs_file_ ( locs_file )
{
  this->Setup( );
}

/**
 * CShowLinks
 * GoInteractive
 */
void CShowLinks::GoInteractive( ) const
{
  while( 1 ) {
    cout << "> " << flush;
    
    // Wait for command.
    String cmnd;
    getline( cin, cmnd );
    if ( ! cin ) {
      cout << endl;
      break;
    }
    
    // "q" (quit).
    if ( cmnd == "q" ) break;
    
    // Parse command string.
    vec<String> tokens;
    Tokenize( cmnd, tokens );
    if ( tokens.size( ) < 4 ) {
      cout << "malformed command\n\n";
      continue;
    }

    int tig1 = tokens[0].Int( );
    int tig2 = tokens[2].Int( );
    bool rct1 = ( tokens[1] == "+" ) ? false : true;
    bool rct2 = ( tokens[3] == "+" ) ? false : true;

    int c1 = rct1 ? - tig1 - 1 : tig1;
    int c2 = rct2 ? - tig2 - 1 : tig2;

    String out_file = "";
    if ( tokens.size( ) >= 5 )
      out_file = tokens[4];

    ofstream fout;
    if ( out_file != "" ) {
      fout.open( out_file.c_str( ) );
      cout << "\nSending output to " << out_file << "\n" << endl;
    }
    ostream &out = ( out_file != "" ) ? fout : * (ostream *) &cout;

    this->ShowLinks( c1, c2, out );
  }
  
}

/**
 * CShowLinks
 * Setup
 * private
 *
 * As usual, -1: unmapped, -2: multiply mapped
 */
void CShowLinks::Setup( )
{
  int n_tigs = bases_.size( );
  
  to_hit_.resize( n_tigs, -1 );
  for (int ii=0; ii<hits_.isize( ); ii++) {
    int cid = hits_[ii].query_id;
    if ( to_hit_[cid] == -1 ) to_hit_[cid] = ii;
    else to_hit_[cid] = -2;
  }

}

/**
 * CShowLinks
 * FindLinks
 * private
 */
void CShowLinks::FindLinks( int c1, int c2, vec<SLinkInfo> &slinks ) const
{
  slinks.clear( );

  // Unpack ids, load locs.
  int tig1 = c1 < 0 ? - c1 - 1 : c1;
  int tig2 = c2 < 0 ? - c2 - 1 : c2;
  bool rct1 = ( c1 < 0 );
  bool rct2 = ( c2 < 0 );
  
  vec<read_loc> locs;
  locs_file_.LoadContig( tig1, locs );
  
  // Contig lengths.
  const int lent1 = bases_[tig1].size( );
  const int lent2 = bases_[tig2].size( );

  // Condensate linking info (loop over all locs).
  for (int loc_id=0; loc_id<locs.isize( ); loc_id++) {
    const read_loc &loc = locs[loc_id];
    const int rclass = loc.ReadClass( );
    const int lib = loc.LibraryId( ); 
    const int dev = ( rclass == 1 ) ? ijumps_.devs_[lib] : iJumps_.devs_[lib];
    if ( (int)loc.PartnerContigId( ) != tig2 ) continue;
      
    // Absolute orientation of reads.
    bool partFw = loc.PartnerFw( );
    int fw1 = ( ( loc.Fw( ) && ( !rct1 ) ) || ( ( !loc.Fw( ) ) && rct1 ) );
    int fw2 = ( ( partFw && ( !rct2 ) ) || ( ( !partFw ) && rct2 ) );
    int sep = 0;
      
    // Case 1: fw-fw, or rc-rc (set dev to 0).
    if ( fw1 == fw2 ) {
      slinks.push_back( SLinkInfo( rclass, lib, !fw1, !fw2, sep, 0 ) );
      continue;
    }

    // Case 2: fw-rc.
    int part_start = loc.PartnerStart( );
    int part_stop = loc.PartnerStop( ); 
    if ( fw1 ) {
      sep = ( rclass == 1 ) ? ijumps_.seps_[lib] : iJumps_.seps_[lib];
      sep += ( rct1 ) ? - loc.Start( ) : - ( lent1 - loc.Stop( ) );
      sep += ( rct2 ) ? - ( lent2 - part_stop ) : - part_start;
      slinks.push_back( SLinkInfo( rclass, lib, !fw1, !fw2, sep, dev ) );
      continue;
    }
      
    // Case 3: rc-fw.
    sep = ( rclass == 1 ) ? - ijumps_.seps_[lib] : - iJumps_.seps_[lib];
    sep += ( rct1 ) ? - loc.Stop( ) : - ( lent1 - loc.Start( ) );
    sep += ( rct2 ) ? - ( lent2 - part_start ) : - part_stop;
    slinks.push_back( SLinkInfo( rclass, lib, !fw1, !fw2, sep, dev ) );
  }

  // Sort links.
  sort( slinks.begin( ), slinks.end( ) );
  
}

/**
 * CShowLinks
 * ShowLinks
 * private
 */
void CShowLinks::ShowLinks( int c1, int c2, ostream &out ) const
{
  vec<SLinkInfo> slinks;
  this->FindLinks( c1, c2, slinks );

  // Unpack ids, load locs.
  int tig1 = c1 < 0 ? - c1 - 1 : c1;
  int tig2 = c2 < 0 ? - c2 - 1 : c2;
  bool rct1 = ( c1 < 0 );
  bool rct2 = ( c2 < 0 );
  
  // Contig lengths.
  const int lent1 = bases_[tig1].size( );
  const int lent2 = bases_[tig2].size( );

  // Generate table.
  vec< vec<String> > table;
  table. reserve( slinks.size( ) );
  for (int ii=0; ii<slinks.isize( ); ii++)
    table.push_back( slinks[ii].ToLine( ) );

  // Print gap on reference.
  while ( hits_.size( ) > 0 ) {
    out << "\nGap on reference: ";
    int hid1 = to_hit_[tig1];
    int hid2 = to_hit_[tig2];
    if ( hid1 < 0 || hid2 < 0 ) {
      out << "na (not uniquely mapped)\n";
      break;
    }
    const look_align &hit1 = hits_[hid1];
    const look_align &hit2 = hits_[hid2];
    int t1 = hit1.TargetId( );
    int t2 = hit2.TargetId( );
    if ( t1 != t2 ) {
      out << "na (on different targets)\n";
      break;
    }
    if ( hit1.Rc1( ) != rct1 || hit2.Rc1( ) != rct2 ) {
      out << "na (wrong orientation)\n";
      break;
    }
    out << hit2.StartOnTarget( ) - hit1.EndOnTarget( ) << "\n";
    break;
  }

  // Print table.
  out << "\n"
      << table.size( )  << " links between "
      << tig1 << ( rct1 ? "[-]" : "[+]" ) << " (" << lent1 << " bases) and "
      << tig2 << ( rct2 ? "[-]" : "[+]" ) << " (" << lent2 << " bases)\n"
      << endl;
  if ( table.size( ) < 1 )
    out << "No links found\n" << endl;
  else {
    String rjust = "";
    for (int ii=0; ii<table[0].isize( ); ii++)
      rjust += "r";
    PrintTabular( out, table, 3, rjust );
    out << endl;
  }
    
}

