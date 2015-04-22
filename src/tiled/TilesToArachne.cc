// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

#include "AnnotatedContig.h"
#include "Basevector.h"
#include "FastIfstream.h"
#include "Qualvector.h"
#include "ReadLocation.h"
#include "String.h"
#include "STLExtensions.h"
#include "system/System.h"
#include "tiled/Build33.h"
#include "tiled/ChrMapper.h"
#include "tiled/PaddedSeq.h"
#include "tiled/ReadsTiling.h"
#include "tiled/TilesToArachne.h"
#include "tiled/TilesContigParser.h"



/*
 * tiles_to_arachne
 * Constructor
 */
tiles_to_arachne::tiles_to_arachne ( const String &tiles_dir,
				     const String &consensus_dir ) :
  tiles_dir_ ( tiles_dir ),
  consensus_dir_ ( consensus_dir ),
  n_reads_ ( -1 ),
  n_bases_ ( -1 )
{
  ForceAssert( IsDirectory( tiles_dir_ ) );
  ForceAssert( IsDirectory( consensus_dir_ ) );
}



/*
 * tiles_to_arachne
 * SaveChrMap
 */
void tiles_to_arachne::SaveChrMap( const String &map_file )
{
  if ( chr_map_.size( ) < 1 )
    this->GenerateChrMap( );
  
  WRITE( map_file, chr_map_ );
}



/*
 * tiles_to_arachne
 * SaveFastbQualb
 */
void tiles_to_arachne::SaveFastbQualb( const String &bases_file,
				       const String &quals_file,
				       const String &pads_file )
{
  // Bases and quals files.
  int n_files = 0;
  vector<String> all_files = AllFiles( consensus_dir_ );
  for (int ii=0; ii<(int)all_files.size( ); ii++)
    if ( all_files[ii].Contains( "chr" ) && 
	 all_files[ii].Contains( "bases" ) )
      n_files++;
  
  // Create a tc_parser.
  tc_parser tile_parser;
  
  if ( n_reads_ < 0 )
    this->CountReadsAndBases( );
  tile_parser.Reserve( n_reads_, n_bases_ );

  // Append bases and quals
  cout << Date( ) << ": loading fasta and quals ("
       << n_files << " pairs)" << endl;
  
  for (int ii=0; ii<n_files; ii++) {
    Dot( cout, ii );
    String sz_chr = ToStandard( ii );
    String in_bases = consensus_dir_ + "/chr" + sz_chr + ".bases";
    String in_quals = consensus_dir_ + "/chr" + sz_chr + ".quals";
    tile_parser.Append( in_bases, in_quals );
  }
  cout << endl;
  
  // Save to disk.
  tile_parser.Write( bases_file, quals_file, pads_file );
  
}



/*
 * tiles_to_arachne
 * SaveLocs
 */
void tiles_to_arachne::SaveLocs( const String &contigs_file,
				 const String &reads_file,
				 const String &pads_file,
				 const String &locs_file )
{
  if ( chr_map_.size( ) < 1 )
    this->GenerateChrMap( );
  
  // Load.
  cout << Date( ) << ": loading read lengths" << endl;
  vec<int> read_len;
  {
    vecbasevector reads;
    reads.ReadAll( reads_file );
    read_len.resize( reads.size( ) );
    for (int ii=0; ii<(int)reads.size( ); ii++)
      read_len[ii] = reads[ii].size( );
  }
  
  cout << Date( ) << ": loading contig lengths" << endl;
  vec<int> cg_len;
  {
    vecbasevector contigs;
    contigs.ReadAll( contigs_file );
    cg_len.resize( contigs.size( ) );
    for (int ii=0; ii<(int)contigs.size( ); ii++)
      cg_len[ii] = contigs[ii].size( );
  }

  cout << Date( ) << ": loading pads informations" << endl;
  READ( pads_file, vec<padded_seq>, pads );

  cout << Date( ) << ": searching for tiled contigs files" << endl;
  vec<String> tiles_files;
  {
    vector<String> all_files = AllFiles( tiles_dir_ );
    for (int ii=0; ii<(int)all_files.size( ); ii++) {
      if ( all_files[ii].Contains( "chr." ) ) {
	if ( all_files[ii].Contains( ".contigs" ) ) {
	  tiles_files.push_back( all_files[ii] );
	}
      }
    }
  }

  // Fill locs.
  cout << Date( ) << ": loading tiled contigs ("
       << tiles_files.size( ) << " files)" << endl;
  
  vec<read_location> locs;
  locs.reserve( read_len.size( ) );

  for (int ii=0; ii<(int)tiles_files.size( ); ii++) {
    Dot( cout, ii );
    
    int chr_id = tiles_files[ii].After( "chr." ).Before( ".contigs" ).Int( );
    const String in_file = tiles_dir_ + "/" + tiles_files[ii];

    int pos_on_chr = 0;
    reads_tiling tiled_cg;
    vec<String> full_names;
    vec<int> begins_on_contig;
    ifstream in_tiles( in_file.c_str( ) );

    while ( 1 ) {
      in_tiles >> tiled_cg;
      if ( !in_tiles )
	break;

      chr_mapper contigizer( chr_id, -1, pos_on_chr, -1 );
      vec<chr_mapper>::const_iterator iter
	= lower_bound( chr_map_.begin( ), chr_map_.end( ), contigizer );
      if ( iter->Chr( ) != chr_id || iter->Pos( ) != pos_on_chr ) {
	ForceAssert( iter->Chr( ) == chr_id );
	ForceAssert( iter->Pos( ) == pos_on_chr );
      }
      int contig_id = iter->Id( );
      const padded_seq &contig_pads = pads[contig_id];
      full_names.clear( );
      begins_on_contig.clear( );
      tiled_cg.ReadsFullNames( full_names );

      for (int jj=0; jj<(int)full_names.size( ); jj++) {
	int read_id;
	orientation orient;
	if ( full_names[jj].Contains( "_rc" ) ) {
	  read_id = full_names[jj].Before( "_rc" ).Int( );
	  orient = ReverseOr;
	}
	else {
	  read_id = full_names[jj].Int( );
	  orient = ForwardOr;
	}
	int read_length = read_len[read_id];
	int contig_length = cg_len[contig_id];
 	int start_on_contig = tiled_cg.BeginOnConsensus( jj );
	int padded_end = start_on_contig;
	int index = 0;
	while ( index< contig_pads.PadsCount( ) &&
		contig_pads.Pad(index) < padded_end ) {
	  start_on_contig--;
	  index++;
	}
	
	read_location newloc( read_id,
			      read_length,
			      contig_id,
			      start_on_contig,
			      orient,
			      contig_length );
	
	locs.push_back( newloc );
      }
      pos_on_chr++;
    }
  }
  cout << endl;

  // Save.
  WriteLocs( locs_file, locs, cg_len.size(), read_len.size() );
  
}



/*
 * tiles_to_arachne
 * SaveSuperInfo
 */
void tiles_to_arachne::SaveSuperInfo( const String &asupers_file )
{
  if ( chr_map_.size( ) < 1 )
    this->GenerateChrMap( );
  
  vec<annotated_supercontig> asupers;
  asupers.reserve( chr_map_.size( ) );

  basevector bv( 0 );
  qualvector qv( 0 );
  vec<semiannotation> semiann( 0 );
  vec<annotation> ann( 0 );
  vec<arachne_contig_placement> pls;
  
  for (int ii=0; ii<(int)chr_map_.size( ); ii++) {
    vec<annotated_contig> cgs;
    vec<int> gaps;
    vec<int> gaps_sd;
    int cg_id = chr_map_[ii].Id( );
    annotated_contig acg( cg_id, bv, qv, False, False, semiann, ann, pls );
    cgs.push_back( acg );
    annotated_supercontig as( cgs, gaps, gaps_sd );
    asupers.push_back( as );
  }

  PipeOstream( out, asupers_file );
  out << asupers.size( ) << "\n";
  for (int ii=0; ii<(int)asupers.size( ); ii++) 
    out << asupers[ii] << "\n";    
}



/*
 * tiles_to_arachne
 * GenerateChrMap
 */
void tiles_to_arachne::GenerateChrMap( )
{
  // Load ancillary files.
  vec<String> ancillary;
  {
    vector<String> all_files = AllFiles( consensus_dir_ );
    for (int ii=0; ii<(int)all_files.size( ); ii++) {
      if ( all_files[ii].Contains( "chr" ) &&
	   all_files[ii].Contains( "ancillary" ) )
	ancillary.push_back( all_files[ii] );
    }
  }

  // Count contigs.
  int n_contigs = 0;
  for (int ii=0; ii<(int)ancillary.size( ); ii++)
    n_contigs += LineCount(  consensus_dir_ + "/" + ancillary[ii] );
  
  // Fill chr_map_.
  chr_map_.resize( n_contigs );
  
  int loc_counter = 0;
  for (int ii=0; ii<(int)ancillary.size( ); ii++) {
    String in_file = consensus_dir_ + "/" + ancillary[ii];
    String st_chr = ancillary[ii].After( "chr" ).Before( ".ancillary" );
    int chr_id = ToBuild33( st_chr );
    
    ifstream in( in_file.c_str( ) );
    
    while ( 1 ) {
      int pos_in_chr;
      int n_reads;
      int len;
      int begin;
      int end;
      
      in >> pos_in_chr >> n_reads >> len >> begin >> end;
      if ( !in )
	break;
      
      chr_map_[loc_counter].Set( chr_id, loc_counter, pos_in_chr, begin );
      loc_counter++;
    }

    in.close( );
  }

  // Sort.
  if ( !is_sorted( chr_map_.begin( ), chr_map_.end( ) ) )
    sort( chr_map_.begin( ), chr_map_.end( ) );

  // Renumber contig ids.
  for (int ii=0; ii<(int)chr_map_.size( ); ii++)
    chr_map_[ii].SetId( ii );
}



/*
 * tiles_to_arachne
 * CountReadsAndBases
 */
void tiles_to_arachne::CountReadsAndBases( )
{
  vec<String> bases_files;

  vector<String> all_files = AllFiles( consensus_dir_ );
  for (int ii=0; ii<(int)all_files.size( ); ii++) {
    if ( all_files[ii].Contains( "chr" ) &&
	 all_files[ii].Contains( ".bases" ) )
      bases_files.push_back( all_files[ii] );
  }

  n_reads_ = 0;
  n_bases_ = 0;

  for (int ii=0; ii<(int)bases_files.size( ); ii++) {
    fast_ifstream in( consensus_dir_ + "/" + bases_files[ii] );
    ForceAssert ( !in.fail( ) );
    String line, gt( ">" );
    while( 1 ) {
      getline( in, line );
      if ( in.fail( ) ) break;
      if ( line.Contains( gt, 0 ) ) ++n_reads_;
      else n_bases_ += line.size( );
    }
  }
}



