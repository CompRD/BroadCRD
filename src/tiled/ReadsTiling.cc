// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#include "Basevector.h"
#include "math/Functions.h"
#include "Qualvector.h"
#include "String.h"
#include "TokenizeString.h"
#include "tiled/ColumnConsensus.h"
#include "tiled/Haploqual.h"
#include "tiled/PaddedSeq.h"
#include "tiled/ReadsTiling.h"



/*
 * reads_tiling
 * Constructor
 */
reads_tiling::reads_tiling( )
{
  this->Clear( );
}



/*
 * reads_tiling
 * SetPointers
 */
void reads_tiling::SetPointers( const basevector *master_bases,
				const vecbasevector *primate_bases,
				const vecqualvector *primate_quals )
{
  master_bases_ = master_bases;
  primate_bases_ = primate_bases;
  primate_quals_ = primate_quals;
}



/*
 * reads_tiling
 * SetFromAligns
 *
 * Set from the given aligns of each read with the master sequence.
 */
void reads_tiling::SetFromAligns( const vec<int> *read_ids,
				  const vec<align> *aligns,
				  const vec<Bool> *isRC,
				  ostream *log )
{
  tile_name_.resize( 0 );
  tile_pad_.resize( 0 );
  
  if ( aligns->size( ) < 1 )
    return;
  
  // Windows begin.
  win_begin_ = (*aligns)[0].pos2( );
  for (int ii=1; ii<(int)aligns->size( ); ii++) {
    int al_begin = (*aligns)[ii].pos2( );
    if ( al_begin < win_begin_ )
      win_begin_ = al_begin;
  }
  
  // Windows end.
  int win_end = 0;
  for (int ii=0; ii<(int)aligns->size( ); ii++) {
    int al_end = (*aligns)[ii].Pos2( );
    if ( al_end > win_end )
      win_end = al_end;
  }

  int extra_space = 0;
  for (int ii=0; ii<(int)aligns->size( ); ii++)
    for (int jj=0; jj<(*aligns)[ii].Nblocks( ); jj++)
      extra_space += (int)fabs( (float) (*aligns)[ii].Gaps( jj ) );
  win_end += extra_space;

  // Tile names.
  tile_name_.resize( read_ids->size( ) );
  for (int ii=0; ii<(int)read_ids->size( ); ii++) {
    tile_name_[ii] = ToString( (*read_ids)[ii] );
    tile_name_[ii] += (*isRC)[ii] ? "_rc" : "   ";
  }
  
  // Initial positioning of reads.
  master_pad_.Set( 0, 0, win_end );
  tile_pad_.resize( read_ids->size( ) );
  for (int ii=0; ii<(int)tile_pad_.size( ); ii++) {
    int begin = (*aligns)[ii].pos2( ) - (*aligns)[ii].pos1( );
    int len = (*aligns)[ii].Pos1( ) - (*aligns)[ii].pos1( );
    tile_pad_[ii].Set( begin, (*aligns)[ii].pos1( ), len );
  }
  
  // Adjust for gaps.
  for (int ii=0; ii<(int)read_ids->size( ); ii++)
    this->Adjust( aligns, ii, log );

}



/*
 * reads_tiling
 * Clear
 */
void reads_tiling::Clear( )
{
  
  master_bases_ = 0;
  primate_bases_ = 0;
  primate_quals_ = 0;

  win_begin_ = 0;
  master_pad_.Clear( );
  tile_name_.resize( 0 );
  tile_pad_.resize( 0 );

  rbases_.resize( 0 );
  rquals_.resize( 0 );

  rectangle_begin_found_ = false;
}



/*
 * reads_tiling
 * Compactify
 *
 * Remove pads-only columns.
 */
void reads_tiling::Compactify( )
{
  // Loop over all master pads.
  for (int ii=0; ii<master_pad_.PadsCount( ); ii++) {
    int pos = master_pad_.Pad( ii );
    
    // Check if this is an all-pads column.
    bool all_pads = true;
    for (int jj=0; jj<(int)tile_pad_.size( ); jj++) {
      int actual_begin = tile_pad_[jj].Begin( ) + tile_pad_[jj].AlBegin( );
      if ( pos < actual_begin )
	continue;
      if ( pos >= actual_begin + tile_pad_[jj].PaddedLength( ) )
	continue;
      if ( !tile_pad_[jj].IsPad( pos ) ) {
	all_pads = false;
	break;
      }
    }
  
    // Remove pads.
    if ( all_pads ) {
      master_pad_.RemovePad( pos );
      for (int jj=0; jj<(int)tile_pad_.size( ); jj++)
	tile_pad_[jj].RemovePad( pos );
    }
  }
}



/*
 * reads_tiling
 * GenerateConsensus
 */
void reads_tiling::GenerateConsensus( bool tag_snp )
{
  this->GenerateRectangle( tag_snp );
}



/*
 * reads_tiling
 * ConsensusBegin
 */
int reads_tiling::ConsensusBegin( ) const
{
  return win_begin_;
}



/*
 * reads_tiling
 * ConsensusLength
 *
 * Only portions of reads aligning master matter for the calculation of
 * the consensus length.
 */
int reads_tiling::ConsensusLength( ) const
{
  int c_len = win_begin_;
  for (int ii=0; ii<(int)tile_pad_.size( ); ii++) {
    int new_end
      = tile_pad_[ii].Begin( )
      + tile_pad_[ii].AlBegin( ) + 
      + tile_pad_[ii].PaddedLength( );
    c_len = ( new_end > c_len ) ? new_end : c_len;
  }
  
  return ( c_len - win_begin_ );
}



/*
 * reads_tiling
 * ReadsCount
 */
int reads_tiling::ReadsCount( ) const
{
  return (int)tile_name_.size( );
}



/*
 * reads_tiling
 * WinBegin
 *
 * The generated rectangle may begin before win_begin_, if reads not aligning
 * the master sequence extend on the left of the consensus.
 */
int reads_tiling::WinBegin( ) const
{
  if ( !rectangle_begin_found_ ) {
    rectangle_begin_ = win_begin_;
    for (int ii=0; ii<(int)tile_pad_.size( ); ii++)
      if ( tile_pad_[ii].Begin( ) < rectangle_begin_ )
	rectangle_begin_ = tile_pad_[ii].Begin( );
    rectangle_begin_found_ = true;
  }
  
  return rectangle_begin_;
}



/*
 * reads_tiling
 * WinEnd
 *
 * Twin to WinBegin
 */
int reads_tiling::WinEnd( ) const
{
  ForceAssert( rbases_.size( ) > 0 );
  return this->WinBegin( ) + rbases_[0].size( );
}



/*
 * reads_tiling
 * BeginOnConsensus
 *
 * This is the beginning of the _read_ on the _extended_ consensus, where:
 * 1) _read_ means that we take the full read into account, and not only the
 * aligned part (in other words it is begin_); and
 * 2) _extended_ means that even though consensus is only built on the
 * aligned portions of reads, in this method we look at the begin with
 * respect to the consensus extended by pasting before and after the real
 * consensus all bits of reads sticking out of it.
 */
int reads_tiling::BeginOnConsensus( int read_index ) const
{
  return + tile_pad_[read_index].Begin( ) - this->WinBegin( );
}



/*
 * reads_tiling
 * ReadsIds
 */
void reads_tiling::ReadsIds( vec<int> &ids ) const
{
  ids.resize( tile_name_.size( ) );
  
  for (int ii=0; ii<(int)tile_name_.size( ); ii++)
    ids[ii]
      = ( tile_name_[ii].Contains( "_" ) )
      ? tile_name_[ii].Before( "_" ).Int( )
      : tile_name_[ii].Int( );
}



/*
 * reads_tiling
 * ReadsFullNames
 */
void reads_tiling::ReadsFullNames( vec<String> &names ) const
{
  names = tile_name_;
}



/*
 * reads_tiling
 * MasterToFasta
 */
void reads_tiling::MasterToFasta( ostream &out, int width ) const
{
  // Check printable rectangles have been generated.
  ForceAssert( rbases_.size( ) > 0 );
  
  // Send to out.
  int col = 1;
  for (int ii=0; ii<(int)rbases_[0].size( ); ii++) {
    if ( col > width ) {
      out << "\n";
      col = 1;
    }    
    out << rbases_[0][ii];
    col++;
  }
}



/*
 * reads_tiling
 * ConsensusBasesToFasta
 */
void reads_tiling::ConsensusBasesToFasta( ostream &out, int width ) const
{
  // Check printable rectangles have been generated.
  ForceAssert( rbases_.size( ) > 0 );
  
  // Send to out.
  int col = 1;
  for (int ii=0; ii<(int)rbases_[1].size( ); ii++) {
    if ( col > width ) {
      out << "\n";
      col = 1;
    }    
    out << rbases_[1][ii];
    col++;
  }
}



/*
 * reads_tiling
 * ConsensusQualsToQual
 */
void reads_tiling::ConsensusQualsToQual( ostream &out, int width ) const
{
  // Check printable rectangles have been generated.
  ForceAssert( rbases_.size( ) > 0 );
  
  // Send to out.
  int col = 1;
  for (int ii=0; ii<(int)rquals_[1].size( ); ii++) {
    if ( col > width ) {
      out << "\n";
      col = 1;
    }    
    out << rquals_[1][ii] << " ";
    col++;
  }
}



/*
 * reads_tiling
 * MasterSeq
 */
const padded_seq &reads_tiling::MasterSeq( ) const
{
  return master_pad_;
}



/*
 * reads_tiling
 * PaddedSeq
 */
const padded_seq &reads_tiling::PaddedSeq( int ii ) const
{
  return tile_pad_[ii];
}



/*
 * reads_tiling
 * SlimPrint
 *
 * Similar to FoldPrint, but it only prints the essentials to save
 * hard disk space.
 *
 * width: folding width.
 */
void reads_tiling::SlimPrint( ostream &out, int width ) const
{
  // Check printable rectangles have been generated.
  ForceAssert( rbases_.size( ) > 0 );
  
  // Constants.
  vec<String> read_name;
  this->GenerateNames( read_name );
  
   // Print.
  int w_begin = this->WinBegin( );
  int w_end = this->WinEnd( );

  int chunk = 0;
  int pos = w_begin;
  while ( pos < w_end ) {
    int end_pos = ( width + pos < w_end ) ? width + pos : w_end;
    
    out << "  chunk " << chunk << ": [" << pos << ", " << end_pos - 1 << ")\n";
    
    // Bases.
    for (int ii=0; ii<2; ii++) {
      out << read_name[ii];
      for (int jj=pos; jj<end_pos; jj++)
	out << rbases_[ii][jj-w_begin];
      out << "\n";
    }
    
    // Quals.
    for (int ii=1; ii<2; ii++) {
      out << read_name[ii];
      for (int jj=pos; jj<end_pos; jj++)
	out << this->PrettyQual( rquals_[ii][jj-w_begin] );
      out << "\n";
    }
    
    // Slide on.
    pos += width;
    chunk++;
  }
  
}



/*
 * reads_tiling
 * FoldPrint
 *
 * quals: if false it will not print quality scores;
 * width: folding width;
 * from: start of subinterval to print;
 * to: end of subinterval to print.
 */
void reads_tiling::FoldPrint( ostream &out,
			      bool quals,
			      int width,
			      int from,
			      int to ) const
{
  // Check printable rectangles have been generated.
  ForceAssert( rbases_.size( ) > 0 );

  // Constants.
  vec<String> read_name;
  this->GenerateNames( read_name );
  int tot_width = (int)rbases_[0].size( );
  int win_begin = this->WinBegin( );
  int win_end = this->WinEnd( );
  
  // Determine beginning and ending for all reads.
  vec< pair<int, int> > reads_win( read_name.size( ), make_pair( 0, 0 ) );
  for (int ii=2; ii<(int)read_name.size( ); ii++) {
    const padded_seq &pad_seq = tile_pad_[ii-2];
    bool found_beginning = false;
    bool found_ending = false;

    // If a guess for the beginning does not work, search all the way through.
    for (int iteration=0; iteration<2; iteration++) {
      int guess_begin = Max( win_begin, pad_seq.Begin( ) - 600 );
      if ( iteration == 1 )
	guess_begin = win_begin;
      for (int jj=guess_begin; jj<tot_width+win_begin; jj++) {
	if ( !found_beginning ) {
	  if ( rbases_[ii][jj-win_begin] != empty_base ) {
	    reads_win[ii].first = jj;
	    found_beginning = true;
	    continue;
	  }
	}
	else {
	  if ( rbases_[ii][jj-win_begin] == empty_base ) {
	    reads_win[ii].second = jj;
	    found_ending = true;
	    break;
	  }
	}
      }
      if ( found_beginning )
	break;
      else
	ForceAssert( 0 == iteration );
    }
    
    if ( !found_ending )
      reads_win[ii].second = tot_width;
  }

  // Print.
  int w_begin = ( from == -1 ) ? win_begin : from;
  int w_end = ( to == -1 ) ?  win_end : to;
  
  int chunk = 0;
  int pos = w_begin;
  while ( pos < w_end ) {
    int end_pos = ( width + pos < w_end ) ? width + pos : w_end;
    
    out << "  chunk " << chunk << ": [" << pos << ", " << end_pos - 1 << ")\n";
    
    // Bases.
    vec<int> print_it( read_name.size( ), 0 );
    for (int ii=0; ii<(int)read_name.size( ); ii++) {
      if ( ii < 2 )
	print_it[ii] = 1;
      else {
	int left = Max( reads_win[ii].first, pos );
	int right = Min( reads_win[ii].second, end_pos );
	if ( left < right )
	  print_it[ii] = 1;
      }
      
      if ( 1 == print_it[ii] ) {
	out << read_name[ii];
	for (int jj=pos; jj<end_pos; jj++)
	  out << rbases_[ii][jj-win_begin];
	out << "\n";
      }
    }

    // Quals.
    if ( quals ) {
      for (int ii=0; ii<(int)read_name.size( ); ii++) {
	if ( 1 == print_it[ii] ) {
	  out << read_name[ii];
	  for (int jj=pos; jj<end_pos; jj++)
	    out << this->PrettyQual( rquals_[ii][jj-win_begin] );
	  out << "\n";
	}
      }
    }

    // Slide on.
    pos += width;
    chunk++;
  }
  
}



/*
 * reads_tiling
 * operator>>
 */
istream &operator>> ( istream &in, reads_tiling &tiles )
{
  // Info line.
  int n_reads;
  in >> tiles.win_begin_
     >> n_reads;

  if ( !in )
    return in;

  // Master.
  in >> tiles.master_pad_;
  
  // Tiles.
  tiles.tile_name_.resize( n_reads );
  tiles.tile_pad_.resize( n_reads );
  for (int ii=0; ii<n_reads; ii++)
    in >> tiles.tile_name_[ii]
       >> tiles.tile_pad_[ii];
  
  return in;
}



/*
 * reads_tiling
 * operator<<
 */
ostream &operator<< ( ostream &out, const reads_tiling &tiles )
{
  // Info line.
  out << tiles.win_begin_ << "\t"
      << tiles.tile_name_.size( ) << "\n";
  
  // Master.
  out << tiles.master_pad_ << "\n";
  
  // Tiles.
  for (int ii=0; ii<(int)tiles.tile_name_.size( ); ii++)
    out << tiles.tile_name_[ii] << "\t"
	<< tiles.tile_pad_[ii] << "\n";

  return out;
}



/*
 * reads_tiling
 * Adjust
 *
 * If you pass a nonzero log, then runtime will go up and a truly massive
 * amount of log will be produced (in other words, use only if debugging).
 */
void reads_tiling::Adjust( const vec<align> *aligns,
			   int read_index,
			   ostream *log )
{
  const align &the_al = (*aligns)[ read_index ];
  int n_blocks = (int)the_al.Nblocks( );

  if ( log ) {
    this->GenerateRectangle( );
    this->FoldPrint( *log, false );
  }

  // Go at the beginning of alignment.
  int beg = tile_pad_[read_index].Begin( ) + (*aligns)[read_index].pos1( );
  int pos = beg;

  // Walk along the alignment.
  for (int ii=0; ii<n_blocks; ii++) {

    // Gap on clone.
    if ( the_al.Gaps( ii ) < 0 ) {
      for (int jj=0; jj<-(int)the_al.Gaps(ii); jj++)
	this->PadClone( pos, read_index );
      
      if ( log ) {
	*log << tile_name_[read_index] << ": "
	     << the_al.Gaps(ii) << " gaps on clone\n";
	this->GenerateRectangle( );
	this->FoldPrint( *log, false );
      }
    }

    // Gap on read.
    if ( the_al.Gaps( ii ) > 0 ) {
      for (int jj=0; jj<(int)the_al.Gaps(ii); jj++)
	this->PadRead( pos, read_index );

      if ( log ) {
	*log << tile_name_[read_index] << ": "
	     << the_al.Gaps(ii) << " gaps on read\n";
	this->GenerateRectangle( );
	this->FoldPrint( *log, false );
      }
    }

    // Length.
    int len = the_al.Lengths( ii );
    
    if ( the_al.Gaps( ii ) < 0 )
      len += -the_al.Gaps( ii );

    while ( len > 0 ) {
      if ( pos >= beg + tile_pad_[read_index].PaddedLength( ) )
	break;
      if ( !tile_pad_[read_index].IsPad( pos ) )
	len--;
      pos++;
    }
    
    if ( log ) {
      *log << tile_name_[read_index] << ": "
	   << the_al.Lengths(ii) << " len\n";
      this->GenerateRectangle( );
      this->FoldPrint( *log, false );
    }
  }
  
}



/*
 * reads_tiling
 * PadClone
 */
void reads_tiling::PadClone( int &pos, int read_index )
{
  if ( tile_pad_[read_index].IsPad( pos ) ) {
    // read_index has bases extending into a master-padded region.
    int loc_pos = pos;
    while ( tile_pad_[read_index].IsPad( loc_pos ) )
      loc_pos++;
    loc_pos--;
    tile_pad_[read_index].RemovePad( loc_pos );
    tile_pad_[read_index].ShiftPadsAfter( loc_pos, +1 );

    // If last base is pad, remove it.
    int last_pos
      = tile_pad_[read_index].Begin( )
      + tile_pad_[read_index].AlBegin( )
      + tile_pad_[read_index].PaddedLength( )
      - 1;
    while ( tile_pad_[read_index].IsPad( last_pos ) ) {
      tile_pad_[read_index].RemovePad( last_pos );
      last_pos--;
    }
  }
  else {
    // Padding master.
    master_pad_.AddPad( pos );
    for (int jj=0; jj<(int)tile_pad_.size( ); jj++) {
      if ( jj != read_index ) {
	// Reads other than read_index must be padded too.
	tile_pad_[jj].AddPad( pos );

	// Adding a pad may push the read into a master padded region.
	int loc_pos
	  = tile_pad_[jj].Begin( )
	  + tile_pad_[jj].AlBegin( )
	  + tile_pad_[jj].PaddedLength( )
	  - 1;
	while ( master_pad_.IsPad( loc_pos ) ) {
	  tile_pad_[jj].AddPad( loc_pos );
	  loc_pos++;
	}
      }
      else {
	// Push pads on read_index mirroring pads on master after pos.
	tile_pad_[jj].ShiftPadsAfter( pos, +1 );
      }
    }
  }
}



/*
 * reads_tiling
 * PadRead
 */
void reads_tiling::PadRead( int &pos, int read_index )
{
  // Move to the first non-padded master base.
  while ( master_pad_.IsPad( pos ) )
    pos++;
  
  // Add the pad.
  tile_pad_[read_index].AddPad( pos );

  // Shift all pads after pos.
  tile_pad_[read_index].ShiftPadsAfter( pos, -1 );
  
  // Adding a pad may push the read into a master padded region.
  int loc_pos
    = tile_pad_[read_index].Begin( )
    + tile_pad_[read_index].AlBegin( )
    + tile_pad_[read_index].PaddedLength( )
    - 1;
  while ( master_pad_.IsPad( loc_pos ) ) {
    tile_pad_[read_index].AddPad( loc_pos );
    loc_pos++;
  }

  // Move on.
  pos++;
}



/*
 * reads_tiling
 * ColumConsensus
 */
void reads_tiling::ColumnConsensus( int column, bool tag_snp )
{
  // Check printable rectangles have been generated.
  ForceAssert( rbases_.size( ) > 0 );

  // Instantiate and fill a column consensus.
  column_consensus consensus;

  for (int ii=2; ii<(int)rbases_.size( ); ii++) {
    char base = rbases_[ii][column];
    if ( empty_base == base ) {
      // Nothing to do.
      continue;
    }
    else if ( gap_base == base ) {
      // A pad.
      int left_qual;
      {
	int loc_pos = column - 1;
	while ( loc_pos > -1 &&
		gap_base == rbases_[ii][loc_pos] )
	  loc_pos--;
	if ( empty_base == loc_pos || loc_pos < 0 ) {
	  // We should never get in here.
	  left_qual = 0;
	}
	left_qual = rquals_[ii][loc_pos];
      }
      int right_qual;
      {
	int loc_pos = column + 1;
	while ( loc_pos < (int)rbases_[ii].size( ) &&
		gap_base == rbases_[ii][loc_pos] )
	  loc_pos++;
	if ( empty_base == loc_pos || loc_pos >= (int)rbases_[ii].size( ) ) {
	  // We should never get in here.
	  right_qual = 0;
	}
	right_qual = rquals_[ii][loc_pos];
      }
      int qual = ( left_qual < right_qual ) ? left_qual : right_qual;
      consensus.Add( base, qual );
      continue;
    }
    else {
      // A non-pad base.
      int qual = rquals_[ii][column];
      consensus.Add( base, qual );
      continue;
    }
  }

  // Find consensus.
  char base;
  int qual;
  consensus.Consensus( base, qual, tag_snp );
  if ( base == gap_base && qual > -1 )
    qual = gap_qual;

  rbases_[1][column] = base;
  rquals_[1][column] = qual;
  return;

}



/*
 * reads_tiling
 * GenerateRectangle
 *
 * Although consensus is built only upon aligned portions of reads, it
 * saves into the rectangle all the bases of the reads. The bases of
 * reads which are not aligned to master are saved in lower case.
 */
void reads_tiling::GenerateRectangle( bool tag_snp )
{
  ForceAssert( master_bases_ && primate_bases_ && primate_quals_ );
  
  // Prepare original fastb's amd qualb's.
  const basevector* mbases = master_bases_;

  vec< basevector > bases_rc( tile_name_.size( ) );
  vec< qualvector > quals_rc( tile_name_.size( ) );
  vec< const basevector* > bases( tile_name_.size( ) );
  vec< const qualvector* > quals( tile_name_.size( ) );
  for (int ii=0; ii<(int)tile_name_.size( ); ii++) {
    int read_id;
    bool isRC;
    if ( tile_name_[ii].Contains( "_" ) ) {
      read_id = tile_name_[ii].Before( "_" ).Int( );
      isRC = true;
    }
    else {
      if ( tile_name_[ii].Contains( " " ) )
	read_id = tile_name_[ii].Before( " " ).Int( );
      else
	read_id = tile_name_[ii].Int( );
      isRC = false;
    }
    
    if ( isRC ) {
      bases_rc[ii] = (*primate_bases_)[read_id];
      quals_rc[ii] = (*primate_quals_)[read_id];
      bases_rc[ii].ReverseComplement( );
      quals_rc[ii].ReverseMe( );
    }
    
    bases[ii] = isRC ? &bases_rc[ii] : &( (*primate_bases_)[read_id] );
    quals[ii] = isRC ? &quals_rc[ii] : &( (*primate_quals_)[read_id] );
  }

  // Resize rectangle.
  int the_begin = win_begin_;
  int the_end = win_begin_;
  for (int ii=0; ii<(int)tile_pad_.size( ); ii++) {
    int rd_begin = tile_pad_[ii].Begin( );
    int rd_end = rd_begin + bases[ii]->size( ) + tile_pad_[ii].PadsCount( );
    the_begin = ( rd_begin < the_begin ) ? rd_begin : the_begin;
    the_end = ( rd_end > the_end ) ? rd_end : the_end;
  }
  int height = (int)tile_name_.size( ) + 2;
  int width = the_end - the_begin;
  
  rbases_.resize( height );
  rquals_.resize( height );
  
  // Master sequence.
  const basevector *mfb = mbases;
  const qualvector *mqb = 0;
  vec<char> &mvecb = rbases_[0];
  vec<int> &mvecq = rquals_[0];
  master_pad_.ResetUnpaddedLength( the_begin+width-master_pad_.PadsCount( ) );
  master_pad_.ToVec( mfb, mqb, mvecb, mvecq, the_begin, the_end );

  rbases_[0].resize( width, empty_base );
  rquals_[0].resize( width, empty_qual );
  
  // Consensus sequence.
  rbases_[1].resize( width, empty_base );
  rquals_[1].resize( width, empty_qual );

  // Tile sequences.
  for (int ii=0; ii<(int)tile_name_.size( ); ii++) {
    const basevector *fb = bases[ii];
    const qualvector *qb = quals[ii];
    vec<char> &vecb = rbases_[ii+2];
    vec<int> &vecq = rquals_[ii+2];

    tile_pad_[ii].ToVec( fb, qb, vecb, vecq, the_begin, the_end );
    
    rbases_[ii+2].resize( width, empty_base );
    rquals_[ii+2].resize( width, empty_qual );
  }
  
  // Build consensus.
  for (int ii=0; ii<(int)rbases_[0].size( ); ii++)
    this->ColumnConsensus( ii, tag_snp );

}



/*
 * reads_tiling
 * GenerateNames
 */
void reads_tiling::GenerateNames( vec<String> &read_name ) const
{
  int max_len_name = 0;
  int n_lines = 2 + (int)tile_name_.size( );
  read_name.resize( n_lines );
  for (int ii=0; ii<n_lines; ii++) {
    if ( 0 == ii )
      read_name[ii] = "clone";
    else if ( 1 == ii )
      read_name[ii] = "consensus";
    else
      read_name[ii] = tile_name_[ii-2];
    
    if ( (int)read_name[ii].size( ) > max_len_name )
      max_len_name = read_name[ii].size( );
  }
  for (int ii=0; ii<n_lines; ii++) {
    while ( (int)read_name[ii].size( ) < 4 + max_len_name )
      read_name[ii] += " ";
  }
}



/*
 * reads_tiling
 * PrettyQual
 */
char reads_tiling::PrettyQual( int n_qual ) const
{
  if ( snp_qual == n_qual )
    return '+';
  
  if ( untrusted_qual == n_qual )
    return '-';

  if ( gap_qual == n_qual )
    return gap_base;

  if ( empty_qual == n_qual )
    return empty_base;
  
  if ( n_qual < 10 )
    return '0';
  
  return ToString( n_qual )[0];
}



/*
 * LoadReadsTiling
 */
void  LoadReadsTilings( const String &file_name,
			vec<reads_tiling> &tilings )
{
  tilings.clear( );

  reads_tiling tiled_contig;
  ifstream in( file_name.c_str( ) );
  while ( 1 ) {
    in >> tiled_contig;
    if ( !in ) break;
    tilings.push_back( tiled_contig );
  }
}


