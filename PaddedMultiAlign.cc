// Copyright (c) 2000-2003 Whitehead Institute of Biomedical Research
// 


#include "PaddedMultiAlign.h"
#include <set>

//
// PaddedMultiAlign methods
//


PaddedMultiAlign::PaddedMultiAlign( int contig_id,
				    vec<unsigned int>& contig_pads,
				    vec<padded_read_location>& reads )
  : contig_id_( contig_id ), 
    contig_pads_( contig_pads ),
    reads_( reads )
{
}

ostream& 
operator<< ( ostream& out, const PaddedMultiAlign& pma )
{
  out.write( (char*) &pma.contig_id_, sizeof(int) );

  unsigned int num_pads = pma.contig_pads_.size();
  out.write( (char*) &num_pads, sizeof(unsigned int) );
  for ( unsigned int i = 0; i < num_pads; ++i )
    out.write( (char*) &pma.contig_pads_[i], sizeof(unsigned int) );
  
  out << pma.reads_;

  return out;
}

istream&
operator>> ( istream& in, PaddedMultiAlign& pma )
{
  in.read( (char*) &pma.contig_id_, sizeof(int) );

  unsigned int num_pads;
  in.read( (char*) &num_pads, sizeof(unsigned int) );
  pma.contig_pads_.resize( num_pads );
  for ( unsigned int i = 0; i < num_pads; ++i )
    in.read( (char*) &pma.contig_pads_[i], sizeof(unsigned int) );
  
  in >> pma.reads_;

  return in;
}  

void 
PaddedMultiAlign::setContigId( int contig_id )
{
  contig_id_ = contig_id;
}

int 
PaddedMultiAlign::getContigId() const
{
  return contig_id_;
}

void 
PaddedMultiAlign::setContigPads( const vec<unsigned int>& contig_pads )
{
  contig_pads_ = contig_pads;
}

vec<unsigned int> 
PaddedMultiAlign::getContigPads() const
{
  return contig_pads_;
}

vec<int>
PaddedMultiAlign::getReadIds() const
{
  vec<int> read_ids;
  read_ids.resize( reads_.size() );
  
  transform( reads_.begin(), reads_.end(),
	     read_ids.begin(),
	     mem_fun_ref( &padded_read_location::getReadNumber ) );
  
  sort( read_ids.begin(), read_ids.end() );

  return read_ids;
}
  
void
PaddedMultiAlign::setReads( const vec<padded_read_location>& reads )
{
  reads_ = reads;
}

const vec<padded_read_location>&
PaddedMultiAlign::getReads() const
{
  return reads_;
}
 
void
PaddedMultiAlign::print( ostream &out, 
                         const basevector &contig_bases,
                         const veccompseq &reads_bases,
                         const vec<int>& left_trims,
                         const vec<int>& right_trims ) const
{
    out << "Contig " << contig_id_ << endl;

    vec<unsigned int> contig_pads( contig_pads_ );
    sort( contig_pads.begin(), contig_pads.end() );
    
    vec<int> read_ids( reads_.size(), -1 );
    vec<int> current_base( reads_.size(), 0 );
    vec<String> sequence( reads_.size() );

    int minStart = 0;
    for ( unsigned int read_idx = 0; read_idx < reads_.size(); ++read_idx )
    {
        const padded_read_location &read_loc = reads_[ read_idx ];
        if ( read_loc.getStartOnContig() < minStart )
          minStart = read_loc.getStartOnContig();

        int read_id = read_loc.getReadNumber();
        read_ids[ read_idx ] = read_id;
        const CompressedSequence &read_bases = reads_bases[ read_id ];
        String read_sequence = read_bases.asString();
        int trimmed_length = read_sequence.size() - left_trims[read_id] - right_trims[read_id];

        read_sequence = read_sequence.substr( left_trims[read_id], trimmed_length ); 

        if ( read_loc.getComplementation() )
            StringReverseComplement( read_sequence, read_sequence );

        sequence[ read_idx ] = read_loc.ConstructPaddedSequence( read_sequence );
    }

    vec<int> columns;
    int golden_column = 0;

    vec<unsigned int>::const_iterator pad_iter = contig_pads.begin();
    int contig_base = minStart;
    int padded_contig_base = minStart;

    while ( contig_base < (int) contig_bases.size() )
    {
        for ( unsigned int read_idx = 0; read_idx < reads_.size(); ++read_idx )
            if ( reads_[ read_idx ].getStartOnContig() == padded_contig_base &&
                 current_base[ read_idx ] == 0 )
            {
                unsigned int col_idx = 0;
                for ( ; col_idx < columns.size(); ++col_idx )
                    if ( columns[ col_idx ] == -2 )
                    {
                        columns[ col_idx ] = read_idx;
                        break;
                    }
                if ( col_idx == columns.size() )
                    columns.push_back( read_idx );
            }

        set<int> golden_cols;
        for ( unsigned int col_idx = 0; col_idx < columns.size(); ++col_idx )
        {
            int read_idx = columns[ col_idx ];
            if ( read_idx < 0 )
                continue;
            int read_base = current_base[ read_idx ];
            if ( contig_base >= 0 &&
                 sequence[ read_idx ][ read_base ] == as_base( contig_bases[ contig_base ] ) )
                golden_cols.insert( col_idx );
        }

        golden_column = -1;
        if ( ! golden_cols.count( golden_column ) &&
             ( pad_iter == contig_pads.end() || 
               (int) *pad_iter != contig_base ) )
            if ( ! golden_cols.empty() )
                golden_column = *(golden_cols.begin());

        out << setw( 7 ) << padded_contig_base << " ";

        if ( pad_iter != contig_pads.end() &&
             (int) *pad_iter == contig_base )
        {
            out << '*';
            ++pad_iter;
        }
        else 
        {
          if ( contig_base >= 0 )
            out << as_base( contig_bases[ contig_base ] );
          ++contig_base;
        }
        out << ' ';

        ++padded_contig_base;

        out << setw( 7 ) << ( golden_column < 0 ? -1 : read_ids[ columns[ golden_column ] ] ) << " ";

        for ( unsigned int col_idx = 0; col_idx < columns.size(); ++col_idx )
        {
            int read_idx = columns[ col_idx ];
            if ( read_idx < 0 )
            {
                columns[ col_idx ] = -2;
                continue;
            }

            int read_base = current_base[ read_idx ];
            if ( read_base < (int) sequence[ read_idx ].size() )
            {
                char base = sequence[ read_idx ][ read_base ];
                if ( (int) col_idx != golden_column )
                    base = tolower( base );
                out << base;
                ++current_base[ read_idx ];
            }
            else
            {
                out << ' ';
                columns[ col_idx ] = -1;
            }
        }
        out << endl;
    }
}
    
    


//
// PaddedMultiAlignBuilder methods
//

void 
PaddedMultiAlignBuilder::addPaddedReadLocation( padded_read_location &rc )
{
  realign_candidates_.push_back( rc ); 
}


vec<padded_read_location>* 
PaddedMultiAlignBuilder::getPaddedReadLocations()
{
  return &realign_candidates_; 
}

void PaddedMultiAlignBuilder::getPaddedMultiAlign( PaddedMultiAlign& pma,
						   const int contig_id,
						   const PaddedBasevector& padded_contig )
{
  pma.setContigId( contig_id );
  pma.setContigPads( padded_contig.getPadPositions() );
  pma.setReads( realign_candidates_ );
}


bool
PaddedMultiAlignBuilder::CheckIfNeedRealignment( const padded_read_location &realign, 
						 const vec<unsigned int> &contig_pads )
{

  unsigned int contig_pads_when_aligned = realign.getNumContigPads();
  unsigned int num_contig_pads = contig_pads.size();

  if ( contig_pads_when_aligned < num_contig_pads )
  {

    vec< pair<unsigned int, bool> > contigPads;
    contigPads.reserve( num_contig_pads );

    //  all pads that were in the contig when the read was aligned are 
    //  "false", new pads are "true"
    for ( unsigned int i=0; i<num_contig_pads; ++i )
    {
      if ( i < realign.getNumContigPads() )
	contigPads.push_back( make_pair( contig_pads[i], False ) );
      else
	contigPads.push_back( make_pair( contig_pads[i], True ) );
    }

    //  sort the list
    sort( contigPads.begin(), contigPads.end() );

    //  convert to "new padded" coordinate system
    for ( unsigned int i=0; i<contigPads.size(); ++i )
      contigPads[i].first += i;

    //  read start and end indicies on the contig
    //    unsigned int startpt = realign.getStartOnContig();
    unsigned int endpt = realign.getEndOnContig();

    vec< pair<unsigned int, bool> >::iterator cp_iter = contigPads.begin();
    while ( cp_iter != contigPads.end() )
    {
      //  if a new pad, check to see if it falls in the read
      if ( cp_iter->second )
      {
	unsigned int pad_index = cp_iter->first;	
// 	if ( pad_index < startpt )
// 	{
// 	  ++startpt;
// 	  ++endpt;
// 	}
// 	else if ( pad_index < endpt ) 
	if ( pad_index<endpt )
	{
	  return True;
	}
      }

      ++cp_iter;
    }
    
    return False;
  }
	
  return False;

}

void
PaddedMultiAlignBuilder::addNewPadsToRead( const PaddedBasevector &padded_contig, 
					   padded_read_location &realign )
			
{

  vec<unsigned int> contig_pads = padded_contig.getPadPositions();
  vec< pair<unsigned int, bool> > contigPads;
  contigPads.reserve( contig_pads.size() );

  //  all pads that were in the contig when the read was aligned 
  //  are "false", new pads are "true"
  for ( unsigned int i=0; i<contig_pads.size(); ++i )
  {
    if ( i < realign.getNumContigPads() )
      contigPads.push_back( make_pair( contig_pads[i], False ) );
    else
      contigPads.push_back( make_pair( contig_pads[i], True ) );
  }

  //  sort the list
  sort( contigPads.begin(), contigPads.end() );

  //  convert to "new padded" coordinate system
  for ( unsigned int i=0; i<contigPads.size(); ++i )
    contigPads[i].first += i;

  //  define the vecs to hold the read pads
  vec<unsigned int> read_pads;

   //  read start and end indicies on the contig
  int startpt = realign.getStartOnContig();
  int endpt = realign.getEndOnContig();
  
  int origstartpt = startpt;

  vec< pair<unsigned int, bool> >::iterator cp_iter = contigPads.begin();
  while ( cp_iter != contigPads.end() )
  {
    //  if a new pad, check to see if it falls in the read
    if ( cp_iter->second )
    {
      int pad_index = cp_iter->first;
      if ( pad_index < startpt )
      {
	++startpt;
	++endpt;
      }
      else if ( pad_index < endpt ) 
      {
	//  the pad index in padded coords
	int padplace = pad_index-startpt;
	int preNum(0);

	read_pads = realign.getReadPads();

	if ( !read_pads.empty() )
	{

	  sort( read_pads.begin(), read_pads.end() );

	  // convert to padded coords
	  for ( unsigned int i=0; i<read_pads.size(); ++i )
	    read_pads[i] += i;
	
	  for ( unsigned int i=0; i<read_pads.size(); ++i )
	  {
	    if ( (int) read_pads[i] < padplace )
	      ++preNum;
	    else
	      break;
	  }

	}
	
	//  add the pad to the read
	realign.addPad( padplace-preNum );
	
	//  added a pad, update the end point
	++endpt;

	read_pads.clear();

      }
      else
	//  past the end of the read
	break;

    }
    
    ++cp_iter;
  
  }

  //  reset the start and end indicies of the read on the contig
  realign.setStartOnContig( startpt );
  realign.setEndOnContig( endpt );
  
}


void
PaddedMultiAlignBuilder::UpdatePads( PaddedBasevector &padded_read, 
				     PaddedBasevector &padded_consensus,
				     vec<char> &padded_consensus_sequence,
				     alignment &the_alignment )
{

  //  unpack the alignment
  int start_position_of_read(0), start_position_on_contig(0), errors(0);
  avector<int> gaps(0), lengths(0);

  the_alignment.Unpack(start_position_of_read, 
		       start_position_on_contig, 
		       errors, gaps, lengths);

  //  if no gaps to add, don't do anything
  if (gaps.length>1)
  {

    //  get the pads
    vec<unsigned int> pre_existing_pads_contig = 
      padded_consensus.getTranslatedPadPositions();

    vec<unsigned int> pre_existing_pads_read = 
      padded_read.getTranslatedPadPositions();

    //  initialize counters 
    int pep_contig_counter(0);
    int pep_read_counter(0);

    unsigned int next_consensus_pad( start_position_on_contig );
    unsigned int next_read_pad( start_position_of_read );

    int pads_added_to_consensus(0);

    for ( unsigned int i=0; i< lengths.length; ++i)
    {
      int consensus_pad_counter = 0;
      int read_pad_counter = 0;

      if ( gaps(i) < 0 )     //  put a pad in the consensus
      {
	//  could be more than one pad
	while ( consensus_pad_counter < -gaps(i) )
	{
	  //  check to find out which pads occur before those to be inserted as
	  //  they affect the insertion point.
	  for (unsigned int j=0; j < pre_existing_pads_contig.size(); ++j )
	    if ( pre_existing_pads_contig[j] < next_consensus_pad )
	      ++pep_contig_counter;

	  padded_consensus.addPadPosition( next_consensus_pad-
					   pep_contig_counter );
	  padded_consensus_sequence.insert( padded_consensus_sequence.begin() + 
					    next_consensus_pad + 
					    pads_added_to_consensus,
					    '*' );
	  
	  ++pads_added_to_consensus;
	  ++consensus_pad_counter;
	  ++next_read_pad;
	  pep_contig_counter = 0;   
	}
      }
      else if ( gaps(i) > 0 )  //  pad the read
      {
	//  could be more than one pad
	while ( read_pad_counter < gaps(i) )
	{

	  //  check to find out which pads occur before those to be inserted as
	  //  they affect the insertion point.
	  for (unsigned int j=0; j< pre_existing_pads_read.size(); ++j)
	    if ( pre_existing_pads_read[j] < next_read_pad )
	      ++pep_read_counter;

	  padded_read.addPadPosition( next_read_pad-pep_read_counter );

	  ++read_pad_counter;
	  ++next_consensus_pad;
	  pep_read_counter=0;

	}
      }
       
      next_consensus_pad += lengths(i);
      next_read_pad += lengths(i) ;
     
    }

//     static String contig_name;

//     vec<char> pc_check = padded_consensus.getSequence();
//     if ( pc_check != padded_consensus_sequence && 
// 	 contig_name != padded_consensus.getReadName() )
//     {
//       contig_name = padded_consensus.getReadName();

//       cout << "Pads were not updated properly on " << contig_name << "." << endl;

//       cout << "Consensus object has pads at: " << endl;
//       for ( int i = 0; i < pc_check.size(); ++i )
// 	if ( pc_check[i] == '*' )
// 	  cout << i << endl;

//       cout << "The consensus object's internal pad positions were: " << endl;
//       copy( pre_existing_pads_contig.begin(), pre_existing_pads_contig.end(), 
//  	    ostream_iterator<int>( cout, "\n" ) );

//       vec<unsigned int> pad_positions = padded_consensus.getPadPositions();
//       sort( pad_positions.begin(), pad_positions.end() );

//       cout << "The new internal pad positions are: " << endl;
//       set_difference( pad_positions.begin(), pad_positions.end(),
// 		      pre_existing_pads_contig.begin(), pre_existing_pads_contig.end(),
// 		      ostream_iterator<int>( cout, "\n" ) );
		      
      
//       cout << "Consensus vector has pads at: " << endl;
//       for ( int i = 0; i < padded_consensus_sequence.size(); ++i )
// 	if ( padded_consensus_sequence[i] == '*' )
// 	  cout << i << endl;

//       cout << "There were " << pads_added_to_consensus << " pads added to consensus." << endl;

//       cout << "Here is the start position: " << start_position_on_contig << endl;

//       cout << "Here are the gaps: " << endl;
//       for ( int i = 0; i < gaps.length; ++i )
// 	cout << gaps(i) << endl;

//       cout << "Here are the lengths: " << endl;
//       for ( int i = 0; i < lengths.length; ++i )
// 	cout << lengths(i) << endl;

//       cin.ignore( 100, '\n' );

//     }
    
  }

}

