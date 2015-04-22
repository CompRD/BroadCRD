// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 




#include "PaddedReadLocation.h"


padded_read_location::padded_read_location() 
  : start_on_contig_(0), 
    end_on_contig_(0), 
    contig_pads_when_aligned_(0), 
    read_number_(0), 
    name_("")
{ }

padded_read_location::padded_read_location( const padded_read_location& RC )
  : start_on_contig_(RC.start_on_contig_), 
    end_on_contig_(RC.end_on_contig_), 
    contig_pads_when_aligned_(RC.contig_pads_when_aligned_), 
    read_number_(RC.read_number_), 
    name_(RC.name_), 
    complementation_(RC.complementation_), 
    pads_(RC.pads_)
{ }

padded_read_location& padded_read_location::operator=(const padded_read_location& rhs )
{
  start_on_contig_ = rhs.start_on_contig_;
  end_on_contig_ = rhs.end_on_contig_;
  contig_pads_when_aligned_ = rhs.contig_pads_when_aligned_;
  read_number_ = rhs.read_number_;
  name_ = rhs.name_;
  complementation_ = rhs.complementation_;
  pads_ = rhs.pads_;
  
  return *this;
}

int padded_read_location::getEndOnContig() const
{  return end_on_contig_; } 

unsigned int padded_read_location::getReadNumber() const
{  return read_number_; }

unsigned int padded_read_location::getNumContigPads() const
{  return contig_pads_when_aligned_; }

String padded_read_location::getReadName() const
{  return name_; }

bool padded_read_location::getComplementation() const
{  return complementation_; }

vec<unsigned int> 
padded_read_location::getReadPads() const
{  return pads_; }

void  padded_read_location::setStartOnContig( const int start_on_contig )
{ start_on_contig_ = start_on_contig; }


void  padded_read_location::setEndOnContig( const int end_on_contig )
{  end_on_contig_ = end_on_contig; }

void  padded_read_location::setReadNumber( const unsigned int read_number )
{  read_number_ = read_number; }

void  padded_read_location::setNumContigPads(  const unsigned int contig_pads_when_aligned )
{  contig_pads_when_aligned_ = contig_pads_when_aligned; }

void  padded_read_location::setReadName( const String read_name )
{  name_ = read_name; }

void  padded_read_location::setComplementation( const bool complementation )
{  complementation_ = complementation; }

void  padded_read_location::setReadPads( const vec<unsigned int>& read_pads )
{  pads_ = read_pads; }

void padded_read_location::addPad( unsigned int pad )
{  pads_.push_back( pad ); }

char padded_read_location::returnRC() const
{  return this->getComplementation() ? 'C' : 'U'; }

PaddedBasevector*  
padded_read_location::ConstructPaddedBasevector( const basevector& the_basevector ) const
{
  PaddedBasevector *pbv = new PaddedBasevector( the_basevector, pads_, name_, complementation_ );

  return pbv;
}

String padded_read_location::ConstructPaddedSequence( const String& the_sequence ) const
{

  unsigned int sequence_size = the_sequence.size();
  vec<unsigned int> the_read_pads = getReadPads();
  unsigned int num_pads = the_read_pads.size();

  String padded_sequence;
  padded_sequence.resize( sequence_size+num_pads );


  unsigned int pad_index(0), sequence_index(0), pads_counter(0), basevector_index(0);


  //  if there are no pads, we just have the original basevector
  if ( num_pads == 0 )
  {  
    for (unsigned int i=0; i < sequence_size; ++i)
      padded_sequence[i] = the_sequence[i];
  }
  else
  {

    //  put the bases from basevector into the sequence, adding pads where 
    //  directed by pad_positions
    vec<unsigned int> sorted_pad_positions = getReadPads();
    sort( sorted_pad_positions.begin(), sorted_pad_positions.end() );
 
    while( pads_counter < num_pads )
    {
	pad_index = sorted_pad_positions[ pads_counter ]; 
	if ( pad_index >= sequence_size )
	  break;

      //  add the bases to the sequence up until the pad
      while( basevector_index < pad_index )
      {
	padded_sequence[ sequence_index ] = the_sequence[ basevector_index ];
	
	++sequence_index;
	++basevector_index;
      }

      //  now add the pad to the sequence
      padded_sequence[ sequence_index ] = '*';
      ++sequence_index; 
      ++pads_counter;

    }

    if ( pad_index < sequence_size )
    {
      //  pick up the stragglers after the last pad
      while( basevector_index < sequence_size )
      {  
	padded_sequence[ sequence_index ] = the_sequence[ basevector_index ];
	++sequence_index;  
	++basevector_index; 
      }

    }
  }

  return padded_sequence;

}


bool operator<(const padded_read_location& lhs, const padded_read_location& rhs)
{
  return (lhs.start_on_contig_ < rhs.start_on_contig_);
}


ostream& operator<< (ostream& out, const padded_read_location& rc )
{
  out.write( (char*) &rc.start_on_contig_, sizeof(int) );
  out.write( (char*) &rc.end_on_contig_, sizeof(int) );
  out.write( (char*) &rc.contig_pads_when_aligned_, sizeof(unsigned int) );
  out.write( (char*) &rc.read_number_, sizeof(unsigned int) );
  out.write( (char*) &rc.complementation_, sizeof(bool) );

  unsigned int num_pads = rc.pads_.size();
  out.write( (char*) &num_pads, sizeof(unsigned int) );
  for ( unsigned int i = 0; i < num_pads; ++i )
    out.write( (char*) &rc.pads_[i], sizeof(unsigned int) );

  return out;
}
 
istream& operator>> (istream& in, padded_read_location& rc )
{
  in.read( (char*) &rc.start_on_contig_, sizeof(int) );
  in.read( (char*) &rc.end_on_contig_, sizeof(int) );
  in.read( (char*) &rc.contig_pads_when_aligned_, sizeof(unsigned int) );
  in.read( (char*) &rc.read_number_, sizeof(unsigned int) );
  in.read( (char*) &rc.complementation_, sizeof(bool) );

  unsigned int num_pads;
  in.read( (char*) &num_pads, sizeof(unsigned int) );
  rc.pads_.resize( num_pads );
  for ( unsigned int i = 0; i < num_pads; ++i )
    in.read( (char*) &rc.pads_[i], sizeof(unsigned int) );

  return in;
}  
