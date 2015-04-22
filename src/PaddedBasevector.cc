// Copyright (c) 2000, 2001 Whitehead Institute of Biomedical Research
// 
// Class encpsulating the idea of a basevector with pads, where
// a pad represents a gap introduced by the alignment of the basevector
// to another basevector (i.e. read vs. read, read vs. consensus)
//
// Used in the creatation of ace files from Arachne assemblies.
//


#include "PaddedBasevector.h"


PaddedBasevector::PaddedBasevector()
{
  read_name_ = "";
  number_pads_=0;
  complementation_ = False;  
}

PaddedBasevector::PaddedBasevector( const PaddedBasevector &padded_sequence ) 
  : the_basevector_( padded_sequence.the_basevector_ ), 
    number_pads_( padded_sequence.number_pads_ ), 
    read_name_(padded_sequence.read_name_), 
    complementation_(padded_sequence.complementation_), 
    pad_positions_(padded_sequence.pad_positions_)
{ }

PaddedBasevector::PaddedBasevector( const basevector &theBasevector, 
				    const vec<unsigned int> &padPositions, 
				    String readName, 
				    bool Complementation )
{
 
  setBasevector( theBasevector );
  setPadPositions( padPositions );
  setReadName( readName );
  setComplementation( Complementation );
  
}

PaddedBasevector &PaddedBasevector::operator=(const PaddedBasevector &rhs_sequence)
{

  the_basevector_  = rhs_sequence.the_basevector_;
  number_pads_ = rhs_sequence.number_pads_;
  read_name_   = rhs_sequence.read_name_;
  complementation_ = rhs_sequence.complementation_;
  pad_positions_ = rhs_sequence.pad_positions_;
  translated_positions_ = rhs_sequence.translated_positions_;
 
  return *this;
}

PaddedBasevector::~PaddedBasevector()
{ }
 

unsigned int PaddedBasevector::PaddedLength() const
{ return ( number_pads_ + the_basevector_.size() ); }

unsigned int PaddedBasevector::NumberOfPads() const
{  return number_pads_; }


vec<unsigned int> 
PaddedBasevector::getPadPositions() const
{ return pad_positions_; }


void
PaddedBasevector::updateTranslatedPadPositions_() const
{ 
  if ( translated_positions_.size() != pad_positions_.size() )
  {
    translated_positions_ = pad_positions_;
    sort( translated_positions_.begin(), translated_positions_.end() );

    for ( unsigned int i = 0; i < translated_positions_.size(); ++i )
      translated_positions_[i] += i;
  }
}

vec<unsigned int> 
PaddedBasevector::getTranslatedPadPositions() const
{ 
  this->updateTranslatedPadPositions_();

  return translated_positions_;
}

int PaddedBasevector::getPadsInRange( int padded_begin, int padded_end ) const
{
  this->updateTranslatedPadPositions_();

  int pads_in_range = 0;

  for ( unsigned int i = 0; i < translated_positions_.size(); ++i )
    if ( (int) translated_positions_[i] >= padded_begin &&
	 (int) translated_positions_[i] < padded_end )
      ++pads_in_range;

  return pads_in_range;
}
      

int PaddedBasevector::getUnpaddedPosition( int padded_position )
{
  int unpadded_position = padded_position;

  this->updateTranslatedPadPositions_();

  for ( unsigned int i = 0; i < translated_positions_.size(); ++i )
    if ( (int) translated_positions_[i] < padded_position )
      --unpadded_position;

  return unpadded_position;
}

int PaddedBasevector::getPaddedPosition( int unpadded_position )
{
  int padded_position = unpadded_position;
  for ( unsigned int i = 0; i < pad_positions_.size(); ++i )
    if ( (int) pad_positions_[i] < unpadded_position )
      ++padded_position;

  return padded_position;
}

  
    

const basevector PaddedBasevector::getBasevector() const
{ return the_basevector_; }

const String PaddedBasevector::getReadName() const
{ return read_name_; }

bool PaddedBasevector::getComplementation() const
{  return complementation_; }

char PaddedBasevector::returnRC() const
{  return this->getComplementation() ? 'C' : 'U'; }

void PaddedBasevector::setPadPositions( const vec<unsigned int> &padPositions )
{
  
  number_pads_ = padPositions.size();
  pad_positions_ = padPositions;
  translated_positions_.clear();

} 

void PaddedBasevector::addPadPosition( unsigned int pad_index )
{  
   pad_positions_.push_back( pad_index ); 
   ++number_pads_; 
}


void PaddedBasevector::setBasevector( const basevector &theBasevector ) 
{ the_basevector_ = theBasevector; }


void PaddedBasevector::setReadName( const String &readName )
{ read_name_ = readName; }

void PaddedBasevector::setComplementation( bool Complementation )
{ complementation_ = Complementation; }  
  


ostream &operator<<(ostream &out, const PaddedBasevector &padded_sequence) 
{
  vec<char> s = padded_sequence.getSequence();
  for (unsigned int j=0; j<s.size(); ++j)
  {
    if ( (j%50) != 0 )         //  ace files seem to like 50 chars/line
     out << s[j];
    else
      out << "\n" << s[j];
  }
  out << endl;
  return out;
}


void PaddedBasevector::DisplaySequence( unsigned int start_position, unsigned int number_bases_to_display ) const
{

  if ( number_bases_to_display == 0 )
    number_bases_to_display = PaddedLength();

  vec<char> s = getSequence(); 
  unsigned int end_pos = start_position+number_bases_to_display;

  if ( end_pos>PaddedLength() )
  {
    cerr << " Requested range exceeds the number of bases in the sequence. "<<endl;
    end_pos = PaddedLength();
  }

  for (unsigned int j=start_position; j<end_pos; ++j)
    cout<< s[j];
  cout<<"\n"<<endl;
  
}

vec<char> PaddedBasevector::getSequence() const
{
  //  return a vector of chars consisting of the bases and pads.  For
  //  use in the alignment code (SmithWatBandedA)

  vec<char> padded_sequence( PaddedLength()  );
  unsigned int pad_index(0), sequence_index(0), pads_counter(0), basevector_index(0);


  //  if there are no pads, we just have the original basevector
  if ( NumberOfPads() == 0 )
  {  
    for (unsigned int i=0; i < the_basevector_.size(); ++i)
      padded_sequence[i] = as_base(the_basevector_[i]);
  }
  else
  {

    //  put the bases from basevector into the sequence, adding pads where 
    //  directed by pad_positions
    vec<unsigned int> sorted_pad_positions = getPadPositions();
    sort( sorted_pad_positions.begin(), sorted_pad_positions.end() );
 
    while( pads_counter < NumberOfPads() )
    {
	pad_index = sorted_pad_positions[ pads_counter ]; 
	if ( pad_index >= the_basevector_.size() )
	  break;

      //  add the bases to the sequence up until the pad
      while( basevector_index < pad_index )
      {
	padded_sequence[ sequence_index ] = as_base(the_basevector_[ basevector_index ]);
	
	++sequence_index;
	++basevector_index;
      }

      //  now add the pad to the sequence
      padded_sequence[ sequence_index ] = '*';
      ++sequence_index; 
      ++pads_counter;

    }

    if ( pad_index < the_basevector_.size() )
    {
      //  pick up the stragglers after the last pad
      while( basevector_index < the_basevector_.size())
      {  
	padded_sequence[ sequence_index ] = as_base(the_basevector_[ basevector_index ]);
	++sequence_index;  
	++basevector_index; 
      }

    }
  }

  return padded_sequence;
 
}
