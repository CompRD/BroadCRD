// Copyright (c) 2000, 2001 Whitehead Institute of Biomedical Research
// 
// Class encpsulating the idea of a basevector with pads, where
// a pad represents a gap introduced by the alignment of the basevector
// to another basevector (i.e. read vs. read, read vs. consensus)
//
// Used in the creatation of ace files from Arachne assemblies.
//


#ifndef PADDEDBASEVECTOR
#define PADDEDBASEVECTOR

#include <algorithm>
#include "Basevector.h"
#include "Vec.h"
#include "ShortVector.h"
#include "String.h"



class PaddedBasevector {
  
  public:

    PaddedBasevector();

    PaddedBasevector( const PaddedBasevector &padded_sequence );

    PaddedBasevector( const basevector &theBasevector, 
		      const vec<unsigned int> &padPositions, 
		      String readName="", 
		      bool Complementation=False );

    PaddedBasevector &operator=(const PaddedBasevector &rhs_sequence);

    ~PaddedBasevector();

    unsigned int NumberOfPads() const;
    unsigned int PaddedLength() const;

  // get pad positions in unpadded space (before which base does a pad occur)
    vec<unsigned int> getPadPositions() const;
  // get pad positions in padded space (where does a pad occur in the padded sequence)
    vec<unsigned int> getTranslatedPadPositions() const;

  // get number of pads in range (range is in padded space)
    int getPadsInRange( int padded_begin, int padded_end ) const;

  // translate from padded to unpadded space
    int getUnpaddedPosition( int padded_position );
  // translate from unpadded to padded space
    int getPaddedPosition( int unpadded_position );

    const basevector getBasevector() const;
    const String getReadName() const;
    bool getComplementation() const;
    char returnRC() const;
    

    void setPadPositions( const vec<unsigned int> &pads );
    void addPadPosition( unsigned int pad_index );

    void setBasevector( const basevector &theBasevector );
    void setReadName( const String & readName );
    void setComplementation( bool Complementation );

    //  unsigned char operator()(const unsigned int &i) const;

    void DisplaySequence( unsigned int start_position=0, 
			  unsigned int number_bases_to_display=0 ) const;

    vec<char> getSequence() const;

    friend ostream &operator<<(ostream &out, const PaddedBasevector &padded_sequence);

  private:

    void updateTranslatedPadPositions_() const;

    basevector the_basevector_;
    unsigned int number_pads_;
    String read_name_;
    bool complementation_;
    vec<unsigned int> pad_positions_;
    mutable vec<unsigned int> translated_positions_;
};

#endif
