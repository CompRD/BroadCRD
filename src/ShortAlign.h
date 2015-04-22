// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

#ifndef SHORT_ALIGN_H
#define SHORT_ALIGN_H

// A short_align consists of numerical indices for two sequences, a flag to 
// indicate if the second sequence is reverse complemented in the alignment,
// an offset (represented as a short), and a score (represented as an unsigned 
// short).  The alignment detail is not part of the data.

#include "CoreTools.h"
#include "feudal/BinaryStreamTraits.h"

class short_align {

     public:

     short_align( ) { }
     short_align( int id1, int id2, Bool rc2, short offset, unsigned short score )
          : id1_(id1), offset_(offset), score_(score)
     {    id2_rc_ = ( !rc2 ? +(id2 + 1) : -(id2 + 1) );    }

     void Set( int id1, int id2, Bool rc2, short offset, unsigned short score )
     {    id1_ = id1;
          offset_ = offset;
          score_ = score;
          id2_rc_ = ( !rc2 ? +(id2 + 1) : -(id2 + 1) );    }

     void SetScore( unsigned short score )
     {    score_ = score;    }

     int Id1( ) const { return id1_; }
     int Id2( ) const { return ( id2_rc_ > 0 ? (id2_rc_ - 1) : (-id2_rc_ - 1) ); }

     Bool Rc2( ) const { return id2_rc_ < 0; }

     short Offset( ) const { return offset_; }
     unsigned short Score( ) const { return score_; }

     void SetId1( int id1 ) { id1_ = id1; }
     void SetId2( int id2 ) { id2_rc_ = id2 + 1; }

     // Note: SetRc2 must be called AFTER SetId2.

     void SetRc2( Bool rc ) { if (rc) id2_rc_ = -id2_rc_; }

     friend Bool operator<( const short_align& a1, const short_align& a2 )
     {    if ( a1.Id1( ) < a2.Id1( ) ) return True;
          if ( a1.Id1( ) > a2.Id1( ) ) return False;
          if ( a1.Id2( ) < a2.Id2( ) ) return True;
          return False;    }

     private:

     int id1_;
     int id2_rc_;
     short offset_;
     unsigned short score_;

};
TRIVIALLY_SERIALIZABLE(short_align);

// BuildIndex assumes that the aligns are ordered by id1.

void BuildIndex( const vec<short_align>& aligns, vec<int>& index, int nseq );

#endif
