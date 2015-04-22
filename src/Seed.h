// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#ifndef SEED_H
#define SEED_H

#include <set>

#include "Alignment.h"
#include "math/Arith.h"
#include "CoreTools.h"
#include "math/Functions.h"
     
// ================================================================================
//
// class unit
//
// ================================================================================

// A unit should be thought of as a pair of reads (with identifiers "left_read" and
// "right_read"), which come from the same chromosome and whose approximate 
// separation "sep" is known.  But units exist in seeds (see below), which are 
// vectors of units, and the unit data also includes the position of the given unit 
// relative to the first unit in the seed.  Specifically, "left_offset" is the 
// position of the starting position of the left_read relative to the left_read of 
// the first unit, and "right_offset" is similarly defined.

// Example:    ----------------->                 <-----------------  first unit
//                --------------->               <-------------       this unit
//
// left_offset = 3, right_offset = 5

class unit {

     public:

     int LeftRead( ) const { return left_read_; }
     int RightRead( ) const { return right_read_; }
     int LeftLength( ) const { return left_length_; }
     int RightLength( ) const { return right_length_; }
     int Sep( ) const { return sep_; }
     int LeftOffset( ) const { return left_offset_; }
     int RightOffset( ) const { return right_offset_; }

     void SetLeftOffset( int o ) { left_offset_ = o; }
     void SetRightOffset( int o ) { right_offset_ = o; }

     int LeftStart( ) const { return left_offset_; }
     int LeftStop( ) const { return left_offset_ + left_length_; }
     int RightStart( ) const { return right_offset_; }
     int RightStop( ) const { return right_offset_ + right_length_; }

     unit( int left_read, int right_read, int left_length, int right_length, 
          int sep ) : left_read_(left_read), right_read_(right_read), 
          left_length_(left_length), right_length_(right_length), sep_(sep), 
          left_offset_(0), right_offset_(0) { }
     unit( ) { }

     void Shift( int left_shift, int right_shift )
     {    left_offset_ += left_shift;
          right_offset_ += right_shift;    }

     void Swap( )
     {    swap( left_read_, right_read_ );
          swap( left_length_, right_length_ );
          swap( left_offset_, right_offset_ );    }

     private:

     int left_read_, right_read_;
     int left_length_, right_length_;
     int sep_;
     int left_offset_, right_offset_;

};

// Given units u1 and u2, SepDiscrep( u1, u2 ) returns the discrepancy between their 
// separations, as a signed integer.  It is assumed that their offsets are related
// to each other, as if u1 and u2 came from the same seed.  If u2.sep seems too big,
// the value returned is positive.

inline int SepDiscrep( const unit& u1, const unit& u2 )
{    return u2.Sep( ) - ( u1.Sep( )  
          - ( u2.LeftOffset( ) - u1.LeftOffset( ) + u2.LeftLength( ) 
               - u1.LeftLength( ) ) 
          - ( u2.RightOffset( ) - u1.RightOffset( ) + u2.RightLength( ) 
               - u1.RightLength( ) ) );    }

// ================================================================================
//
// class seed
//
// ================================================================================

// A seed is a collection of units, organized as a vector, such that the overlap
// graph of the left reads is connected, the overlap graph of the right reads
// is connected, and the unit separations are reasonably consistent.

// A seed can be initialized from a single unit, which should have left_offset = 0 
// and right_offset = 0.

// To try to add a unit u to a seed, use:
//
//   Merge( u, left, right, a_left, a_right ),
//
// where (*this)[left].idL aligns with u.idL via a_left, and
// (*this)[right].idR aligns with u.idR via a_right.  Note that the  
// left_offset and right_offset members of u are set as part of the merger process.

// Merge tests to see if the unit separations are reasonably consistent.  If so, it
// appends the unit to the seed and returns true.

// There is a different Merge, wherein a given seed absorbs another seed, assuming
// of course that the two seeds share at least one read.

class seed : public vec<unit> {

     public:
     
     seed( ) { }
     seed( const unit& u ) { push_back(u); }

     int LeftMember( int read ) const
     {    for ( unsigned int i = 0; i < size( ); i++ )
               if ( (*this)[i].LeftRead( ) == read ) return i;
          return -1;    }

     int RightMember( int read ) const
     {    for ( unsigned int i = 0; i < size( ); i++ )
               if ( (*this)[i].RightRead( ) == read ) return i;
          return -1;    }

     int LeftmostMember( ) const
     {    int best_start = (*this)[0].LeftOffset( ), best_read = 0;
          for ( unsigned int j = 1; j < size( ); j++ )
               if ( (*this)[j].LeftOffset( ) < best_start )
               {    best_read = j;
                    best_start = (*this)[j].LeftOffset( );    }
          return best_read;    }

     int RightmostMember( ) const
     {    int best_start = (*this)[0].RightOffset( ), best_read = 0;
          for ( unsigned int j = 1; j < size( ); j++ )
               if ( (*this)[j].RightOffset( ) < best_start )
               {    best_read = j;
                    best_start = (*this)[j].RightOffset( );    }
          return best_read;    }

     Bool Merge( unit u, int left, int right, const alignment& a_left,
          const alignment& a_right, const vec<alignment_plus>& all_aligns, 
          vec<alignment_plus>& extra_aligns, vec< vec<int> >& all_aligns_index, 
          const vecbasevector& EE, const vecqualvector& Q, 
          int MaxSeparationDiscrepancy, ostream& out,
          set< pair<int, int> >& no_alignment, Float max_score_for_compatibility );

     seed Swap( ) const
     {    seed s(*this);
          for ( unsigned int i = 0; i < size( ); i++ )
               s[i].Swap( );
          return s;    }

     Bool Merge( const seed& s2 );

     int LeftContigLength( ) const
     {    int left_start = (*this)[0].LeftStart( ), 
               left_stop = (*this)[0].LeftStop( );
          for ( unsigned int i = 1; i < size( ); i++ )
          {    left_start = Min( left_start, (*this)[i].LeftStart( ) );
               left_stop = Max( left_stop, (*this)[i].LeftStop( ) );    }
          return left_stop - left_start;    }

     int RightContigLength( ) const
     {    int right_start = (*this)[0].RightStart( ), 
               right_stop = (*this)[0].RightStop( );
          for ( unsigned int i = 1; i < size( ); i++ )
          {    right_start = Min( right_start, (*this)[i].RightStart( ) );
               right_stop = Max( right_stop, (*this)[i].RightStop( ) );    }
          return right_stop - right_start;    }
     
     friend ostream& operator<<( ostream& o, seed& s )
     {    o << "SEED with " << s.size( ) << " reads, having contig lengths "
               << s.LeftContigLength( ) << "/" << s.RightContigLength( ) << ":\n";
          for ( unsigned int j = 0; j < s.size( ); j++ )
               o << s[j].LeftRead( ) << " @" << s[j].LeftOffset( ) << " ";
          o << "##";
          for ( unsigned int j = 0; j < s.size( ); j++ )
               o << " " << s[j].RightRead( ) << " @" << s[j].RightOffset( );
          return o << "\n";    }
     
};

#endif
