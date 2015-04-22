/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Corrector.h"

/**
 * corrector
 * Constructor
 */
corrector::corrector( ) :
  id_ ( -1 ),
  pos_ ( -1 ),
  old_base_ ( empty_base ),
  new_base0_ ( empty_base ),
  new_base1_ ( empty_base ),
  old_qual_ ( empty_qual ),
  new_qual0_ ( empty_qual ),
  new_qual1_ ( empty_qual )
{ }

/**
 * corrector
 * Constructor
 */
corrector::corrector( int id,
		      int pos,
		      char oldB,
		      char newB,
		      int oldQ,
		      int newQ ) :
  id_ ( id ),
  pos_ ( pos ),
  old_base_ ( oldB ),
  new_base0_ ( newB ),
  new_base1_ ( empty_base ),
  old_qual_ ( oldQ ),
  new_qual0_ ( newQ ),
  new_qual1_ ( empty_qual )
{ }

/**
 * corrector
 * Constructor
 */
corrector::corrector( int id,
		      int pos,
		      char oldB,
		      char newB_0,
		      char newB_1,
		      int oldQ,
		      int newQ_0,
		      int newQ_1) :
  id_ ( id ),
  pos_ ( pos ),
  old_base_ ( oldB ),
  new_base0_ ( newB_0 ),
  new_base1_ ( newB_1 ),
  old_qual_ ( oldQ ),
  new_qual0_ ( newQ_0 ),
  new_qual1_ ( newQ_1 )
{ }

/**
 * corrector
 * OldBase
 *
 * Return a string representing the old base (the old base only!) as in
 * "(G,32)".
 */
String corrector::OldBase( ) const
{
  String base;
  base = old_base_;
  String qual = ToString( old_qual_ );
  String str = "(" + base + "," + qual + ")";

  return str;
}

/**
 * corrector
 * NewBase
 *
 * Return a string representing the new base (the new base only!) as a
 * SNP, as "(G,50)", or "(G,50) (*,50)".
 */
String corrector::NewBase( ) const
{
  String base0;
  base0 = new_base0_;
  String qual0 = ToString( new_qual0_ );
  String str0 = "(" + base0 + "," + qual0 + ")";

  String base1 = "";
  String qual1 = "";
  String str1 = "";
  if ( !IsEmpty( new_base1_ ) ) {
    base1 = new_base1_;
    qual1 = ToString( new_qual1_ );
    str1 = " (" + base1 + "," + qual1 + ")";
  }

  return str0 + str1;
}

/**
 * corrector
 * BinaryOut
 */
void corrector::BinaryOut( std::ostream &out ) const
{
  out.write( (char *) &( id_ ), sizeof( int ) );
  out.write( (char *) &( pos_ ), sizeof( int ) );
  out.write( (char *) &( old_base_ ), sizeof( char ) );
  out.write( (char *) &( new_base0_ ), sizeof( char ) );
  out.write( (char *) &( new_base1_ ), sizeof( char ) );
  out.write( (char *) &( old_qual_ ), sizeof( int ) );
  out.write( (char *) &( new_qual0_ ), sizeof( int ) );
  out.write( (char *) &( new_qual1_ ), sizeof( int ) );
}

/**
 * corrector
 * BinaryIn
 */
void corrector::BinaryIn ( std::istream &in )
{
  in.read( (char *) &( id_ ), sizeof( int ) );
  in.read( (char *) &( pos_ ), sizeof( int ) );
  in.read( (char *) &( old_base_ ), sizeof( char ) );
  in.read( (char *) &( new_base0_ ), sizeof( char ) );
  in.read( (char *) &( new_base1_ ), sizeof( char ) );
  in.read( (char *) &( old_qual_ ), sizeof( int ) );
  in.read( (char *) &( new_qual0_ ), sizeof( int ) );
  in.read( (char *) &( new_qual1_ ), sizeof( int ) );
}

/**
 * corrector
 * operator<<
 */
std::ostream &operator<< ( std::ostream &out, const corrector &crc )
{
  out << crc.id_ << "\t"
      << crc.pos_ << "\t"
      << crc.old_base_ << "\t"
      << crc.new_base0_ << "\t"
      << crc.new_base1_ << "\t"
      << crc.old_qual_ << "\t"
      << crc.new_qual0_ << "\t"
      << crc.new_qual1_ << "\n";

  return out;
}

/**
 * corrector
 * operator>>
 */
std::istream &operator>> ( std::istream &in, corrector &crc )
{
  in >> crc.id_
     >> crc.pos_
     >> crc.old_base_
     >> crc.new_base0_
     >> crc.new_base1_
     >> crc.old_qual_
     >> crc.new_qual0_
     >> crc.new_qual1_;
  
  return in;
}

