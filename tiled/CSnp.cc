/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "String.h"
#include "math/Functions.h"
#include "system/Assert.h"
#include "tiled/CharBaseUtil.h"
#include "tiled/CSnp.h"

/**
 * class CSnp
 * Constructor
 */
CSnp::CSnp( ) :
  b1_( empty_base ),
  q1_( empty_qual ),
  b2_( empty_base ),
  q2_( empty_qual )
{ }
  
/**
 * class CSnp
 * Constructor
 */
CSnp::CSnp( char b1, int q1, char b2, int q2 ) :
  b1_( b1 ),
  q1_( q1 ),
  b2_( b2 ),
  q2_ ( q2 )
{
  this->Standardize( );
}

/**
 * class CSnp
 * Set
 */
void CSnp::Set( char b1, int q1, char b2, int q2 )
{
  b1_ = b1;
  q1_ = q1;
  b2_ = b2;
  q2_ = q2;

  this->Standardize( );
}


#ifdef SKIP
/**
 * class CSnp
 * Unpack
 */
void CSnp::Unpack( char snp_base, int snp_qual )
{
  ForceAssert( ! IsEmpty( snp_base ) );
  ToUpperCase( snp_base );
  
  // Not a SNP.
  if ( IsGap( snp_base ) || IsConventionalBase( snp_base ) ) {
    this->Set( snp_base, snp_qual );
    return;
  }

  // SNP qual.
  int q2 = ( snp_qual & 0x000000FF );
  int q1 = ( snp_qual & 0x0000FF00 ) >> 8;
  
  // SNP base.
  if ( snp_base == snp_AC_base ) this->Set( 'A', q1, 'C', q2 );
  else if ( snp_base == snp_AG_base ) this->Set( 'A', q1, 'G', q2 );
  else if ( snp_base == snp_AT_base ) this->Set( 'A', q1, 'T', q2 );
  else if ( snp_base == snp_Apad_base ) this->Set( 'A', q1, gap_base, q2 );
  else if ( snp_base == snp_CG_base ) this->Set( 'C', q1, 'G', q2 );
  else if ( snp_base == snp_CT_base ) this->Set( 'C', q1, 'T', q2 );
  else if ( snp_base == snp_Cpad_base ) this->Set( 'C', q1, gap_base, q2 );
  else if ( snp_base == snp_GT_base ) this->Set( 'G', q1, 'T', q2 );
  else if ( snp_base == snp_Gpad_base ) this->Set( 'G', q1, gap_base, q2 );
  else if ( snp_base == snp_Tpad_base ) this->Set( 'T', q1, gap_base, q2 );
  else ForceAssert( 1 == 0 );
}
#endif

/**
 * class CSnp
 * Pack
 */
void CSnp::Pack( char &snp_base, int &snp_qual, const int *cap ) const
{
  // Not a SNP.
  if ( IsEmpty( b2_ ) ) {
    snp_base = b1_;
    snp_qual = cap ? Min( *cap, q1_ ) : q1_;
    return;
  }
  
  // Put pads (if any) in b2. Sort bases.
  char b1 = b1_;
  char b2 = b2_;
  int q1 = cap ? Min( *cap, q1_ ) : q1_;
  int q2 = cap ? Min( *cap, q2_ ) : q2_;
  if ( IsGap( b1 ) ) {
    swap( b1, b2 );
    swap( q1, q2 );
  }
  else if ( ( !IsGap( b2 ) ) && ( b2 < b1 ) ) {
    swap( b1, b2 );
    swap( q1, q2 );
  }
  
  // SNP base.
  if ( b1 == 'A' ) {
    if ( b2 == 'C' ) snp_base = snp_AC_base;
    else if ( b2 == 'G' ) snp_base = snp_AG_base;
    else if ( b2 == 'T' ) snp_base = snp_AT_base;
    else if ( IsGap( b2 ) ) snp_base = snp_Apad_base;
    else ForceAssert( 1 == 0 );
  }
  else if ( b1 == 'C' ) {
    if ( b2 == 'G' ) snp_base = snp_CG_base;
    else if ( b2 == 'T' ) snp_base = snp_CT_base;
    else if ( IsGap( b2 ) ) snp_base = snp_Cpad_base;
    else ForceAssert( 1 == 0 );
  }
  else if ( b1 == 'G' ) {
    if ( b2 == 'T' ) snp_base = snp_GT_base;
    else if ( IsGap( b2 ) ) snp_base = snp_Gpad_base;
    else ForceAssert( 1 == 0 );
  }
  else if ( b1 == 'T' ) {
    if ( IsGap( b2 ) ) snp_base = snp_Tpad_base;
    else ForceAssert( 1 == 0 );
  }
  else ForceAssert( 1 == 0 );

  // SNP qual.
  snp_qual = ( q1 << 8 ) | q2;
}

/**
 * class CSnp
 * IsEmptySnp
 */
bool CSnp::IsEmptySnp( ) const
{
  return ( IsEmpty( b1_ ) && IsEmpty( b2_ ) );
}

/**
 * class CSnp
 * Base
 */
void CSnp::Base( char &base, int &qual, bool second ) const
{
  base = second ? b2_ : b1_;
  qual = second ? q2_ : q1_;
}

/**
 * class CSnp
 * PrintBase
 */
String CSnp::PrintBase( bool second ) const
{
  String base;
  base = ( second ? b2_ : b1_ );
  String qual = ToString( second ? q2_ : q1_ );

  return "(" + base + "," + qual + ")";
}

/**
 * class CSnp
 * PrintSnp
 */
String CSnp::PrintSnp( ) const
{
  String base1;
  base1 = b1_;
  String qual1 = ToString( q1_ );
  String str1 = "(" + base1 + "," + qual1 + ")";

  String base2 = "";
  String qual2 = "";
  String str2 = "";
  if ( !IsEmpty( b2_ ) ) {
    base2 = b2_;
    qual2 = ToString( q2_ );
    str2 = " (" + base2 + "," + qual2 + ")";
  }

  return str1 + str2;
}

/**
 * class CSnp
 * operator==
 */
bool operator== ( const CSnp &left, const CSnp &right )
{
  if ( left.b1_ == right.b1_ && left.b2_ == right.b2_ ) return true;
  if ( left.b1_ == right.b2_ && left.b2_ == right.b1_ ) return true;
  return false;
}

/**
 * class CSnp
 * operator!=
 */
bool operator!= ( const CSnp &left, const CSnp &right )
{
  return ( ! ( left == right ) );
}

/**
 * class CSnp
 * operator+
 */
CSnp operator+ ( const CSnp &left, const CSnp &right )
{
  ForceAssert( left == right );
  
  char b1Left = left.b1_;
  char b2Left = left.b2_;
  int q1Left = left.q1_;
  int q2Left = left.q2_;
  
  char b1Right = right.b1_;
  char b2Right = right.b2_;
  int q1Right = right.q1_;
  int q2Right = right.q2_;
  
  if ( b1Left != b1Right ) {
    swap( b1Right, b2Right );
    swap( q1Right, q2Right );
  }
  ForceAssert( b1Left == b1Right && b2Left == b2Right );

  int q1Sum = q1Left + q1Right;
  int q2Sum = q2Left + q2Right;
  
  CSnp result( b1Left, q1Sum, b2Left, q2Sum );
  return result;
}

/**
 * class CSnp
 * Standardize
 */
void CSnp::Standardize( )
{
  // Empty case first.
  if ( IsEmpty( b1_ ) && IsEmpty( b2_ ) ) {
    ForceAssert ( q1_ == empty_qual );
    ForceAssert( q2_ == empty_qual );
    return;
  }

  // Cap quality scores (2^15).
  q1_ = Min( q1_, 32768 );
  q2_ = Min( q2_, 32768 );
  
  // Allow for gap_quals only if (b2,q2) is empty.
  ForceAssert( IsEmpty( b2_ ) || q1_ >= 0 );

  // Either a regular base, or a true SNP.
  ForceAssert( IsEmpty( b2_ ) || ( b1_ != b2_ ) );
  
  // Force upper case on bases.
  if ( !IsEmpty( b1_ ) ) ToUpperCase( b1_ );
  if ( !IsEmpty( b2_ ) ) ToUpperCase( b2_ );

  // Store in (b1,q1) the base with highest qual.
  if ( ( !IsEmpty( b2_ ) ) && ( q1_ < q2_ ) ) { 
    swap( b1_, b2_ );
    swap( q1_, q2_ );
  }
}

