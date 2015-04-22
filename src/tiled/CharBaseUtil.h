/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// Utilities to lower/upper case of bases, needed by Tiling and PaddedSeq.

#ifndef CHAR_BASE_UTIL_H
#define CHAR_BASE_UTIL_H

#include "dna/Bases.h"

/**
 * Hard coded base types and quals.
 */
const char snp_AC_base      = 'M';
const char snp_AG_base      = 'R';
const char snp_AT_base      = 'W';
const char snp_Apad_base    = 'E';
const char snp_CG_base      = 'S';
const char snp_CT_base      = 'Y';
const char snp_Cpad_base    = 'F';
const char snp_GT_base      = 'K';
const char snp_Gpad_base    = 'J';
const char snp_Tpad_base    = 'P';
const char gap_base         = '*';
const char lowqual_gap_base = '#';
const char empty_base       = '.';
const int gap_qual          = -1;
const int empty_qual        = -2;
const int snp_qual          = -3;
const int untrusted_qual    = -4;
const int fin_qual          = 50;    // Default value for finished qual.

/**
 * IsEmpty
 */
inline bool IsEmpty( char base ) {
  return ( base == empty_base );
}

/**
 * IsGap
 */
inline bool IsGap( char base ) {
  return ( base == gap_base || base == lowqual_gap_base );
}

/**
 * IsLowerCase
 */
inline bool IsLowerCase( char base ) {
  return ( ( 96 < base && base < 123 ) || base == lowqual_gap_base );
}

/**
 * IsUpperCase
 */
inline bool IsUpperCase( char base ) {
  return ( ( 64 < base && base < 91 ) || base == gap_base );
}

/**
 * ToLowerCase
 * On error set base to 'x'.
 */
inline void ToLowerCase( char &base ) {
  if ( IsLowerCase( base ) ) ; // Already lower case, do nothing.
  else if ( 64 < base && base < 91 ) base += 32;
  else if ( base == gap_base ) base = lowqual_gap_base;
  else base = 'x';
}

/**
 * ToUpperCase
 * On error set base to 'x'.
 */
inline void ToUpperCase( char &base ) {
  if ( IsUpperCase( base ) ) ; // Already upper case, do nothing.
  else if ( 96 < base && base < 123 ) base += -32;
  else if ( base == lowqual_gap_base ) base = gap_base;
  else base = 'x';
}

/**
 * IsConventionalBase
 * A conventional base is one of {a,c,g,t,A,C,G,T}.
 */
inline bool IsConventionalBase( char base ) {
    return Base::isBase(base);
}

#endif
