// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology
//

#include "tiled/CKDiff.h"



/*
 * ck_diff
 * Constructor
 */
ck_diff::ck_diff( ) :
  n_total ( 0 ),
  n_pads_k ( 0 ),
  n_pads_c ( 0 ),
  n_confirm_c ( 0 ),
  n_confirm_k ( 0 ),
  n_confirm_both ( 0 ),
  n_confirm_none ( 0 )
{ }



/*
 * ck_diff
 * TotalDiff
 *
 * Total number of differences.
 */
int ck_diff::TotalDiff( ) const
{
  return ( n_confirm_c + n_confirm_k + n_confirm_both + n_confirm_none );
}



/*
 * ck_diff
 * PrintCompact
 */
void ck_diff::PrintCompact( ostream &out, bool oneline ) const
{
  String str_c
    = n_confirm_c > 0 ? " confirm_c=" + ToString( n_confirm_c ) : "";
  String str_k
    = n_confirm_k > 0 ? " confirm_k=" + ToString( n_confirm_k ) : "";
  String str_both
    = n_confirm_both > 0 ? " confirm_both=" + ToString( n_confirm_both ) : "";
  String str_none
    = n_confirm_none > 0 ? " confirm_neither=" + ToString( n_confirm_none ) : "";
  String appendix = oneline ? "" : "\n";

  out << "diff=" << this->TotalDiff( )
      << " (out of " << n_total
      << " checked):"
      << str_c
      << str_k
      << str_both
      << str_none
      << appendix;
}



/*
 * ck_diff
 * Print
 */
void ck_diff::Print( ostream &out, bool oneline, bool newline ) const
{
  String separator = ( oneline ) ? " " : "\n";
  String ending = ( newline ) ? "\n" : "";
  
  out << "total=" << n_total << separator
      << "diff=" << this->TotalDiff( ) << separator
      << "pads_k=" << n_pads_k << separator
      << "pads_c=" << n_pads_c << separator
      << "confirm_c=" << n_confirm_c << separator
      << "confirm_k=" << n_confirm_k << separator
      << "confirm_both=" << n_confirm_both << separator
      << "confirm_neither=" << n_confirm_none << ending;
}



/*
 * ck_diff
 * operator+=
 */
ck_diff &ck_diff::operator+= ( const ck_diff &addendum )
{
  n_total += addendum.n_total;
  n_pads_k += addendum.n_pads_k;
  n_pads_c += addendum.n_pads_c;
  n_confirm_c += addendum.n_confirm_c;
  n_confirm_k += addendum.n_confirm_k;
  n_confirm_both += addendum.n_confirm_both;
  n_confirm_none += addendum.n_confirm_none;

  return *this;
}
