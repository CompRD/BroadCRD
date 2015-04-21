/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "text/RangedRegEx.h"
#include "Vec.h"

String Digit( int d ) {    
  static String digits[10] = { String("0"), 
			       String("1"), 
			       String("2"), 
			       String("3"),
			       String("4"), 
			       String("5"), 
			       String("6"), 
			       String("7"), 
			       String("8"), 
			       String("9") };
  return digits[d];
}

String DigitRange( int d1, int d2 ) {    
  ForceAssert( d1 <= d2 && d1 >= 0 && d2 <= 9 );
  if ( d1 == d2 ) 
    return Digit(d1);
  else if ( d1 + 1 == d2 )
    return String("[") + Digit(d1) + Digit(d2) + "]";
  else 
    return String("[") + Digit(d1) + "-" + Digit(d2) + "]";    
}

String Dup( String s, int n ) {    
  String t;
  for ( int i = 0; i < n; i++ )
    t += s;
  return t;
}

// RangeToRegex accepts as input two positive integers, from <= to.
// It produces as output an extended regular expression which matches any
// integer in the range [from, to].

String RangeToRegex( int from, int to ) {    
  ForceAssert( from >= 1 && from <= to );

  vec<int> from_digits, to_digits;

  int x = from;
  while(x != 0) {    
    from_digits.insert( from_digits.begin( ), x % 10 );
    x = x / 10;
  }

  x = to;
  while(x != 0) {
    to_digits.insert( to_digits.begin( ), x % 10 );
    x = x / 10;
  }

  int from_places = from_digits.size( ), to_places = to_digits.size( );

  vec<String> subexpressions;

  // If from_places = to_places, generate pattern for [from, to].

  if ( from_places == to_places ) {

    // Let agreeing_digits = the number of digits of agreement.  Put in prefix.
    // e.g. from = 1234, to = 1256, generate prefix 12
    String prefix;
    int agreeing_digits;
    for ( agreeing_digits = 0; agreeing_digits < to_places; agreeing_digits++ ) {
      if ( from_digits[agreeing_digits] == to_digits[agreeing_digits] ) 
	prefix += Digit( to_digits[agreeing_digits] );
      else 
	break;    
    }

    // If the number of agreeing digits == number of digits in to 
    // (which == number of digits in from), then from and to are 
    // the same (i.e. we're done)
    if ( agreeing_digits == to_places ) 
      return prefix;

    // If they agree up to the next to the last digit, the desired expression
    // is the prefix and the range from the last digit of from to the last 
    // digit of to.
    // e.g. from = 1234, to = 1238, return 123[4-8]
    if ( agreeing_digits == to_places-1 ) 
      return prefix + DigitRange( from_digits[agreeing_digits], 
				  to_digits[agreeing_digits] );

    // Generate pattern to raise first nonagreeing digit.
    // Example: if from = 1234, to = 1256, generate pattern for 1234 - 1239.
    // Example: if from = 1072, to = 1358, generate pattern for 1072 - 1099.

    for ( int i = to_places-1; i > agreeing_digits; i-- ) {
      String subex = prefix + Digit(from_digits[agreeing_digits]);
      for ( int j = agreeing_digits + 1; j < i; j++ )
	subex += Digit( from_digits[j] );
      int begin = (i == to_places-1) ? from_digits[i] : from_digits[i] + 1;
      if ( begin > 9 ) 
	continue;
      subex += DigitRange( begin, 9 );
      for ( int j = i + 1; j < to_places; j++ )
	subex += DigitRange(0, 9);
      subexpressions.push_back(subex);    
    }

    // Generate pattern for the rest.
    // Example: if from = 1234, to = 1256, generate pattern for 1240 - 1256.
    // Example: if from = 1072, to = 1358, generate pattern for 1100 - 1358.

    bool diff = false;
    for ( int i = agreeing_digits; i < to_places; i++ ) {
      diff = false;
      if ( ( i == agreeing_digits && to_digits[i] > from_digits[i] + 1 ) ||
	   ( i > agreeing_digits && to_digits[i] > 0 ) )
	diff = true;
      else
	continue;
      String subex = prefix;
      for ( int j = agreeing_digits; j < i; j++ )
	subex += Digit( to_digits[j] );
      int first = (i == agreeing_digits) ? from_digits[i] + 1 : 0;
      int last = (i < to_places - 1) ? to_digits[i] - 1 : to_digits[i];
      if ( first < last ) 
	subex += DigitRange( first, last );
      else 
	subex += Digit(first);
      for ( int j = i + 1; j < to_places; j++ )
	subex += DigitRange(0, 9);
      if ( subexpressions.size( ) == 0 || subex != subexpressions.back( ) ) 
	subexpressions.push_back(subex);    
    }
    
    if ( !diff ) {
      String subex = prefix;
      for ( int j = agreeing_digits; j < to_places; j++ )
	subex += Digit( to_digits[j] );
      subexpressions.push_back(subex);    
    }    
  }

  // If from_places < to_places, generate pattern for [from, 10^from_places - 1].

  if ( from_places < to_places ) {
    for ( int i = from_places-1; i >= 0; i-- ) {
      String subex;
      for ( int j = 0; j < i; j++ )
	subex += Digit( from_digits[j] );
      int begin = (i == from_places-1) ? from_digits[i] : from_digits[i] + 1;
      if ( begin == 10 ) 
	continue;
      // int begin = from_digits[i];
      subex += DigitRange( begin, 9 );
      for ( int j = i + 1; j < from_places; j++ )
	subex += DigitRange(0, 9);
      subexpressions.push_back(subex);    
    } 
  }

  // If from_places < to_places-1, generate pattern for 
  // [10^from_places, 10^(to_places-1) - 1].

  for ( int i = from_places+1; i <= to_places-1; i++ )
    subexpressions.push_back( DigitRange(1, 9) + Dup( DigitRange(0, 9), i-1 ) );

  // If from_places < to_places, generate pattern for [10^(to_places-1), to].

  if ( from_places < to_places ) {
    bool diff = false;
    for ( int i = 0; i < to_places; i++ ) {
      diff = false;
      if ( ( i == 0 && to_digits[i] > 1 ) ||
	   ( i > 0 && to_digits[i] > 0 ) )
	diff = true;
      else
	continue;
      String subex;
      for ( int j = 0; j < i; j++ )
	subex += Digit( to_digits[j] );
      int first = (i == 0) ? 1 : 0;
      int last = (i < to_places - 1) ? to_digits[i] - 1 : to_digits[i];
      if ( first < last ) 
	subex += DigitRange( first, last );
      else 
	subex += Digit(first);
      for ( int j = i + 1; j < to_places; j++ )
	subex += DigitRange(0, 9);
      if ( subexpressions.size( ) == 0 || subex != subexpressions.back( ) ) 
	subexpressions.push_back(subex);    
    }
    if ( !diff ) {
      String subex;
      for ( int j = 0; j < to_places; j++ )
	subex += Digit( to_digits[j] );
      subexpressions.push_back(subex); 
    }
  }

  String answer;
  for ( unsigned int i = 0; i < subexpressions.size( ); i++ ) {
    if ( i > 0 ) answer += "|";
    answer += "(" + subexpressions[i] + ")";    
  }
  return answer;    
}


void SubstituteIntegerRanges( String& expression ) {
  // find '[#'
  for ( unsigned int range_start = 0; 
	range_start < expression.size() - 1 ;
	range_start++ ) {
    if ( expression[range_start] == '[' &&
	 expression[range_start+1] == '#' ) {
      // when we find a '[#', look for a matching ']'
      unsigned int range_end;
      for ( range_end = range_start + 1; 
	    range_end < expression.size(); 
	    range_end++ )
	if ( expression[range_end] == ']' ) break;
      // if we've exited the above loop before the end of expression, we found a ']'
      if ( range_end < expression.size() ) {
	unsigned int range_length = range_end - range_start + 1;
	String range = expression.substr(range_start + 2, range_length - 3);
	if ( range.Contains( "-" ) ) {
	  String from = range.Before( "-" );
	  String to = range.After( "-" );
	  unsigned int place;
	  for ( place = 0; place < from.size( ); place++ )
	    if ( !isdigit( from[place] ) ) 
	      break;
	  if ( place < from.size( ) ) 
	    break;
	  for ( place = 0; place < to.size( ); place++ )
	    if ( !isdigit( to[place] ) ) 
	      break;
	  if ( place < to.size( ) ) 
	    break;
	  if ( from.Int( ) > to.Int( ) ) 
	    break;
	  if ( from.Int( ) <= 0 ) 
	    break;
	  String regex = "(" + RangeToRegex( from.Int( ), to.Int( ) ) + ")";
	  expression.replace( range_start, range_length, regex.c_str( ) );
	  range_start += regex.size() - 1; // move range_start to end of new regex
	}
      }
    }
  }
}

