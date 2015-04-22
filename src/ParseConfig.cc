// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


// This file contains the code which parses a project config file and evaluates
// read names using it.  See ParseConfigFile and TranslateReadName, below.
// Invoked by FetchAndTrim2.cc.

// See also ReadTypes.h.

#include <ctype.h>
#include <regex.h>

#include <algorithm>
#include <map>

#include "system/Assert.h"
#include "math/Functions.h"
#include "ParseConfig.h"
#include "ReadTypes.h"
#include "String.h"
#include "system/Types.h"
#include "Vec.h"

#include "system/System.h"

namespace
{

class markid {
     public:
     String var;
     String val;
     int pos;
     markid( String var_arg,
	     String val_arg,
	     int pos_arg )
       : var(var_arg),
	 val(val_arg),
	 pos(pos_arg)
     { }
};

String digits[10] = { String("0"), String("1"), String("2"), String("3"),
     String("4"), String("5"), String("6"), String("7"), String("8"), String("9") };

String Digit( int d )
{    return digits[d];    }

String DigitRange( int d1, int d2 )
{    ForceAssert( d1 <= d2 && d1 >= 0 && d2 <= 9 );
     if ( d1 == d2 ) return Digit(d1);
     else return String("[") + Digit(d1) + "-" + Digit(d2) + "]";    }

String Dup( String s, int n )
{    String t;
     for ( int i = 0; i < n; i++ )
          t += s;
     return t;    }

// RangeToRegex accepts as input two positive integers, from <= to.
// It produces as output an extended regular expression which matches any
// integer in the range [from, to].

String RangeToRegex( int from, int to )
{    ForceAssert( from >= 1 && from <= to );
     vec<String> pats;

     vec<int> from_digits, to_digits;
     int x = from;
     while(1)
     {    from_digits.insert( from_digits.begin( ), x % 10 );
          x = x / 10;
          if ( x == 0 ) break;    }
     x = to;
     while(1)
     {    to_digits.insert( to_digits.begin( ), x % 10 );
          x = x / 10;
          if ( x == 0 ) break;    }
     int f = from_digits.size( ), t = to_digits.size( );

     // If f = t, generate pattern for [from, to].

     if ( f == t )
     {    bool diff = false;

          // Let ti = the number of digits of agreement.  Put in prefix.

          String prefix;
          int ti;
          for ( ti = 0; ti < t; ti++ )
          {    if ( from_digits[ti] == to_digits[ti] ) 
                    prefix += Digit( to_digits[ti] );
               else break;    }

          if ( ti == t ) return prefix;

          if ( ti == t-1 ) 
               return prefix + DigitRange( from_digits[ti], to_digits[ti] );

          // Generate pattern to raise first nonagreeing digit.
          // Example: if from = 1234, to = 1256, generate pattern for 1234 - 1239.
          // Example: if from = 1072, to = 1358, generate pattern for 1072 - 1099.

          for ( int i = t-1; i > ti; i-- )
          {    String s = prefix + Digit(from_digits[ti]);
               for ( int j = ti + 1; j < i; j++ )
                    s += Digit( from_digits[j] );
               int begin = (i == t-1) ? from_digits[i] : from_digits[i] + 1;
               if ( begin > 9 ) continue;
               s += DigitRange( begin, 9 );
               for ( int j = i + 1; j < t; j++ )
                    s += DigitRange(0, 9);
               pats.push_back(s);    }

          // Generate pattern for the rest.
          // Example: if from = 1234, to = 1256, generate pattern for 1240 - 1256.
          // Example: if from = 1072, to = 1358, generate pattern for 1100 - 1358.

          for ( int i = ti; i < t; i++ )
          {    diff = false;
               String s = prefix;
               if ( i == ti && to_digits[i] > from_digits[i] + 1 ) diff = true;
               if ( i > ti && to_digits[i] > 0 ) diff = true;
               if ( !diff ) continue;
               for ( int j = ti; j < i; j++ )
                    s += Digit( to_digits[j] );
               int first = (i == ti) ? from_digits[i] + 1 : 0;
               int last = (i < t - 1) ? to_digits[i] - 1 : to_digits[i];
               if ( first < last ) s += DigitRange( first, last );
               else s += Digit(first);
               for ( int j = i + 1; j < t; j++ )
                    s += DigitRange(0, 9);
               if ( pats.size( ) == 0 || s != pats.back( ) ) pats.push_back(s);    }
          if ( !diff )
          {    String s = prefix;
               for ( int j = ti; j < t; j++ )
                    s += Digit( to_digits[j] );
               pats.push_back(s);    }    }

     // If f < t, generate pattern for [from, 10^f - 1].

     if ( f < t )
     {    for ( int i = f-1; i >= 0; i-- )
          {    String s;
               for ( int j = 0; j < i; j++ )
                    s += Digit( from_digits[j] );
               int begin = (i == f-1) ? from_digits[i] : from_digits[i] + 1;
               if ( begin == 10 ) continue;
               // int begin = from_digits[i];
               s += DigitRange( begin, 9 );
               for ( int j = i + 1; j < f; j++ )
                    s += DigitRange(0, 9);
               pats.push_back(s);    }    }

     // If f < t-1, generate pattern for [10^f, 10^(t-1) - 1].

     for ( int i = f+1; i <= t-1; i++ )
          pats.push_back( DigitRange(1, 9) + Dup( DigitRange(0, 9), i-1 ) );

     // If f < t, generate pattern for [10^(t-1), to].

     if ( f < t )
     {    bool diff = false;
          for ( int i = 0; i < t; i++ )
          {    diff = false;
               String s;
               if ( i == 0 && to_digits[i] > 1 ) diff = true;
               if ( i > 0 && to_digits[i] > 0 ) diff = true;
               if ( !diff ) continue;
               for ( int j = 0; j < i; j++ )
                    s += Digit( to_digits[j] );
               int first = (i == 0) ? 1 : 0;
               int last = (i < t - 1) ? to_digits[i] - 1 : to_digits[i];
               if ( first < last ) s += DigitRange( first, last );
               else s += Digit(first);
               for ( int j = i + 1; j < t; j++ )
                    s += DigitRange(0, 9);
               if ( pats.size( ) == 0 || s != pats.back( ) ) pats.push_back(s);    }
          if ( !diff )
          {    String s;
               for ( int j = 0; j < t; j++ )
                    s += Digit( to_digits[j] );
               pats.push_back(s);    }    }

     String answer;
     for ( unsigned int i = 0; i < pats.size( ); i++ )
     {    if ( i > 0 ) answer += "|";
          answer += "(" + pats[i] + ")";    }
     return answer;    }

void SubstituteIntegerRanges( String& expression )
{    restart:
     for ( unsigned int j = 0; j < expression.size( ); j++ )
     {    if ( expression[j] == '<' )
          {    unsigned int k;
               for ( k = j + 1; k < expression.size( ); k++ )
                    if ( expression[k] == '>' ) break;
               if ( k < expression.size( ) )
               {    String range = expression.substr(j + 1, k - j - 1);
                    if ( range.Contains( "," ) )
                    {    String from = range.Before( "," );
                         String to = range.After( "," );
                         unsigned int l;
                         for ( l = 0; l < from.size( ); l++ )
                              if ( !isdigit( from[l] ) ) break;
                         if ( l < from.size( ) ) break;
                         for ( l = 0; l < to.size( ); l++ )
                              if ( !isdigit( to[l] ) ) break;
                         if ( l < to.size( ) ) break;
                         if ( from.Int( ) > to.Int( ) ) break;
                         if ( from.Int( ) == 0 ) break;
                         String r 
                              = "(" + RangeToRegex( from.Int( ), to.Int( ) ) + ")";
                         expression.replace( j, k - j + 1, r.c_str( ) );
                         goto restart;    }    }    }    }    }

// Definition.  A *marked regular expression* is an extended regular expression,
// with the following modifications.  First, The bracket symbols { and } may not
// appear in the expression, except as part of a *mark*, which we now describe.

// A mark has the form 
//
//     {variable1=value1,variable2=value2,...,variablen=valuen}
//
// where each variablei is an alphanumeric string (allowing _), and each
// valuei is an alphanumeric string (allowing _), or "@".

// A mark must always be followed immediately by a right parenthesis ")".
// There must be a matching left parenthesis somewhere before the mark.  What
// is in between the two parentheses is called the *mark site*.

// UnpackMarkedRegex( String marked_exp, String& regexp, vec<markid>& marklist )
//
// takes a marked regular expression marked_exp, and separates the marks from
// the regular expression.  The marks are placed in marklist (which shows their
// positions in marked_exp: pos is the number of left parentheses which precede
// the left parenthesis of the mark site), and the regular expression itself is
// placed in regexp.

// Later: return an actual regular expression?

void UnpackMarkedRegexError( String marked_exp, int i )
{    cout << "At ";
     int head = 3;
     int j;
     int start = Max( 0, i - 20 ), stop = Min( (int) marked_exp.size( ), i + 20 );
     if ( start > 0 )
     {    cout << "...";
          head += 3;    }
     for ( j = start; j < stop; j++ )
          cout << marked_exp[j];
     if ( stop < (int) marked_exp.size( ) ) cout << "...";
     cout << "\n";
     for ( int k = 0; k < head + i - start; k++ )
          cout << " ";
     cout << "^\n\n";
     exit(1);    }

void UnpackMarkedRegex( String marked_exp, regex_t& regexp, vec<markid>& marklist )
{    
     // Look for quoted expressions.  These are temporarily converted into
     // %0, %1, etc.

     vec<String> quotes;
     for ( int i = 0; i < (int) marked_exp.size( ); i++ )
          if ( marked_exp[i] == '%' )
          {    cout << "Illegal character.\n";
               UnpackMarkedRegexError( marked_exp, i );    }
     for ( int i = 0; i < (int) marked_exp.size( ); i++ )
     {    if ( marked_exp[i] == '"' )
          {    int j;
               for ( j = i + 1; j < (int) marked_exp.size( ); j++ )
                    if ( marked_exp[j] == '"' )
                    {    String quote = marked_exp.substr( i+1, j-i-1 );
                         String head = marked_exp.substr( 0, i ) 
                              + "%" + ToString( quotes.size( ) );
                         marked_exp = head
                              + marked_exp.substr( j+1, marked_exp.size( ) - j - 1 );
                         quotes.push_back(quote);
                         i = head.size( ) - 1;
                         break;    }    }    }

     String exp = "";
     for ( int i = 0; i < (int) marked_exp.size( ); i++ )
     {    if ( marked_exp[i] != '{' )
          {    exp += marked_exp[i];
               continue;    }
          int j;
          for ( j = i + 1; j < (int) marked_exp.size( ); j++ )
               if ( marked_exp[j] == '}' ) break;
          if ( j == (int) marked_exp.size( ) )
          {    cout << "\nUnmatched left curly bracket ({) detected.\n";
               UnpackMarkedRegexError( marked_exp, i );    }
          if ( j == (int) marked_exp.size( ) - 1 || marked_exp[j+1] != ')' )
          {    cout << "\nRight curly bracket (}) not followed by right "
                    << "parenthesis.\n";
               UnpackMarkedRegexError( marked_exp, j );    }

          //     { at pos i     ...     } at pos j

          // Find the matching left paren.

          int k, paren_count = 1;
          for ( k = i - 1; k >= 0; k-- )
          {    if ( marked_exp[k] == '(' ) --paren_count;
               else if ( marked_exp[k] == ')' ) ++paren_count;
               if ( paren_count == 0 ) break;    }

          if ( k < 0 )
          {    cout << "\nMatching right parenthesis not found.\n";
               UnpackMarkedRegexError( marked_exp, j+1 );    }

          // Find number of left parens prior to the matching left paren.

          int left_paren_count = 0;
          for ( int l = 0; l < k; l++ )
               if ( marked_exp[l] == '(' ) ++left_paren_count;

          // Parse the mark.  First unpack, breaking at commas.

          vec<String> varvals;
          char delim = ',';
          String s = marked_exp.substr( i+1, j-i-1 );
          int count = 0;
          if ( s.size( ) > 0 )
          {    if ( s[ s.size( ) - 1 ] == delim )
               {    cout << "Unexpected use of comma.\n";
                    UnpackMarkedRegexError( marked_exp, i + count );    }
               String part;
               while( s.size( ) != 0 )
               {    if ( s[0] == delim )
                    {    s.erase( 0, 1 );
                         ++count;
                         varvals.push_back(part);
                         part = "";    }
                    else
                    {    part += s[0];
                         s.erase( 0, 1 );
                         ++count;    }    }
               varvals.push_back(part);    }

          for ( unsigned int l = 0; l < varvals.size( ); l++ )
          {    int var_start = i + 1;
               for ( unsigned int m = 0; m < l; m++ )
                    var_start += varvals[m].size( ) + 1;
               if ( !varvals[l].Contains( "=" ) )
               {    cout << "Missing equals sign.\n";
                    int place = i + 1;
                    for ( unsigned int m = 0; m < l+1; m++ )
                         place += varvals[l].size( ) + 1;
                    UnpackMarkedRegexError( marked_exp, 
                         var_start + varvals[l].size( ) );    }
               String var = varvals[l].Before( "=" ), val = varvals[l].After( "=" );
               if ( var.size( ) == 0 ) 
               {    cout << "Variable not found.\n";
                    UnpackMarkedRegexError( marked_exp, var_start );    }
               int val_start = var_start + var.size( ) + 1;
               if ( val.size( ) == 0 )
               {    cout << "Value not found.\n";
                    UnpackMarkedRegexError( marked_exp, val_start );    }
               for ( unsigned int m = 0; m < var.size( ); m++ )
                    if ( !isalnum( var[m] ) && var[m] != '_' )
                    {    cout << "Illegal character.\n";
                         UnpackMarkedRegexError( marked_exp, var_start + m );    }
               if ( val != "@" )
                    for ( unsigned int m = 0; m < val.size( ); m++ )
                         if ( !isalnum( val[m] ) && val[m] != '_' && val[m] != '%' )
                         {    cout << "Illegal character.\n";
                              UnpackMarkedRegexError( marked_exp, val_start + m );  }
               if ( val.size( ) > 0 && val[0] == '"' 
                    && val[ val.size( ) - 1 ] == '"' )
               {    val.erase( 0, 1 );
                    val.erase( val.size( ) - 1, 1 );    }
               marklist.push_back( markid( var, val, left_paren_count ) );    }

          i = j;    }

     String exp_plus = "^(" + exp + ")$";
     if ( regcomp( &regexp, exp_plus.c_str( ), REG_EXTENDED ) != 0 )
     {    cout << "Failed to translate " << exp << " into regular expression.\n";
          exit(1);    }

     for ( unsigned int i = 0; i < marklist.size( ); i++ )
     {    if ( marklist[i].val.Contains( "%", 0 ) )
          {    marklist[i].val 
                    = quotes[ marklist[i].val.After( "%" ).Int( ) ];    }    }    }

// Bool MatchMarkedRegex( String s, const regex_t& regexp, 
//                        const vec<markid>& marklist, int maxpos,
//                        vec<String>& variables, vec<String>& values )
//
// [where regexp and marklist come from UnpackMarkedRegex]
// 
// tries to match s to regexp, completely.  If it fails, return False.  Otherwise, 
// return True, and set variables and values, as follows.  For each matching mark 
// site, and each variablei=valuei part of it, append variablei to variables and 
// valuei to values, unless valuei = @, in which case append to values the substring
// of s which matched the mark site.

// A variable may not be set to a value more than once.

// Example.  Let m = (t(abcd{x=1,y=@})|(rs{z=cc}){w=5}) 
//
// Then m matches tabcd, yielding x=1, y=abcd, w=5
// and  m matches trs,   yielding z=cc, w=5.

Bool MatchMarkedRegex( String s, const regex_t& regexp, const vec<markid>& marklist, 
     int maxpos, vec<String>& variables, vec<String>& values )
{
     int nmatches = maxpos + 3;
     size_t nmatch(nmatches);
     vec<regmatch_t> pmatch(nmatches);

     if ( regexec( &regexp, s.c_str( ), nmatch, &pmatch[0], 0 ) == 0 )
     {    for ( int j = 0; j < nmatches; j++ )
          {    int start = pmatch[j].rm_so, stop = pmatch[j].rm_eo;
               if ( start >= 0 )
               {    for ( unsigned int i = 0; i < marklist.size( ); i++ )
                    {    const markid& m = marklist[i];
                         if ( m.pos == j - 2 )
                         {    variables.push_back( m.var );
                              if ( m.val != "@" ) values.push_back( m.val );
                              else values.push_back( s.substr( 
                                   start, stop - start ) );    }    }    }    }
          for ( unsigned int i = 0; i < variables.size( ); i++ )
               for ( unsigned int j = i + 1; j < variables.size( ); j++ )
                    if ( variables[i] == variables[j] )
                    {    cout << "The variable " << variables[i] << " was assigned "
                              << "a value twice, when matching to " << s << ".\n";
                         exit(1);    }

          return True;    }

     return False;    }

// WhiteSpaceFreeExceptQuotes: remove white space, except for white space between
// double quotes " ... ".

String WhiteSpaceFreeExceptQuotes( String s )
{    String answer;
     Bool in_quote = False;
     for ( unsigned int i = 0; i < s.size( ); i++ )
     {    if (in_quote) 
          {    answer += s[i];
               if ( s[i] == '"' ) in_quote = False;    }
          else
          {    if ( !isspace( s[i] ) ) answer += s[i];
               if ( s[i] == '"' ) in_quote = True;    }    }
     return answer;    }

// Perform macro substitutions in expression.  Then squeeze out white
// space and substitute for integer ranges.

void SubstituteMacrosEtc( String& expression, map<String,String>& macros)
{
  // We exploit the ordering of map keys here to prevent short macros from 
  // clobbering longer ones.  For example, say we have a macro "$alpha" that
  // is replaced by "[a-z]" and a macro "$alphanum" that is replaced by 
  // "[a-z0-9]".  If we traversed forward, any instance of "$alphanum" would
  // be turned into "[a-z]num" by $alpha before $alphanum ever got substituted.
  // By traversing the map in reverse order, we guarantee that these longer-
  // named macros will be checked first. (The $ prevents something like $num
  // causing a problem for $alphanum.)

  map<String,String>::reverse_iterator macro_iter;
  for ( macro_iter = macros.rbegin(); 
	macro_iter != macros.rend(); 
	macro_iter++ ) {
    const String& macro = macro_iter->first;
    int pos;
    while( (pos = expression.Position( macro )) >= 0 ) {
      const String& replacement = macro_iter->second;
      expression.replace( (unsigned int) pos, 
			  macro.size(), 
			  replacement.c_str( ) );    
    }    
  }
  expression = WhiteSpaceFreeExceptQuotes(expression);
  SubstituteIntegerRanges(expression);
}

void InvalidStatement( String& s )
{    cout << "The statement\n" << s << "\nin reads.config is invalid.\n";
     exit(1);    }

class string_int {
     public:
     String s;
     int i;
     string_int( const String& s_arg, int i_arg )
       : s(s_arg), i(i_arg)
     { }

     friend bool operator<( const string_int& si1, const string_int& si2 )
     {    return si1.s < si2.s;    }
};

}

// TranslateReadName takes as input a list of statements, as would be generated by
// ParseConfigFile, and a read name.  The first time TranslateReadName is called, it
// parses the statements.  This creates matched vectors global_vars and global_vals,
// as well as internal structures which provide a reduced representation of the
// statements.  Each time TranslateReadName is called, it evaluates the read name,
// by comparing with the internal structures.  It then generates:
//
//    * a boolean "unmatched", which tells if the read did not match any include
//      statement (presumably an error condition);
//    * a boolean "excluded", which tells if the read is to be ignored;
//    * matched vectors "variables" and "values".

void TranslateReadName( vec<String>& statements, String read_name,
     vec<String>& variables, vec<String>& values, Bool& excluded, Bool& unmatched,
     vec<String>& global_vars, vec<String>& global_vals )
{    
     static bool first_call(true);
     static vec<regex_t> includes, excludes, unpairs;
     static vec< vec<markid> > marklists;
     static vec<int> maxpos;
     static vec<String> exclude_heads;
     static vec< pair< vec<String>, vec<int> > > exheads;

     if (first_call)
     {    
          first_call = false;
	  
	  // we use a map here because it keeps its keys in order
          map<String,String> macros;
	  // this is exploited in SubstituteMacrosEtc() above

          for ( unsigned int i = 0; i < statements.size( ); i++ )
          {    String st = statements[i];
               DeleteTrailingWhiteSpace(st);
               DeleteLeadingWhiteSpace(st);

               // Check for and delete trailing semicolon.
               // Convert white space to blanks.

               ForceAssert( st.Contains( ";", -1 ) );
               st = st.Before( ";" );
               for ( unsigned int j = 0; j < st.size( ); j++ )
                    if ( isspace( st[j] ) ) st[j] = ' ';

               if ( st.Contains( "fix", 0 ) )
               {    st = st.After( "fix" );
                    if ( !st.Contains( "=" ) )
                    {    cout << "I can't make sense out of " <<
                              statements[i] << ".\n";
                         exit(1);    }
                    String variable = WhiteSpaceFree( st.Before( "=" ) );
                    String value = WhiteSpaceFree( st.After( "=" ) );
                    global_vars.push_back(variable);
                    global_vals.push_back(value);    }

               else if ( st.Contains( "exclude ", 0 ) )
               {    String stx = st.After( "exclude " );
                    SubstituteMacrosEtc( stx, macros );

                    // For each excluded class of reads, we extract its head:
                    // the part of the regular expression which is just an ordinary
                    // string, which has to match if the entire regular expression
                    // is going to match.  By sorting these heads (later), we are 
                    // able to efficiently perform the comparison between large 
                    // numbers of read names and large numbers of excludes.

                    String head;
                    for ( unsigned int j = 0; j < stx.size( ); j++ )
                    {    if ( !isalnum( stx[j] ) ) break;
                         head += stx[j];    }
                    exclude_heads.push_back(head);

                    stx = "^(" + stx + ")$";
                    regex_t r;
                    if ( regcomp( &r, stx.c_str( ), REG_NOSUB|REG_EXTENDED ) != 0 )
                    {    cout << "Unable to compile regular expression in:\n";
                         cout << st << "\n";
                         InvalidStatement( statements[i] );    }
                    excludes.push_back(r);    }

               else if ( st.Contains( "include ", 0 ) )
               {    String stx = st.After( "include " );
                    SubstituteMacrosEtc( stx, macros );
                    stx = "(" + stx + ")";
                    regex_t regexp;
                    vec<markid> marklist;
                    UnpackMarkedRegex( stx, regexp, marklist );
                    includes.push_back(regexp);
                    marklists.push_back(marklist);    
                    int mp = 0;
                    for ( unsigned int j = 0; j < marklist.size( ); j++ )
                         mp = Max( mp, marklist[j].pos );
                    maxpos.push_back(mp);    }

               else if ( st.Contains( "unpair ", 0 ) )
               {    String stx = st.After( "unpair " );
                    SubstituteMacrosEtc( stx, macros );
                    stx = "^(" + stx + ")$";
                    regex_t r;
                    if ( regcomp( &r, stx.c_str( ), REG_NOSUB|REG_EXTENDED ) != 0 )
                    {    cout << "Unable to compile regular expression in:\n";
                         cout << st << "\n";
                         InvalidStatement( statements[i] );    }
                    unpairs.push_back(r);    }

               else if ( st.Contains( "define ", 0 ) )
               {    st = st.After( "define " );
                    int eq = st.Position( "=" );
                    if ( eq <= 0 ) InvalidStatement( statements[i] );
                    String macro_name = WhiteSpaceFree( st.Before( "=" ) );
                    String expression = st.After( "=" ) + " ";
                    SubstituteMacrosEtc( expression, macros );
		    expression = "(" + expression + ")";
                    macros.insert(make_pair("$"+macro_name,expression)); }
               else InvalidStatement( statements[i] );    }

          // Now we set up the structure for efficiently determining (later) which
          // of the exclude_head's are suffixes of a given read name.  Although
          // suffix trees would provide a cleaner solution, I doubt that much time
          // would be saved.
     
          // exhead[i] = list of pairs (j, k), where j is an exclude_head of length
          // i, and k is the index in excludes of the corresponding regular 
          // expression.  The list exhead[i] is sorted lexicographically by j.

          for ( unsigned int sz = 0; ; sz++ )
          {    vec<string_int> exhead;
               Bool done = True;
               for ( unsigned int j = 0; j < exclude_heads.size( ); j++ )
               {    if ( exclude_heads[j].size( ) >= sz ) done = False;
                    if ( exclude_heads[j].size( ) == sz ) 
                         exhead.push_back( string_int( exclude_heads[j], j ) );    }
               if (done) break;
               sort( exhead.begin( ), exhead.end( ) );
               vec<String> ex_s;
               vec<int> ex_i;
               for ( unsigned int j = 0; j < exhead.size( ); j++ )
               {    ex_s.push_back( exhead[j].s );
                    ex_i.push_back( exhead[j].i );    }
               exheads.push_back( make_pair( ex_s, ex_i ) );    }    }

     // Test to see if a given read is excluded.

     unmatched = False;
     {    size_t nmatch(0);
          regmatch_t* pmatch(0);
          for ( unsigned int i = 0; i < exheads.size( ); i++ )
          {    if ( read_name.size( ) < i ) break;
               vec<String>& suf = exheads[i].first;
               if ( suf.size( ) == 0 ) continue;
               vec<int>& loc = exheads[i].second;
               String read_head = read_name.substr( 0, i );
               int pos = BinPosition( suf, read_head );
               if ( pos >= 0 )
               {    int pos0, pos1;
                    for ( pos0 = pos - 1; pos0 >= 0; pos0-- )
                         if ( suf[pos0] != read_head ) break;
                    ++pos0;
                    for ( pos1 = pos + 1; pos1 < (int) suf.size( ); pos1++ )
                         if ( suf[pos1] != read_head ) break;
                    for ( int j = pos0; j < pos1; j++ )
                    {    if ( regexec( &excludes[ loc[j] ], read_name.c_str( ), 
                              nmatch, &pmatch[0], 0 ) == 0 )
                         {    excluded = True;
                              return;    }    }    }    }    }

     Bool unpaired = False;
     {    size_t nmatch(0);
          regmatch_t* pmatch(0);
          for ( unsigned int i = 0; i < unpairs.size( ); i++ )
          {    if ( regexec( &unpairs[i], read_name.c_str( ), nmatch, 
                    &pmatch[0], 0 ) == 0 )
               {    unpaired = True;
                    break;    }    }    }
     excluded = False;
     variables.resize(0);
     values.resize(0);
     for ( unsigned int i = 0; i < includes.size( ); i++ )
     {    if ( MatchMarkedRegex( 
               read_name, includes[i], marklists[i], maxpos[i], variables, values ) )
          {    
               if (unpaired)
               {    for ( unsigned int j = 0; j < variables.size( ); j++ )
                         if ( variables[j] == "Orientation" )
                         {    values[j] = "none";
                              break;    }    }

               // for ( unsigned int i = 0; i < variables.size( ); i++ )
               //      cout << variables[i] << " = " << values[i] << "\n";
               return;    }    }

     cout << "No match found for read name " << read_name << "." << endl;
     unmatched = True;    }

// ParseConfigFile: Read in the configuration file, doing minimal processing,
// yielding as output a list of statements.

void ParseConfigFile( istream& in, vec<String>& statements )
{    String splus;
     while(1)
     {    String s;
          getline( in, s );
          if ( !in ) break;
          if ( s.size( ) >= 2 && s[0] == '/' && s[1] == '/' ) continue;
          for ( unsigned int i = 0; i < s.size( ); i++ )
          {    splus += s[i];
               if ( s[i] == ';' )
               {    statements.push_back(splus);
                    splus.resize(0);    }    }
          splus += " ";    }
     for ( unsigned int i = 0; i < splus.size( ); i++ )
     {    if ( !isspace( splus[i] ) )
          {    cout << "The read config file has some stuff at its end "
                    << "which is not followed by a semicolon:\n" << splus << "\n";
               exit(1);    }    }    }
