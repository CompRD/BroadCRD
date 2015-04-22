// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 



// PaddedSmithWat( S, T, offset, bandwidth )
//
// Let S and T be PaddedBasevectors.  Find the best (lowest) scoring alignment of S
// with T, relative to the following rules:
//
// (a) a mismatch scores +1.0;
// (b) each blank in a gap scores +1.5;
// (c) gaps on either end of either read do not count.

// Only search for alignments with the given offset, plus or minus the given
// bandwidth.  The meaning of "offset" is explained by:
//
// example: offset = 2
// S -----------------
// T   ---------------------

// example: offset = -2
// S   ---------------------
// T -----------------

// Return the score of the alignment, and by reference, the alignment itself.

#undef VERBOSE_SW

#include "Alignment.h"
#include "Basevector.h"
#include "PaddedBasevector.h"
#include "math/Functions.h"
#include "PackAlign.h"
#include "ShortVector.h"
#include "PaddedSmithWat.h"
#include "Vec.h"

#include "system/System.h"


float PaddedSmithWat( const vec<char>& s, const vec<char>& t, 
                      int offset, int bandwidth, alignment& a, 
                      const int mismatch_penalty, const int gap_penalty,
                      const float divider )
{
     /*
     if ( offset > (int) S.PaddedLength( ) || offset < - (int) T.PaddedLength( ) )
     {    cout << "Warning: PaddedSmithWat passed nonsense arguments:" << endl;
          PRINT2( S.PaddedLength( ), T.PaddedLength( ) );
          PRINT2( offset, bandwidth );    }
     */

     // If arguments are nonsense, reform them.

     if ( offset > (int) s.size() ) offset = (int) s.size();
     if ( offset < - (int) t.size() ) offset = - (int) t.size();

     int left = offset - bandwidth, right = offset + bandwidth;

     const int pad_gap_penalty = 0;

     const int Infinity = 1000000000;

     int n = (int) s.size(), N = (int) t.size();

     int best_score = Infinity;
     int best_j = 0, best_i = 0;
     avector<int> x(n+1);
     int istart = 0, istop = 0;

     for ( int i = 0; i <= n; i++ )
     {  
       x(i) = 0;
       if ( !( left <= i && i <= right ) ) 
	 x(i) = Infinity;    }

     int jstart = Max( 0, -right-1 ), jstop = Min( N-1, n-left );
     
     // [i,j] to be stored at from[ j - jstart ][ i - left - j + 1 ]
     vec< vec<unsigned char> > from( jstop - jstart + 1, 
				     vec<unsigned char>( right - left + 3 ,'n') );

#ifdef VERBOSE_SW
     cout << " left: " << left 
	  << " right: " << right
	  << " jstart: " << jstart
	  << " jstop: " << jstop << endl;

     istart = Max( 0, left + jstart - 1 );
     istop = Min( n, right + jstop - 1 + 1 );
     cout << "  ";
     for ( int i = istart; i < istop; ++i )
       cout << s[i];
     cout << endl;
#endif

     for ( int j = jstart; j <= jstop; ++j )
     {    
#ifdef VERBOSE_SW
        cout << t[j] << " ";
#endif
       x(0) = 0;
       if ( !( left <= -j && -j <= right ) ) 
	 x(0) = Infinity;
       int lastx = 0;
       if ( !( left <= -(j-1) && -(j-1) <= right ) ) 
	 lastx = Infinity;
       const char* sp = &s[0];

       istart = Max( 0, left + j - 1 );
       istop = Min( n - 1, right + j + 1 );
       int* xp = x.x + istart;

       if ( istart > 0 )
       {
	 lastx = *xp;
	 *xp = Infinity;    
       }

#ifdef VERBOSE_SW
       for ( int i = 0; i < istart; ++i)
	 cout << " ";
       
       cout << endl;
       cout << "  j: " << setw(2) << j 
	    << "  istart: " << istart 
	    << "  istop: " << istop;

       if ( lastx == Infinity )
	 cout << "  lastx: I";
       else
	 cout << "  lastx: " << lastx;

       cout << "  x:";
       for ( int ii = 0; ii < x.length; ++ii )
	 if ( x(ii) != Infinity )
	   cout << setw(3) << x(ii);
	 else
	   cout << "  I";
#endif
  
       #define SWMIN(I)                        \
   /* \
       cout << endl;                           \
       if ( a < Infinity )                     \
	    cout << " a: " << setw(2) << a;    \
       else                                    \
	    cout << " a:  I";                  \
       if ( b < Infinity )                     \
            cout << " b: " << setw(2) << b;    \
       else                                    \
	    cout << " b:  I";                  \
       if ( c < Infinity )                     \
	    cout << " c: " << setw(2) << c;    \
       else                                    \
            cout << " c:  I";                  \
    */ \
         lastx = *xp;                          \
	 if ( a <= b )                         \
	 {				       \
	   if ( a <= c ) {		       \
	     from[j-jstart][I-left-j+1] = 'a'; \
             *xp = a;                          \
           }                                   \
	   else {			       \
	     from[j-jstart][I-left-j+1] = 'c'; \
	     *xp = c;                          \
	   }                                   \
	 }                   		       \
	 else                                  \
	 {				       \
	   if ( b <= c ) { 		       \
	     from[j-jstart][I-left-j+1] = 'b'; \
             *xp = b;                          \
           }                                   \
	   else				       \
           {                                   \
	     from[j-jstart][I-left-j+1] = 'c'; \
             *xp = c;                          \
           }                                   \
	 }                                     \
//        cout << from[j-jstart][I-left-j+1];

   /* \
       cout << "  x:";                             \
       for ( int ii = 0; ii < x.length; ++ii )     \
	 if ( x(ii) != Infinity )                  \
           if ( x.x + ii == xp )                   \
	     cout << ">" << setw(2) << x(ii);      \
           else                                    \
	     cout << setw(3) << x(ii);             \
	 else                                      \
	   cout << "  I";                          \
       cout << " from: ";                          \
       for ( int ii = 0; ii < from[j-jstart].size(); ++ii ) \
	 cout << from[j-jstart][ii] << " ";        \
   */ 


       #define SWCORE(J)                                             \
	 {                                                           \
	   if ( istart <= istop )                                    \
	   {                                                         \
	     int a = lastx + mismatch_penalty * (sp[istart] != J);   \
	     int b = *xp + gap_penalty;                              \
	     ++xp;                                                   \
	     int c = *xp + gap_penalty;                              \
	     if ( !( left <= (istart-1) - (j-1) &&                   \
		     (istart-1) - (j-1) <= right ) )                 \
	       a = Infinity;                                         \
	     if ( !( istart - (j-1) <= right ) )                     \
	       b = Infinity;                                         \
	     if ( !( left <= (istart-1) - j ) )                      \
	       c = Infinity;                                         \
	     SWMIN(istart);                                          \
	   }                                                         \
	   if ( istart + 1 <= istop )                                \
	   {                                                         \
	     int a = lastx + mismatch_penalty * (sp[istart+1] != J); \
	     int b = *xp + gap_penalty;                              \
	     ++xp;                                                   \
	     int c = *xp + gap_penalty;                              \
	     if ( !( left <= istart - (j-1) &&                       \
		     istart - (j-1) <= right ) )                     \
	       a = Infinity;                                         \
	     if ( !( istart+1 - (j-1) <= right ) )                   \
	       b = Infinity;                                         \
	     if ( !( left <= istart - j ) )                          \
	       c = Infinity;                                         \
	     SWMIN(istart+1);                                        \
	   }                                                         \
	   for ( int i = istart + 2; i <= istop - 2; i++ )           \
	   {                                                         \
	     int a = lastx + mismatch_penalty * (sp[i] != J);        \
	     int b = *xp + gap_penalty;                              \
	     ++xp;                                                   \
	     int c = *xp + gap_penalty;                              \
	     SWMIN(i);                                               \
	   }                                                         \
	   if ( istop - 1 >= istart + 2 )                            \
	   {                                                         \
	     int a = lastx + mismatch_penalty * (sp[istop-1] != J);  \
	     int b = *xp + gap_penalty;                              \
	     ++xp;                                                   \
	     int c = *xp + gap_penalty;                              \
	     if ( !( istop - j <= right ) )                          \
	       b = Infinity;                                         \
	     SWMIN(istop-1);                                         \
	   }                                                         \
	   if ( istop >= istart + 2 )                                \
	   {                                                         \
	     int a = lastx + mismatch_penalty * (sp[istop] != J);    \
	     int b = *xp + gap_penalty;                              \
	     ++xp;                                                   \
	     int c = *xp + gap_penalty;                              \
	     if ( !( istop - j <= right ) )                          \
	       a = Infinity;                                         \
	     if ( !( istop - (j-1) <= right ) )                      \
	       b = Infinity;                                         \
	     SWMIN(istop);                                           \
	   }                                                         \
	 }

       #define PADSWCORE(J)                                          \
	 {                                                           \
	   if ( istart <= istop )                                    \
	   {                                                         \
	     int a = lastx + gap_penalty * (sp[istart] != J);        \
	     int b = *xp + gap_penalty;                              \
	     ++xp;                                                   \
	     int c = *xp + pad_gap_penalty;                          \
	     if ( !( left <= (istart-1) - (j-1) &&                   \
		     (istart-1) - (j-1) <= right ) )                 \
	       a = Infinity;                                         \
	     if ( !( istart - (j-1) <= right ) )                     \
	       b = Infinity;                                         \
	     if ( !( left <= (istart-1) - j ) )                      \
	       c = Infinity;                                         \
	     SWMIN(istart);                                          \
	   }                                                         \
	   if ( istart + 1 <= istop )                                \
	   {                                                         \
	     int a = lastx + gap_penalty * (sp[istart+1] != J);      \
	     int b = *xp + gap_penalty;                              \
	     ++xp;                                                   \
	     int c = *xp + pad_gap_penalty;                          \
	     if ( !( left <= istart - (j-1) &&                       \
		     istart - (j-1) <= right ) )                     \
	       a = Infinity;                                         \
	     if ( !( istart+1 - (j-1) <= right ) )                   \
	       b = Infinity;                                         \
	     if ( !( left <= istart - j ) )                          \
	       c = Infinity;                                         \
	     SWMIN(istart+1);                                        \
	   }                                                         \
	   for ( int i = istart + 2; i <= istop - 2; i++ )           \
	   {                                                         \
	     int a = lastx + gap_penalty * (sp[istart+1] != J);      \
	     int b = *xp + gap_penalty;                              \
	     ++xp;                                                   \
	     int c = *xp + pad_gap_penalty;                          \
	     SWMIN(i);                                               \
	   }                                                         \
	   if ( istop - 1 >= istart + 2 )                            \
	   {                                                         \
	     int a = lastx + gap_penalty * (sp[istart+1] != J);      \
	     int b = *xp + gap_penalty;                              \
	     ++xp;                                                   \
	     int c = *xp + pad_gap_penalty;                          \
	     if ( !( istop - j <= right ) )                          \
	       b = Infinity;                                         \
	     SWMIN(istop-1);                                         \
	   }                                                         \
	   if ( istop >= istart + 2 )                                \
	   {                                                         \
	     int a = lastx + gap_penalty * (sp[istart+1] != J);      \
	     int b = *xp + gap_penalty;                              \
	     ++xp;                                                   \
	     int c = *xp + pad_gap_penalty;                          \
	     if ( !( istop - j <= right ) )                          \
	       a = Infinity;                                         \
	     if ( !( istop - (j-1) <= right ) )                      \
	       b = Infinity;                                         \
	     SWMIN(istop);                                           \
	   }                                                         \
	 }

       switch ( t[j] ) {
       case 'A':
	 SWCORE('A');
	 break;
       case 'C':
	 SWCORE('C');
	 break;
       case 'G':
	 SWCORE('G');
	 break;
       case 'T':
	 SWCORE('T');
	 break;
       case '*':
	 PADSWCORE('*');
	 break;
       default:
	 FatalErr( "t[j] = " << t[j] << ", not one of [ACGT*]" << endl );
       }

       if ( istop < n - 1 )
       {    
	 ++xp;
	 *xp = Infinity;    
       }
       
       if ( istop == n - 1 ) 
       {
	 if ( x(n) < best_score ) 
	 {
	   best_j = j;
	   best_i = n-1;
	   best_score = x(n);  
	 }
       }

#ifdef VERBOSE_SW
       cout << endl;

       cout << " -> ";
       cout << "x:";
       for ( int ii = 0; ii < x.length; ++ii )
	 if ( x(ii) != Infinity )
	   cout << setw(3) << x(ii);
	 else
	   cout << "  I";

       cout << " from: ";
       for ( int ii = 0; ii < from[j-jstart].size(); ++ii )
	 cout << from[j-jstart][ii] << " ";

       cout << " best_i: " << best_i;
       cout << " best_j: " << best_i;
       cout << " best_score: " << best_score;
       
       cout << endl;
#endif
     }
     
     for ( int i = istart; i <= istop; i++ )
     {
       if ( x(i) < best_score )
       {    
	 // cout << "type 2 set\n"; // XXX
	 // PRINT3( i, istart, istop ); // XXX
	 // PRINT( x[i] ); // XXX
	 best_j = jstop;
	 best_i = i - 1;     // NEW NEW NEW NEW NEW NEW!!!!!
	 // best_i = i;    
	 // PRINT2( best_i, best_j ); // XXX
       }
       best_score = Min( x(i), best_score );    
     }
     
     ForceAssertLt( best_i, (int) s.size() );
     ForceAssertLt( best_j, (int) t.size() );
     
     int j = best_j, i = best_i;
     int lcount = 0, g1count = 0, g2count = 0, last_length = 0;
     avector<int> gaps(0), lengths(0);
     
     // int a_count = 0; // XXX
     while(1)
     {    
       if ( j < 0 || i < 0 || j-jstart < 0 || i-left-j+1 < 0 ||
	    i-left-j+1 >= right-left+3 ) 
	 break;

       /*
	 char letter = from[j-jstart][i-left-j+1]; // XXX
	 if ( letter == 'a' ) ++a_count; // XXX
	 else  // XXX
	 {
	 if ( a_count > 0 ) cout << " a^" << a_count << " "; // XXX
	 a_count = 0; // XXX
	 cout << letter;    } // XXX
       */
       
       if ( from[j-jstart][i-left-j+1] == 'a' )
       {
	 if ( g1count > 0 )
	 {
	   gaps.Prepend( g1count );
	   lengths.Prepend( last_length );
	   g1count = 0;    }
	 if ( g2count > 0 )
	 {
	   gaps.Prepend( -g2count );
	   lengths.Prepend( last_length );
	   g2count = 0;    
	 }
	 if ( i == 0 || j == 0 ) 
	   break; // NEW NEW NEW !!!
	 ++lcount;
	 --i;
	 --j;
       }
       else if ( from[j-jstart][i-left-j+1] == 'b' )  // gap on long sequence
       {
	 if ( lcount > 0 )
	 {
	   last_length = lcount;
	   lcount = 0;
	 }
	 //ForceAssertEq( g1count, 0 );
	 ++g2count;
	 --i;
       }
       else if ( from[j-jstart][i-left-j+1] == 'c' )  // gap on short sequence
       {
	 if ( lcount > 0 )
	 {
	   last_length = lcount;
	   lcount = 0;    
	 }
	 //ForceAssertEq( g2count, 0 );
	 ++g1count;
	 if ( j == 0 ) 
	   break; // NEW NEW NEW !!!
	 --j;
       }
       else
	 break;
     }
     // if ( i < 1 ) break;
     // ForceAssert( j >= 0 );    }
     // if ( a_count > 0 ) cout << " a^" << a_count << " "; // XXX
     // cout << "\n"; // XXX
     
     // ForceAssert( g1count == 0 );
     // ForceAssert( g2count == 0 );
     // gaps.Prepend(0);

     if ( g1count != 0 ) 
       gaps.Prepend( g1count );
     else if ( g2count != 0 ) 
       gaps.Prepend( -g2count );
     else
       gaps.Prepend(0);

     lengths.Prepend( lcount + 1 );
     // lengths.push_back(lcount);

     int pos1 = i, pos2 = j;
     if ( pos1 < 0 || pos2 < 0 )
     { 
       ++pos1;
       ++pos2;    
     }
     int errors = best_score;
     a = alignment( pos1, pos2, errors, gaps, lengths );
 
     if ( pos1 < 0 || pos2 < 0 || 
	  a.Pos1( ) > (int) s.size() ||
          a.Pos2( ) > (int) t.size() )
     {    
          // See no evil...

          /*
          cout << "Warning: SmithWatBandedA produced a garbage alignment, "
               << "presumably because of a bug.\n";
          cout << "Replacing the garbage alignment by an alignment which has "
               << "one-base overlap.\n";
          PRINT3( pos1, a.Pos1( ), S.Length( ) );
          PRINT3( pos2, a.Pos2( ), T.Length( ) );
          int pos1x, pos2x, errors;
          avector<int> gapsx, lengthsx;
          a.Unpack( pos1x, pos2x, errors, gapsx, lengthsx );
          cout << "gaps/lengths:";
          for ( unsigned int i = 0; i < lengthsx.length; i++ )
               cout << " " << gapsx(i) << "/" << lengthsx(i);
          cout << endl;
          PRINT2(offset, bandwidth);
          if ( S.Length( ) < 1000 ) PRINT( S.ToString( ) );
          if ( T.Length( ) < 1000 ) PRINT( T.ToString( ) );
          */

       ForceAssertGt( (int) s.size(), 0 );
       ForceAssertGt( (int) t.size(), 0 );
       Bool base_matches = ( s.front() == t.back() );

       avector<int> gapsz(1), lengthsz(1);
       gapsz(0) = 0;
       lengthsz(0) = 1;
       a.Set( 0, t.size() - 1, (base_matches ? 0 : 1), gapsz, lengthsz );

       if (base_matches) 
	 return 0;
       else
	 return float(mismatch_penalty)/divider;    
     }

     return float(best_score)/divider;
}
