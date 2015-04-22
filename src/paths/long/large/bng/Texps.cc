///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Go through the DISCOVAR lines, generating a data structure texps, that
// we will use for mapping the lines to the BNG maps.
//
// This code also prints out a human-readable representation of the 
// translation of each line to "map space".
//
// Example:
//
// 40630:
// {{11426,2888,7665,382,540,2717,1527,375,2627,723}}
// 
// This is the simplest case: the line is unambiguously expanded as a list
// of restrictions-site-free intervals.  Note that the first and last invervals
// are incomplete and thus special.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "CoreTools.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/Lines.h"
#include "paths/long/large/bng/Texps.h"

void Texps( const HyperBasevectorX& hb, const vec<vec<vec<vec<int>>>>& lines,
     const String& cut, const int min_line, const int max_ignored_indel,
     vec<vec<vec<int>>>& texps, const Bool verbose, ostream& out )
{
     String rcut;
     StringReverseComplement( cut, rcut );
     int K = hb.K( );
     vec<int> lens;
     GetLineLengths( hb, lines, lens );

     cout << Date( ) << ": start traversing lines" << endl;
     int N;
     for ( N = 0; N < lines.isize( ); N++ )
          if ( lens[N] < min_line ) break;
     texps.clear( );
     texps.resize(N);
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int li = 0; li < N; li++ )
     {    const vec<vec<vec<int>>>& L = lines[li];
          vec<vec<vec<int>>> TAGS;

          // Go through the given line L.

          for ( int j = 0; j < L.isize( ); j++ )
          {    
               // Go through the given part of the line.  Its 'tags' are the
               // unique sequences of integers representing restriction-site-free
               // intervals.

               vec<vec<int>> tags;
               for ( int l = 0; l < L[j].isize( ); l++ )
               {    const vec<int>& p = L[j][l];
                    vec<int> tag;
                    if ( p.nonempty( ) )
                    {    String s = hb.Cat(p).ToString( );
                         if ( j != L.isize( ) - 1 )
                              s.resize( s.isize( ) - K + 1 + cut.isize( ) );
                         int count = 0;
                         for ( int j = 0; j < s.isize( ); j++ )
                         {    if ( j >= cut.isize( ) 
                                   && s.Contains( cut, j - cut.isize( ) ) )
                              {    tag.push_back(count);
                                   count = 0;    }
                              else if ( s.Contains( rcut, j ) )
                              {    tag.push_back(count);
                                   count = 0;    }
                              else count++;    }
                         tag.push_back(count);    }
                    tags.push_back(tag);    }
               UniqueSort(tags);

               // Ignore gaps.

               if ( tags.solo( ) && tags[0].empty( ) ) continue;

               // Handle the 'matrix' case.  This is where all the tag sequences
               // agree except for small indels.  In that case we take the average.

               Bool all_same = True;
               for ( int j = 1; j < tags.isize( ); j++ )
                    if ( tags[j].size( ) != tags[j-1].size( ) ) all_same = False;
               if (all_same) // transpose
               {    vec<vec<int>> t( tags[0].size( ) );
                    for ( int j = 0; j < tags.isize( ); j++ )
                    for ( int k = 0; k < tags[j].isize( ); k++ )
                         t[k].push_back( tags[j][k] );
                    vec<int> T( t.isize( ) );
                    Bool ok = True;
                    for ( int l = 0; l < t.isize( ); l++ )
                    {    Sort( t[l] );
                         // Squish small indels.
                         if ( t[l].back( ) - t[l].front( ) <= max_ignored_indel )
                         {    int m = ( t[l].back( ) + t[l].front( ) ) / 2;
                              T[l] = {m};    }
                         else ok = False;    }
                    if (ok)
                    {    TAGS.push_back( {T} );
                         continue;    }    }
          
               /*
               Bool all_solo = True;
               for ( int j = 0; j < tags.isize( ); j++ )
                    if ( !tags[j].solo( ) ) all_solo = False;
               if ( all_solo && tags.nonempty( )
                    && tags.back( )[0] - tags.front( )[0] <= max_ignored_indel )
               {    vec<vec<int>> t = {{ ( tags.front( ) + tags.back( ) ) / 2 }};
                    tags = t;    }
               */

               // Save.

               TAGS.push_back(tags);    }

          // Simplify.

          vec<Bool> to_delete( TAGS.size( ), False );
          for ( int j = 1; j < TAGS.isize( ); j++ )
          {    if ( TAGS[j-1].solo( ) && TAGS[j-1][0].solo( )
                    && TAGS[j].solo( ) && TAGS[j][0].solo( ) )
               {    to_delete[j-1] = True;
                    TAGS[j][0][0] += TAGS[j-1][0][0] - cut.isize( );    }    }
          EraseIf( TAGS, to_delete );

          vec<Bool> to_delete2( TAGS.size( ), False );
          for ( int j = 1; j < TAGS.isize( ); j++ )
          {    if ( TAGS[j].solo( ) && TAGS[j-1].solo( )
                    && TAGS[j][0].nonempty( ) && TAGS[j-1][0].nonempty( ) )
               {    to_delete2[j-1] = True;
                    vec<int> t = TAGS[j-1][0];
                    t.back( ) += TAGS[j][0][0] - cut.isize( );
                    for ( int m = 1; m < TAGS[j][0].isize( ); m++ )
                         t.push_back( TAGS[j][0][m] );
                    TAGS[j][0] = t;
                    /*
                    vec<Bool> del( TAGS[j][0].size( ), False );
                    for ( int l = 0; l < TAGS[j][0].isize( ); l++ )
                         if ( TAGS[j][0][l] < 0 ) del[l] = True;
                    EraseIf( TAGS[j][0], del );    
                    */
                         }    }
          EraseIf( TAGS, to_delete2 );

          // Get rid of -1 values.  Patching over a bug.

          for ( int i = 0; i < TAGS.isize( ); i++ )
          for ( int l = 0; l < TAGS[i].isize( ); l++ )
          for ( int m = 0; m < TAGS[i][l].isize( ); m++ )
               if ( TAGS[i][l][m] < 0 ) TAGS[i][l][m] = 0;

          // Print.

          if (verbose)
          {    
               #pragma omp critical
               {    out << "\n" << li << ":\n";
                    for ( int j = 0; j < TAGS.isize( ); j++ )
                    {    const vec<vec<int>>& tags = TAGS[j];
                         out << "{";
                         for ( int l = 0; l < tags.isize( ); l++ )
                         {    if ( l > 0 ) out << ",";
                              out << "{" << printSeq( tags[l] ) << "}";    }
                         out << "}\n";    }    }    }

          // Expand TAGS.

          vec<vec<int>> texp(1);
          for ( int j = 0; j < TAGS.isize( ); j++ )
          {    vec<vec<int>> texp_new;
               for ( int i = 0; i < texp.isize( ); i++ )
               for ( int l = 0; l < TAGS[j].isize( ); l++ )
               {    vec<int> t = texp[i];
                    t.append( TAGS[j][l] );
                    texp_new.push_back(t);    }
               texp = texp_new;    }
          // if (verbose) PRINT( texp.size( ) );

          // Delete end values.  They are incomplete intervals!

          for ( int j = 0; j < texp.isize( ); j++ )
          {    vec<int> x;
               for ( int l = 1; l < texp[j].isize( ) - 1; l++ )
                    x.push_back( texp[j][l] );
               texp[j] = x;    }
          /*
          #pragma omp critical
          {    cout << "texps[" << li << "] = {";
               for ( int i = 0; i < texp.isize( ); i++ )
               {    if ( i > 0 ) cout << ",";
                    cout << "{" << printSeq( texp[i] ) << "}";    }
               cout << "}" << endl;    }
          */
          texps[li] = texp;    }

     cout << "\n" << Date( ) << ": done" << endl;    }
