///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "Bitvector.h"
#include "CoreTools.h"
#include "math/HoInterval.h"
#include "paths/reporting/PerfAlign.h"
#include "paths/reporting/Sog.h"

sog::sog( const basevector& s, const bitvector& g, const vec<perf_align>& aligns )
          : aligns_(aligns)
{    s_.resize( s.size( ) );
     for ( size_t i = 0; i < s.size( ); i++ )
     {    s_[i] = as_base( s[i] );
          if ( g[i] ) s_[i] = 'N';    }    }

void sog::Print( const int id ) const
{    if ( aligns_.empty( ) ) cout << id << "[" << BaseCount( ) << "] free\n";
     vec< vec<String> > rows;
     for ( int l = 0; l < AlignsCount( ); l++ )
     {    const perf_align& a = Align(l);
          vec<String> row;
          row.push_back( 
               ToString( a.Id1( ) ) + "[" + ToString( BaseCount( ) ) + "]",
               ToString( a.pos1( ) ) + "-" + ToString( a.Pos1( ) ),
               ToString( a.Id2( ) ),
               ToString( a.pos2( ) ) + "-" + ToString( a.Pos2( ) ),
               ( a.Rc( ) ? "-" : "+" ), ToString( a.Len( ) ) );
          rows.push_back(row);    }
     PrintTabular( cout, rows, 2, "rrrrrr" );    }

void sog::PrintFasta( ostream& out, const int start, const int stop, 
     const String title ) const
{    if ( title != "" ) out << ">" << title << "\n";
     int count = 0;
     for ( int j = start; j < stop; j++ )
     {    if ( count++ == 80 ) 
          {    out << "\n";
               count = 1;    }
          out << Base(j);    }
     out << "\n";    }

vec<ho_interval> sog::Tigs( ) const
{    vec<ho_interval> tigs;
     for ( int i = 0; i < BaseCount( ); i++ )
     {    if ( Gap(i) ) continue;
          int j;
          for ( j = i + 1; j < BaseCount( ); j++ )
               if ( Gap(j) ) break;
          tigs.push( i, j );
          i = j;    }
     return tigs;    }

void sog::Reverse( )
{    s_.ReverseMe( );
     for ( int i = 0; i < BaseCount( ); i++ )
          SetBase( i, CompBase(i) );
     for ( int i = 0; i < AlignsCount( ); i++ )
     {    int len = Align(i).Pos1( ) - Align(i).pos1( );
          Align(i).pos1( ) = BaseCount( ) - Align(i).Pos1( );
          Align(i).Pos1( ) = Align(i).pos1( ) + len;
          Align(i).Rc( ) = !Align(i).Rc( );    }
     Sort(aligns_);    }

void sog::Validate( ) const
{    for ( int i = 0; i < AlignsCount( ); i++ )
     {    const perf_align& a = Align(i);
          int id2 = a.Id2( ), off = a.Offset( );
          Bool fail = False;
          if ( a.Fw( ) )
          {    for ( int j = a.pos1( ); j < a.Pos1( ); j++ )
                    if ( Base(j) != GenomeBase( id2, j-off ) ) fail = True;    }
          else
          {    for ( int j = a.pos1( ); j < a.Pos1( ); j++ )
               {    if ( CompBase(j) != 
                         GenomeBase( id2, a.Pos2( ) - ( j- a.pos1( ) ) - 1 ) )
                    {    fail = True;    }    }    }
          if (fail)
          {    cout << "\nsog validation failed" << endl;
               PRINT3( a.Id1( ), a.Id2( ), int( a.Fw( ) ) );
               PRINT4( a.pos1( ), a.Pos1( ), a.pos2( ), a.Pos2( ) );
               cout << "\nassembly sequence:\n";
               int count = 0;
               for ( int j = a.pos1( ); j < a.Pos1( ); j++ )
               {    if ( count++ % 80 == 0 ) cout << "\n";
                    cout << Base(j);    }
               cout << "\n\nreference sequence:\n";
               count = 0;
               vec<int> diffs;
               if ( a.Fw( ) )
               {    for ( int j = a.pos2( ); j < a.Pos2( ); j++ )
                    {    if ( count++ % 80 == 0 ) cout << "\n";
                              if ( GenomeBase( id2, j ) != Base( j + a.Offset( ) ) )
                              diffs.push_back( j - a.pos2( ) );
                         cout << GenomeBase( id2, j );    }    }
               else
               {    for ( int j = a.Pos2( ) - 1; j >= a.pos2( ); j-- )
                    {    if ( count++ % 80 == 0 ) cout << "\n";
                         if ( Comp( GenomeBase( id2, j ) )
                              != Base( a.Pos2( ) - 1 - j + a.pos1( ) ) )
                         {    diffs.push_back( a.Pos2( ) - 1 - j );    }
                         cout << Comp( GenomeBase( id2, j ) );    }    }
               cout << "\n\ndifferences at ";
               for ( int j = 0; j < diffs.isize( ); j++ )
               {    if ( j > 0 ) cout << ", ";
                    if ( j == 5 ) 
                    {    cout << "...";
                         break;    }
                    cout << diffs[j];    }
               cout << endl;
               TracebackThisProcess( );    }    }    }

void sog::Replace( const int start, const int stop, const vec<char>& replacement,
     const int first_align_index, const perf_align& new_align )
{    ForceAssertLe( start, stop );
     ForceAssertLt( first_align_index + 1, AlignsCount( ) );

     // Replace s_.[start,stop) by replacement.

     int new_size = start + replacement.isize( ) + BaseCount( ) - stop;
     int old_size = s_.size( );
     if ( new_size > old_size ) s_.resize(new_size);
     if ( replacement.isize( ) != stop - start )
          memmove( &s_[start] + replacement.isize( ), &s_[stop], old_size - stop );
     for ( int i = 0; i < replacement.isize( ); i++ )
          s_[ start + i ] = replacement[i];
     if ( new_size < s_.isize( ) ) s_.resize(new_size);

     // Change alignments.

     const perf_align& a1 = Align(first_align_index);
     const perf_align& a2 = Align(first_align_index+1);
     ForceAssertGe( start, a1.pos1( ) ); ForceAssertLe( start, a1.Pos1( ) );
     ForceAssertGe( a2.pos1( ), stop ); ForceAssertLe( stop, a2.Pos1( ) );
     ForceAssert( Comparable( a1, a2 ) );
     for ( int i = 0; i < AlignsCount( ); i++ )
     {    if ( i == first_align_index || i == first_align_index + 1 ) continue;
          const perf_align& a = Align(i);
          if ( !( a.Pos1( ) <= start || a.pos1( ) >= stop ) )
          {    PRINT2( a.Id1( ), a.Id2( ) );
               PRINT4( a.pos1( ), a.Pos1( ), a.pos2( ), a.Pos2( ) );    
               PRINT2( a1.Id1( ), a1.Id2( ) );
               PRINT4( a1.pos1( ), a1.Pos1( ), a1.pos2( ), a1.Pos2( ) );    
               PRINT2( a2.Id1( ), a2.Id2( ) );
               PRINT4( a2.pos1( ), a2.Pos1( ), a2.pos2( ), a2.Pos2( ) );
               PRINT2( start, stop );    }
          ForceAssert( a.Pos1( ) <= start || a.pos1( ) >= stop );    }
     vec<perf_align> new_aligns;
     new_aligns.reserve( AlignsCount( ) - 1 );
     for ( int i = 0; i < AlignsCount( ); i++ )
     {    if ( i == first_align_index ) new_aligns.push_back(new_align);
          else if ( i == first_align_index + 1 );
          else if ( Align(i).Pos1( ) <= start ) 
               new_aligns.push_back( Align(i) );
          else if ( Align(i).pos1( ) >= stop )
          {    perf_align a = Align(i);
               a.pos1( ) = a.pos1( ) + replacement.isize( ) - ( stop - start );
               a.Pos1( ) = a.Pos1( ) + replacement.isize( ) - ( stop - start );
               new_aligns.push_back(a);    }    }
     aligns_ = new_aligns;    }

Bool sog::ReplaceOK( const int start, const int stop, const int first_align_index )
{    Bool verbose = False;
     const perf_align& a1 = Align(first_align_index);
     const perf_align& a2 = Align(first_align_index+1);
     if (verbose)
     {    cout << "\nentering ReplaceOK" << endl;
          PRINT2( start, stop ); PRINT(a1); PRINT(a2);    }
     if ( !( start <= stop ) ) return False;
     if ( !( start >= a1.pos1( ) ) ) return False;
     if ( !( start <= a1.Pos1( ) ) ) return False;
     if ( !( a2.pos1( ) >= stop ) ) return False;
     if ( !( stop <= a2.Pos1( ) ) ) return False;
     if (verbose) cout << "passes initial tests" << endl;
     for ( int i = 0; i < AlignsCount( ); i++ )
     {    if ( i == first_align_index || i == first_align_index + 1 ) continue;
          const perf_align& a = Align(i);
          if ( !( a.Pos1( ) <= start || a.pos1( ) >= stop ) ) 
          {    if (verbose) cout << "conflict with " << a << endl;
               return False;    }    }
     return True;    }

vecbasevector* sog::genome_;
vecbitvector* sog::genome_amb_;
