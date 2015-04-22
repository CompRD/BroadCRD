/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "bias/DinukeBias.h"
#include "lookup/LookAlign.h"
#include <sstream>

void DinukeBias( const vecbasevector& reads, const vecbasevector& ref,
     const vecbasevector& rref, const vec<look_align>& aligns )
{
     vec<int> in_read(4);
     vec< vec<int> > dinuke_all(4), dinuke_reads(4), dinuke_at_break(4);
     for ( int i = 0; i < 4; i++ )
     {    dinuke_all[i].resize( 4, 0 );
          dinuke_reads[i].resize( 4, 0 );
          dinuke_at_break[i].resize( 4, 0 );    }
     for ( int j = 0; j < aligns.isize( ); j++ )
     {    const look_align& la = aligns[j];
          int id = la.query_id;
          int id2 = la.target_id;
          align a = la.a;
          const basevector& rd1 = reads[la.query_id];
          for ( unsigned int k = 0; k < rd1.size( ); k++ )
               ++in_read[ rd1[k] ];
          const basevector& rd2 = ( !la.rc1 ? ref[id2] : rref[id2] );
          if (la.rc1) a.ReverseThis( rd1.size( ), rd2.size( ) );
          if ( a.pos2( ) > 0 && !la.rc1 ) 
          {    ++dinuke_at_break[ rd2[ a.pos2( ) - 1 ] ][ rd2[ a.pos2( ) ] ];
               for ( int j = a.pos2( ) + 1; j < a.Pos2( ); j++ )
                    ++dinuke_reads[ rd2[ j - 1 ] ][ rd2[ j ] ];    }    }

     for ( size_t j = 0; j < ref.size( ); j++ )
     {    for ( unsigned int i = 1; i < ref[j].size( ); i++ )
          {    ++dinuke_all[ ref[j][i-1] ][ ref[j][i] ];
               ++dinuke_all[ rref[j][i-1] ][ rref[j][i] ];    }    }
     int dinuke_all_total = 0, dinuke_reads_total = 0, dinuke_at_break_total = 0;
     for ( int i = 0; i < 4; i++ )
     {    for ( int j = 0; j < 4; j++ )
          {    dinuke_all_total += dinuke_all[i][j];
               dinuke_reads_total += dinuke_reads[i][j];
               dinuke_at_break_total += dinuke_at_break[i][j];    }    }
     vec< vec<String> > rows;
     vec<String> row1, row2;
     row1.push_back( "", "at", "in", "in" );
     row2.push_back( "", "break", "read", "reference" );
     rows.push_back( row1, row2 );
     for ( int i = 0; i < 4; i++ )
     {    for ( int j = 0; j < 4; j++ )
          {    vec<String> row;
               String d(2);
               d[0] = as_base(i);
               d[1] = as_base(j);
               row.push_back(d);
               ostringstream p1, p2, p3;
               p1 << PERCENT_RATIOB( 5, 2, dinuke_at_break[i][j], 
                    dinuke_at_break_total );
               p2 << PERCENT_RATIOB( 5, 2, dinuke_reads[i][j], 
                    dinuke_reads_total );
               p3 << PERCENT_RATIOB( 5, 2, dinuke_all[i][j], 
                    dinuke_all_total );
               row.push_back( p1.str( ), p2.str( ), p3.str( ) );
               rows.push_back(row);    }    }
     PrintTabular( cout, rows, 2, "rrrr" );
     cout << "\n";
     rows.clear( ), row1.clear( ), row2.clear( );
     row1.push_back( "", "before", "after", "in",   "in" );
     row2.push_back( "", "break",  "break", "read", "reference" );
     rows.push_back( row1, row2 );
     for ( int i = 0; i < 4; i++ )
     {    vec<String> row;
          String d(1);
          d[0] = as_base(i);
          row.push_back(d);
          ostringstream p1, p2, p3, p4;
          p1 << PERCENT_RATIOB( 5, 2, dinuke_at_break[i][0] + dinuke_at_break[i][1]
               + dinuke_at_break[i][2] + dinuke_at_break[i][3],
               dinuke_at_break_total );
          p2 << PERCENT_RATIOB( 5, 2, dinuke_at_break[0][i] + dinuke_at_break[1][i]
               + dinuke_at_break[2][i] + dinuke_at_break[3][i],
               dinuke_at_break_total );
          p3 << PERCENT_RATIOB( 5, 2, in_read[i], Sum(in_read) );
          p4 << PERCENT_RATIOB( 5, 2, dinuke_all[i][0] + dinuke_all[i][1]
               + dinuke_all[i][2] + dinuke_all[i][3],
               dinuke_all_total );
          row.push_back( p1.str( ), p2.str( ), p3.str( ), p4.str( ) );
          rows.push_back(row);    }
     PrintTabular( cout, rows, 2, "lrrr" );
     cout << "\n";    }
