// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#include "Alignment.h"
#include "math/Arith.h"
#include "CoreTools.h"
#include "pairwise_aligners/RecoverAlignment.h"
#include "Seed.h"
#include "Set.h"

Bool seed::Merge( unit u, int left, int right, const alignment& a_left,
     const alignment& a_right, const vec<alignment_plus>& all_aligns, 
     vec<alignment_plus>& extra_aligns, vec< vec<int> >& all_aligns_index, 
     const vecbasevector& EE, const vecqualvector& Q, 
     int MaxSeparationDiscrepancy, ostream& out,
     set< pair<int, int> >& no_alignment, Float max_score_for_compatibility )
{    unit& uL = (*this)[left], uR = (*this)[right];
     // unit& u0 = (*this)[0];
     u.SetLeftOffset( uL.LeftOffset( ) + a_left.pos1( ) - a_left.pos2( ) );
     u.SetRightOffset( uR.RightOffset( ) + a_right.pos1( ) - a_right.pos2( ) );
     vec<int> discrep( size( ) );

     for ( unsigned int i = 0; i < size( ); i++ )
          discrep[i] = SepDiscrep( u, (*this)[i] );
     sort( discrep.begin( ), discrep.end( ) );

     if ( Max( Abs( discrep[ (2*size( ))/3 ] ), Abs( discrep[ size( )/3 ] ) )
          > MaxSeparationDiscrepancy ) return False;

     for ( unsigned int j1 = 0; j1 < size( ); j1++ )
     {    unit& u1 = (*this)[j1];
          unit& u2 = u;
          int id1, id2, offset, length1, length2, overlap;
          for ( int pass = 1; pass <= 2; pass++ )
          {    id1 = (pass == 1) ? u1.LeftRead( ) : u1.RightRead( );
               id2 = (pass == 1) ? u2.LeftRead( ) : u2.RightRead( );
               offset = (pass == 1) ? u2.LeftOffset( ) - u1.LeftOffset( ) :
                    u2.RightOffset( ) - u1.RightOffset( );
               length1 = (pass == 1) ? u1.LeftLength( ) : u1.RightLength( );
               length2 = (pass == 1) ? u2.LeftLength( ) : u2.RightLength( );
               if ( offset >= 0 ) overlap = Min(length1 - offset, length2);
               else overlap = Min( length2 + offset, length1 );

               if ( overlap >= 100 )
               {    bool overlap_found = false;
                    int reported_offset = 0;
                    Bool if_rc = False;
                    Float score = -1.0;
                    for (unsigned int j = 0; j < all_aligns_index[id1].size( ); j++)
                    {    int entry = all_aligns_index[id1][j];
                         const alignment_plus& a = entry >= 0
                              ? all_aligns[entry] : extra_aligns[-entry-1];
                         ForceAssert( a.Id1( ) == id1 ); // XXX AAA
                         if ( a.Id2( ) != id2 ) continue;
                         if_rc = a.Rc2( );
                         reported_offset = a.a.pos1( ) - a.a.pos2( );
                         score = a.score;
                         overlap_found = true;
                         break;    }
                    if ( overlap_found && score > max_score_for_compatibility )
                         return False;
                    if ( overlap_found &&
                         ( Abs(offset - reported_offset) > 30 || if_rc ) )
                         return False;

                    /*
                    if ( !overlap_found || Abs(offset - reported_offset) > 30 )
                    {    PRINT(id1);
                         PRINT(id2);
                         PRINT(overlap);
                         PRINT( int(overlap_found) );
                         PRINT(offset);
                         PRINT(reported_offset);    }    
                    */

                    /*
                    Make sure this code works first....
                    if ( !overlap_found )
                    {    Bool alignment_found = false;
                         alignment_plus ap;
                         if ( RecoverAlignment( id1, id2, EE, Q, False, 
                              ap, 12, offset - 30, offset + 30 ) )
                         {    alignment_found = True;
                              out << "pushing back alignment between " // XXX
                                   << id1 << " and " << id2 << "\n"; // XXX
                              extra_aligns.push_back(ap);
                              all_aligns_index[id1].push_back( 
                                   -extra_aligns.size( ) );
                              ap.Swap( EE[ap.Id1( )].Length( ),
                                   EE[ap.Id2( )].Length( ) );
                              extra_aligns.push_back(ap);
                              all_aligns_index[id2].push_back( 
                                   -extra_aligns.size( ) );    }
                         if ( !alignment_found ) return False;   }   }   }   }
                    */

                    if ( !overlap_found )
                    {    if ( Member( no_alignment, make_pair(id1, id2) ) )
                              return False;
                         Bool alignment_found = false;
                         alignment_plus ap;
                         if ( RecoverAlignment( id1, id2, EE, Q, False, 
                              ap, 12 ) )
                         {    int new_offset = ap.a.pos1( ) - ap.a.pos2( );
                              if ( Abs( offset - new_offset ) <= 30 &&
                                   !ap.Rc2( ) )
                                   alignment_found = True;
                              extra_aligns.push_back(ap);
                              all_aligns_index[id1].push_back( 
                                   -extra_aligns.size( ) );
                              ap.Swap( EE[ap.Id1( )].size( ),
                                   EE[ap.Id2( )].size( ) );
                              extra_aligns.push_back(ap);
                              all_aligns_index[id2].push_back( 
                                   -extra_aligns.size( ) );    }
                         else no_alignment.insert( make_pair( id1, id2 ) );
                         if ( !alignment_found ) return False;   }   }   }   }

     push_back(u); 
     return True;    }

Bool seed::Merge( const seed& s2 )
{    seed& s1 = *this;
     seed s(s1);
     for ( unsigned int j = 0; j < s2.size( ); j++ )
     {    int i = s1.LeftMember( s2[j].LeftRead( ) );
          if ( i >= 0 )
          {    for ( unsigned int k = 0; k < s2.size( ); k++ )
               {    if ( s1.LeftMember( s2[k].LeftRead( ) ) >= 0 ) continue;
                    unit u( s2[k] );
                    u.Shift( s1[i].LeftOffset( ) - s2[j].LeftOffset( ),
                         s1[i].RightOffset( ) - s2[j].RightOffset( ) );
                    s.push_back(u);    }
               *this = s;
               return True;    }
          i = s1.LeftMember( s2[j].RightRead( ) );
          if ( i >= 0 )
          {    seed s2s = s2.Swap( );
               for ( unsigned int k = 0; k < s2s.size( ); k++ )
               {    if ( s1.LeftMember( s2s[k].LeftRead( ) ) >= 0 ) continue;
                    unit u( s2s[k] );
                    u.Shift( s1[i].LeftOffset( ) - s2s[j].LeftOffset( ),
                         s1[i].RightOffset( ) - s2s[j].RightOffset( ) );
                    s.push_back(u);    }
               *this = s;
               return True;    }    }
     return False;    }
