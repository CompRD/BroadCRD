///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Form iterative consensus.  Extremely slow, and not necessarily good.  Over and
// over, find best edit of first sequence and make it.

#include "Basevector.h"
#include "CoreTools.h"
#include "paths/long/magic/BasicScore.h"
#include "paths/long/magic/IterCon.h"

basevector IterCon( const vecbasevector& segs )
{
     basevector b = segs[0];

     int errs = 0;
     for ( int j = 0; j < (int) segs.size( ); j++ )
          errs += BasicScore( b, segs[j] );
     // PRINT(errs);

     for ( int pass = 1; pass <= 100; pass++ )
     {    Bool improved = False;
          {    
               int best_errs2_ins = errs;
               int best_p_ins = -1;
               int best_ins = -1;

               int best_errs2_sub = errs;
               int best_p_sub = -1;
               int best_add = -1;

               int best_errs2_del = errs;
               int best_p_del = -1;
     
               for ( int p = 0; p < b.isize( ); p++ )
               {    
                    // Try insertion.

                    for ( int ins = 0; ins < 4; ins++ )
                    {    basevector s;
                         for ( int j = 0; j < b.isize( ); j++ )
                         {    s.push_back( b[j] );
                              if ( j == p ) s.push_back(ins);    }
                         int errs2 = 0;
                         for ( int j = 0; j < (int) segs.size( ); j++ )
                              errs2 += BasicScore( s, segs[j] );
                         if ( errs2 < best_errs2_ins )
                         {    best_errs2_ins = errs2;
                              best_p_ins = p;
                              best_ins = ins;    }    }
     
                    // Try deletion.

                    basevector s;
                    for ( int j = 0; j < b.isize( ); j++ )
                         if ( j != p ) s.push_back( b[j] );
                    int errs2 = 0;
                    for ( int j = 0; j < (int) segs.size( ); j++ )
                         errs2 += BasicScore( s, segs[j] );
                    if ( errs2 < best_errs2_del )
                    {    best_errs2_del = errs2;
                         best_p_del = p;    }
               
                    // Try substitution.
     
                    for ( int add = 1; add <= 3; add++ )
                    {    basevector s( b );
                         s.Set( p, ( s[p] + add ) % 4 );
                         int errs2 = 0;
                         for ( int j = 0; j < (int) segs.size( ); j++ )
                              errs2 += BasicScore( s, segs[j] );
                         if ( errs2 < best_errs2_sub )
                         {    best_errs2_sub = errs2;
                              best_p_sub = p;
                              best_add = add;    }    }    }
     
               if ( best_errs2_ins < errs && best_errs2_ins <= best_errs2_sub 
                    && best_errs2_ins <= best_errs2_del )
               {    int errs2 = best_errs2_ins;
                    int p = best_p_ins;
                    int ins = best_ins;
                    basevector s;
                    for ( int j = 0; j < b.isize( ); j++ )
                    {    s.push_back( b[j] );
                         if ( j == p ) s.push_back(ins);    }
                    b = s;
                    errs = errs2;
                    // PRINT2(pass,errs);    
                    improved = True;    }

               else if ( best_errs2_sub < errs && best_errs2_sub <= best_errs2_ins
                    && best_errs2_sub <= best_errs2_del )
               {    int errs2 = best_errs2_sub;
                    int p = best_p_sub;
                    int add = best_add;
                    basevector s( b );
                    s.Set( p, ( s[p] + add ) % 4 );
                    b = s;
                    errs = errs2;
                    // PRINT2(pass,errs);    
                    improved = True;    }
     
               else if ( best_errs2_del < errs && best_errs2_del <= best_errs2_ins
               && best_errs2_del <= best_errs2_sub )
               {    int errs2 = best_errs2_del;
                    int p = best_p_del;
                    basevector s;
                    for ( int j = 0; j < b.isize( ); j++ )
                         if ( j != p ) s.push_back( b[j] );
                    b = s;
                    errs = errs2;
                    // PRINT2(pass,errs);    
                    improved = True;    }    }
     
          if ( !improved ) break;    }

     return b;    }
