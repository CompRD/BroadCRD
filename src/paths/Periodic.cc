///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "paths/Periodic.h"

Bool IsPeriodic( const basevector& x, const int top_period )
{    for ( int p = 1; p <= top_period; p++ )
     {    Bool fail = False;
          for ( int j = 0; j < p; j++ )
          {    if (fail) break;
               for ( int k = j + p; k < x.isize( ); k += p )
               {    if ( x[k] != x[j] )
                    {    fail = True;
                         break;    }    }    }
          if ( !fail ) return True;    }
     return False;    }
