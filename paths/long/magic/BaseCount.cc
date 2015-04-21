///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Compute the "adjusted correct base count", the number of correct aligned bases
// at a given position, minus 1/3 the number of incorrect aligned bases.

#include "Basevector.h"
#include "MainTools.h"
#include "lookup/LookAlign.h"
#include "paths/long/SupportedHyperBasevector.h"

// QueryLookupTable K=12 MM=12 MC=0.005 SEQS=Jan30.2014.3/all.fastb L=scs.lookup 
//      SMITH_WAT=True FW_ONLY=True PARSEABLE=True > Jan30.2014.3/all.aligns

int main( )
{    RunTime( );

     SupportedHyperBasevector shb;
     BinaryReader::readFile( "/wga/scr4/macro/scs.shbv", &shb );

     vec<look_align> aligns;
     LoadLookAligns( "/wga/scr4/macro/Jan30.2014.3/all.aligns", aligns );

     vec< vec<int> > cov1( shb.EdgeObjectCount( ) ), cov2( shb.EdgeObjectCount( ) );
     for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
     {    cov1[e].resize( shb.EdgeLengthBases(e), 0 );
          cov2[e].resize( shb.EdgeLengthBases(e), 0 );    }

     vecbasevector magic;
     magic.ReadAll( "/wga/scr4/macro/Jan30.2014.3/all.fastb" );

     for ( int i = 0; i < aligns.isize( ); i++ )
     {    const look_align& la = aligns[i];
          basevector b = magic[ la.query_id ];
          int e = la.target_id;
          const basevector& rd2 = shb.EdgeObject(e);
          if ( la.Rc1( ) ) b.ReverseComplement( );
          align a = la.a;
          int p1 = a.pos1( ), p2 = a.pos2( );
          for ( int j = 0; j < a.Nblocks( ); j++ ) 
          {    if ( a.Gaps(j) > 0 ) p2 += a.Gaps(j);
               if ( a.Gaps(j) < 0 ) p1 -= a.Gaps(j);
               for ( int x = 0; x < a.Lengths(j); x++ ) 
               {    if ( b[p1] == rd2[p2] ) cov1[e][p2]++;
                    else cov2[e][p2]++;
                     ++p1; ++p2;    }    }    }

     for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
     {    cout << "\n";
          for ( int j = 0; j < cov1[e].isize( ); j++ )
          {    cout << e << "." << j << " " << cov1[e][j] - cov2[e][j]/3.0 
                    << "\n";    }    }    }
