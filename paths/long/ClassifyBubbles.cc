///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// ClassifyBubbles.  Find the bubbles in an assembly and determine what fraction
// are SNPs.  Should be more general and should be part of something more general.
//
// Added some useful assembly stats.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "paths/long/SupportedHyperBasevector.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(IN_HEAD, "looks for IN_HEAD.shbv");
     EndCommandArguments;

     SupportedHyperBasevector shb;
     BinaryReader::readFile( IN_HEAD + ".shbv", &shb );

     int snp = 0, other = 0;
     #pragma omp parallel for
     for ( int v = 0; v < shb.N( ); v++ )
     {    if ( shb.To(v).size( ) != 1 || shb.From(v).size( ) != 2 ) continue;
          if ( shb.From(v)[0] != shb.From(v)[1] ) continue;
          int w = shb.From(v)[0];
          if ( shb.To(w).size( ) != 2 || shb.From(w).size( ) != 1 ) continue;
          int x = shb.To(v)[0], y = shb.From(w)[0];
          basevector b1 = shb.EdgeObjectByIndexFrom( v, 0 );
          basevector b2 = shb.EdgeObjectByIndexFrom( v, 1 );
          if ( b1.size( ) > b2.size( ) ) swap( b1, b2 );
          int best_loc;
          alignment al;
          SmithWatFree( b1, b2, best_loc, al );
          align a(al);
          vec<int> mgg = a.MutationsGap1Gap2( b1, b2 );
          #pragma omp critical
          {    if ( a.pos1( ) == 0 && a.Pos1( ) == b1.isize( )
                    && a.pos2( ) == 0 && a.Pos2( ) == b2.isize( ) 
                    && mgg[0] == 1 && Sum(mgg) == 1 ) 
               {    snp++;    }
               else other++;    }    }
     cout << "bubbles: ";
     PRINT2( snp, other );
     cout << PERCENT_RATIO( 3, snp, snp+other ) << " of bubbles are SNPs" << endl;

     cout << ToStringAddCommas( shb.EdgeObjectCount( ) ) << " edges and " 
          << shb.NComponents( ) << " components" << endl;
     vec<int> len;
     for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
          len.push_back( shb.EdgeLengthKmers(e) );
     Sort(len);
     cout << "N50 edge = " << ToStringAddCommas( N50(len) ) << endl;
     int64_t kmers = 0;
     for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
          if ( shb.Inv(e) <= e ) kmers += shb.EdgeLengthKmers(e);
     cout << "total kmers = " << ToStringAddCommas(kmers) << endl;    }

