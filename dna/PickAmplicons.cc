/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// PickAmplicons.  Pick nonoverlapping amplicons with given properties from a 
// given genome.  Print the primers needed to amplify.
//
// Caveat: this code may run forever.

#include "Basevector.h"
#include "MainTools.h"
#include "dna/DNAHybridization.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "random/Random.h"

Bool InGCRange( const basevector& b, const int start, const int stop,
     const int low_GC, const int high_GC )
{    int gc = 0;
     for ( int j = start; j < stop; j++ )
          if ( IsGC( b[j] ) ) ++gc;
     if ( double(gc)/double(stop-start) < double(low_GC) / 100.0 ) return False;
     if ( double(gc)/double(stop-start) > double(high_GC) / 100.0 ) return False;
     return True;    }

Bool GetPrimer( const basevector& amplicon, const int low_primer,
     const int high_primer, const int low_melt, const int high_melt, basevector& p )
{    vec< pair< double, int > > X;
     for ( int l = low_primer; l <= high_primer; l++ )
     {    p.SetToSubOf( amplicon, 0, l );
          double T = Tm_NearestNeighbor(p);
          if ( T < low_melt || T > high_melt ) continue;
          X.push( Abs( T - double(low_melt+high_melt)/2.0 ), l );    }
     Sort(X);
     if ( X.empty( ) ) return False;
     p.SetToSubOf( amplicon, 0, X[0].second );
     return True;    }

int main( int argc, char *argv[] )
{
     RunTime( );
     
     BeginCommandArguments;
     CommandArgument_String_Doc(GENOME, "fastb file for genome");
     CommandArgument_Int_Doc(LENGTH, "amplicon length");
     CommandArgument_String_Doc(GC, "GC range, e.g. 50-55");
     CommandArgument_Int_Doc(N, "number of amplicons");
     CommandArgument_String_Doc(PRIMER_LENGTH,
          "generate primers in this length range");
     CommandArgument_String_OrDefault_Doc(TITLE, "",
          "title prefix for individual fasta headers");
     CommandArgument_String_Doc(MELT, "primer melting temperature range" );
     CommandArgument_Int_OrDefault_Doc(LEFT_TRIM, 0, 
        "if provided, for each chosen amplicon (length LENGTH), also consider "
        "the amplicon obtain by removing LEFT_TRIM bases from its left; we require "
        "that this new amplicon also have GC content in the given range; with this "
        "option you get three primers per amplicon");
     CommandArgument_String_OrDefault_Doc(TAIL_LEFT, "", "tail for left primers");
     CommandArgument_String_OrDefault_Doc(TAIL_RIGHT, "", "tail for right primers");
     EndCommandArguments;

     if ( TITLE != "" ) TITLE += " ";

     int low_GC = GC.Before( "-" ).Int( ), high_GC = GC.After( "-" ).Int( );
     int low_primer = PRIMER_LENGTH.Before( "-" ).Int( ); 
     int high_primer = PRIMER_LENGTH.After( "-" ).Int( );
     int low_melt = MELT.Before( "-" ).Int( ), high_melt = MELT.After( "-" ).Int( );

     vecbasevector genome(GENOME);
     vec< vec<ho_interval> > amplicons( genome.size( ) );
     for ( int i = 0; i < N; i++ )
     {    size_t tig;
          int xstart = randomx( ) % genome.sumSizes();
          for ( tig = 0; tig < genome.size( ); tig++ )
          {    xstart -= genome[tig].isize( );
               if ( xstart < 0 ) break;    }
          const basevector& g = genome[tig];
          while(1)
          {    int gc = 0;
               int start = randomx( ) % ( g.isize( ) - LENGTH );
               int stop = start + LENGTH;
               if ( !InGCRange( g, start, stop, low_GC, high_GC ) ) continue;
               int start2 = start + LEFT_TRIM;
               if ( !InGCRange( g, start2, stop, low_GC, high_GC ) ) continue;
               basevector b, p;
               b.SetToSubOf( genome[tig], start, stop - start );
               if ( !GetPrimer(b, low_primer, high_primer, low_melt, high_melt, p ) )
                    continue;
               b.ReverseComplement( );
               GetPrimer( b, low_primer, high_primer, low_melt, high_melt, p );
               if ( !GetPrimer(b, low_primer, high_primer, low_melt, high_melt, p ) )
                    continue;
               b.SetToSubOf( genome[tig], start2, stop - start2 );
               if ( !GetPrimer(b, low_primer, high_primer, low_melt, high_melt, p ) )
                    continue;
               b.ReverseComplement( );
               GetPrimer( b, low_primer, high_primer, low_melt, high_melt, p );
               if ( !GetPrimer(b, low_primer, high_primer, low_melt, high_melt, p ) )
                    continue;
               ho_interval h( start, stop );
               if ( Overlap( h, amplicons[tig] ) > 0 ) continue;
               cout << "found amplicon " << i+1 << " of " << N << endl;
               amplicons[tig].push_back(h);    
               break;    }    }

     int count = 0;
     for ( size_t tig = 0; tig < genome.size( ); tig++ )
     {    for ( int j = 0; j < amplicons[tig].isize( ); j++ )
          {    int start = amplicons[tig][j].Start( ); 
               int stop = amplicons[tig][j].Stop( );
               basevector b, p;
               b.SetToSubOf( genome[tig], start, stop - start );
               GetPrimer( b, low_primer, high_primer, low_melt, high_melt, p );
               double T = Tm_NearestNeighbor(p);
               p = Cat( basevector(TAIL_LEFT), p );
               p.Print( cout, TITLE + "amplicon " + ToString(++count) 
                    + ", primer 1, Tm = " + ToString( T, 1 ) );
               if ( LEFT_TRIM > 0 )
               {    int start2 = start + LEFT_TRIM;
                    b.SetToSubOf( genome[tig], start2, stop - start2 );
                    GetPrimer( b, low_primer, high_primer, low_melt, high_melt, p );
                    T = Tm_NearestNeighbor(p);
                    p = Cat( basevector(TAIL_LEFT), p );
                    p.Print( cout, TITLE + "amplicon " + ToString(count) 
                         + ", primer 1', Tm = " + ToString( T, 1 ) );    }
               b.ReverseComplement( );
               GetPrimer( b, low_primer, high_primer, low_melt, high_melt, p );
               T = Tm_NearestNeighbor(p);
               p = Cat( basevector(TAIL_RIGHT), p );
               p.Print( cout, TITLE + "amplicon " + ToString(count) 
                    + ", primer 2, Tm = " + ToString( T, 1 ) );    }    }    }
