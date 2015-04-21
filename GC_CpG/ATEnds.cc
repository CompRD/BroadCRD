// Copyright (c) 2000, 2001 Whitehead Institute for Biomedical Research

// ATEnds: compute AT-content of first and last 1kb of supercontigs >= 100kb.
// (If first or last contig is < 1kb, just compute content of it.)

#include <iomanip>

#include "Alignment.h"
#include "AnnotatedContig.h"
#include "Basevector.h"
#include "math/Functions.h"
#include "MainTools.h"
#include "ReadLocation.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String(SUBDIR);
     CommandArgument_UnsignedInt_OrDefault(MIN_SUPER, 100000);
     CommandArgument_UnsignedInt_OrDefault(WINDOW, 1000);
     CommandArgument_UnsignedInt_OrDefault(MIN_AT_PERCENT, 70);
     EndCommandArguments;

     int window_size = WINDOW;

     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/" + SUBDIR;
     String source_dir = run_dir + 
          ( IsDirectory( run_dir + "/preFindSeeds" ) ? "/preFindSeeds" : "" );

     // Set up standard data structures.

     READ( sub_dir + "/mergedcontigs.answer", vec<annotated_supercontig>, asupers );
     vecbasevector tigs;
     tigs.ReadAll( sub_dir + "/mergedcontigs.fastb" );

     vec<float> AT;
     for ( int i = 0; i < (int) asupers.size( ); i++ )
     {    annotated_supercontig& as = asupers[i];
          longlong total = 0;
          int n = as.NumContigs( );
          for ( int j = 0; j < n; j++ )
               total += tigs[ as.Contig(j).ID( ) ].size( );
          if ( total < (int) MIN_SUPER ) continue;
          float frontAT, backAT;
          {    int bases = 0, ATbases = 0;
               int m = as.Contig(0).ID( );
               int w = Min( WINDOW, tigs[m].size( ) );
               for ( int j = 0; j < w; j++ )
               {    ++bases;
                    unsigned char b = tigs[m][j];
                    if ( as_base(b) == 'A' || as_base(b) == 'T' ) ++ATbases;    }
               frontAT = float(ATbases)/float(bases);
               AT.push_back(frontAT);    
               basevector b;
               b.SetToSubOf( tigs[m], 0, w );
               b.Print( cout, "super " + ToString(i) + " begin, AT = "
                    + ToString(ATbases) + " of " + ToString(bases) );    }
          {    int bases = 0, ATbases = 0;
               int m = as.Contig(n-1).ID( );
               int N = tigs[m].size( );
               int w = Min( (int) WINDOW, N );
               for ( int j = (int) N - 1; j >= N - w; j-- )
               {    ++bases;
                    unsigned char b = tigs[m][j];
                    if ( as_base(b) == 'A' || as_base(b) == 'T' ) ++ATbases;    }
               backAT = float(ATbases)/float(bases);
               AT.push_back(backAT);
               basevector b;
               b.SetToSubOf( tigs[m], N - w, w );
               b.Print( cout, "super " + ToString(i) + " end, AT = "
                    + ToString(ATbases) + " of " + ToString(bases) );    }    }

     int count = 0;
     for ( int i = 0; i < (int) AT.size( ); i++ )
          if ( AT[i] >= float(MIN_AT_PERCENT)/100.0 ) ++count;
     cout << setprecision(3) << 100.0 * float(count) / float( AT.size( ) )
          << "% of begin and end supercontig windows have " << " >= "
          << MIN_AT_PERCENT << "% AT\n";

     // Compute fraction of AT-rich windows which are in 5kb ends of 
     // supercontigs.

     int AT_rich_windows = 0, end_AT_rich_windows = 0;
     for ( int i = 0; i < (int) asupers.size( ); i++ )
     {    annotated_supercontig& as = asupers[i];
          longlong total = 0;
          int n = as.NumContigs( );
          for ( int j = 0; j < n; j++ )
               total += tigs[ as.Contig(j).ID( ) ].size( );
          if ( total < (int) MIN_SUPER ) continue;
          int pos = 0;
          for ( int j = 0; j < n; j++ )
          {    int m = as.Contig(j).ID( );
               for ( int k = 0; k < (int) tigs[m].size( ) - window_size;
                    k += window_size/10 )
               {    static basevector b;
                    b.SetToSubOf( tigs[m], k, window_size );
                    int AT = 0;
                    for ( int l = 0; l < window_size; l++ )
                    {    unsigned char c = b[l];
                         if ( as_base(c) == 'A' || as_base(c) == 'T' ) ++AT;    }
                    if ( float(AT)/float(window_size) 
                         >= float(MIN_AT_PERCENT)/100.0 )
                    {    ++AT_rich_windows;
                         if ( pos + k + window_size <= 5000
                              || pos + k >= (int) tigs[m].size( ) - 5000 )
                              ++end_AT_rich_windows;    }    }    }    }
     cout << end_AT_rich_windows << " of " << AT_rich_windows
          << " (" << setprecision(3) << 100.0 * float(end_AT_rich_windows)
          / float(AT_rich_windows) << "%) are in terminal 5kb\n";    }
