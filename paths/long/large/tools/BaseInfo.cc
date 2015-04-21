///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// BaseInfo.  Extract bases from an assembly.

#include "Basevector.h"
#include "MainTools.h"
#include "TokenizeString.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault_Doc(DIR, ".", "looks for DIR/a.{fastb}");
     CommandArgument_String_Doc(E, "semicolon-separated list of arguments, each "
          "defining a path in the assembly, where each path is specified by a "
          "comma-separated list of segment identifiers, and each segment "
          "identifier has either the form e (an edge identifier), or e.start-stop, "
          "specifying a zero-based half open interval on edge e");
     EndCommandArguments;

     int K;
     if ( IsRegularFile( DIR + "/a.k" ) )
     {    Ifstream( in, DIR + "/a.k" );
          in >> K;    }
     else K = 200; // FOR BACKWARD COMPATIBILITY, NOTE HARDCODED!!!!!!!!!!

     vec<int> ids;
     vec<String> seq;
     Tokenize( E, ';', seq );
     for ( int i = 0; i < seq.isize( ); i++ )
     {    vec<String> x;
          Tokenize( seq[i], ',', x );
          for ( int j = 0; j < x.isize( ); j++ )
          {    if ( x[j].IsInt( ) ) ids.push_back( x[j].Int( ) );
               else if ( x[j].Contains( "." ) && x[j].Before( "." ).IsInt( ) )
               {    ids.push_back( x[j].Before( "." ).Int( ) );    }
               else
               {    cout << "I can't parse " << x[j] << "." << endl;
                    Scram(1);    }    }    }
     UniqueSort(ids);
     vecbasevector tigs;
     tigs.Read( DIR + "/a.fastb", ids );
     vec<basevector> X( seq.size( ) );
     for ( int i = 0; i < seq.isize( ); i++ )
     {    vec<String> x;
          Tokenize( seq[i], ',', x );
          int n = x.size( );
          if ( n == 0 )
          {    cout << "That doesn't make sense." << endl;
               Scram(1);    }
          vec< triple<int,int,int> > segs(n);
          for ( int j = 0; j < n; j++ )
          {    int e;
               if ( !x[j].Contains( "." ) ) e = x[j].Int( );
               else e = x[j].Before( "." ).Int( );
               int p = BinPosition( ids, e ), start, stop;
               if ( !x[j].Contains( "." ) )
               {    start = 0, stop = tigs[p].size( );    }
               else
               {    if ( !x[j].Contains( "-" )
                         || !x[j].Between( ".", "-" ).IsInt( )
                         || !x[j].After( "-" ).IsInt( ) )
                    {    cout << "I can't parse " << x[j] << "." << endl;
                         Scram(1);    }
                    start = x[j].Between( ".", "-" ).Int( );
                    stop = x[j].After( "-" ).Int( );
                    if ( !(start >= 0 && start < stop && stop <= tigs[p].isize( )) )
                    {    cout << "The start/stop values specified in " << x[j]
                              << " don't make sense." << endl;
                         Scram(1);    }    }
               segs[j] = make_triple( p, start, stop );    }
          vec<basevector> B(n);
          for ( int j = 0; j < n; j++ )
          {    int p = segs[j].first, start = segs[j].second, stop = segs[j].third;
               B[j] = basevector( tigs[p], start, stop - start );
               if ( j > 0 )
               {    if ( start != 0 
                         || segs[j-1].third != tigs[ segs[j-1].first ].isize( ) )
                    {    cout << x[j-1] << "," << x[j]
                              << " doesn't make sense." << endl;
                         Scram(1);    }    }    }
          basevector A = B[0];
          for ( int j = 1; j < n; j++ )
          {    if ( A.isize( ) < K || B[j].isize( ) < K )
               {    cout << "Problem converting " << x[i] << "." << endl;
                    Scram(1);    }
               basevector left( A, A.isize( ) - (K-1), K-1 );
               basevector right( B[j], 0, K-1 );
               if ( left != right )
               {    cout << "It looks like " << x[i] << " is inconsistent "
                         << "with the graph structure." << endl;
                    Scram(1);    }
               A.resize( A.isize( ) - (K-1) );
               A.append( B[j] );    }
          X[i] = A;    }
     for ( int i = 0; i < X.isize( ); i++ )
     {    X[i].Print( cout, seq[i] + " (bases=" + ToString( X[i].isize( ) )
               + ",kmers=" + ToString( X[i].isize( ) - K + 1 ) + ")" );    }    }
