///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// NEED TO LOWER CASE HOMOPOLYMERS.

// Decode a pseudo-path file, defined as a sequence of lines of the following types:
//
// 1. edge id
// 2. edge id[-n] (delete n bases on right)
// 3. [-n]edge id (delete n bases on left)
// 4. base sequence (allowing ambiguous bases and lower case).
//
// These sequences are concatenated except 1/1 or 1/2 or 3/1 or 3/2 (all effectively
// two edge ids abutting), where we delete K-1 bases and concatenate.
//
// Anything # or after is ignored.
// Blank lines are ignored.
//
// Homopolymers of length >= 30 are converted to lower case.

#include "FastIfstream.h"
#include "MainTools.h"
#include "paths/HyperBasevector.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault_Doc(DIR, ".", "looks for DIR/a.{hbv}");
     CommandArgument_String_Doc(IN, "input file");
     EndCommandArguments;

     HyperBasevector hb;
     BinaryReader::readFile( DIR + "/a.hbv", &hb );
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     int K = hb.K( );

     fast_ifstream in(IN);
     String line;
     String b;
     Bool at_edge = False;
     int last_e = -1;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          if ( line.Contains( "#" ) ) line = line.Before( "#" );
          line.GlobalReplaceBy( " ", "" );
          if ( line.size( ) == 0 ) continue;
          if ( line.IsInt( ) )
          {    int e = line.Int( );
               if (at_edge) 
               {    b.resize( b.isize( ) - (K-1) );
                    if ( to_right[last_e] != to_left[e] )
                    {    cout << e << " doesn't follow" << endl;
                         Scram(1);    }    }
               b += hb.O(e).ToString( );
               last_e = e;
               at_edge = True;    }
          else if ( line.Contains( "[-", 0 ) )
          {    ForceAssert( line.Contains( "]" ) );
               int d = line.Between( "[-", "]" ).Int( );
               int e = line.After( "]" ).Int( );
               String s = hb.O(e).ToString( );
               b += hb.O(e).ToString( ).substr( d, s.isize( ) - d );
               last_e = e;
               at_edge = True;    }
          else if ( line.Contains( "[-" ) )
          {    ForceAssert( line.Contains( "]", -1 ) );
               int e = line.Before( "[" ).Int( );
               int d = line.Between( "[-", "]" ).Int( );
               String s = hb.O(e).ToString( );
               if (at_edge) 
               {    b.resize( b.isize( ) - (K-1) );
                    if ( to_right[last_e] != to_left[e] )
                    {    cout << e << " doesn't follow" << endl;
                         Scram(1);    }    }
               b += hb.O(e).ToString( ).substr( 0, s.isize( ) - d );
               at_edge = False;    }
          else
          {    for ( int i = 0; i < line.isize( ); i++ )
               {    if ( !isalpha( line[i] ) )
                    {    PRINT(line);
                         ForceAssert( isalpha( line[i] ) );    }    }
               b += line;
               at_edge = False;    }    }

     const int min_poly = 30;
     for ( int i = 0; i < b.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < b.isize( ); j++ )
               if ( b[j] != b[i] ) break;
          if ( j - i >= min_poly && b[i] != 'N' )
          {    for ( int k = i; k < j; k++ )
                    b[k] = tolower(b[k]);    }
          i = j - 1;    }

     cout << ">seq" << endl;
     for ( int j = 0; j < b.isize( ); j++ )
     {    if ( j > 0 && j % 80 == 0 ) cout << "\n";
          cout << b[j];    }
     cout << endl;    }
