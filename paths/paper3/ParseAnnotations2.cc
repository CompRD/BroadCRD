///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Parse segmental duplication annotations.

#include "FastIfstream.h"
#include "MainTools.h"
#include "TokenizeString.h"

int main( )
{    
     RunTime( );

     for ( int pass = 1; pass <= 2; pass++ )
     {    String dir = ( pass == 1 ? "/wga/scr1/ALLPATHS/M.musculus"
               : "/wga/scr1/ALLPATHS/H.sapiens.NA12878" );
          String headers = ( pass == 1 ? "mm8" : "other_ref/build18" );
          fast_ifstream hin( dir + "/" + headers + ".headers" );
          String line;
          vec<String> ids, tokens;
          while(1)
          {    getline( hin, line );
               if ( hin.fail( ) ) break;
               ids.push_back( line.After( ">" ) );    }
          String afile = ( pass == 1 ? "WGAC.tab" :
               "oo.weild10kb.join.all.cull.xwparse.txt" );
          fast_ifstream rin( dir + "/annotation/" + afile );
          getline( rin, line );
          vec<char> tab;
          tab.push_back( '\t' );
          Ofstream( out, dir + "/annotation/annotations_segdup" );
          vec< triple<int,int,int> > dups;
          while(1)
          {    getline( rin, line );
               if ( rin.fail( ) ) break;
               TokenizeStrictly( line, tab, tokens );
               int ntokens = ( pass == 1 ? 10 : 8 );
               ForceAssertGe( tokens.isize( ), ntokens );
               if ( tokens[0].Contains( "random" ) ) continue;
               int chr1 = Position( ids, tokens[0] );
               if ( chr1 < 0 ) PRINT( tokens[0] );
               ForceAssertGe( chr1, 0 );
               int start1 = tokens[1].Int( ) - 1, stop1 = tokens[2].Int( );
               int ind2 = ( pass == 1 ? 6 : 4 );
               if ( tokens[ind2].Contains( "random" ) ) continue;
               int chr2 = Position( ids, tokens[ind2] );
               if ( chr2 < 0 ) PRINT( tokens[ind2] );
               ForceAssertGe( chr2, 0 );
               if ( tokens[ind2+1].Int( ) > tokens[ind2+2].Int( ) ) 
                    swap( tokens[ind2+1], tokens[ind2+2] );
               int start2 = tokens[ind2+1].Int( ) - 1, stop2 = tokens[ind2+2].Int( );
               dups.push( chr1, start1, stop1 );
               dups.push( chr2, start2, stop2 );    }
          UniqueSort(dups);
          for ( int i = 0; i < dups.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < dups.isize( ); j++ )
               {    if ( dups[j].first != dups[i].first ) break;
                    if ( dups[j].second > dups[j-1].third ) break;    }
               out << dups[i].first << " " << dups[i].second << " "
                    << dups[j-1].third << "\n";
               i = j - 1;    }    }    }
