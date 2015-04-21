/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// ParsePathsToLocs.  Parse output of PathsToLocs, run with options
// UNIPATH=True UNIPATHS=unipaths SHOW_PLACEMENTS=True.  Ugly.

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "graph/Digraph.h"
#include "paths/GetNexts.h"
#include "paths/UnibaseUtils.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(IN);
     CommandArgument_String(UNIBASES);
     EndCommandArguments;

     // Set up to read.

     fast_ifstream in(IN);
     String line;

     // Find K.

     int K = -1;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) 
          {    cout << "Couldn't find K." << endl;
               exit(1);    }
          if ( line.Contains( " K=" ) )
          {    K = line.Between( " K=", " " ).Int( );
               break;    }    }

     // Get number of unipaths.

     int nuni = -1;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) 
          {    cout << "Couldn't find unipaths count." << endl;
               exit(1);    }
          if ( line.Contains( " unipaths, " ) )
          {    nuni = line.Before( " unipaths, " ).Int( );
               break;    }    }

     // Load unibases and compute graph structure.

     vecbasevector unibases(UNIBASES);
     vec<int> to_rc;
     UnibaseInvolution( unibases, to_rc );
     vec< vec<int> > before(nuni), after;
     GetNexts( K, unibases, after );
     for ( int u = 0; u < nuni; u++ )
     {    before[u] = after[ to_rc[u] ];
          for ( int j = 0; j < before[u].isize( ); j++ )
               before[u][j] = to_rc[ before[u][j] ];
          Sort( after[u] );
          Sort( before[u] );    }
     digraph A( after, before );

     // Get placements.

     vec<String> placements;
     vec< vec<int> > index(nuni);
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) 
          {    cout << "Couldn't find adjacency graph." << endl;
               exit(1);    }
          if ( line.Contains( "Adjacency", 0 ) ) break;
          if ( line.Contains( "path", 0 ) ) 
          {    int id = line.Between( " ", " " ).Int( );
               index[id].push_back( placements.size( ) );
               placements.push_back(line);    }    }

     // Build blocks.

     String place;
     vec<String> block;
     vec< vec<String> > blocks;
     int last_tig = -1, last_end = -1;
     for ( int i = 0; i < placements.isize( ); i++ )
     {    const String& p = placements[i];
          place = p.Between( "--> ", " " );
          int tig = place.Before( "." ).Int( );
          int start = place.Between( ".", "-" ).Int( );
          int end = place.After( "-" ).Int( );
          if ( block.empty( ) );
          else if ( tig == last_tig && last_end - start == K - 1 );
          else
          {    blocks.push_back(block);
               block.clear( );    }
          block.push_back(p);
          last_tig = tig, last_end = end;    }
     if ( block.nonempty( ) ) blocks.push_back(block);

     // Try to extend and merge blocks.

     vec< vec<String> > blocks2;
     for ( int i = 0; i < blocks.isize( ); i++ )
     {    if ( blocks[i].front( ).Contains( "path " )
               && blocks[i].front( ).After( "path " ).Contains( " " )
               && !blocks[i].front( ).Contains( " or " ) )
          {    String b = blocks[i].front( );
               b.GlobalReplaceBy( "  ", " " );
               int id0 = b.Between( " ", " " ).Int( );
               int npushes = 0;
               grow_back: ++npushes;
               if ( npushes <= 10 )
               {    if ( A.To(id0).empty( ) ) blocks[i].push_front( "stop" );
                    else if ( A.To(id0).solo( ) )
                    {    id0 = A.To(id0)[0];
                         blocks[i].push_front( "path " + ToString(id0) );
                         goto grow_back;    }
                    else if ( A.To(id0).size( ) == 2 )
                    {    int x1 = A.To(id0)[0], x2 = A.To(id0)[1];
                         blocks[i].push_front( 
                              "path " + ToString(x1) + " or " + ToString(x2) );
                         if ( A.To(x1) == A.To(x2) && A.To(x1).solo( ) )
                         {    id0 = x1;
                              goto grow_back;    }    }    }    }
          if ( i < blocks.isize( ) - 1 )
          {    String b1 = blocks[i].back( ), b2 = blocks[i+1].front( );
               b1.GlobalReplaceBy( "  ", " " );
               b2.GlobalReplaceBy( "  ", " " );
               int id1 = b1.Between( " ", " " ).Int( );
               int id2 = b2.Between( " ", " " ).Int( );
               vec<String> between;
               Bool bridged = False;
               while(1)
               {    if ( A.From(id1).solo( ) )
                    {    solo:
                         id1 = A.From(id1)[0];
                         if ( id1 == id2 )
                         {    bridged = True;
                              break;    }
                         else 
                         {    between.push_back( "path " + ToString(id1) );    
                              if ( between.size( ) > 10 ) break;    }    }
                    else if ( A.From(id1).size( ) == 2 )
                    {    int x1 = A.From(id1)[0], x2 = A.From(id1)[1];
                         if ( x1 == id2 || x2 == id2 )
                         {    bridged = True;
                              break;    }
                         between.push_back( 
                              "path " + ToString(x1) + " or " + ToString(x2) );
                         if ( between.size( ) > 10 ) break;
                         if ( A.From(x1) == A.From(x2) && A.From(x1).solo( ) )
                         {    id1 = x1;
                              goto solo;    }
                         else break;    }
                    else if ( A.From(id1).empty( ) )
                    {    between.push_back( "stop" );
                         break;    }
                    else break;    }
               if (bridged)
               {    blocks[i].append(between);
                    blocks[i].append( blocks[i+1] );
                    blocks[i+1] = blocks[i];    }
               else
               {    blocks[i].append(between);
                    blocks2.push_back( blocks[i] );    }    }
          else blocks2.push_back( blocks[i] );    }

     // Print blocks.
     
     for ( int i = 0; i < blocks2.isize( ); i++ )
     {    cout << "\n[" << i+1 << "]\n";
          for ( int j = 0; j < blocks2[i].isize( ); j++ )
               cout << blocks2[i][j] << "\n";    }    }
