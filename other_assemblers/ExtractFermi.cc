///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Extract a fermi assembly.

#include "FastIfstream.h"
#include "CoreTools.h"
#include "other_assemblers/ConnectionAssembly.h"
#include "other_assemblers/ExtractFermi.h"

void ExtractFermi( const String& IN_DIR, connection_assembly& A )
{
     // Check for file.

     if ( !IsDirectory(IN_DIR) )
     {    FatalErr("Can't find input directory.");    }
     String fn = IN_DIR + "/fmdef.p2.mag.gz";
     if ( !IsRegularFile(fn) )
     {    FatalErr("Failed to find " << fn << ".");    }
     String fn0 = IN_DIR + "/fmdef.p2.mag";

     // Parse sequences.

     vecbasevector bases;
     String line;
     PipeIstream( xin, fn0 );
     while(1)
     {    getline( xin, line );
          if ( xin.fail( ) ) break;
          getline( xin, line );
          ForceAssert( !xin.fail( ) );
          basevector b( line.size( ) );
          for ( int j = 0; j < b.isize( ); j++ )
               b.Set( j, as_char( line[j] ) );
          getline( xin, line );
          ForceAssert( !xin.fail( ) );
          getline( xin, line );
          ForceAssert( !xin.fail( ) );
          bases.push_back_reserve(b);    }
     xin_pipe.close( );

     // Parse connections.

     fast_pipe_ifstream in( "zcat " + fn + " | grep @" );
     map< int64_t, pair<int,int> > verts;
     vec< vec< pair<int64_t,int> > > lefts, rights;
     int nv = 0;
     Bool verbose = False;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          if (verbose) cout << line << endl;
          ForceAssert( line.Contains( ":" ) );
          ForceAssert( line.After( ":" ).Contains( "\t" ) );
          String v1 = line.Between( "@", ":" ), v2 = line.Between( ":", "\t" );
          ForceAssert( v1.IsInt( ) ), ForceAssert( v2.IsInt( ) );
          int64_t n1 = v1.Int( ), n2 = v2.Int( );
          ForceAssert( verts.find(n1) == verts.end( ) );
          ForceAssert( verts.find(n2) == verts.end( ) );
          verts[n1] = make_pair(nv,0), verts[n2] = make_pair(nv,1);
          nv++;
          ForceAssert( line.After( "\t" ).Contains( "\t" ) );
          String post = line.After( "\t" ).After( "\t" );
          ForceAssert( post.Contains( "\t" ) );
          String l1 = post.Before( "\t" ), l2 = post.After( "\t" );
          if ( l1 == "." ) l1 = "";
          if ( l2 == "." ) l2 = "";
          vec< pair<int64_t,int> > x1, x2;
          while( l1 != "" )
          {    ForceAssert( l1.Contains( ";" ) );
               String s = l1.Before( ";" );
               l1 = l1.After( ";" );
               ForceAssert( s.Contains( "," ) );
               ForceAssert( s.Before( "," ).IsInt( ) );
               ForceAssert( s.After( "," ).IsInt( ) );
               x1.push( s.Before( "," ).Int( ), s.After( "," ).Int( ) );    }
          while( l2 != "" )
          {    ForceAssert( l2.Contains( ";" ) );
               String s = l2.Before( ";" );
               l2 = l2.After( ";" );
               ForceAssert( s.Contains( "," ) );
               ForceAssert( s.Before( "," ).IsInt( ) );
               ForceAssert( s.After( "," ).IsInt( ) );
               x2.push( s.Before( "," ).Int( ), s.After( "," ).Int( ) );    }
          lefts.push_back(x1), rights.push_back(x2);    }

     // Build graph.

     vec< vec<int> > from(nv), to(nv), from_edge_obj(nv), to_edge_obj(nv);
     vec<connection> edges;
     vec< pair< pair<int,int>, triple<int,int,int> > > X;
     int min_overlap = 1000000000;
     for ( int id1 = 0; id1 < nv; id1++ )
     {    for ( int j = 0; j < lefts[id1].isize( ); j++ )
          {    int z = lefts[id1][j].first;
               ForceAssert( verts.find(z) != verts.end( ) );
               int id2 = verts[z].first, pos2 = verts[z].second, pos1 = 0;
               int over = lefts[id1][j].second;
               X.push( make_pair(id1,id2), make_triple(pos1,pos2,over) );    }
          for ( int j = 0; j < rights[id1].isize( ); j++ )
          {    int z = rights[id1][j].first;
               ForceAssert( verts.find(z) != verts.end( ) );
               int id2 = verts[z].first, pos2 = verts[z].second, pos1 = 1;
               int over = rights[id1][j].second;
               X.push( make_pair(id1,id2), make_triple(pos1,pos2,over) );    }    }
     UniqueSort(X);
     for ( int i = 0; i < X.isize( ); i++ )
     {    int id1 = X[i].first.first, id2 = X[i].first.second;
          int pos1 = X[i].second.first, pos2 = X[i].second.second;
          int over = X[i].second.third;
          min_overlap = Min( min_overlap, over );
          basevector b1 = bases[id1], b2 = bases[id2];
          Bool rc1 = ( pos1 != 1 ), rc2 = ( pos2 != 0 );
          if (rc1) b1.ReverseComplement( );
          if (rc2) b2.ReverseComplement( );
          ForceAssert( b1.Overlap( b2, over ) );
          from[id1].push_back(id2), to[id2].push_back(id1);
          from_edge_obj[id1].push_back( edges.size( ) );
          to_edge_obj[id2].push_back( edges.size( ) );
          edges.push( rc1, rc2, over );    }    
     A.bases = bases;
     A.G.Initialize( from, to, vec<connection_assembly::node_t>(bases.size()), edges, to_edge_obj, from_edge_obj );    }
