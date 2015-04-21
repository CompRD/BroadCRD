/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "Qualvector.h"
#include "cancer/AssembleSNP.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/Unipath.h"

class call {
     public:
     unsigned char base;
     unsigned char qual;
     int read_id;
     int read_pos;
     call( const unsigned char base, const unsigned char qual, 
          const int read_id, const int read_pos ) : base(base), qual(qual),
          read_id(read_id), read_pos(read_pos) { }
};

Bool AssembleSNP( const int maxread, const vecbasevector& bases, 
     const vecqualvector& quals, const vec<look_align>& aligns, 
     const vecbasevector& ref, const int tig, const int pos, const char altbase, 
     const Bool verbose )
{
     // Verbose logging.

     if (verbose) 
     {    cout << "\nentering AssembleSNP\n";
          PRINT3( tig, pos, as_base(altbase) );
          cout << "\nreads:\n";
          int rcount = 0;
          for ( int i = 0; i < aligns.isize( ); i++ )
          {    const look_align& la = aligns[i];
               if ( la.target_id != tig || pos < la.pos2( ) || pos >= la.Pos2( ) ) 
                    continue;
               basevector rd = bases[ la.query_id ];
               if ( la.Rc1( ) ) rd.ReverseComplement( );
               cout << "[" << ++rcount << "] " << rd.ToString( ) << "\n";    }    }

     // Find the reads that cross this locus.  If there are a huge number,
     // discard most of them.

     vec<look_align> alignsx;
     for ( int i = 0; i < aligns.isize( ); i++ )
     {    const look_align& la = aligns[i];
          int p2 = la.a.pos2( );
          if ( la.target_id != tig || pos < p2 || pos >= la.a.Pos2( ) ) continue;
          alignsx.push_back( aligns[i] );    }
     const int max_aligns = 500;
     if ( alignsx.isize( ) > max_aligns )
     {    if (verbose)
          {    cout << "reducing " << aligns.size( ) << " aligns to "
                    << max_aligns << endl;    }
          vec<look_align> alignsx2;
          int incr = alignsx.isize( ) / max_aligns;
          for ( int i = 0; i < alignsx.isize( ); i += incr )
          {    alignsx2.push_back( alignsx[i] );
               if ( alignsx2.isize( ) == max_aligns ) break;    }
          alignsx = alignsx2;    }

     // Build a cheap multiple alignment.  Mark indel positions.  

     vec< vec<call> > calls( 2*maxread + 1 );
     vec<int> indel( 2*maxread + 1, 0 ), indel_loc( 2*maxread + 1, -1 );
     basevector rd;
     qualvector q;
     for ( int i = 0; i < alignsx.isize( ); i++ )
     {    const look_align& la = alignsx[i];
          int p1 = la.pos1( ), p2 = la.a.pos2( ), id = la.query_id;
          rd = bases[id];
          q = quals[id];
          if ( la.Rc1( ) ) 
          {    rd.ReverseComplement( );
               q.ReverseMe( );    }
          p2 = p2 - pos + maxread;
          for ( int j = 0; j < la.a.Nblocks( ); j++ ) 
          {    if ( la.a.Gaps(j) > 0 )
               {    if ( p2 >= 0 && p2 < calls.isize( ) )
                    {    indel[p2]++;
                         indel_loc[p2] = id;    }
                    p2 += la.a.Gaps(j);
                    if ( p2 >= 0 && p2 < calls.isize( ) )
                    {    indel[p2]++;
                         indel_loc[p2] = id;    }    }
               if ( la.a.Gaps(j) < 0 ) 
               {    if ( p2 >= 0 && p2 < calls.isize( ) )
                    {    indel[p2]++;
                         indel_loc[p2] = id;    }
                    p1 -= la.a.Gaps(j);    }
               for ( int x = 0; x < la.a.Lengths(j); x++ ) 
               {    if ( p2 >= 0 && p2 < calls.isize( ) )
                         calls[p2].push( rd[p1], q[p1], id, p1 );
                    ++p1; ++p2;    }    }    }
     if (verbose)
     {    for ( int i = 0; i < calls.isize( ); i++ )
          {    for ( int j = 0; j < calls[i].isize( ); j++ )
                    cout << as_base( calls[i][j].base );
               cout << endl;    }    }

     // Find reads containing an unsupported indel.

     vec<Bool> unsupported( bases.size( ), False );
     for ( int i = 0; i < indel.isize( ); i++ )
          if ( indel[i] == 1 ) unsupported[ indel_loc[i] ] = True;

     // Define good positions.  The statistical approach here is primitive.  We
     // require that no indel on any read was detected at the position, that the
     // runner-up base has not more than 10% of the quality score points of the
     // winning base, and that the runner-up base has at most 40 quality score
     // points.  Also there must be at least one read.

     vec<Bool> good( 2*maxread + 1, False );
     basevector goodb( 2*maxread + 1 );
     for ( int i = 0; i < good.isize( ); i++ )
     {    if ( indel[i] > 1 ) continue;
          vec<int> count( 4, 0 ), ids( 4, vec<int>::IDENTITY );
          for ( int j = 0; j < calls[i].isize( ); j++ )
               count[ calls[i][j].base ] += calls[i][j].qual;
          ReverseSortSync( count, ids );
          goodb.Set( i, ids[0] );
          if ( count[0] > 0 && count[1] <= 40 && 10 * count[1] <= count[0] ) 
               good[i] = True;    }
     if (verbose)
     {    for ( int i = 0; i < good.isize( ); i++ )
          {    if ( good[i] ) cout << "1";
               else cout << "0";    }
          cout << endl;
          for ( int i = 0; i < good.isize( ); i++ )
               cout << as_base( goodb[i] );
          cout << endl;    }
          
     // Find flanking K-mers consisting of good positions.

     const int K = 20;
     int jl, jr;
     Bool lgoodness(False), rgoodness(False);
     for ( jl = 1; jl <= maxread - K; jl++ )
     {    lgoodness = True;
          for ( int u = 0; u < K; u++ )
               if ( !good[ maxread - jl - u ] ) lgoodness = False;
          if (lgoodness) break;    }
     for ( jr = 1; jr <= maxread - K; jr++ )
     {    rgoodness = True;
          for ( int u = 0; u < K; u++ )
               if ( !good[ maxread + jr + u ] ) rgoodness = False;
          if (rgoodness) break;    }
     if (verbose) PRINT3( maxread, int(lgoodness), int(rgoodness) );
     if ( !lgoodness || !rgoodness ) return False;
     int left = maxread - jl - (K-1), right = maxread + jr + K;
     basevector leftkmer, rightkmer;
     leftkmer.SetToSubOf( goodb, left, K );
     rightkmer.SetToSubOf( goodb, right - K, K );
     vecbasevector ends(2);
     ends[0] = leftkmer, ends[1] = rightkmer;
     if (verbose)
     {    cout << "leftkmer = " << leftkmer.ToString( ) << endl; 
          cout << "rightkmer = " << rightkmer.ToString( ) << endl;    }

     // Find read kmers between the left and right kmers.  For now we just use
     // the entire reads.  Add the interval on the reference.

     vecbasevector readlets;
     for ( int i = 0; i < alignsx.isize( ); i++ )
     {    const look_align& la = alignsx[i];
          int id = la.query_id;
          if ( unsupported[id] ) continue;
          rd = bases[id];
          if ( la.Rc1( ) ) rd.ReverseComplement( );
          int p1 = la.a.pos1( ), p2 = la.a.pos2( );
          if ( la.target_id != tig || pos < p2 || pos >= la.a.Pos2( ) ) continue;
          readlets.push_back_reserve(rd);    }
     basevector reflet;
     reflet.SetToSubOf( ref[tig], left + pos - maxread, right - left );
     // PRINT( reflet.ToString( ) ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     // PRINT( maxread - left ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     readlets.push_back_reserve(reflet);

     // Build the unipath graph.

     vecKmerPath paths, pathsrc, unipaths;
     vec<tagged_rpint> pathsdb, unipathsdb;
     int nb = readlets.size( );
     vecbasevector basesplus(readlets);
     basesplus.Append(ends);
     ReadsToPathsCoreY( basesplus, K, paths );
     CreateDatabase( paths, pathsdb );
     KmerBaseBroker kbb( K, paths, pathsrc, pathsdb, basesplus );
     Unipath( paths, pathsrc, pathsdb, unipaths, unipathsdb );
     if (verbose) PRINT( unipaths.size( ) );
     digraph A;
     BuildUnipathAdjacencyGraph( paths, pathsrc, pathsdb, unipaths, unipathsdb, A );
     HyperKmerPath h;
     BuildUnipathAdjacencyHyperKmerPath( K, A, unipaths, h );
     if (verbose)
     {    for ( size_t i = 0; i < unipaths.size( ); i++ )
               cout << "unipath " << i << " = "
                    << kbb.Seq( unipaths[i] ).ToString( ) << endl;    }

     // Find the ends.

     longlong x1 = paths[nb].FirstSegment( ).Start( );
     longlong x2 = paths[nb+1].LastSegment( ).Stop( );
     static vec<longlong> con1, con2;
     Contains( unipathsdb, x1, con1 ), Contains( unipathsdb, x2, con2 );
     ForceAssert( con1.nonempty( ) ); ForceAssert( con2.nonempty( ) );
     int u1 = unipathsdb[ con1[0] ].ReadId( ), u2 = unipathsdb[ con2[0] ].ReadId( );

     // Define end trim amounts.

     const tagged_rpint &t1 = unipathsdb[ con1[0] ], &t2 = unipathsdb[ con2[0] ];
     const KmerPath &p1 = unipaths[ t1.ReadId( ) ], &p2 = unipaths[ t2.ReadId( ) ];
     int left_trim = x1 - t1.Start( );
     for ( int j = 0; j < t1.PathPos( ); j++ )
          left_trim += p1.Length(j);
     int right_trim = t2.Stop( ) - x2;
     for ( int j = t2.PathPos( ) + 1; j < p2.NSegments( ); j++ )
          right_trim += p2.Length(j);

     // Print adjacencies.

     if (verbose)
     {    cout << "\nAdjacencies:\n";
          for ( int v = 0; v < h.N( ); v++ )
          {    for ( int j1 = 0; j1 < h.To(v).isize( ); j1++ )
               {    int e1 = h.EdgeObjectIndexByIndexTo( v, j1 );
                    for ( int j2 = 0; j2 < h.From(v).isize( ); j2++ )
                    {    int e2 = h.EdgeObjectIndexByIndexFrom( v, j2 );
                         PRINT2( e1, e2 );    }    }    }    }

     // Find everything between u1 and u2.

     vec<int> to_left, to_right;
     h.ToLeft(to_left), h.ToRight(to_right);
     vec<int> from_u1_end, to_u2_end;
     h.GetSuccessors1( to_right[u1], from_u1_end );
     if (verbose)
     {    cout << "from_u1_end:";
          for ( int j = 0; j < from_u1_end.isize( ); j++ )
               cout << " " << from_u1_end[j];
          cout << endl;    }
     if ( !BinMember( from_u1_end, to_left[u2] ) ) return False;
     h.GetPredecessors1( to_left[u2], to_u2_end );
     if (verbose)
     {    cout << "to_u2_end:";
          for ( int j = 0; j < to_u2_end.isize( ); j++ )
               cout << " " << to_u2_end[j];
          cout << endl;    }
     vec<int> between = Intersection( from_u1_end, to_u2_end );

     // Find all paths from u1 to u2.  Fail if cycle encountered or two paths
     // have different lengths.  Also fail if altbase does not occur on a path.

     vec< vec<int> > allpaths0;
     vec<int> U1;
     U1.push_back(u1);
     allpaths0.push_back(U1);
     vec<int> lengths;
     Bool have_alt = False;
     while( allpaths0.nonempty( ) )
     {    vec<int> p = allpaths0.back( );
          allpaths0.resize( allpaths0.size( ) - 1 );
          int v = to_right[ p.back( ) ];
          for ( int j = 0; j < h.From(v).isize( ); j++ )
          {    int w = h.From(v)[j];
               int e = h.EdgeObjectIndexByIndexFrom( v, j );
               if ( e != u2 && !BinMember( between, w ) ) continue;
               if (verbose)
               {    if ( Member( p, e ) )
                    {    for ( int l = 0; l < p.isize( ); l++ )
                         {    if ( l > 0 ) cout << " ";
                              cout << p[l];    }
                         cout << " " << e << endl;    }    }
               if ( Member( p, e ) ) return False;
               vec<int> q(p);
               q.push_back(e);
               allpaths0.push_back(q);    }
          if ( p.back( ) == u2 ) 
          {    basevector b = kbb.Seq( h.EdgeObject( p[0] ) );
               for ( int j = 1; j < p.isize( ); j++ )
               {    b.resize( b.size( ) - K + 1 );
                    b = Cat( b, kbb.Seq( h.EdgeObject( p[j] ) ) );    }
               if (verbose)
               {    cout << "allpath " << ":";
                    for ( int j = 0; j < p.isize( ); j++ )
                         cout << " " << p[j];
                    cout << ", len = " << b.size( ) << "\n";
                    cout << b.ToString( ) << "\n";    }
               b.SetToSubOf( b, left_trim, b.isize( ) - left_trim - right_trim );
               if ( lengths.nonempty( ) && b.isize( ) != lengths[0] )
               {    if (verbose)
                    {    cout << "see different lengths: " << lengths[0]
                              << " " << b.size( ) << endl;    }
                    return False;    }
               if ( lengths.empty( ) ) lengths.push_back( b.size( ) );
               if ( maxread-left >= 0 && maxread-left < b.isize( ) 
                    && b[maxread-left] == altbase ) 
               {    have_alt = True;    }    }    }
     UniqueSort(lengths);
     if (verbose) PRINT( lengths.size( ) );
     if ( !lengths.solo( ) ) return False;
     if (verbose) PRINT( int(have_alt) );
     return have_alt;    }
