///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Follow up from Snorgle.  This substitutes for PHASE=3.
//
// Let K = 100.  First find the closures whose first K-mer occurs uniquely in the 
// assembly, in the appropriate sense, as follows.
//
// This is done by restricting each assembly edge v --e--> w to a subinterval.  The 
// default is to use [0,|e|-(KK-K)), where KK is "big K" (nominally 200).  However 
// if exactly one edge enters v, then we change the start to KK-K, so long as this 
// keeps the interval length positive.  Finally if e is the only edge ending in w, 
// then we change the stop to |e|.
//
// These are start closures.  Similarly define stop closures.  Next define a link 
// from one closure to another if there is a proper overlap between them of 
// length >= K.  Then look for paths from start closures to stop closures, via a 
// sequence of overlaps, giving up if too many paths are seen.  These paths define 
// patches.  
// 
// Then canonicalize the patch ends.  To do this, starting at the left end of the
// patch, we match up the patch to the graph, continuing until it no longer matches,
// and mark this point.  We do the same thing on the other end.  If we do not have
// KK bases on the left end that match, we traverse the graph backwards uniquely
// to extend.  (To continue.)

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "paths/HyperBasevector.h"
#include "paths/RemodelGapTools.h"

class connector {

     public:

     connector( ) { }
     int e1, e2;
     int p1, p2;
     basevector seq;
     vec< pair<int,int> > rid;

     friend Bool operator<( const connector& c1, const connector& c2 )
     {    if ( c1.e1 < c2.e1 ) return True;
          if ( c1.e1 > c2.e1 ) return False;
          if ( c1.e2 < c2.e2 ) return True;
          if ( c1.e2 > c2.e2 ) return False;
          if ( c1.p1 < c2.p1 ) return True;
          if ( c1.p1 > c2.p1 ) return False;
          if ( c1.p2 < c2.p2 ) return True;
          if ( c1.p2 > c2.p2 ) return False;
          if ( c1.seq < c2.seq ) return True;
          if ( c1.seq > c2.seq ) return False;
          if ( c1.rid < c2.rid ) return True;
          return False;    }

};

int main( )
{    RunTime( );
     double clock = WallClockTime( );

     // Directory path.

     // String dir = "/wga/scr4/jaffe/GapToy/51400.newchem";
     String dir = "/wga/scr4/jaffe/GapToy/helico1/16";

     // Logging.

     Bool print_initial_links = True;
     Bool print_closure_paths = False;

     // Edges to trace.

     // vec<int> trace = { 18572790,18572791,18572792 };
     vec<int> trace;

     // Load assembly.

     cout << Date( ) << ": load genome" << endl;
     const int K = 100;
     vecbasevector genome( dir + "/a.200/a.fastb" );
     HyperBasevectorX hb;
     BinaryReader::readFile( dir + "/a.200/a.hbx", &hb );

     // Build patches based on unique matches between hanging ends.

     vecbasevector patches;
     vec< triple<int,int,int> > joins;
     {    vec< triple<kmer<K>,int,int> > kmers_plus;
          cout << Date( ) << ": building lookup table" << endl;
          MakeKmerLookup0( genome, kmers_plus );
          cout << Date( ) << ": looking for joins" << endl;
          for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
          {    int64_t j;
               for ( j = i + 1; j < kmers_plus.jsize( ); j++ )
                    if ( kmers_plus[j].first != kmers_plus[i].first ) break;
               if ( j - i == 2 )
               {    int e1 = kmers_plus[i].second, pos1 = kmers_plus[i].third;
                    int e2 = kmers_plus[i+1].second, pos2 = kmers_plus[i+1].third;
                    if ( e1 != e2 )
                    {    int offset = pos1 - pos2;
                         Bool mismatch = False;
                         for ( int p1 = 0; p1 < genome[e1].isize( ); p1++ )
                         {    int p2 = p1 - offset;
                              if ( p2 < 0 ) continue;
                              if ( p2 >= genome[e2].isize( ) ) break;
                              if ( genome[e1][p1] != genome[e2][p2] )
                              {    mismatch = True;
                                   break;    }    }
                         if ( !mismatch ) joins.push( e1, e2, offset );    }    }
               i = j - 1;    }    }
     ParallelUniqueSort(joins);
     for ( int64_t i = 0; i < joins.jsize( ); i++ )
     {    int e1 = joins[i].first, e2 = joins[i].second, offset = joins[i].third;
          int low = Min( 0, offset ); 
          int high = Max( genome[e1].isize( ), offset + genome[e2].isize( ) );
          basevector b( high - low );
          for ( int i = 0; i < genome[e1].isize( ); i++ )
               b.Set( i-low, genome[e1][i] );
          for ( int i = 0; i < genome[e2].isize( ); i++ )
               b.Set( i-low+offset, genome[e2][i] );
          patches.push_back(b);    }

     // Load closures.

     cout << Date( ) << ": loading closures" << endl;
     vecbasevector clo( dir + "/closures.fastb" );
     cout << Date( ) << ": there are " << ToStringAddCommas( clo.size( ) )
          << " closures" << endl;
     vecbasevector clo2( dir + "/closures2.fastb" );
     vec<int64_t> clo_ids;
     BinaryReader::readFile( dir + "/closures.ids", &clo_ids );
     vec<int64_t> clo2_ids(clo_ids);
     clo2_ids.append(clo_ids);
     vec<int> clo2_n( clo2.size( ), 0 );
     for ( int64_t i = 1; i < (int64_t) clo2.size( ); i++ )
     {    if ( clo2_ids[i] == clo2_ids[i-1] ) clo2_n[i] = clo2_n[i-1] + 1;
          else clo2_n[i] = 0;    }

     // Fetch closures having a unique genomic start/stop.  

     cout << Date( ) << ": fetching closures having unique start/stop" << endl;
     vec< pair<int,int> > astart, astop;
     BinaryReader::readFile( dir + "/closures.astart", &astart );
     BinaryReader::readFile( dir + "/closures.astop", &astop );

     // Fetch overlaps between closures.  Note that we're loading overlaps between
     // some closures that have just been killed, which doesn't really make sense.

     cout << Date( ) << ": looking for overlaps" << endl;
     vec< vec< pair<int,int> > > O;
     BinaryReader::readFile( dir + "/closures.overlaps", &O );

     // Look for paths.

     cout << Date( ) << ": finding paths" << endl;
     vec<connector> con;
     vec<vec< pair<int,int> >> paths, paths2;
     vec<int> stops;
     for ( int i = 0; i < (int) clo2.size( ); i++ )
     {    if ( astart[i].first < 0 ) continue;
          if ( clo2[i].size( ) == 0 ) continue; // needed?
          paths.clear( ), paths2.clear( );
          vec< pair<int,int> > p = { make_pair( i, 0 ) };
          paths.push_back(p);
          const int max_paths = 100;
          const int max_iterations = 10;
          int iter = 0;
          if ( Member( trace, astart[i].first ) ) 
          {    cout << "\n";
               PRINT2( i, astart[i].first );    }
          stops.clear( );
          basevector b;
          while(1)
          {    if ( ++iter > max_iterations ) break;
               if ( paths.empty( ) ) break;
               stops.clear( );
               for ( int j = 0; j < paths.isize( ); j++ )
               {    if ( astop[ paths[j].back( ).first ].first >= 0 ) 
                         stops.push_back(j);    }
               if ( stops.nonempty( ) )
               {    for ( int s = 0; s < stops.isize( ); s++ )
                    {    int j = stops[s];
                         int e = paths[j].back( ).first;
                         connector c;
                         c.e1 = astart[i].first, c.e2 = astop[e].first;
                         c.p1 = astart[i].second, c.p2 = astop[e].second;
                         c.rid = paths[ stops[s] ];
                         b = clo2[ paths[j].front( ).first ];
                         for ( int l = 1; l < paths[j].isize( ); l++ )
                         {    int rid = paths[j][l].first;
                              int rpos = paths[j][l].second;
                              if ( rpos > b.isize( ) ) // SHOULD ASSERT!
                              {    b = basevector( "ACGT" ); // code for failure
                                   break;    }
                              int ridm = paths[j][l-1].first;
                              b.resize( 
                                   b.isize( ) - ( clo2[ridm].isize( ) - rpos ) );
                              b = Cat( b, clo2[rid] );    }
                         c.seq = b;
                         con.push_back(c);    
                         if (print_closure_paths)
                         {    cout << "\ncon " << con.size( ) << ", closure " << i 
                                   << "(" << c.e1
                                   << "." << c.p1 << ") --..--> closure " << e
                                   << "(" << c.e2 << "." << c.p2 << ")\n";    }    }
                    break;    }

               if ( Member( trace, astart[i].first ) )
               {    cout << "\niter = " << iter << ", paths =\n";
                    for ( int j = 0; j < paths.isize( ); j++ )
                    {    const vec< pair<int,int> >& p = paths[j];
                         vec<int> x;
                         for ( int k = 0; k < p.isize( ); k++ )
                              x.push_back( p[k].first );
                         cout << "[" << j+1 << "] " << x << "\n";    }    }

               if ( paths.isize( ) > max_paths ) break;
               paths2.clear( );
               {    for ( int j = 0; j < paths.isize( ); j++ )
                    {    vec< pair<int,int> > p = paths[j];
                         int e = p.back( ).first;
                         for ( int l = 0; l < O[e].isize( ); l++ )
                         {    vec< pair<int,int> > q = p;
                              if ( clo2[ O[e][l].first ].size( ) == 0 ) 
                                   continue; // needed?
                              q.push_back( O[e][l] );
                              paths2.push_back(q);    }   }   }
               paths = paths2;    }    }

     // Summarize.

     cout << Date( ) << ": sorting" << endl;
     ParallelSort(con);
     cout << Date( ) << ": done" << endl;
     if (print_initial_links) cout << "\n";
     vec< pair<basevector,int64_t> > cans;
     for ( int64_t i = 0; i < con.jsize( ); i++ )
     {    int e1 = con[i].e1, e2 = con[i].e2;
          int64_t j;
          for ( j = i + 1; j < con.jsize( ); j++ )
               if ( con[j].e1 != e1 || con[j].e2 != e2 ) break;
          if ( j - i == 1 ) continue;
          // if ( e1 >= 200000 ) break; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          if (print_initial_links)
          {    cout << "\ncon " << i << ".." << j-1 << ", " << j - i 
                    << " links from " << e1 << " to " << e2 << endl;    }
          for ( int64_t k = i; k < j; k++ )
          {    const connector& c = con[k];
               int p1 = c.p1, p2 = c.p2;
               if (print_initial_links)
               {    Bool show_seq = True;
                    cout << p1 << ".";
                    if (show_seq) cout<< c.seq.ToString( );
                    cout << "." << p2 << " (";
                    for ( int l = 0; l < c.rid.isize( ); l++ )
                    {    if ( l > 0 ) cout << ",";
                         int id = c.rid[l].first, pos = c.rid[l].second;
                         Bool fw = ( id < (int) clo.size( ) );
                         cout << ( fw ? "+" : "-" ) << clo2_ids[id]/2 << "." 
                              << clo2_n[id] << "." << pos;    }
                    cout << ")\n";    }

               // Start to canonicalize closure.  Shorten or extend the closure 
               // sequence on both ends so that there is exactly one hb.K()-mer
               // [  (hb.K()+1)-mer??  ]
               // on each end that matches the assembly.  This could fail.  First 
               // find the path forward through the graph, until the closure no 
               // longer matches the graph.  Then doing the same thing from the 
               // other side.  Finally, as needed, extend c on both ends.

               basevector d = c.seq;

               // First set pos1 to the first mismatch position.

               int e = e1, p = p1-1;
               int pos1;
               for ( pos1 = 0; pos1 < d.isize( ); pos1++ ) // could start at K-1
               {    ++p;
                    if ( p < hb.Bases(e) ) 
                    {    if ( hb.O(e)[p] != d[pos1] ) break;    }
                    else
                    {    Bool match = False;
                         int v = hb.ToRight(e);
                         for ( int r = 0; r < (int) hb.From(v).size( ); r++ )
                         {    int en = hb.IFrom(v,r);
                              if ( hb.O(en)[hb.K( )-1] == d[pos1] )
                              {    e = en;
                                   p = hb.K( )-1;
                                   match = True;
                                   break;    }    }
                         if ( !match ) break;    }    }

               // Now set pos2 to the first match position.

               int pos2;
               e = e2, p = p2;
               for ( pos2 = d.isize( ) - 1; pos2 >= 0; pos2-- ) // could...
               {    --p;
                    if ( p >= 0 )
                    {    if ( hb.O(e)[p] != d[pos2] ) break;    }
                    else
                    {    Bool match = False;
                         int v = hb.ToLeft(e);
                         for ( int r = 0; r < (int) hb.To(v).size( ); r++ )
                         {    int en = hb.ITo(v,r);
                              if ( hb.O(en)[ hb.Bases(en) - hb.K( ) ] == d[pos2] )
                              {    e = en;
                                   p = hb.Bases(en) - hb.K( );
                                   match = True;
                                   break;    }    }
                         if ( !match ) break;    }    }
               pos2++; // now pos2 is the first match position

               // Back up until we have a hb.K( ) + 1 base match (if possible).
               // Note that this is not exactly what we actually do!

               e = e1;
               if ( pos1 < hb.K( ) + 1 && p1 > 0 )
               {    int delta = Min( p1, hb.K( ) + 1 - pos1 );
                    pos1 += delta;
                    basevector dnew( hb.O(e), p1 - delta, delta );
                    dnew.append(d);
                    d = dnew;    }
               while( pos1 < hb.K( ) + 1 )
               {    int v = hb.ToLeft(e);
                    if ( hb.To(v).size( ) != 1 ) break;
                    e = hb.ITo(v,0);
                    pos1 += hb.Kmers(e);
                    basevector dnew( hb.Kmers(e) + d.isize( ) );
                    for ( int l = 0; l < d.isize( ); l++ )
                         dnew.Set( hb.Kmers(e) + l, d[l] );
                    for ( int l = 0; l < hb.Kmers(e); l++ )
                         dnew.Set( l, hb.O(e)[l] );
                    d = dnew;    }

               // Go the other way.

               e = e2;
               if ( d.isize( ) - pos2 < hb.K( ) && p2 < hb.Bases(e2) )
               {    int add 
                         = Min( hb.K( ) - ( d.isize( ) - pos2 ), hb.Bases(e2) - p2 );
                    basevector a( hb.O(e2), p2, add );
                    d.append(a);    }
               while( d.isize( ) - pos2 < hb.K( ) + 1 )
               {    int v = hb.ToRight(e);
                    if ( hb.From(v).size( ) != 1 ) break;
                    e = hb.IFrom(v,0);
                    int nd = d.size( );
                    d.resize( nd + hb.Kmers(e) );
                    for ( int l = 0; l < hb.Kmers(e); l++ )
                         d.Set( nd + l, hb.O(e)[ hb.K( ) - 1 + l ] );    }

               // Trim.

               if ( pos1 < hb.K( ) || d.isize( ) - pos2 < hb.K( ) )
               {    if (print_initial_links) 
                    {    cout << "can't be canonicalized" << endl;
                         PRINT4( pos1, pos2, c.seq.size( ), d.size( ) );    }    }
               else
               {    int ltrim = Max( 0, pos1 - hb.K( ) - 1 );
                    int rtrim = Max( 0, d.isize( ) - pos2 - hb.K( ) - 1 );
                    if ( d.isize( ) - ltrim - rtrim < hb.K( ) )
                    {    if (print_initial_links)
                              cout << "too much trim" << endl;    }
                    else
                    {    basevector dnew( d.isize( ) - ltrim - rtrim );
                         for ( int l = ltrim; l < d.isize( ) - rtrim; l++ )
                              dnew.Set( l-ltrim, d[l] );
                         d = dnew;    
                         if (print_initial_links)
                              cout << "canonicalized: " << d.ToString( ) << endl;    
                         cans.push( d, k );    }    }    }
          i = j - 1;    }
     if (print_initial_links) cout << "\n";

     // Sort out canonical closures.

     cout << Date( ) << ": sorting " << cans.size( ) 
          << " canonical closures" << endl;
     ParallelSort(cans);
     cout << Date( ) << ": select canonical closures" << endl;
     for ( int64_t i = 0; i < cans.jsize( ); i++ )
     {    int64_t j;
          const basevector& c = cans[i].first;
          for ( j = i + 1; j < cans.jsize( ); j++ )
               if ( cans[j].first != c ) break;
          patches.push_back(c);
          cout << "\nclosure: " << c.ToString( ) << endl;
          for ( int64_t z = i; z < j; z++ )
          {    int64_t k = cans[z].second;
               const connector& c = con[k];
               int e1 = c.e1, e2 = c.e2, p1 = c.p1, p2 = c.p2;
               cout << e1 << "." << p1 << " --> " << e2 << "." << p2 << " (";
               for ( int l = 0; l < c.rid.isize( ); l++ )
               {    if ( l > 0 ) cout << ",";
                    int id = c.rid[l].first, pos = c.rid[l].second;
                    Bool fw = ( id < (int) clo.size( ) );
                    cout << ( fw ? "+" : "-" ) << clo2_ids[id] << "." 
                         << clo2_n[id] << "." << pos;    }
               cout << ")\n";    }
          i = j - 1;    }

     // Write patches.

     cout << Date( ) << ": writing patches" << endl;
     patches.WriteAll( dir + "/new_stuff.exp" );
     cout << "\n" << Date( ) << ": done, time used = " << TimeSince(clock) << endl;
     Scram(0);    }
