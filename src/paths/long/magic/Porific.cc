///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: dependency QueryLookupTable
// MakeDepend: dependency MakeLookupTable

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/SupportedHyperBasevector.h"

int main( )
{
     RunTime( );

     SupportedHyperBasevector shb;

     BinaryReader::readFile( "/wga/scr4/macro/ecoli/scs.shbv", &shb );
     vecbasevector bases( "/wga/scr4/macro/ecoli/ecoli.Feb7.1.fastb" );
     bases.ReadAll( "/wga/scr4/macro/ecoli/ecoli.Jan30.fastb", True );
     bases.ReadAll( "/wga/scr4/macro/ecoli/ecoli_plasmid_reads.fastb", True );
     bases.ReadAll( "/wga/scr4/macro/ecoli/mindthegaps_full.fastb", True );

     /*
     // r48688:LongProto SAMPLE=scardovia READS=#picard TMP=tmp.xxx OUT_INT_HEAD=...
     BinaryReader::readFile( "/wga/scr4/macro/scardo/scardo.shbv", &shb );
     vecbasevector bases( "/wga/scr4/macro/scardo/scardo.Feb7.1.fastb" );
     */

     vecbasevector all(bases);
     vecbasevector basesrc(bases);
     for ( int id = 0; id < (int) basesrc.size( ); id++ )
          basesrc[id].ReverseComplement( );
     all.Append(basesrc);
     for ( int e = 0; e < shb.EdgeObjectCount( ); e++ )
          all.push_back( shb.EdgeObject(e) );

     const int K = 20;
     vec< triple<kmer<K>,int,int> > kmers_plus;
     vec< vec< triple<int,int,int> > > hits( 2 * bases.size( ) );
     MakeKmerLookup0( all, kmers_plus );
     for ( int i = 0; i < (int) kmers_plus.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < (int) kmers_plus.isize( ); j++ )
               if ( kmers_plus[j].first != kmers_plus[i].first ) break;
          vec<int> g;
          for ( int k = i; k < j; k++ )
               if ( kmers_plus[k].second >= 2 * (int) bases.size( ) ) g.push_back(k);
          if ( g.size( ) == 1 )
          {    int ke = g[0];
               int e = kmers_plus[ke].second - ( 2 * bases.size( ) );
               int epos = kmers_plus[ke].third;
               for ( int k = i; k < j; k++ )
               {    if ( k == g[0] ) continue;
                    int id = kmers_plus[k].second;
                    int rpos = kmers_plus[k].third;
                    int offset = rpos - epos;
                    hits[id].push( e, rpos, offset );    }    }
          i = j - 1;    }

     vec<int> L;
     for ( int i = 0; i < shb.EdgeObjectCount( ); i++ )
          L.push_back( shb.EdgeLengthKmers(i) );
     digraphE<int> G( shb, L );
     vec<int> to_left, to_right;
     shb.ToLeft(to_left), shb.ToRight(to_right);

     vec< pair< vec<int>, int > > allpaths;
     #pragma omp parallel for
     for ( int id = 0; id < 2 * (int) bases.size( ); id++ )
     {    ostringstream out;
          vec< pair< vec<int>, int > > allpaths1, allpaths2;
          vec< triple<int,int,int> > x = hits[id];
          Sort(x);
          vec< triple<int,int,int> > places;
          vec<int> places_count;
          vec< pair<int,int> > lefts, rights; // (rpos,epos) at ends
          for ( int i = 0; i < x.isize( ); i++ )
          {    int e = x[i].first;
               int j;
               for ( j = i + 1; j < x.isize( ); j++ )
                    if ( x[j].first != e ) break;
               int rpos1 = x[i].second, epos1 = x[i].second - x[i].third;
               int rpos2 = x[j-1].second, epos2 = x[j-1].second - x[j-1].third;
               places.push( x[i].third, x[j-1].third + shb.EdgeLengthKmers(e), e );
               places_count.push(j-i);
               lefts.push( rpos1, epos1 );
               rights.push( rpos2, epos2 );
               i = j - 1;    }
          SortSync( places, places_count, lefts, rights );

          const int max_over = 1000;
          const int max_del = 5;
          const int min_mult = 1;
          for ( int pass = 1; pass <= 3; pass++ )
          {    vec<Bool> to_delete( places.size( ), False );
               for ( int i = 1; i < places.isize( ); i++ )
               {    if ( IntervalOverlap( places[i-1].first, places[i-1].second,
                         places[i].first, places[i].second ) > max_over )
                    {    if ( places_count[i-1] <= max_del
                              && places_count[i] >= min_mult * places_count[i-1] )
                         {    to_delete[i-1] = True;    }
                         if ( places_count[i] <= max_del
                              && places_count[i-1] >= min_mult * places_count[i] )
                         {    to_delete[i] = True;    }    }    }
               EraseIf( places, to_delete ), EraseIf( places_count, to_delete );
               EraseIf( lefts, to_delete ), EraseIf( rights, to_delete );    }

          if ( places.solo( ) ) continue;

          out << "\nhits for read " << id << endl;
          const int ladd    = 50;
          const double lmul = 1.5;
          const int max_paths = 100;
          allpaths2.push( vec<int>( ), 1 );
          for ( int i = 0; i < places.isize( ); i++ )
          {    if ( i > 0 )
               {    int e1 = places[i-1].third, e2 = places[i].third;
                    int v = to_right[e1], w = to_left[e2];
                    int d = Max( 0, places[i].first - places[i-1].second );
                    int L1 = int( floor( double(d) / lmul ) ) - ladd;
                    int L2 = int( ceil( double(d) * lmul ) ) + ladd;
                    vec< vec<int> > paths;
                    G.AllPathsLengthRange( v, w, L1, L2, to_right, paths );
                    if ( paths.isize( ) > max_paths ) paths.clear( );
                    int rpos1 = rights[i-1].first, epos1 = rights[i-1].second;
                    int rpos2 = lefts[i].first, epos2 = lefts[i].second;
                    if ( paths.size( ) == 1 )
                    {    for ( int l = 0; l < paths[0].isize( ); l++ )
                              out << "[edge " << paths[0][l] << "]\n";
                         for ( int m = 0; m < allpaths2.isize( ); m++ )
                              allpaths2[m].first.append( paths[0] );    }
                    else if ( paths.empty( ) ) 
                    {    out << "see 0 paths\n";
                         allpaths1.append(allpaths2);
                         allpaths2.clear( );
                         allpaths2.push( vec<int>( ), 1 );    }
                    else 
                    {    vec<int> errs( paths.size( ) );
                         for ( int j = 0; j < paths.isize( ); j++ )
                         {    int rstart = rpos1 + K, rstop = rpos2;
                              int lext = 0, rext = 0;
                              if ( rstart > rstop )
                              {    lext = (rstart - rstop)/2;
                                   rext = rstart - rstop - lext;    }
                              rstart -= lext;
                              rstop += rext;
                              int estart = epos1 + K - lext, estop = epos2 + rext;
                              int l1 = shb.EdgeLengthBases(e1) - estart, l2 = estop;
                              if ( l1 < shb.K( ) )
                              {    estart -= shb.K( ) - l1;
                                   rstart -= shb.K( ) - l1;    }
                              if ( l2 < shb.K( ) )
                              {    estop += shb.K( ) - l2;
                                   rstop += shb.K( ) - l2;    }
                              int rlen = rstop - rstart;
                              int elen = shb.EdgeLengthBases(e1) - estart;
                              for ( int l = 0; l < paths[j].isize( ); l++ )
                                   elen += shb.EdgeLengthKmers( paths[j][l] );
                              elen += estop - ( shb.K( ) - 1 );
                              if ( elen < rlen )
                              {    estart -= (rlen-elen)/2;
                                   estop += (rlen-elen) - (rlen-elen)/2;    }
                              basevector r( all[id], rstart, rstop - rstart );
                              basevector e( shb.EdgeObject(e1), estart,
                                   shb.EdgeLengthBases(e1) - estart );
                              for ( int l = 0; l < paths[j].isize( ); l++ )
                              {    e.resize( e.isize( ) - ( shb.K( ) - 1 ) );
                                   e = Cat( e, shb.EdgeObject( paths[j][l] ) );    }
                              e = Cat( e, basevector( shb.EdgeObject(e2), 
                                   shb.K( ) - 1, estop - ( shb.K( ) - 1 ) ) );
                              int best_loc;
                              alignment a;
                              errs[j] = SmithWatFree( r, e, best_loc, a );    }
                         SortSync( errs, paths );
                         out << "see " << paths.size( ) << " paths" << endl;    
                         for ( int j = 0; j < paths.isize( ); j++ )
                         {    if ( j > 0 && errs[j] > errs[j-1] )
                              {    paths.resize(j);
                                   break;    }    }
               
                         vec< pair< vec<int>, int > > allpaths2x;
                         for ( int j = 0; j < paths.isize( ); j++ )
                         {    out << "path = " << printSeq( paths[j] )
                                   << ", errs = " << errs[j] << "\n";    
                              for ( int l = 0; l < allpaths2.isize( ); l++ )
                              {    vec<int> p = allpaths2[l].first;
                                   p.append( paths[j] );
                                   allpaths2x.push( p, allpaths2.size( ) 
                                        * paths.size( ) );    }    }
                         allpaths2 = allpaths2x;   }    }

               int ncov = IntervalOverlap( 0, all[id].isize( ),
                    places[i].first, places[i].second );
               double cov = double(places_count[i]) / double(ncov);
               for ( int l = 0; l < allpaths2.isize( ); l++ )
                    allpaths2[l].first.push_back( places[i].third );
               out << "edge " << places[i].third << ": from " << places[i].first
                    << " to " << places[i].second << " [" << places_count[i] 
                    << ",cov=" << cov << "]" << "\n";    }
          #pragma omp critical
          {    cout << out.str( );    
               allpaths.append(allpaths1);
               allpaths.append(allpaths2);    }    }

     Sort(allpaths);
     cout << "\nPATHS:\n";
     int count = 0;
     for ( int i = 0; i < allpaths.isize( ); i++ )
     {    if ( allpaths[i].first.size( ) <= 1 ) continue;
          int j;
          for ( j = i + 1; j < allpaths.isize( ); j++ )
               if ( allpaths[j].first != allpaths[i].first ) break;
          double weight = 0;
          for ( int k = i; k < j; k++ )
               weight += 1 / double( allpaths[k].second );
          cout << "[" << ++count << "] " << printSeq( allpaths[i].first )
               << " (" << weight << ")\n";
          i = j - 1;    }    }
