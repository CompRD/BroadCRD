///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Basevector.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "paths/GetHomes.h"
#include "paths/UnibaseUtils.h"
#include "paths/Uniseq.h"
#include "util/SearchFastb2Core.h"

class vepath {

     public:

     vepath( ) { }
     vepath( const Bool vfirst, const vec< pair<int,int> >& x )
          : vfirst(vfirst), x(x) { }

     Bool vfirst;
     vec< pair<int,int> > x;

     friend ostream& operator<<( ostream& out, const vepath& p )
     {    out << "(";
          for ( int j = 0; j < p.x.isize( ); j++ )
          {    if ( j > 0 ) out << ",";
               Bool vert = ( (j%2 == 0) ^ !p.vfirst );
               if (vert) out << "v" << ToString( p.x[j].first );
               else 
               {    out << "e" << ToString( p.x[j].first ) << "."
                         << ToString( p.x[j].second );    }    }
          return out << ")";    }

};

void GetHomes( const String run_dir, const int K2, const vecbasevector& unibases2, 
     snark& S, const Bool VERBOSE )
{
     // Heuristics.

     const int max_dist = 10000;
     const int min_dist = 0;
     const int len = 20;
     const int K = 20;

     // Set up data.

     uniseq dummy;
     dummy.SetUnibases(unibases2);
     vec<int> to_rc;
     UnibaseInvolution( unibases2, to_rc );
     S.SetUnibases(unibases2);
     S.SetToRc(to_rc);

     // Make vector of everything in the snark.

     Mkdir777( run_dir + "/GetHomes" );
     String F1 = run_dir + "/GetHomes/jumps.fastb"; 
     String F2 = run_dir + "/GetHomes/all.fastb";
     vecbasevector all;
     vec<uniseq> uniseqs;
     vec< pair<int,int> > where;
     for ( int v = 0; v < S.VertN( ); v++ )
     {    all.push_back_reserve( S.Vert(v).Bases( ) );
          uniseqs.push_back( S.Vert(v) );
          where.push( v, 0 );    }
     for ( int e = 0; e < S.EdgeN( ); e++ )
     {    const gapster& g = S.Edge(e);
          for ( int j = 0; j < g.ClosureCount( ); j++ )
          {    all.push_back_reserve( g.Closure(j).Bases( ) );
               uniseqs.push_back( g.Closure(j) );
               where.push( e, j );    }    }
     all.WriteAll(F2);

     // Get jumps.

     PairsManager jpairs( run_dir + "/jump_reads_filt.pairs" );
     vecbasevector jbases( run_dir + "/jump_reads_filt.fastb" );
     int64_t nreads = jbases.size( );
     for ( size_t id = 0; id < jbases.size( ); id++ )
     {    ForceAssertGe( jbases[id].isize( ), len );
          jbases[id].SetToSubOf( jbases[id], jbases[id].isize( ) - len, len );    }
     jbases.WriteAll(F1);
     jbases.clear( );
     jbases.ReadAll( run_dir + "/jump_reads_filt.fastb" );

     // Align the reads.

     vec< triple<int64_t,int64_t,int> > aligns;
     SearchFastb2( F1, F2, K, &aligns, 0, -1, 0.90, False );

     // Screen alignments.

     ParallelSort(aligns);
     vec<Bool> aligns_to_delete0( aligns.size( ), False );
     vec<int> readlen( jbases.size( ), -1 );
     for ( int i = 0; i < aligns.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < aligns.isize( ); j++ )
               if ( aligns[j].first != aligns[i].first ) break;
          vec<int> ext(j-i, 0), ids( j-i, vec<int>::IDENTITY );
          for ( int k = i; k < j; k++ )
          {    const basevector& r = jbases[ aligns[k].first ];
               const basevector& a = all[ aligns[k].second ];
               int pos = aligns[k].third;
               Bool fw = ( pos >= 0 );
               if ( pos < 0 ) pos = -pos-1;

               // pos = position at which last K bases of r starts aligning to a

               if (fw)
               {    for ( int m = r.isize( ) - K - 1; m >= 0; m-- )
                    {    if ( pos + m - ( r.isize( ) - K ) < 0 ) break;
                         if ( r[m] != a[ pos + m - ( r.isize( ) - K ) ] ) break;
                         ext[k-i]++;    }    }
               else
               {    for ( int m = r.isize( ) - K - 1; m >= 0; m-- )
                    {    if ( pos + K + ( r.isize( ) - K - 1 ) - m >= a.isize( ) )
                              break;
                         if ( 3 - r[m] != a[ pos + K + ( r.isize( ) - K - 1 ) - m ] )
                              break;
                         ext[k-i]++;    }    }    }
          ReverseSortSync( ext, ids );
          readlen[ aligns[i].first ] = K + ext[0];
          for ( int k = 0; k < ext.isize( ); k++ )
          {    if ( ext[k] < ext[0] )
               {    for ( int r = k; r < ext.isize( ); r++ )
                         aligns_to_delete0[ i + ids[r] ] = True;
                    break;    }    }
          i = j - 1;    }
     EraseIf( aligns, aligns_to_delete0 );

     // Remove alignments to edges that could be on vertices instead.

     vec<Bool> aligns_to_delete( aligns.size( ), False );
     for ( int i = 0; i < aligns.isize( ); i++ )
     {    if ( aligns[i].second < S.VertN( ) ) continue;
          const uniseq& u = uniseqs[ aligns[i].second ];
          int pos = aligns[i].third;
          Bool fw = ( pos >= 0 );
          if ( pos < 0 ) pos = -pos-1;
          int Pos = pos;
          int ext = readlen[ aligns[i].first ] - K;
          if (fw) pos -= ext;
          else Pos += ext;
          int head_bases = unibases2[ u.U( ).front( ) ].size( );
          int tail_bases = unibases2[ u.U( ).back( ) ].size( );
          int all_bases = all[ aligns[i].second ].size( );
          if ( Pos <= head_bases || pos >= all_bases - tail_bases )
               aligns_to_delete[i] = True;    }
     EraseIf( aligns, aligns_to_delete );

     // Index the alignments.

     vec< vec<int> > aligns_index(nreads);
     for ( int i = 0; i < aligns.isize( ); i++ )
          aligns_index[ aligns[i].first ].push_back(i);

     // Get graph info.

     vec<int> to_left, to_right;
     S.G( ).ToLeft(to_left), S.G( ).ToRight(to_right);

     // Go through the pairs.
     
     longlong ndone = 0;;
     DPRINT( jpairs.nPairs( ) );
     vec< vec< vec< pair<int,int> > > > Rall;
     #pragma omp parallel for
     for ( size_t pid = 0; pid < jpairs.nPairs( ); pid++ )
     {    
           #pragma omp critical 
            { ndone++;
       
	      if ( ndone % (jpairs.nPairs()/100 +1 ) == 0 ) 
		Dot( cout, 100.0 * ndone/(double)jpairs.nPairs() );
	    }
	  //DPRINT3( pid, Sum(jpairs_done), jpairs.nPairs() );
	   // Dot( cout, 100 * Sum(jpairs_done)/jpairs.nPairs() );
          int64_t id1 = jpairs.ID1(pid), id2 = jpairs.ID2(pid);

          // Ignore pairs that do not have both ends placed.

          if ( aligns_index[id1].empty( ) || aligns_index[id2].empty( ) ) continue;
	  if ( aligns_index[id1].size() * aligns_index[id2].size() > 100 ) continue;

          // First test for valid placement on a vertex.

          /*
          Bool on_vertex = False;
          for ( int j1 = 0; j1 < aligns_index[id1].isize( ); j1++ )
          {    if (on_vertex) break;
               for ( int j2 = 0; j2 < aligns_index[id2].isize( ); j2++ )
               {    int p1 = aligns_index[id1][j1], p2 = aligns_index[id2][j2];
                    if ( aligns[p1].second != aligns[p2].second ) continue;
                    if ( aligns[p1].second >= S.VertN( ) ) continue;
                    if ( aligns[p1].third < 0 ) swap( p1, p2 );
                    if ( aligns[p1].third < 0 || aligns[p2].third >= 0 ) continue;
                    int m = aligns[p1].second;
                    int pos1 = aligns[p1].third, pos2 = -aligns[p2].third-1;
                    int Pos1 = pos1 + jbases[ aligns[p1].first ].isize( );
                    int sep = pos2 - Pos1;
                    if ( sep < 0 || sep > max_dist ) continue;
                    on_vertex = True;
                    break;    }    }
          if (on_vertex) continue;
          */

   
	  // Now explore other placements.

          Bool found = False;
          vec< vec< pair<int,int> > > R;
          for ( int j1 = 0; j1 < aligns_index[id1].isize( ); j1++ )
          for ( int j2 = 0; j2 < aligns_index[id2].isize( ); j2++ )
          {    int p1 = aligns_index[id1][j1], p2 = aligns_index[id2][j2];
               if ( aligns[p1].third < 0 ) swap( p1, p2 );
               if ( aligns[p1].third < 0 || aligns[p2].third >= 0 ) continue;
               int m1 = aligns[p1].second, m2 = aligns[p2].second;
               int pos1 = aligns[p1].third, pos2 = -aligns[p2].third-1;
	       int Pos1 = pos1 + jbases[ aligns[p1].first ].isize( );

               // Check for placements on same sequence that are too far apart.

               if ( m1 == m2 )
               {    int sep = pos2 - Pos1;
                    if ( sep < 0 || sep > max_dist ) continue;    }

               // Define vertices.

               int v1 = ( m1 < S.VertN( ) ? m1 : to_right[ where[m1].first ] );
               int v2 = ( m2 < S.VertN( ) ? m2 : to_left[ where[m2].first ] );

               // Find all plausible paths.

               vec<vepath> partials, fulls;
               vepath init;
               if ( m1 < S.VertN( ) ) init.vfirst = True;
               else
               {    init.vfirst = False;
                    init.x.push_back( where[m1] );    }
               init.x.push( v1, 0 );
               partials.push_back(init);
               while( partials.nonempty( ) && partials.size() < 1000 )
               {    vepath p = partials.back( );
                    partials.pop_back( );
                    if ( p.x.back( ).first == v2 )
                    {    vepath f(p);
                         if ( m2 >= S.VertN( ) ) f.x.push_back( where[m2] );
                         fulls.push_back(f);    }
                    int v = p.x.back( ).first;
                    for ( int j = 0; j < S.From(v).isize( ); j++ )
                    {    int w = S.From(v)[j];
                         const gapster& g = S.G( ).EdgeObjectByIndexFrom( v, j );
                         int ei = S.G( ).EdgeObjectIndexByIndexFrom( v, j );
                         for ( int m = 0; m < g.ClosureCount( ); m++ )
                         {    vepath q(p);
                              q.x.push( ei, m );
                              q.x.push( w, 0 );

			      int len   = 0; 
			      int start = 1;
			      int end   = q.x.isize(); 
			      // if second read lands on ei or w then 'end' position should change
			      if ( where[m2].first == w )       end = q.x.isize() -1;
			      else if ( where[m2].first == ei ) end = q.x.isize() -2;
			      
			      for ( int l = start; l < end; l += 1 )
                              {    Bool vert = ( (l%2 == 0) ^ !q.vfirst );
				   const uniseq& u
				     = ( vert ? S.Vert( q.x[l].first )
					 : S.Edge( q.x[l].first ).Closure(q.x[l].second) );
				   len += u.Len( ) - unibases2[ u.U(0) ].isize( );   
				   if ( l == end -1 )
				     len -= unibases2[ u.U( u.N()-1 ) ].isize();
			      }
                              
			      if ( len <= max_dist ) 
                                   partials.push_back(q);    }    }    }

	       if ( partials.nonempty() ) continue;


               // Screen again on separation.

               vec<Bool> fulls_to_delete( fulls.isize( ), False );
               vec<int> seps;
               for ( int j = 0; j < fulls.isize( ); j++ )
               {    const vepath& p = fulls[j];
                    int Pos1 = pos1 + jbases[ aligns[p1].first ].isize( );
                    int sep;
                    if ( m1 == m2 ) sep = pos2 - pos1;
                    else
                    {    sep = all[m1].isize( ) - Pos1;
                         for ( int l = 1; l < p.x.isize( ); l++ )
                         {    Bool vert = ( (l%2 == 0) ^ !p.vfirst );
                              const uniseq& u
                                   = ( vert ? S.Vert( p.x[l].first )
                                   : S.Edge( p.x[l].first ).Closure(p.x[l].second) );
                              sep += u.Len( ) - unibases2[ u.U(0) ].isize( );
                              if ( l == p.x.isize( ) - 1 )
                                   sep -= u.Len( ) - pos2;    }    }
                    if ( sep < min_dist || sep > max_dist ) 
                         fulls_to_delete[j] = True;
                    seps.push_back(sep);    }
               EraseIf( fulls, fulls_to_delete);
               if ( fulls.empty( ) ) continue;

               // Extract essential info.

               for ( int j = 0; j < fulls.isize( ); j++ )
               {    const vepath& p = fulls[j];
                    vec< pair<int,int> > E;
                    int start = ( p.vfirst ? 1 : 0 );
                    for ( int j = start; j < p.x.isize( ); j += 2 )
                    {    pair<int,int> e = p.x[j];
                         if ( S.Edge(e.first).ClosureCount( ) > 1 ) 
                              E.push_back(e);    }
                    UniqueSort(E);
                    if ( E.nonempty( ) ) R.push_back(E);    }

               // Print first.

               ostringstream pout;
               pout << "{";
               for ( int j = 0; j < fulls.isize( ); j++ )
               {    if ( j > 0 ) pout << ",";
                    pout << fulls[j];
                    pout << ":" << seps[j];    }
               pout << "}";
               String paths = pout.str( );
               found = True;    }

          if ( R.empty( ) ) continue;
          UniqueSort(R);

          // Test for no all of one edge, stupid case.

          int e = R[0][0].first;
          Bool all_e = True;
          for ( int j = 1; j < R.isize( ); j++ )
          {    if ( R[j].size( ) > 1 ) all_e = False;
               if ( R[j][0].first != e ) all_e = False;    }
          if (all_e)
          {    if ( R.isize( ) == S.Edge(e).ClosureCount( ) ) continue;    }

          // Save.

          #pragma omp critical
          {    Rall.push_back(R);    }

          // Print.

          if (VERBOSE)
          {
               #pragma omp critical
               {    cout <<  "\n";
                    PRINT2( id1, id2 );
                    for ( int j = 0; j < R.isize( ); j++ )
                    {    cout << "[" << j+1 << "]";
                         for ( int l = 0; l < R[j].isize( ); l++ )
                              cout << " e" << R[j][l].first << "." << R[j][l].second;
                         cout << "\n";    }    }    }    }

     if (VERBOSE) cout << "\n";

     // Cleaning.  If Rall[i] contains rows 
     // {ej.1,ei.r}, {ej.2,ei.r}, ..., {ej.n,ei.r}, where n = |ej|, then these rows
     // may be replaced by the single row {ei.r}.

     for ( int i = 0; i < Rall.isize( ); i++ )
     {    vec< vec< pair<int,int> > >& R = Rall[i];
          for ( int e = 0; e < S.EdgeN( ); e++ )
          {    int n = S.Edge(e).ClosureCount( );
               if ( n == 0 ) continue;
               vec< vec< pair<int,int> > > partners(n);
               for ( int j = 0; j < R.isize( ); j++ )
               {    if ( R[j].size( ) != 2 ) continue;
                    if ( R[j][0].first != e && R[j][1].first != e ) continue;
                    if ( R[j][0].first == e && R[j][1].first == e ) continue;
                    if ( R[j][0].first == e )
                         partners[ R[j][0].second ].push_back( R[j][1] );
                    if ( R[j][1].first == e )
                         partners[ R[j][1].second ].push_back( R[j][0] );    }
               for ( int j = 0; j < n; j++ )
                    UniqueSort( partners[j] );
               vec< pair<int,int> > pi;
               Intersection( partners, pi );
               vec<Bool> to_delete( R.size( ), False );
               for ( int j = 0; j < R.isize( ); j++ )
               {    if ( R[j].size( ) != 2 ) continue;
                    if ( R[j][0].first != e && R[j][1].first != e ) continue;
                    if ( R[j][0].first == e && R[j][1].first == e ) continue;
                    if ( R[j][0].first == e && BinMember( pi, R[j][1] ) )
                         to_delete[j] = True;
                    if ( R[j][1].first == e && BinMember( pi, R[j][0] ) )
                         to_delete[j] = True;    }
               EraseIf( R, to_delete );
               for ( int j = 0; j < pi.isize( ); j++ )
               {    vec< pair<int,int> > x;
                    x.push_back( pi[j] );
                    R.push_back(x);    }
               UniqueSort(R);    }    }

     // Cleaning.  If Rall[i] contains rows {ej.1}, {ej.2}, ..., {ej.n}, where
     // n = |ej|, then Rall[i] is always satisfied, and may be deleted.

     vec<Bool> to_delete( Rall.size( ), False );
     for ( int i = 0; i < Rall.isize( ); i++ )
     {    const vec< vec< pair<int,int> > >& R = Rall[i];
          Bool vacuous = False;
          for ( int e = 0; e < S.EdgeN( ); e++ )
          {    int n = S.Edge(e).ClosureCount( );
               if ( n == 0 ) continue;
               vec<Bool> have( n, False );
               for ( int l = 0; l < R.isize( ); l++ )
               {    if ( R[l].solo( ) && R[l][0].first == e )
                         have[ R[l][0].second ] = True;    }
               if ( Sum(have) == n ) vacuous = True;    }
          if (vacuous) to_delete[i] = True;    }
     cout << "deleting " << Sum(to_delete) << " vacuous cases" << endl;
     EraseIf( Rall, to_delete );

     // Cleaning.  Check for subsets.

     for ( int i = 0; i < Rall.isize( ); i++ )
     {    vec< vec< pair<int,int> > >& R = Rall[i];
          UniqueSort(R);
          vec<Bool> to_delete( R.size( ), False );
          for ( int j1 = 0; j1 < R.isize( ); j1++ )
          for ( int j2 = 0; j2 < R.isize( ); j2++ )
               if ( j2 != j1 && BinSubset( R[j2], R[j1] ) ) to_delete[j1] = True;
          EraseIf( R, to_delete );    }

     // Catalog results.

     Sort(Rall);
     vec< vec< vec< pair<int,int> > > > RR;
     vec<int> RR_mult;
     if (VERBOSE) cout << "\nRESULTS:\n";
     for ( int i = 0; i < Rall.isize( ); i++ )
     {    int j = Rall.NextDiff(i);
          const vec< vec< pair<int,int> > >& R = Rall[i];
          RR.push_back(R);
          RR_mult.push_back(j-i);
          if (VERBOSE)
          {    cout <<  "\ncount = " << j - i << "\n";
               for ( int j = 0; j < R.isize( ); j++ )
               {    cout << "[" << j+1 << "]";
                    for ( int l = 0; l < R[j].isize( ); l++ )
                         cout << " e" << R[j][l].first << "." << R[j][l].second;
                    cout << "\n";    }    }
          i = j - 1;    }

     // Do some primitive voting.

     for ( int e = 0; e < S.EdgeN( ); e++ )
     {    int n = S.Edge(e).ClosureCount( );
          if ( n <= 1 ) continue;
          const int min_count = 5;
          const int max_count = 1;
          vec<int> count1( n, 0 ), count2( n, 0 );
          for ( int i = 0; i < RR.isize( ); i++ )
          {    vec< pair<int,int> > I;
               Intersection( RR[i], I );
               Bool found = False;
               for ( int j = 0; j < I.isize( ); j++ )
               {    if ( I[j].first == e ) count1[ I[j].second ] += RR_mult[i];
                    found = True;    }
               if (found) continue;
               vec<Bool> hit( n, False );
               for ( int j = 0; j < RR[i].isize( ); j++ )
               {    for ( int k = 0; k < RR[i][j].isize( ); k++ )
                    {    if ( RR[i][j][k].first == e )
                              hit[ RR[i][j][k].second ] = True;    }    }
               for ( int j = 0; j < n; j++ )
                    if ( hit[j] ) count2[j] += RR_mult[i];    }
          if (VERBOSE)
          {    cout << "\nEDGE " << e << "\n";
               for ( int j = 0; j < n; j++ )
               {    cout << "(" << j << ")  " << count1[j] << "  " 
                         << count2[j] << "\n";    }    }    
          int best_count = Max(count1), best = -1, second_best_count = 0;
          for ( int j = 0; j < n; j++ )
               if ( count1[j] == best_count ) best = j;
          for ( int j = 0; j < n; j++ )
          {    if ( j != best ) 
                    second_best_count = Max( second_best_count, count1[j] );    }
          if ( best_count >= min_count && second_best_count <= max_count
               && Sum(count2) == 0 )
          {    vec<Bool> to_delete( n, False );
               for ( int j = 0; j < n; j++ )
               {    if ( count1[j] <= max_count ) 
                    {    to_delete[j] = True;
                         if (VERBOSE)
                         {    cout << "deleting closure " << j << " from edge "
                                   << e << endl;    }    }    }
               S.Gmutable( ).EdgeObjectMutable(e).RemoveSomeClosures(to_delete);
                    }    }

     // Clean up.

     S.SwallowSimpleGaps( );
     S.BringOutTheDead( );    }
