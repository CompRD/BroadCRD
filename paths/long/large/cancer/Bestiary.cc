///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Bestiary.  Classify tumor-only edges.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "FastIfstream.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "PrintAlignment.h"
#include "Set.h"
#include "kmers/KmerRecord.h"
#include "math/Functions.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/Lines.h"
#include "polymorphism/Edit.h"
#include "random/Random.h"
#include "random/Shuffle.h"

void PrintSomaticIds( const String out_dir, const String id, const vec<Bool>& s,
     const vec<int>& inv )
{    Ofstream( out, out_dir + "/somatic" + id );
     for ( int e = 0; e < s.isize( ); e++ )
          if ( s[e] ) out << e << " " << inv[e] << endl;    }

int main( )
{    RunTime( );

     // Define hardcoded paths.

     String work_dir = "/wga/scr4/jaffe/GapToy/49875.HCC1143+BL";
     String DIR = work_dir + "/a.fin";
     String fin_dir = DIR;
     String out_dir = "calls";
     Mkdir777(out_dir);

     // Load assembly.

     HyperBasevector hb;
     BinaryReader::readFile( DIR + "/a.s.hbv", &hb );
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     vec<int> inv;
     BinaryReader::readFile( DIR + "/a.s.inv", &inv );
     vec<vec<vec<vec<int>>>> lines;
     BinaryReader::readFile( DIR + "/a.s.lines", &lines );
     vec<int> tol, lens;
     GetTol( hb, lines, tol );
     GetLineLengths( hb, lines, lens );
     vec< vec< pair<int,int> > > hits;
     BinaryReader::readFile( fin_dir + "/a.aligns", &hits );

     // Load counts.

     vec<String> subsam_names;
     vec<vec<int>> count;
     if ( IsRegularFile( DIR + "/a.counts" ) )
     {    Ifstream( in, DIR + "/a.counts" );
          String line;
          getline( in, line );
          ForceAssert( !in.fail( ) );
          vec<String> x;
          Tokenize( line, x );
          for ( int j = 1; j < x.isize( ); j++ )
               subsam_names.push_back( x[j] );
          int ns = subsam_names.size( );
          count.resize( ns, vec<int>( hb.E( ) ) );
          for ( int e = 0; e < hb.E( ); e++ )
          {    getline( in, line );
               Tokenize( line, x );
               for ( int j = 0; j < ns; j++ )
                    count[j][e] = x[j+1].Int( );    }    }

// =================================================================================

     // Find somatic (T >> N) edges, excluding half of them.

     vec<Bool> somatic( hb.E( ), False );
     for ( int e = 0; e < hb.E( ); e++ )
     {    if ( count[0][e] >= 10 && count[1][e] == 0 && e <= inv[e] )
               somatic[e] = True;    }
     Ofstream( outall, out_dir + "/somatic" );
     for ( int e = 0; e < hb.E( ); e++ )
          if ( somatic[e] ) outall << e << " " << inv[e] << endl;

// =================================================================================

     // Define class 1 = cellular somatic edges.  These are somatic edges that 
     // occur in some cell of some line, for which the cell contains no other edge
     // having normal count zero, and for which there is only one path containing 
     // the somatic edge.

     vec<Bool> cellular_somatic( hb.E( ), False );
     for ( int i = 0; i < lines.isize( ); i++ )
     {    const vec<vec<vec<int>>>& L = lines[i];
          for ( int j = 1; j < L.isize( ); j += 2 )
          {    const vec<vec<int>>& x = L[j];

               // Find the somatic edges.

               vec<int> es;
               for ( int m = 0; m < x.isize( ); m++ )
               for ( int n = 0; n < x[m].isize( ); n++ )
               {    int e = x[m][n];
                    if ( somatic[e] ) es.push_back(e);    }
               UniqueSort(es);

               // Traverse the somatic edges.

               for ( int si = 0; si < es.isize( ); si++ )
               {    int e = es[si];

                    // Must have just one path for the somatic edge.

                    vec<int> ps;
                    for ( int m = 0; m < x.isize( ); m++ )
                         if ( Member( x[m], e ) ) ps.push_back(m);
                    if ( !ps.solo( ) ) continue;

                    // Must have a path through edges having normal support.

                    vec<int> nx;
                    for ( int m = 0; m < x.isize( ); m++ )
                    {    Bool normal = True;
                         for ( int n = 0; n < x[m].isize( ); n++ )
                         {    int e = x[m][n];
                              if ( count[1][e] == 0 ) normal = False;    }
                         if (normal) nx.push_back(m);    }
                    if ( nx.nonempty( ) ) cellular_somatic[e] = True;    }    }    }

     // Print class 1 edges.

     PrintSomaticIds( out_dir, "1", cellular_somatic, inv );

     // Generate info for class 1 edges.

     vec<uint64_t> shuffle;
     Shuffle64( lines.size( ), shuffle );
     vec< pair<int,String> > results;
     #pragma omp parallel for
     for ( int z = 0; z < lines.isize( ); z++ )
     {    
          int i = shuffle[z];

          const vec<vec<vec<int>>>& L = lines[i];
          for ( int j = 1; j < L.isize( ); j += 2 )
          {    const vec<vec<int>>& x = L[j];

               // Find the cellular somatic edges.

               vec<int> es;
               for ( int m = 0; m < x.isize( ); m++ )
               for ( int n = 0; n < x[m].isize( ); n++ )
               {    int e = x[m][n];
                    if ( cellular_somatic[e] ) es.push_back(e);    }
               UniqueSort(es);

               // Traverse the cellular somatic edges.

               for ( int si = 0; si < es.isize( ); si++ )
               {    int e = es[si];

                    // Find the path for the somatic edge, and the normal paths.

                    vec<int> ps, nx;
                    for ( int m = 0; m < x.isize( ); m++ )
                         if ( Member( x[m], e ) ) ps.push_back(m);
                    for ( int m = 0; m < x.isize( ); m++ )
                    {    Bool normal = True;
                         for ( int n = 0; n < x[m].isize( ); n++ )
                         {    int e = x[m][n];
                              if ( count[1][e] == 0 ) normal = False;    }
                         if (normal) nx.push_back(m);    }

                    // Basevectorize the paths.

                    int m = ps[0];
                    basevector tpath = hb.Cat( x[m] );
                    vec<basevector> npaths;
                    for ( int l = 0; l < nx.isize( ); l++ )
                         npaths.push_back( hb.Cat( x[ nx[l] ] ) );

                    // Align normal paths to tumor path.

                    int np = npaths.size( );
                    vec<vec<pair<int,edit0>>> edits(np);
                    vec<align> aligns(np);
                    for ( int l = 0; l < np; l++ )
                    {    alignment al;
                         SmithWatAffine( npaths[l], tpath, al );
                         align a = al;
                         aligns[l] = a;
                         edits[l] = AlignToEdits( a, npaths[l], tpath );    }

                    // Delete inferior paths.

                    for ( int l1 = 0; l1 < np; l1++ )
                    {    Bool best = True;
                         for ( int l2 = 0; l2 < np; l2++ )
                         {    if ( !Subset( edits[l1], edits[l2] ) ) 
                              {    best = False;
                                   break;    }    }
                         if (best)
                         {    if ( l1 > 0 )
                              {    npaths[0] = npaths[l1];
                                   edits[0] = edits[l1];    }
                              npaths.resize(1);
                              edits.resize(1);
                              aligns.resize(1);
                              break;    }    }
                    np = npaths.size( );

                    // Compute reverse alignments.

                    vec<vec<pair<int,edit0>>> editsr(np);
                    for ( int l = 0; l < np; l++ )
                    {    align a = aligns[l];
                         a.Flip( );
                         editsr[l] = AlignToEdits( a, tpath, npaths[l] );    }

                    // Print results.

                    vec<String> edx(np), edy(np);
                    for ( int l = 0; l < np; l++ )
                    {    vec<String> ed;
                         for ( int r = 0; r < edits[l].isize( ); r++ )
                         {    const edit0& d = edits[l][r].second;
                              String x = ToString( edits[l][r].first ) + ".";
                              if ( d.etype == INSERTION ) x += "ins." + d.seq;
                              if ( d.etype == DELETION ) x += "del." + ToString(d.n);
                              if ( d.etype == SUBSTITUTION ) x += "sub." + d.seq;
                              ed.push_back(x);    }
                         ostringstream o;
                         o << printSeq(ed);
                         edx[l] = o.str( );    }
                    for ( int l = 0; l < np; l++ )
                    {    vec<String> ed;
                         for ( int r = 0; r < editsr[l].isize( ); r++ )
                         {    const edit0& d = editsr[l][r].second;
                              String x = ToString( editsr[l][r].first ) + ".";
                              if ( d.etype == INSERTION ) x += "ins." + d.seq;
                              if ( d.etype == DELETION ) x += "del." + ToString(d.n);
                              if ( d.etype == SUBSTITUTION ) x += "sub." + d.seq;
                              ed.push_back(x);    }
                         ostringstream o;
                         o << printSeq(ed);
                         edy[l] = o.str( );    }
                    ostringstream out;
                    tpath.Print( out, "tumor." + ToString(e) 
                         + " normal_count=" + ToString( npaths.size( ) ) );
                    for ( int l = 0; l < npaths.isize( ); l++ )
                    {    npaths[l].Print( out, "normal." + ToString(e) + "." 
                              + ToString(l+1) + " changes_to_tumor=" + edx[l] 
                              + " changes_to_normal=" + edy[l] );    }    
                    #pragma omp critical
                    {    results.push( e, out.str( ) );    }    }    }    }

     // Print results.

     ParallelSort(results);
     Ofstream( out1b, out_dir + "/class1.txt" );
     for ( int i = 0; i < results.isize( ); i++ )
          out1b << results[i].second;
     PRINT( Sum(somatic) );
     PRINT( Sum(cellular_somatic) );
     cout << "cellular frac = "
          << PERCENT_RATIO( 3, Sum(cellular_somatic), Sum(somatic) ) << endl;

// =================================================================================

     // Define class 1b.

     const int max_path = 10;
     const int max_paths = 10000;
     vec<int> one( hb.E( ), 1 );
     digraphE<int> G( hb, one );
     vec< pair<int,String> > results_1b;
     vec<Bool> somatic1b( hb.E( ), False );
     #pragma omp parallel for
     for ( int e = 0; e < hb.E( ); e++ )
     {    if ( !somatic[e] ) continue;
          int v = to_left[e], w = to_right[e];
          int x = v, y = w;
          vec<int> ys = {y}, xs = {x};
          while( hb.To(y).solo( ) )
          {    if ( hb.From(y).size( ) != 1 ) break;
               y = hb.From(y)[0];
               if ( Member( ys, y ) ) break;
               ys.push_back(y);    }
          while( hb.From(x).solo( ) )
          {    if ( hb.To(x).size( ) != 1 ) break;
               x = hb.To(x)[0];
               if ( Member( xs, x ) ) break;
               xs.push_back(x);    }
          if ( hb.To(y).solo( ) || hb.From(x).solo( ) ) continue;
          vec<vec<int>> paths;
          Bool OK = G.AllPathsLengthRange( 
               v, w, 1, max_path, to_right, paths, max_paths );
          if ( !OK ) continue;
          int somatics = 0;
          for ( int i = 0; i < paths.isize( ); i++ )
          {    for ( int j = 0; j < paths[i].isize( ); j++ )
                    if ( count[1][ paths[i][j] ] == 0 ) somatics++;    }
          if ( paths.solo( ) || somatics > 1 ) continue;
          somatic1b[e] = True;

          vec<int> ps;
          for ( int m = 0; m < paths.isize( ); m++ )
               if ( Member( paths[m], e ) ) ps.push_back(m);

          int m = ps[0];
          basevector tpath = hb.Cat( paths[m] );
          vec<basevector> npaths;
          for ( int l = 0; l < paths.isize( ); l++ )
               if ( l != m ) npaths.push_back( hb.Cat( paths[l] ) );

          // Align normal paths to tumor path.

          int np = npaths.size( );
          vec<vec<pair<int,edit0>>> edits(np);
          vec<align> aligns(np);
          for ( int l = 0; l < np; l++ )
          {    alignment al;
               SmithWatAffine( npaths[l], tpath, al );
               align a = al;
               aligns[l] = a;
               edits[l] = AlignToEdits( a, npaths[l], tpath );    }

          // Delete inferior paths.

          for ( int l1 = 0; l1 < np; l1++ )
          {    Bool best = True;
               for ( int l2 = 0; l2 < np; l2++ )
               {    if ( !Subset( edits[l1], edits[l2] ) ) 
                    {    best = False;
                         break;    }    }
               if (best)
               {    if ( l1 > 0 )
                    {    npaths[0] = npaths[l1];
                         edits[0] = edits[l1];    }
                    npaths.resize(1);
                    edits.resize(1);
                    aligns.resize(1);
                    break;    }    }
          np = npaths.size( );

          // Compute reverse alignments.

          vec<vec<pair<int,edit0>>> editsr(np);
          for ( int l = 0; l < np; l++ )
          {    align a = aligns[l];
               a.Flip( );
               editsr[l] = AlignToEdits( a, tpath, npaths[l] );    }

          // Print results.

          vec<String> edx(np), edy(np);
          for ( int l = 0; l < np; l++ )
          {    vec<String> ed;
               for ( int r = 0; r < edits[l].isize( ); r++ )
               {    const edit0& d = edits[l][r].second;
                    String x = ToString( edits[l][r].first ) + ".";
                    if ( d.etype == INSERTION ) x += "ins." + d.seq;
                    if ( d.etype == DELETION ) x += "del." + ToString(d.n);
                    if ( d.etype == SUBSTITUTION ) x += "sub." + d.seq;
                    ed.push_back(x);    }
               ostringstream o;
               o << printSeq(ed);
               edx[l] = o.str( );    }
          for ( int l = 0; l < np; l++ )
          {    vec<String> ed;
               for ( int r = 0; r < editsr[l].isize( ); r++ )
               {    const edit0& d = editsr[l][r].second;
                    String x = ToString( editsr[l][r].first ) + ".";
                    if ( d.etype == INSERTION ) x += "ins." + d.seq;
                    if ( d.etype == DELETION ) x += "del." + ToString(d.n);
                    if ( d.etype == SUBSTITUTION ) x += "sub." + d.seq;
                    ed.push_back(x);    }
               ostringstream o;
               o << printSeq(ed);
               edy[l] = o.str( );    }
          ostringstream out;
          tpath.Print( out, "tumor." + ToString(e) 
               + " normal_count=" + ToString( npaths.size( ) ) );
          for ( int l = 0; l < npaths.isize( ); l++ )
          {    npaths[l].Print( out, "normal." + ToString(e) + "." 
                    + ToString(l+1) + " changes_to_tumor=" + edx[l] 
                    + " changes_to_normal=" + edy[l] );    }    
          #pragma omp critical
          {    results_1b.push( e, out.str( ) );    }    }

     // Print results.


     PrintSomaticIds( out_dir, "1b", somatic1b, inv );
     Ofstream( outx, out_dir + "/class1b.txt" );
     ParallelSort(results_1b);
     for ( int i = 0; i < results_1b.isize( ); i++ )
          outx << results_1b[i].second;
     PRINT( Sum(somatic) );
     PRINT( Sum(somatic1b) );
     cout << "somatic1b frac = "
          << PERCENT_RATIO( 3, Sum(somatic1b), Sum(somatic) ) << endl;

// =================================================================================

     // Define class 2.  Junk edges.

     fast_ifstream in( "/wga/dev/jaffe/BroadCRD/gloops" );
     vec<int> class2;
     String line;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          int e = line.Int( );
          if ( somatic[e] )
               cout << e << " " << inv[e] << endl;
          else if ( somatic[ inv[e] ] )
               cout << inv[e] << " " << e << endl;    }

// =================================================================================

     // Define class 3, decomposable edges.

     // Load remaining edges.

     vec<int> rem;
     fast_ifstream rin( out_dir + "/somatic_rem" );
     while(1)
     {    String line;
          getline( rin, line );
          if ( rin.fail( ) ) break;
          rem.push_back( line.Before( " " ).Int( ) );    }

     Bool use_special = False;

     const int L = 60;
     const int K = hb.K( );
     vec< triple<kmer<L>,int,int> > kmers_plus;
     vecbasevector edges( hb.E( ) );
     vec<int> special // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          = { 1968, 703184, 763876, 4752021, 4760960, 6508319, 6513391, // XXXXXXXXX
          7054649, 8076495, 8204258, 1613977, 1613978, // XXXXXXXXXXXXXXXXXXXXXXXXXX
          9640341, 6584678, 8918689, 9356639 }; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     special.push_back( 2759211,7325690,7866236,6228386,6228387,6228388, // XXXXXXXX
          8687610,9613618,2801113,6678651 ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     special.push_back( 1472378,1472379,7606973 ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     special.push_back( 1378146,1378154,1687546,7180023,7584781 ); // XXXXXXXXXXXXXX
     special.push_back( 436708,436709,7806742,574521,7931194,8112285 ); // XXXXXXXXX
     special.push_back( 8700500,4820655,8407469,8420668,9181468 ); // XXXXXXXXXXXXXX
     special.push_back( 8089374,8116879,8714574,490071,6454496,8093118 ); // XXXXXXX
     special.push_back( 8719350,8706446,5465770 ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     special.push_back( 230935,7466134,7489316,7998978,7695636 ); // XXXXXXXXXXXXXXX
     special.push_back( 41481,8733453,6433685 ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     // Define edge trims.

     vec<int> ltrim( hb.E( ), 0 ), rtrim( hb.E( ), 0 );
     for ( int e = 0; e < hb.E( ); e++ )
     {    int v = to_left[e], w = to_right[e];
          if ( hb.From(w).solo( ) ) rtrim[e] = K-L;
          if ( hb.To(v).solo( ) ) 
          {    if ( K-L + rtrim[e] <= hb.EdgeLengthBases(e) )
                    ltrim[e] = K-L;    }    }

     // Create lookup table.

     cout << Date( ) << ": making lookup table" << endl;
     vec<int64_t> starts;
     starts.push_back(0);
     for ( int e = 0; e < hb.E( ); e++ )
     {    const basevector& u = hb.EdgeObject(e);

          if (use_special) // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          {    if ( !Member( special, e ) ) // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               {    starts.push_back( starts.back( ) ); // XXXXXXXXXXXXXXXXXXXXXXXXX
                    continue;    }    } // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

          if ( count[1][e] == 0 ) 
          {    starts.push_back( starts.back( ) );
               continue;    }
          int nu = u.isize( ) - ltrim[e] - rtrim[e];
          starts.push_back( starts.back( ) + Max( 0, nu - L + 1 ) );    }
     kmers_plus.resize( starts.back( ) );
     #pragma omp parallel for
     for ( int e = 0; e < hb.E( ); e++ )
     {
          if ( use_special && !Member( special, e ) ) continue; // XXXXXXXXXXXXXXXXX

          if ( count[1][e] == 0 ) continue;
          const basevector& u = hb.EdgeObject(e);
          kmer<L> x;
          for ( int j = ltrim[e]; j <= u.isize( ) - L - rtrim[e]; j++ )
          {    int64_t r = starts[e] + j - ltrim[e];
               x.SetToSubOf( u, j ); 
               kmers_plus[r].first = x;
               kmers_plus[r].second = e; 
               kmers_plus[r].third = j;    }    }
     ParallelSort(kmers_plus);

     // Go through tumor-only edges.

     for ( int ir = 0; ir < rem.isize( ); ir++ )
     {    
          // Find matches.

          int e = rem[ir];
          int v = to_left[e], w = to_right[e];
          vec< triple<int,int,int> > P;

          if (use_special) // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          {    if ( e != 1970 && e != 8204258 && e != 4752020 && e != 763876 // XXXX
                    && e != 2759211 && e != 8706446 && e != 41481 ) // XXXXXXXXXXXXX
               {    continue;     }    } // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

          const basevector& E = hb.EdgeObject(e);
          Bool fail = False;
          for ( int j = 0; j <= E.isize( ) - L; j++ )
          {    kmer<L> x;
               x.SetToSubOf( E, j );
               int64_t low = LowerBound1( kmers_plus, x );
               int64_t high = UpperBound1( kmers_plus, x );

               vec<int64_t> candidates;
               Bool local = False;
               for ( int64_t m = low; m < high; m++ )
               {    candidates.push_back(m);
                    int f = kmers_plus[m].second;
                    if ( to_left[f] == v || to_right[f] == v ) local = True;
                    if ( to_left[f] == w || to_right[f] == w ) local = True;    }
               if (local)
               {    vec<Bool> to_delete( candidates.size( ), False );
                    for ( int64_t m = low; m < high; m++ )
                    {    Bool locm = False;
                         int f = kmers_plus[m].second;
                         if ( to_left[f] == v || to_right[f] == v ) locm = True;
                         if ( to_left[f] == w || to_right[f] == w ) locm = True;
                         if ( !locm ) to_delete[m-low] = True;    }
                    EraseIf( candidates, to_delete );     }

               if ( candidates.solo( ) )
               {    int64_t m = candidates[0];
                    P.push( j, kmers_plus[m].second, 
                         kmers_plus[m].third );    }    }

          // Collate matches.

          vec< quad<int,int,int,int> > Q; // (tstart, len, eid, estart)
          for ( int i = 0; i < P.isize( ); i++ )
          {    int l;
               for ( l = i + 1; l < P.isize( ); l++ )
               {    if ( P[l].second != P[i].second ) break;
                    int tshift = P[l].first - P[l-1].first;
                    int eshift = P[l].third - P[l-1].third;
                    if ( tshift != eshift || tshift > L ) break;    }
               Q.push( P[i].first, P[l-1].first - P[i].first + L,
                    P[i].second, P[i].third );
               i = l - 1;    }
          Bool cov = True;
          if ( Q.empty( ) || Q.front( ).first > 0 ) cov = False;
          else if ( Q.back( ).first + Q.back( ).second < E.isize( ) ) cov = False;
          if (cov)
          {    cout << "\n" << "[" << ir+1 << "]\n" << e << " " << inv[e] << ":\n";
               for ( int i = 0; i < Q.isize( ); i++ )
               {    if ( i > 0 )
                    {    int start1 = Q[i-1].first, start2 = Q[i].first; 
                         int len1 = Q[i-1].second, len2 = Q[i].second;
                         int stop1 = Q[i-1].first + len1, stop2 = Q[i].first + len2;
                         int e1 = Q[i-1].third, e2 = Q[i].third;
                         int gap = start2 - stop1;
                         int estart1 = Q[i-1].fourth, estop1 = Q[i-1].fourth + len1;
                         int estart2 = Q[i].fourth, estop2 = Q[i].fourth + len2;

                         if ( gap == -L+1 && to_right[e1] == to_left[e2]
                              && ( hb.EdgeLengthBases(e1) - estop1 ) + estart2 
                              == K - L );
                         else if ( gap <= 0 ) cout << "(" << gap << ")\n";    
                         else
                         {    cout << "(";
                              for ( int z = stop1; z < start2; z++ )
                                   cout << as_base( E[z] );
                              cout << ")\n";    }    }
                    cout << Q[i].third << "." << Q[i].fourth << "-"
                         << Q[i].fourth + Q[i].second << " [" << Q[i].second
                         << "]\n";    }    }    }

// =================================================================================

     // Class 4.  Try to exhibit residual somatic edges as substitutions.  For
     // a given somatic edge e: v --> w, consider paths going forward from v and
     // backwards from w.  Amongst these, find the paths that are at least as long
     // as e, and which define a substitution with e.  If there is only one such
     // substitution, we treat this as the normal sequence from which e originated
     // via a single somatic mutation.

     vec<int> res;
     fast_ifstream lin( out_dir + "/leftovers2" );
     while(1)
     {    String line;
          getline( lin, line );
          if ( lin.fail( ) ) break;
          res.push_back( line.Before( " " ).Int( ) );    }
     vec<Bool> somatic4( hb.E( ), False );
     for ( int ri = 0; ri < res.isize( ); ri++ )
     {    int e = res[ri];
          int v = to_left[e], w = to_right[e], ne = hb.EdgeLengthBases(e);
          vec<vec<int>> paths, paths_fin;
          for ( int j = 0; j < hb.From(v).isize( ); j++ )
          {    int f = hb.IFrom( v, j );
               if ( f == e ) continue;
               if ( count[1][f] == 0 ) continue;
               vec<int> x = {f};
               paths.push_back(x);    }
          while( paths.nonempty( ) )
          {    vec<vec<int>> paths_new;
               for ( int i = 0; i < paths.isize( ); i++ )
               {    const vec<int>& x = paths[i];
                    int n = hb.K( ) - 1;
                    for ( int j = 0; j < x.isize( ); j++ )
                         n += hb.EdgeLengthKmers( x[j] );
                    if ( n >= ne ) paths_fin.push_back(x);
                    else
                    {    int t = to_right[ x.back( ) ];
                         for ( int j = 0; j < hb.From(t).isize( ); j++ )
                         {    int f = hb.IFrom( t, j );
                              if ( hb.EdgeLengthBases(f) == 0 ) continue;
                              if ( count[1][f] == 0 ) continue;
                              vec<int> y(x);
                              y.push_back(f);
                              basevector b = hb.Cat(y);
                              int nb = b.size( );
                              int nsubs = 0;
                              for ( int i = hb.K( ) - 1; i < Min( ne, nb ); i++ )
                              {    if ( hb.EdgeObject(e)[i] != b[i] ) nsubs++;
                                   if ( nsubs > 1 ) break;    }
                              if ( nsubs <= 1 )
                                   paths_new.push_back(y);    }    }    }
               paths = paths_new;    }
          vec< pair<int,char> > subs1;
          for ( int j = 0; j < paths_fin.isize( ); j++ )
          {    basevector b = hb.Cat( paths_fin[j] );
               vec< pair<int,char> > subs;
               for ( int i = hb.K( ) - 1; i < ne; i++ )
               {    if ( hb.EdgeObject(e)[i] != b[i] ) subs.push( i, b[i] );
                    if ( subs.size( ) > 1 ) break;    }
               if ( subs.solo( ) ) subs1.push_back( subs[0] );    }
          paths_fin.clear( );
          for ( int j = 0; j < hb.To(w).isize( ); j++ )
          {    int f = hb.ITo( w, j );
               if ( f == e ) continue;
               if ( count[1][f] == 0 ) continue;
               vec<int> x = {f};
               paths.push_back(x);    }
          while( paths.nonempty( ) )
          {    vec<vec<int>> paths_new;
               for ( int i = 0; i < paths.isize( ); i++ )
               {    const vec<int>& x = paths[i];
                    int n = hb.K( ) - 1;
                    for ( int j = 0; j < x.isize( ); j++ )
                         n += hb.EdgeLengthKmers( x[j] );
                    if ( n >= ne ) paths_fin.push_back(x);
                    else
                    {    int t = to_left[ x.front( ) ];
                         for ( int j = 0; j < hb.To(t).isize( ); j++ )
                         {    int f = hb.ITo( t, j );
                              if ( hb.EdgeLengthBases(f) == 0 ) continue;
                              if ( count[1][f] == 0 ) continue;
                              vec<int> y = {f};
                              y.append(x);
                              basevector b = hb.Cat(y);
                              int nb = b.size( );
                              int nsubs = 0;
                              for ( int i = hb.K( ); i <= Min( ne, nb ); i++ )
                              {    if ( hb.EdgeObject(e)[ne-i] != b[nb-i] ) nsubs++;
                                   if ( nsubs > 1 ) break;    }
                              if ( nsubs <= 1 )
                                   paths_new.push_back(y);    }    }    }
               paths = paths_new;    }
          for ( int j = 0; j < paths_fin.isize( ); j++ )
          {    basevector b = hb.Cat( paths_fin[j] );
               int nb = b.size( );
               vec< pair<int,char> > subs;
               for ( int i = hb.K( ); i <= ne; i++ )
               {    if ( hb.EdgeObject(e)[ne-i] != b[nb-i] ) 
                         subs.push( ne-i, b[nb-i] );
                    if ( subs.size( ) > 1 ) break;    }
               if ( subs.solo( ) ) subs1.push_back( subs[0] );    }
          UniqueSort(subs1);
          if ( subs1.solo( ) )
          {    cout << e << " " << inv[e] << ": " << subs1[0].first << "." 
                    << as_base( subs1[0].second ) << endl;
               somatic4[e] = True;    }    }
     PrintSomaticIds( out_dir, "4", somatic4, inv );

// =================================================================================

     // Class 5.  Exclude immunoglobulin region that rearranges in normal, making
     // it impossible to know whether somatic changes occurred in the normal or the
     // tumor.

     vec<Bool> somatic5( hb.E( ), False );
     const int old_start = 106212273;
     const int old_stop  = 106518806;
     const int new_start = 105745936;
     const int new_stop  = 106062557;
     vec<int> bad_lines;
     for ( int e = 0; e < hits.isize( ); e++ )
     for ( int j = 0; j < hits[e].isize( ); j++ )
     {    if ( hits[e][j].first == 13
               && IntervalOverlap( old_start, old_stop,
               hits[e][j].second, hits[e][j].second + hb.EdgeLengthBases(e) ) > 0 )
          {    bad_lines.push_back( tol[e] );    }    }
     UniqueSort(bad_lines);
     for ( int i = 0; i < res.isize( ); i++ )
          if ( BinMember( bad_lines, tol[ res[i] ] ) ) somatic5[ res[i] ] = True;
     PrintSomaticIds( out_dir, "5", somatic5, inv );    }
