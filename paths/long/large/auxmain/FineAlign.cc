///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Starting from seed placements of assembly edges on a reference sequence, try to
// align all the edges.

// NA12878 total edges = 10,592,170; time = 25.9 minutes (crd29), 8.9% unaligned

// 50919:
//
// LINE = 34882;
// LINE = 34880;
// LINE = 7187; // 200 kb
// LINE = 606; // 500 kb --> 9.23 seconds
// LINE = 13; // 1 Mb --> 23.9 seconds
// LINE = 10; // 1 Mb --> 57.8 seconds (current)

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "PrintAlignment.h"
#include "math/Functions.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "paths/HyperBasevector.h"

class align_data {

     public:

     align_data( ) { }
     align_data( const int g, const align& a, const int offset, const int bandwidth,
          const int start, const int spread )
          : g(g), a(a), offset(offset), bandwidth(bandwidth), start(start), 
          spread(spread) { }
     int g;
     align a;
     int offset;
     int bandwidth;
     int start;
     int spread;

};

Bool Clean( const align& a, const basevector& b1, const basevector& b2,
     const int score )
{    if ( score > 100 ) return False;
     vec<ho_interval> perf1;
     a.PerfectIntervals1( b1, b2, perf1 );
     const int min_flank = 60;
     return ( perf1.nonempty( ) && perf1[0].Start( ) == 0
          && perf1[0].Length( ) >= min_flank && perf1.back( ).Length( ) >= min_flank
          && perf1.back( ).Stop( ) == b1.isize( ) );    }

// Align a path.  This assumes each edge in the path is uniquely aligned to the
// same genome contig.  CURRENTLY: doesn't align, just returns offset and
// bandwidth.

void AlignPath( const vec<int>& p, const HyperBasevectorX& hb, 
     const vecbasevector& genome, const vec<vec<align_data>>& aligns,
     int& offset, int& bandwidth )
{    vec<int> lows, highs;
     for ( int i = 0; i < p.isize( ); i++ )
     {    int e = p[i];
          int offset = aligns[e][0].offset;
          for ( int j = 0; j < i; j++ )
               offset -= hb.Kmers( p[j] );
          int bandwidth = aligns[e][0].bandwidth;
          lows.push_back( offset - bandwidth );
          highs.push_back( offset + bandwidth );    }
     int low = Min(lows), high = Max(highs);
     offset = low + (high-low)/2;
     bandwidth = Max( offset - low, high - offset );    }

// Try to decompose an alignment of a path.  Whether we succeed depends in part
// on how stringently we interpret the definition of an align.

Bool DecomposeAlign( const vec<int>& p, const align& a,
     const HyperBasevectorX& hb, const basevector& G, vec<align>& d )
{    d.clear( );
     int start1 = 0;
     for ( int i = 0; i < p.isize( ); i++ )
     {    int e = p[i];
          int stop1 = start1 + hb.Bases(e);
          int pos1 = 0, pos2 = a.pos2( ), align_pos2 = 0;
          vec<int> gaps, lengths;
          for ( int j = 0; j < a.Nblocks( ); j++ )
          {
               if ( a.Gaps(j) > 0 ) // gap on first sequence
               {    if ( pos1 > start1 ) gaps.push_back( a.Gaps(j) );
                    pos2 += a.Gaps(j);    }

               if ( a.Gaps(j) < 0 ) // gap on second sequence
               {    if ( pos1 >= start1 ) 
                    {    if ( pos1 - a.Gaps(j) >= stop1 ) return False;
                         if ( gaps.empty( ) ) return False;
                         gaps.push_back( a.Gaps(j) );    }
                    else if ( pos1 - a.Gaps(j) > start1 ) return False;
                    pos1 -= a.Gaps(j);    }

               if ( pos1 + a.Lengths(j) > start1 )
               {    if ( gaps.empty( ) ) 
                    {    gaps.push_back(0);
                         align_pos2 = pos2;
                         if ( pos1 < start1 ) align_pos2 += start1 - pos1;    }

                    int l = a.Lengths(j);
                    if ( pos1 < start1 ) l -= ( start1 - pos1 );
                    if ( pos1 + a.Lengths(j) > stop1 )
                         l -= pos1 + a.Lengths(j) - stop1;

                    lengths.push_back(l);    }
               pos1 += a.Lengths(j), pos2 += a.Lengths(j);

               if ( pos1 >= stop1 ) break;    }
     
          ForceAssertEq( gaps.size( ), lengths.size( ) );
          avector<int> agaps( gaps.size( ) ), alengths( lengths.size( ) );
          for ( int j = 0; j < gaps.isize( ); j++ )
          {    agaps(j) = gaps[j], alengths(j) = lengths[j];    }
          d.push( 0, align_pos2, agaps, alengths );
          start1 += hb.Kmers(e);    }
     return True;    }

Bool Consistent( const align& a1, const align& a2, const int K )
{    int count1 = 0, count2 = 0;
     vec<String> s1, s2;
     for ( int m = a1.Nblocks( ) - 1; m >= 0; m-- )
     {    int g = a1.Gaps(m), l = a1.Lengths(m);
          if ( count1 + l > K - 1 ) l = K - 1 - count1;
          s1.push_back( "l=" + ToString(l) );
          count1 += l;
          if ( count1 == K - 1 ) break;
          if ( g < 0 && count1 - g > K - 1 ) g = -( K - 1 - count1 );
          if ( g != 0 ) s1.push_back( "g=" + ToString(g) );
          if ( g < 0 ) count1 -= g;
          if ( count1 == K - 1 ) break;    }
     s1.ReverseMe( );
     for ( int m = 0; m < a2.Nblocks( ); m++ )
     {    int g = a2.Gaps(m), l = a2.Lengths(m);
          if ( g < 0 && count2 - g > K - 1 ) g = -( K - 1 - count2 );
          if ( g != 0 ) s2.push_back( "g=" + ToString(g) );
          if ( g < 0 ) count2 -= g;
          if ( count2 == K - 1 ) break;
          if ( count2 + l > K - 1 ) l = K - 1 - count2;
          s2.push_back( "l=" + ToString(l) );
          count2 += l;
          if ( count2 == K - 1 ) break;    }

     int n = 0;
     for ( int i = 0; i < s1.isize( ); i++ )
          if ( s1[i].Contains( "l=", 0 ) ) n += s1[i].After( "l=" ).Int( );
     else
     {    int g = s1[i].After( "g=" ).Int( );
          if ( g > 0 ) n += g;    }

     // cout << "s1 = " << printSeq(s1) << endl;
     // cout << "s2 = " << printSeq(s2) << endl;
     // PRINT3( a1.Pos2( ) - a2.pos2( ), n, K - 1 );

     if ( s1 != s2 ) return False;
     if ( a1.Pos2( ) - a2.pos2( ) != n ) return False;
     // if ( a1.Pos2( ) - a2.pos2( ) != K - 1 ) return False;
     // cout << "Passes!\n";
     return True;    }

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault_Doc(DIR, ".",
          "looks for DIR/a.{fastb,inv,paths,paths.inv}");
     CommandArgument_Int_OrDefault_Doc(LINE, -1,
          "if specified, process only this line");
     CommandArgument_Bool_OrDefault(ALTREF, False);
     EndCommandArguments;

     // Heuristics.

     const int max_depth = 20;
     const int depth_add = 2;
     const int max_nhood = 1000;
     const int max_offset_diff = 200;
     const int bw_fudge = 3;
     const int score_mult = 10;
     const int max_bandwidth = 1000;

     // Load assembly.

     HyperBasevectorX hb;
     BinaryReader::readFile( DIR + "/a.hbx", &hb );
     vec<int> inv;
     BinaryReader::readFile( DIR + "/a.inv", &inv );
     vec< vec< pair<int,int> > > hits;
     String suffix = ( ALTREF ? "_alt" : "" );
     BinaryReader::readFile( DIR + "/a.aligns" + suffix, &hits );
     vec<vec<vec<vec<int>>>> lines;
     BinaryReader::readFile( DIR + "/a.lines", &lines );

     // Load genome.

     vecbasevector genome;
     if (ALTREF) genome.ReadAll( DIR + "/../genome.fastb_alt" );
     else genome.ReadAll( DIR + "/../genome.fastb" );
     int ng = genome.size( );

     // Create mirror rc hits.

     {    vecbasevector genome_rc(genome);
          for ( int g = 0; g < ng; g++ )
               genome_rc[g].ReverseComplement( );
          genome.Append(genome_rc);    }
     for ( int e = 0; e < hb.E( ); e++ )
     {    int re = inv[e];
          for ( int i = 0; i < hits[e].isize( ); i++ )
          {    int g = hits[e][i].first, gpos = hits[e][i].second;
               if ( g >= ng ) continue;
               int rgpos = genome[g].isize( ) - ( gpos + hb.Bases(e) );
               hits[re].push( g + ng, rgpos );    }    }

     // Make a list of edges to use.

     PRINT( hb.E( ) );
     vec<int> es;
     if ( LINE >= 0 )
     {    for ( int j = 0; j < lines[LINE].isize( ); j++ )
          for ( int m = 0; m < lines[LINE][j].isize( ); m++ )
          for ( int n = 0; n < lines[LINE][j][m].isize( ); n++ )
               es.push_back( lines[LINE][j][m][n] );
          Sort(es);    }
     else
     {    for ( int e = 0; e < hb.E( ); e++ )
               es.push_back(e);    }

     // Propagate locs.

     double clock = WallClockTime( );
     vec< vec<align_data> > aligns( hb.E( ) );
     vec< vec<int> > score( hb.E( ) );
     vec<String> reports( es.size( ) );
     #pragma omp parallel for schedule(dynamic, 100)
     for ( int ie = 0; ie < es.isize( ); ie++ )
     {    int e = es[ie];
          if ( hb.Bases(e) == 0 ) continue;
          int re = inv[e];
          if ( re < e && BinMember( es, re ) ) continue;
          double eclock = WallClockTime( );
          ostringstream out;

          vec< pair<int,int> > locs = hits[e];
          if ( locs.empty( ) )
          {    vec< pair<int,int> > nhood;
               nhood.push( e, 0 );
               int md = max_depth;
               for ( int d = 1; d <= md; d++ )
               {    vec< pair<int,int> > nhood2;
                    for ( int i = 0; i < nhood.isize( ); i++ )
                    {    int f = nhood[i].first, d = nhood[i].second;
                         nhood2.push( f, d );
                         int v = hb.ToLeft(f), w = hb.ToRight(f);
                         for ( int j = 0; j < (int) hb.From(v).size( ); j++ )
                              nhood2.push( hb.IFrom(v,j), d );
                         for ( int j = 0; j < (int) hb.To(v).size( ); j++ )
                         {    int g = hb.ITo( v, j );
                              nhood2.push( g, d - hb.Kmers(g) );    }
                         for ( int j = 0; j < (int) hb.From(w).size( ); j++ )
                              nhood2.push( hb.IFrom(w,j), d + hb.Kmers(f) );
                         for ( int j = 0; j < (int) hb.To(w).size( ); j++ )
                         {    int g = hb.ITo( w, j );
                              nhood2.push( 
                                   g, d + hb.Kmers(f) - hb.Kmers(g) );    }    }
                    UniqueSort(nhood2);
                    if ( nhood2.isize( ) > max_nhood ) break;
                    nhood = nhood2;
                    int old_locs_size = locs.size( );
                    locs.clear( );
                    for ( int i = 0; i < nhood.isize( ); i++ )
                    {    int f = nhood[i].first, d = nhood[i].second;
                         for ( int r = 0; r < hits[f].isize( ); r++ )
                         {    int g = hits[f][r].first, gpos = hits[f][r].second - d;
                              locs.push( g, gpos );    }    }
                    UniqueSort(locs);
                    if ( old_locs_size == 0 && locs.nonempty( ) ) 
                         md = d + depth_add;    }    }

          // Group the locs.

          Sort(locs);
          int count = 0;
          for ( int i = 0; i < locs.isize( ); i++ )
          {    int j, g = locs[i].first;
               for ( j = i + 1; j < locs.isize( ); j++ )
               {    if ( locs[j].first != g ) break;
                    if ( locs[j].second - locs[j-1].second > max_offset_diff )
                         break;    }
               int start = locs[i].second;
               int spread = locs[j-1].second - start;

               // Calculate offset and bandwidth, then align.

               int offset = locs[i].second + (spread+1)/2;
               int bandwidth 
                    = Max( offset - start, locs[j-1].second - offset ) + bw_fudge;
               align a;
               int errors;
               int err = SmithWatAffineBanded( hb.EdgeObject(e), genome[g], -offset,
                    bandwidth, a, errors );
               aligns[e].push( g, a, offset, bandwidth, start, spread );
               score[e].push_back(err);
               i = j - 1;    }

          // Delete inferior alignments, defines as those whose score is 
          // significantly greater than the best score.

          SortSync( score[e], aligns[e] );
          for ( int i = 1; i < score[e].isize( ); i++ )
          {    if ( score[e][i] > score_mult * ( score[e][0] + 1 ) )
               {    aligns[e].resize(i), score[e].resize(i);
                    break;    }    }

          // Copy to rc.

          if ( re != e && BinMember( es, re ) )
          {    score[re] = score[e], aligns[re] = aligns[e];
               for ( int i = 0; i < aligns[e].isize( ); i++ )
               {    int g = aligns[e][i].g;
                    if ( g >= ng ) aligns[re][i].g -= ng;
                    else aligns[re][i].g += ng;
                    int ne = hb.Bases(e);
                    aligns[re][i].offset = genome[g].isize( ) 
                         - ( aligns[e][i].offset + ne );
                    aligns[re][i].start = genome[g].isize( ) 
                         - ( aligns[e][i].start + aligns[e][i].spread + ne );
                    aligns[re][i].a.ReverseThis( ne, genome[g].size( ) );    }    }

          // Print alignments.

          out << "\nplacements of edge " << e;
          if ( aligns[e].solo( ) && Clean( aligns[e][0].a, 
               hb.EdgeObject(e), genome[ aligns[e][0].g ], score[e][0] ) )
          {    int g = aligns[e][0].g, start = aligns[e][0].start; 
               int spread = aligns[e][0].spread;
               out << " --> " << g << "." << start << "+" << spread;
               if ( score[e][0] == 0 ) out << " (perfect)\n";
               else out << " (score=" << score[e][0] << ")\n";    }
          else
          {    out << endl;
               for ( int i = 0; i < aligns[e].isize( ); i++ )
               {    int g = aligns[e][i].g, start = aligns[e][i].start; 
                    int spread = aligns[e][i].spread;
                    const align& a = aligns[e][i].a;
                    int errs = score[e][i];
                    Bool clean = Clean( a, hb.EdgeObject(e), genome[g], errs );
                    out << "[" << i+1 << "] " << g << "." << start << "+" << spread;
                    if ( errs == 0 ) out << " (perfect)\n";
                    else 
                    {    out << " (score=" << errs << ")";
                         if ( errs > 400 ) out << " Looks like garbage." << endl;
                         else if ( !clean )
                         {    PrintVisualAlignment( True, out, 
                                   hb.EdgeObject(e), genome[g], a );    }    
                         else out << "\n";    }    }    }
          if ( WallClockTime( ) - eclock > 2 )
               out << "time for edge " << e << " = " << TimeSince(eclock) << endl;
          reports[ie] = out.str( );    }
     cout << "\n" << TimeSince(clock) << " used in main loop" << endl;
     for ( int i = 0; i < reports.isize( ); i++ )
          cout << reports[i];

     // List unaligned edges.

     vec<int> unaligned;
     for ( int e = 0; e < hb.E( ); e++ )
     {    if ( !BinMember( es, e ) ) continue;
          if ( aligns[e].empty( ) ) unaligned.push_back(e);    }
     cout << "\n" << unaligned.size( ) << " unaligned edges" << " (" 
          << PERCENT_RATIO( 3, unaligned.size( ), es.size( ) ) << ")" << endl;
     if ( unaligned.nonempty( ) && unaligned.size( ) < 100 )
          cout << "unaligned: " << printSeq(unaligned) << "\n" << endl;

     // Find pairs of inconsistent alignments.

     vec< pair<int,int> > bad_pairs;
     for ( int v = 0; v < hb.N( ); v++ )
     for ( int i = 0; i < (int) hb.To(v).size( ); i++ )
     for ( int j = 0; j < (int) hb.From(v).size( ); j++ )
     {    int e1 = hb.ITo(v,i), e2 = hb.IFrom(v,j);
          if ( !BinMember( es, e1 ) ||  !BinMember( es, e2 ) ) continue;

          // NOT RIGHT!
          if ( !aligns[e1].solo( ) || !aligns[e2].solo( ) ) continue;

          int g1 = aligns[e1][0].g, g2 = aligns[e2][0].g;
          const align &a1 = aligns[e1][0].a, &a2 = aligns[e2][0].a;
          // cout << "\nTesting consistency...\n";
          // PRINT2( e1, e2 );
          if ( g1 != g2 || !Consistent( a1, a2, hb.K( ) ) )
               bad_pairs.push( e1, e2 );    }
     cout << "found " << bad_pairs.size( ) << " bad pairs" << endl;

     // Align the inconsistent pairs.

     vec<String> reports2( bad_pairs.size( ) );
     vec< vec< pair<int,align> > > aligns2( hb.E( ) );
     #pragma omp parallel for schedule(dynamic, 50)
     for ( int i = 0; i < bad_pairs.isize( ); i++ )
     {    int e1 = bad_pairs[i].first, e2 = bad_pairs[i].second;
          int re1 = inv[e1], re2 = inv[e2];
          if ( make_pair( re2, re1 ) < make_pair( e1, e2 ) 
               && BinMember( es, re1 ) && BinMember( es, re2 ) )
          {    continue;    }
          ostringstream out;
          vec<int> p = {e1, e2};
          int offset, bandwidth;
          AlignPath( p, hb, genome, aligns, offset, bandwidth );
          int g = aligns[e1][0].g;
          out << "\nalignment of path " << printSeq(p) << " to " << g << endl;
          basevector b = hb.Cat(p);
          align a;
          int errors;
          if ( bandwidth > max_bandwidth )
          {    PRINT_TO( out, bandwidth );
               out << "Bandwidth too large, giving up." << endl;    }
          else
          {    int errx = SmithWatAffineBanded( 
                    b, genome[g], -offset, bandwidth, a, errors );
               out << "score = " << errx << endl;
               if ( errx > 400 ) out << "Looks like garbage." << endl;
               else 
               {    if ( !Clean( a, b, genome[g], errx ) )
                    {    PrintVisualAlignment( True, out, b, 
                              genome[g], a );    }    
                    /*
                    for ( int m = 0; m < a.Nblocks( ); m++ ) // XXXXXXXXXXXXXXXXXXXX
                         PRINT3_TO( out, m, a.Gaps(m), a.Lengths(m) ); // XXXXXXXXXX
                    */

                    vec<align> d;
                    Bool ok = DecomposeAlign( p, a, hb, genome[g], d );
                    if ( !ok ) 
                    {    out << "failed to decompose" << endl;
                         // for ( int j = 0; j < p.isize( ); j++ )
                         //      PRINT2( j, hb.Bases( p[j] ) );    
                              }
                    else
                    {    
                         #pragma omp critical
                         {    for ( int j = 0; j < p.isize( ); j++ )
                              {    aligns2[ p[j] ].push( g, d[j] );    }    }
                         out << "Decomposition:\n";
                         for ( int j = 0; j < p.isize( ); j++ )
                         {    out << "\nalign of " << p[j] << endl;
                              /*
                              PRINT_TO( out, hb.Bases( p[j] ) );
                              PRINT2_TO( out, d[j].pos1( ), d[j].pos2( ) );
                              for ( int m = 0; m < d[j].Nblocks( ); m++ )
                              {    PRINT3_TO( out, m, d[j].Gaps(m), 
                                        d[j].Lengths(m) );    }
                              */
                              PrintVisualAlignment( True, out, 
                                   hb.EdgeObject( p[j] ), 
                                   genome[g], d[j] );    }    }    }
     
                    }
          reports2[i] = out.str( );    }
     for ( int i = 0; i < reports2.isize( ); i++ )
          cout << reports2[i];

     // Harvest new aligns. (NOT YET!)

     for ( int e = 0; e < hb.E( ); e++ )
     {    UniqueSort( aligns2[e] );
          // if ( aligns2[e].nonempty( ) ) PRINT2( e, aligns2[e].size( ) ); // XXXXX
          if ( aligns2[e].size( ) > 1 ) aligns2[e].resize(0);    }

     // Find pairs of inconsistent alignments (again).

     bad_pairs.clear( );
     for ( int v = 0; v < hb.N( ); v++ )
     for ( int i = 0; i < (int) hb.To(v).size( ); i++ )
     for ( int j = 0; j < (int) hb.From(v).size( ); j++ )
     {    int e1 = hb.ITo(v,i), e2 = hb.IFrom(v,j);
          if ( !BinMember( es, e1 ) ||  !BinMember( es, e2 ) ) continue;

          // NOT RIGHT!
          if ( !aligns[e1].solo( ) || !aligns[e2].solo( ) ) continue;

          int g1, g2;
          align a1, a2;

          if ( aligns2[e1].nonempty( ) )
          {    g1 = aligns2[e1][0].first;
               a1 = aligns2[e1][0].second;    }
          else 
          {    g1 = aligns[e1][0].g;
               a1 = aligns[e1][0].a;    }

          if ( aligns2[e2].nonempty( ) )
          {    g2 = aligns2[e2][0].first;
               a2 = aligns2[e2][0].second;    }
          else 
          {    g2 = aligns[e2][0].g;
               a2 = aligns[e2][0].a;    }

          if ( g1 != g2 || !Consistent( a1, a2, hb.K( ) ) )
          {    // cout << "bad: " << e1 << "," << e2 << endl;
               bad_pairs.push( e1, e2 );    }    }
     cout << "found " << bad_pairs.size( ) << " bad pairs" << endl;

     /*
     for ( int j = 0; j < bad_pairs.isize( ); j++ )
     {    int e1 = bad_pairs[j].first, e2 = bad_pairs[j].second;
          // cout << "\nTesting consistency...\n";
          // PRINT2( e1, e2 );
          if ( aligns2[e1].nonempty( ) && aligns2[e2].nonempty( )
               && aligns2[e1][0].first == aligns2[e2][0].first
               && Consistent( 
                    aligns2[e1][0].second, aligns2[e2][0].second, hb.K( ) ) )
          {    cout << "Now " << e1 << " and " << e2 << " are consistent." 
                    << endl;    }    }
     */

     // Write output.

     if ( LINE < 0 )
     {    cout << "\n" << Date( ) << ": writing output" << endl;
          vec< vec< pair<int,align> > > alignsx( hb.E( ) );
          for ( int e = 0; e < hb.E( ); e++ )
          {    alignsx[e].resize( aligns[e].size( ) );
               for ( int i = 0; i < alignsx[e].isize( ); i++ )
               {    alignsx[e][i].first = aligns[e][i].g;
                    alignsx[e][i].second = aligns[e][i].a;    }    }
          BinaryWriter::writeFile( DIR + "/a.alignsx" + suffix, alignsx );
          BinaryWriter::writeFile( DIR + "/a.bad_pairs" + suffix, bad_pairs );    }

     // Done.

     cout << Date( ) << ": done, " << TimeSince(clock) << " used" << endl;
     Scram(0);    }
