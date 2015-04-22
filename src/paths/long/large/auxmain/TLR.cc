///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// NEXT STEP: figure out how to make this faster.

// TLR (ThreadLongReads).  Thread long reads through a DISCOVAR de novo assembly 
// graph.  Experimental.
// 1. Make kmer lookup tables for the assembly using K = 32 and K2 = 16.
// 2. Scan the reads to find K-mer matches in which the assembly K-mer is unique.
// 3. If a read hits a given line, add all its K2-mer matches to that line,
//    so long as the K2-mer is unique in the line, and its offset relative to the
//    line is within max_offset_diff = 100 of some K-mer hit, and it does not 
//    overlap a K-mer hit.
// 4. If a read hits a given line along a single edge, and the spread of the hits
//    is less than min_spread = 100, delete the hits.
// 5. Sort the hits for a given read by their position on the read.
// 6. Whenever two consecutive hits are to edges ---->v and w---> on the same line,
//    find all paths between v and w that are consistent with the distance predicted
//    by the read, up to parameters lfudge and lfudge2.
// ...
//
// nreads = 60000 --> 37.1 Mb, 14.7 minutes = 0.40 min/Mb

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "FastIfstream.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "PrintAlignment.h"
#include "math/HoInterval.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "paths/HyperBasevector.h"
#include "paths/RemodelGapTools.h"
#include "paths/long/large/Lines.h"
#include "random/Shuffle.h"

void Dump( vec<int>& upath, ostringstream& out )
{    out << "# unique path " << printSeq(upath) << endl;
     upath.clear( );    }

void CloseGap( const HyperBasevector& hb, const digraphE<int>& ihb,
     const vec<int>& to_right, const basevector& R, const int rstart, 
     const int rstop, const int eid0, const int eid, const int e0stop, 
     const int estart, vec<int>& upath, const double lfudge, const int lfudge2,
     const int max_paths, const int v, const int w, ostringstream& out )
{
     // Find paths of the right length.

     double clock = WallClockTime( );
     vec<vec<int>> paths;
     int expect = rstop - rstart - hb.Kmers(eid0) + e0stop - estart;
     int L1 = int( floor( ( 1.0 - lfudge ) * expect ) ) - lfudge2;
     int L2 = int( ceil( ( 1.0 + lfudge ) * expect ) ) + lfudge2;
     Bool ok = ihb.AllPathsLengthRange( v, w, L1, L2, to_right, paths, max_paths );
     out << TimeSince(clock) << " used finding paths" << endl;
     if ( !ok )
     {    out << "paths exploded" << endl;
          Dump(upath,out);    }
     else if ( paths.empty( ) )
     {    out << "no paths found, expect = " << expect
               << ", rstop - rstart = " << rstop - rstart << endl;
          Dump(upath,out);    }
     else
     {    basevector r( R, rstart, rstop - rstart );

          // Precheck for paths that would lead to out-of-range condition in
          // scoring.  This is presumably indicative of a bug somewhere else that
          // should be fixed.

          vec<Bool> bad( paths.size( ), False );
          for ( int v = 0; v < paths.isize( ); v++ )
          {    const vec<int>& p = paths[v];
               vec<int> q = {eid0};
               q.append(p);
               q.push_back(eid);
               basevector b = hb.Cat(q);
               if ( e0stop < 0 || e0stop >= b.isize( )
                    || b.isize( ) - e0stop - ( hb.Bases(eid) - estart ) < 0
                    || b.isize( ) - ( hb.Bases(eid) - estart ) > b.isize( ) )
               {    bad[v] = True;    }    }
          EraseIf( paths, bad );

          // Score paths.

          int np = paths.size( );
          vec<int> score(np);
          vec<align> aligns(np);
          vec<basevector> c(np);
          double clock2 = WallClockTime( );
          for ( int v = 0; v < np; v++ )
          {    const vec<int>& p = paths[v];
               vec<int> q = {eid0};
               q.append(p);
               q.push_back(eid);
               basevector b = hb.Cat(q);
               c[v] = basevector( b, e0stop, 
                    b.isize( ) - e0stop - ( hb.Bases(eid) - estart ) );
               alignment al;

               score[v] = SmithWatAffine( r, c[v], al );

               /*
               int nerrors;
               int offset = ***
               int bandwidth = ***
               score[v] = 
                    SmithWatAffineBanded( r, c[v], offset, bandwidth, al, nerrors );
               */

               aligns[v] = al;    }
          out << TimeSince(clock2) << " used scoring " << np
               << " paths on a read segment of length " << r.size( ) << endl;

          // Keep only best paths.

          vec<int> ids( np, vec<int>::IDENTITY );
          SortSync( score, aligns, ids, c );
          for ( int v = 1; v < np; v++ )
          {    if ( score[v] > score[v-1] )
               {    np = v;
                    break;    }    }
          if ( np == 1 ) upath.append( paths[ ids[0] ] );
          else Dump(upath,out);

          // Print the winners.

          for ( int v = 0; v < np; v++ )
          {    const vec<int>& p = paths[ ids[v] ];
               out << "[" << v+1 << "] = " << printSeq(p) << endl;
               int plen = 0;
               for ( int j = 0; j < p.isize( ); j++ )
                    plen += hb.Kmers( p[j] );
               double delta = double(plen-expect) / expect;
               int errs = ActualErrors( r, c[v], aligns[v], 1, 1 );
               double err_rate = double(errs) / r.size( );
               out << "path length = " << plen << ", expect = " << expect 
                    << ", off by " << PERCENT_RATIO(3, plen-expect, expect)
                    << ", error rate = " << PERCENT_RATIO( 3, errs, r.isize( ) )
                    << "\n";
               if ( err_rate > 0.25 )
               {    PrintVisualAlignment( True, out, r, 
                         c[v], aligns[v] );    }    }    }    }

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_Bool_OrDefault(CACHE, False);
     CommandArgument_Int_OrDefault(CACHE_ID, 1 );
     EndCommandArguments;

     double allclock = WallClockTime( );

     // Load assembly.

     cout << Date( ) << ": loading assembly" << endl;
     String work_dir = "/wga/scr4/jaffe/GapToy/51400.newchem";
     HyperBasevector hb;
     BinaryReader::readFile( work_dir + "/a.final/a.hbv", &hb );
     vec<vec<vec<vec<int>>>> lines;
     BinaryReader::readFile( work_dir + "/a.final/a.lines", &lines );

     // Set up ancillary data structures for the assembly.

     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     vec<int> tol, tol2( hb.E( ), -1 );
     GetTol( hb, lines, tol );
     for ( int l = 0; l < lines.isize( ); l++ )
     for ( int i = 0; i < lines[l].isize( ); i++ )
     for ( int j = 0; j < lines[l][i].isize( ); j++ )
     for ( int k = 0; k < lines[l][i][j].isize( ); k++ )
          tol2[ lines[l][i][j][k] ] = i;
     vec<int> lens( hb.E( ) );
     for ( int e = 0; e < hb.E( ); e++ )
          lens[e] = hb.Kmers(e);
     digraphE<int> ihb( hb, lens );

     // Define edge starts on lines.

     cout << Date( ) << ": defining edge starts on lines" << endl;
     vec<int> lstart( hb.E( ), 0 );
     for ( int l = 0; l < lines.isize( ); l++ )
     {    int pos = 0;
          for ( int i = 0; i < lines[l].isize( ); i++ )
          {    vec<int> plens;
               vec< pair<int,int> > places;
               for ( int j = 0; j < lines[l][i].isize( ); j++ )
               {    int plen = 0;
                    for ( int k = 0; k < lines[l][i][j].isize( ); k++ )
                    {    int e = lines[l][i][j][k];
                         places.push( e, pos + plen );
                         plen += hb.Kmers(e);    }
                    plens.push_back(plen);    }
               Sort(plens), Sort(places);
               for ( int j = 0; j < places.isize( ); j++ )
               {    int k;
                    vec<int> p = { places[j].second };
                    for ( k = j + 1; k < places.isize( ); k++ )
                    {    if ( places[k].first != places[j].first ) break;
                         p.push_back( places[k].second );    }
                    Sort(p);
                    if ( p.nonempty( ) ) lstart[ places[j].first ] = Median(p);
                    j = k - 1;    }
               if ( plens.nonempty( ) ) pos += Median(plens);    }    }

     // Get some PacBio reads and toss the shorter ones.

     cout << Date( ) << ": loading pb reads" << endl;
     vecbasevector pb;
     vec<basevector> pbx;
     const int nreads = 20000;
     const int min_read = 4000;
     fast_pipe_ifstream in( "samtools view "
          "/wga/scr4/human_data/NA12878/NA12878.pacbio.bwa-mem.20131224.bam "
          "| head -" + ToString(nreads) + " | Col 10" );
     String line;
     int64_t nbases = 0;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          if ( line.isize( ) >= min_read ) 
          {    pbx.push_back( basevector(line) );
               nbases += line.size( );    }    }
     Shuffle( pbx.begin( ), pbx.end( ) );
     for ( int i = 0; i < pbx.isize( ); i++ )
          pb.push_back( pbx[i] );
     cout << Date( ) << ": found " << pb.size( ) << " reads = " 
          << ToStringAddCommas(nbases) << " bases" << endl;

     // Heuristics.

     const int K = 32;
     const int K2 = 16;
     const int min_spread = 100;
     const int max_paths = 1000;
     double lfudge = 0.2;
     const int lfudge2 = 100;
     const int max_offset_diff = 100;

     vec< quad<int,int,int,ho_interval> > hitsp;
     if (CACHE) BinaryReader::readFile( "hitsp." + ToString(CACHE_ID), &hitsp );
     else
     {    
          // Build a genome lookup table using K2.

          cout << Date( ) << ": building using K2" << endl;
          vecbasevector edges2;
          for ( int e = 0; e < hb.E( ); e++ )
          {    basevector b = hb.EdgeObject(e);
               if ( b.size( ) > 0 ) b.resize( b.isize( ) - hb.K( ) + K2 );
               edges2.push_back(b);    }
          vec< triple<kmer<K2>,int,int> > kmers_plus2;
          cout << Date( ) << ": looking up kmers" << endl;
          MakeKmerLookup0( edges2, kmers_plus2 );

          // Process with K.

          cout << Date( ) << ": proceeding with K" << endl;
          vecbasevector edges;
          for ( int e = 0; e < hb.E( ); e++ )
          {    basevector b = hb.EdgeObject(e);
               if ( b.size( ) > 0 ) b.resize( b.isize( ) - hb.K( ) + K );
               edges.push_back(b);    }

          // Build all.
     
          vecbasevector all(edges);
          // all.Append(pb);

          // Build kmers.
     
          cout << Date( ) << ": building kmers" << endl;
          vec< triple<kmer<K>,int,int> > kmers_plus;
          MakeKmerLookup0( all, kmers_plus );

          // Find matches using K.

          cout << Date( ) << ": finding matches" << endl;
          vec< triple<int,int,int> > hits;
          const int64_t batches = 100;
          #pragma omp parallel for schedule(dynamic, 1)
          for ( int64_t bi = 0; bi < batches; bi++ )
          {    kmer<K> x;
               vec< triple<int,int,int> > hitsi;
               vec< quad<int,int,int,ho_interval> > hitspi;
               for ( int64_t rid = ( bi * (int64_t) pb.size( ) ) / batches; 
                    rid < ( (bi+1) * (int64_t) pb.size( ) ) / batches; rid++ )
               {    for ( int rpos = 0; rpos <= pb[rid].isize( ) - K; rpos++ )
                    {    x.SetToSubOf( pb[rid], rpos );
                         int64_t low = LowerBound1( kmers_plus, x );
                         int64_t high = UpperBound1( kmers_plus, x );
                         if ( high - low == 1 )
                         {    int offset = kmers_plus[low].third - rpos;
                              int eid = kmers_plus[low].second;
                              hitsi.push( rid, eid, offset );    
                              hitspi.push( rid, eid, offset, 
                                   ho_interval( rpos, rpos + K ) );    }    }    }
               #pragma omp critical
               {    hits.append(hitsi);
                    hitsp.append(hitspi);    }    }

          // Sort matches.

          cout << Date( ) << ": sorting matches" << endl;
          ParallelUniqueSort(hits);
          ParallelSort(hitsp);

          // Find line ids.

          vec<vec<int>> lids( pb.size( ) ); // line ids
          for ( int64_t i = 0; i < hits.isize( ); i++ )
          {    int rid = hits[i].first, eid = hits[i].second, offset = hits[i].third;
               lids[rid].push_back( tol[eid] );    }
          for ( int64_t i = 0; i < (int64_t) pb.size( ); i++ )
               UniqueSort( lids[i] );

          // Index hits by reads.

          vec<int64_t> starts0( pb.size( ) + 1, -1 );
          for ( int64_t i = hitsp.jsize( ) - 1; i >= 0; i-- )
               starts0[ hitsp[i].first ] = i;
          starts0[ pb.size( ) ] = hitsp.size( );
          for ( int i = starts0.isize( ) - 1; i >= 0; i-- )
               if ( starts0[i] < 0 ) starts0[i] = starts0[i+1];

          // Look for K2 hits.

          cout << Date( ) << ": looking for K2 hits" << endl;
          #pragma omp parallel for schedule(dynamic, 1)
          for ( int64_t bi = 0; bi < batches; bi++ )
          {    vec< quad<int,int,int,ho_interval> > hitsrp;
               for ( int64_t rid = ( bi * (int64_t) pb.size( ) ) / batches; 
                    rid < ( (bi+1) * (int64_t) pb.size( ) ) / batches; rid++ )
               {
                    kmer<K2> x;
                    vec< triple<int,int,int> > xhits;
                    for ( int rpos = 0; rpos <= pb[rid].isize( ) - K2; rpos++ )
                    {    x.SetToSubOf( pb[rid], rpos );
                         int64_t low = LowerBound1( kmers_plus2, x );
                         int64_t high = UpperBound1( kmers_plus2, x );
          
                         // if ( high - low != 1 ) continue;
     
                         xhits.clear( );
                         for ( int64_t l = low; l < high; l++ )
                         {    int eid = kmers_plus2[l].second;
                              if ( !BinMember( lids[rid], tol[eid] ) ) continue;
                              int offset = kmers_plus2[l].third - rpos;
                              xhits.push( tol[eid], eid, offset );    }
                         Sort(xhits);
                         for ( int m1 = 0; m1 < xhits.isize( ); m1++ )
                         {    int m2;
                              for ( m2 = m1 + 1; m2 < xhits.isize( ); m2++ )
                                   if ( xhits[m2].first != xhits[m1].first ) break;
                              if ( m2 - m1 == 1 )
                              {    int eid = xhits[m1].second; 
                                   int offset = xhits[m1].third;
                                   int loffset = offset + lstart[eid];
          
                                   Bool conflict = False;
                                   for ( int64_t i = starts0[rid]; 
                                        i < starts0[rid+1]; i++ )
                                   {    int eid2 = hitsp[i].second;
                                        if ( tol[eid2] != tol[eid] ) continue;
                                        int loffset2 = hitsp[i].third + lstart[eid2];
                                        const ho_interval& h2 = hitsp[i].fourth;
                                        if ( IntervalOverlap( rpos, rpos + K2, 
                                             h2.Start( ), h2.Stop( ) ) > 0 )
                                        {    conflict = True;
                                             break;    }    }
     
                                   if ( !conflict )
                                   {    Bool match = False;
                                        for ( int64_t i = starts0[rid]; 
                                             i < starts0[rid+1]; i++ )
                                        {    int eid2 = hitsp[i].second;
                                             if ( tol[eid2] != tol[eid] ) continue;
                                             int loffset2 
                                                  = hitsp[i].third + lstart[eid2];
                                             const ho_interval& h2 = hitsp[i].fourth;
                                             if ( Abs( loffset - loffset2 ) 
                                                  <= max_offset_diff )
                                             {    match = True;
                                                  break;    }    }
                                        if (match)
                                        {    hitsrp.push( rid, eid, offset, 
                                                  ho_interval( 
                                                  rpos, rpos + K2 ) );    }    }    }
                              m1 = m2 - 1;    }    }    }
               #pragma omp critical
               {    hitsp.append(hitsrp);    }    }

          // Sort hitsp.

          cout << Date( ) << ": sorting hitsp" << endl;
          ParallelSort(hitsp);

          // Remove hits showing too little spread.  Preliminary, as it only looks
          // at spread along single edges.

          vec<Bool> to_delete( hitsp.size( ), False );
          for ( int64_t i = 0; i < hitsp.jsize( ); i++ )
          {    int rid = hitsp[i].first, eid = hitsp[i].second; 
               int offset = hitsp[i].third;
               int64_t j;
               for ( j = i + 1; j < hitsp.jsize( ); j++ )
                    if ( hitsp[j].first != rid ) break;
          
               // Looking at hits for read rid.
     
               vec< pair<int,int> > ls;
               vec<int> ls1;
               for ( int64_t k = i; k < j; k++ )
               {    int e = hitsp[k].second;
                    ls.push( tol[e], e );    }
               UniqueSort(ls);
               for ( int r = 0; r < ls.isize( ); r++ )
               {    int s;
                    for ( s = r + 1; s < ls.isize( ); s++ )
                         if ( ls[s].first != ls[r].first ) break;
                    if ( s - r == 1 ) ls1.push_back( ls[r].first );
                    r = s - 1;    }
               for ( int64_t k = i; k < j; k++ )
               {    int64_t m;
                    for ( m = k + 1; m < k; m++ )
                         if ( hitsp[m].second != hitsp[k].second ) break;
                    int l = tol[ hitsp[k].second ];
                    if ( BinMember( ls1, l ) )
                    {    vec<ho_interval> h;
                         for ( int64_t n = k; n < m; n++ )
                              h.push_back( hitsp[n].fourth );
                         int x = TotalCovered(h);
                         if ( x < min_spread )
                         {    for ( int64_t n = k; n < m; n++ )
                                   to_delete[n] = True;    }    }
                    k = m - 1;    }
               i = j - 1;   }
          EraseIf( hitsp, to_delete );
          BinaryWriter::writeFile( "hitsp." + ToString(CACHE_ID), hitsp );    }

     // Index hits by reads.

     vec<int64_t> starts( pb.size( ) + 1, -1 );
     for ( int64_t i = hitsp.jsize( ) - 1; i >= 0; i-- )
          starts[ hitsp[i].first ] = i;
     starts[ pb.size( ) ] = hitsp.size( );
     for ( int i = starts.isize( ) - 1; i >= 0; i-- )
          if ( starts[i] < 0 ) starts[i] = starts[i+1];
               
     // Print matches.

     cout << Date( ) << ": printing matches" << endl;
     const int batches2 = 200;
     double pclock = WallClockTime( );
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int64_t bi = 0; bi < batches2; bi++ )
     {    double clock = WallClockTime( );
          ostringstream out;
          for ( int64_t rid = ( bi * (int64_t) pb.size( ) ) / batches2; 
               rid < ( (bi+1) * (int64_t) pb.size( ) ) / batches2; rid++ )
          {    
               // Looking at one read.

               out << "\nrid = " << rid << endl;
               vec< quad<ho_interval,int,int,int> > locs; // (H,eid,offset,n)

               for ( int64_t r = starts[rid]; r < starts[rid+1]; r++ )
               {    int eid = hitsp[r].second, offset = hitsp[r].third;
                    int64_t j;
                    for ( j = r + 1; j < hitsp.jsize( ); j++ )
                    {    if ( hitsp[j].second != hitsp[r].second ) break;
                         if ( hitsp[j].third != hitsp[r].third ) break;    }
                    vec<ho_interval> h;
                    for ( int64_t l = r; l < j; l++ )
                         h.push_back( hitsp[l].fourth );
                    int low = h[0].Start( ), high = h[0].Stop( );
                    for ( int l = 1; l < h.isize( ); l++ )
                    {    low = Min( low, h[l].Start( ) ); 
                         high = Max( high, h[l].Stop( ) );    }
                    ho_interval H( low, high );
                    int n = TotalCovered(h) - K2 + 1;
                    locs.push( H, eid, offset, n );
                    r = j - 1;    }
               Sort(locs);
               for ( int u = 0; u < locs.isize( ); u++ )
               {    
                    // Look for hits to this line.

                    int uu;
                    for ( uu = u + 1; uu < locs.isize( ); uu++ )
                    {    if ( tol[ locs[uu].second ] != tol[ locs[u].second ] ) 
                              break;    }
                    out << "\nexploring hits to line " << tol[ locs[u].second ] 
                         << "\n" << endl;
                    vec<int> upath;

                    for ( int mu = u; mu < uu; mu++ )
                    {    const ho_interval& H = locs[mu].first;
                         int eid = locs[mu].second, offset = locs[mu].third; 
                         int n = locs[mu].fourth, line = tol[eid];
                         if ( mu > u )
                         {    int eid0 = locs[mu-1].second; 
                              int offset0 = locs[mu-1].third;
                              const ho_interval& H0 = locs[mu-1].first;
                              int v = to_right[eid0], w = to_left[eid];
                              int rstart = H0.Stop( ), rstop = H.Start( );
                              int e0stop = rstart + offset0, estart = rstop + offset;
                              if ( eid0 != eid && v != w && rstart >= rstop )
                              {    out << "-- huh, see read overlap of " 
                                        << rstart - rstop << endl;
                                   PRINT2_TO( out, rstart, rstop );
                                   int left = Max( 0, rstop - 32 );
                                   int right = Min( pb[rid].isize( ), rstart + 32 );
                                   basevector b( pb[rid], left, right - left );
                                   b.Print( out, "read." + ToString(left) + "-"
                                        + ToString(right) );
                                   Dump(upath,out);     }
                              if ( eid0 != eid && v != w && rstart < rstop )
                              {    CloseGap( hb, ihb, to_right, pb[rid], rstart, 
                                        rstop, eid0, eid, e0stop, estart, upath, 
                                        lfudge, lfudge2, max_paths, 
                                        v, w, out );    }    }

                         // Print the main hit.

                         out << "e = " << eid << ", offset = " << offset 
                              << ", n" + ToString(K2) << " = " << n << ", rpos = " 
                              << H << "[" << H.Length( ) << "]" << endl;
                         if ( upath.empty( ) || eid != upath.back( ) )
                              upath.push_back(eid);    }
                    Dump(upath,out);
                    u = uu - 1;    }    }
          out << TimeSince(clock) << " used on batch " << bi << endl;
          #pragma omp critical
          {    cout << out.str( );    }    }
  
     // Done.

     cout << TimeSince(pclock) << " spent in parallel loop" << endl;
     cout << Date( ) << ": done, total time used = " << TimeSince(allclock)
          << endl;    }
