///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Flag regions of an assembly that appear to be weak.  Consider the reads that
// are placed on the assembly.  These reads are used in two ways to define weak
// regions:
// (i)  An aligned 20-mer on a read is 'anchored' if its quality score sum goes 
//      up enough when its position is perturbed.  Imagine trimming reads to their
//      extremal anchor points.  Then they should cover the edge, with overlaps.
// (ii) Positions whose naive base call are not the base assigned to the edge.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "feudal/PQVec.h"
#include "math/HoInterval.h"
#include "paths/HyperBasevector.h"
#include "paths/long/FriendAligns.h"
#include "paths/long/ReadPath.h"
#include "paths/long/ReadStack.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_Int_OrDefault(VERBOSITY, 0);
     CommandArgument_Int_OrDefault(E, -1);
     EndCommandArguments;

     // Heuristics.

     const int L = 20;
     const int radius = 5;
     const double min_ratio = 5.0;
     const int min_add = 10;
     const int min_qsum = 2;
     const double min_qsum_ratio = 1.5;
     const int gap_close = 20;

     // Define directories.

     String work_dir = "/wga/scr4/jaffe/GapToy/plasmo";
     String fin_dir = work_dir + "/a.final";

     // Load data.

     cout << Date( ) << ": loading data" << endl;
     HyperBasevector hb;
     BinaryReader::readFile( fin_dir + "/a.hbv", &hb );
     vec<int> inv;
     BinaryReader::readFile( fin_dir + "/a.inv", &inv );
     ReadPathVec paths;
     paths.ReadAll( fin_dir + "/a.paths" );
     VecULongVec paths_index;
     invert( paths, paths_index, hb.EdgeObjectCount( ) );
     vecbasevector bases( work_dir  + "/data/frag_reads_orig.fastb" );
     VecPQVec quals( work_dir + "/data/frag_reads_orig.qualp" );

     // Set up output.

     vecbitvector flagged( hb.EdgeObjectCount( ) );
     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
          flagged[e].resize( hb.EdgeLengthBases(e), False );

     // Declare edge of interest.

     // negative controls
     // vec<int> es = {76646,79734,83885,77902,80514,77264,81901,66972,82768,76785};
     // vec<int> es = {80426,82294,76134,83941,78298,78123,80180,79507,79571,80023};

     // positive controls
     // vec<int> es = {78358,84156,80555,74569,71690,51278,71682,80812,78988,71692};

     double clock = WallClockTime( );
     vec< pair<int,String> > reports;
     cout << Date( ) << ": starting main loop" << endl;
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
     {    if ( E >= 0 && e != E ) continue;

          ostringstream out;
          int ne = hb.EdgeLengthBases(e);
          out << "\n";
          PRINT3_TO( out, e, ne, inv[e] );

          // Find support.
     
          vec< triple<int64_t,int,Bool> > locs;
          for ( int pass = 1; pass <= 2; pass++ )
          {    const int f = ( pass == 1 ? e : inv[e] );
               for ( int i = 0; i < (int) paths_index[f].size( ); i++ )
               {    int64_t id = paths_index[f][i];
                    const ReadPath& p = paths[id];
                    for ( int j = 0; j < (int) p.size( ); j++ )
                    {    if ( p[j] == f )
                         {    int pos = p.getOffset( );
                              for ( int l = 0; l < j; l++ )
                                   pos -= hb.Kmers( p[l] );
                              locs.push( id, pos, pass == 1 );    }    }    }    }

          // Check for no coverage at all.

          if ( locs.empty( ) )
          {    if ( ne > 0 ) out << "flagged: 0-" << ne << endl;
               #pragma omp critical
               {    reports.push( e, out.str( ) );    }
                    continue;    }

          // Build stack.
     
          readstack stack( locs.size( ), ne );
          for ( int i = 0; i < locs.isize( ); i++ )
          {    int64_t id = locs[i].first;
               qualvector q;
               quals[id].unpack(&q);
               int pos = locs[i].second;
               Bool fw = locs[i].third;
               stack.SetLen( i, bases[id].size( ) );
               if (fw)
               {    stack.SetOffset( i, pos );
                    for ( int j = 0; j < bases[id].isize( ); j++ )
                    {    int p = pos + j;
                         if ( p >= 0 && p < ne )
                         {    stack.SetBase( i, p, bases[id][j] );
                              stack.SetQual( i, p, q[j] );    }    }    }
               else 
               {    stack.SetOffset( i, ne - pos - bases[id].isize() );
                    for ( int j = 0; j < bases[id].isize( ); j++ )
                    {    int p = pos + j;
                         if ( p >= 0 && p < ne )
                         {    stack.SetBase( i, ne - p - 1, 3 - bases[id][j] );
                              stack.SetQual( 
                                   i, ne - p - 1, q[j] );    }    }    }    }

          // Test anchoring.  

          vec<ho_interval> cov;
          for ( int j = 0; j < stack.Rows( ); j++ )
          {    ho_interval x( stack.Offset(j), stack.Offset(j) + stack.Len(j) );
               if ( Subset( x, cov ) ) continue;
               vec<int> anchors;
               // effectively executing two passes:
               // for ( int i = 0; i <= ne - L; i++ )
               // for ( int i = ne - L; i >= 0; i-- )
               for ( int xpass = 1; xpass <= 2; xpass++ )
               {    int i = ( xpass == 1 ? 0 : ne - L );
                    while(1)
                    {    if ( i > ne - L || i < 0 ) break;
                         Bool def = True;
                         for ( int l = 0; l < L; l++ )
                         {    if ( !stack.Def( j, i + l ) )
                              {    def = False;
                                   break;    }    }
                         if (def)
                         {    int baseline = 0, alt = 1000000000;
                              for ( int r = -radius; r <= +radius; r++ )
                              {    if ( i + r < 0 || i + r + L > ne ) continue;
                                   int qsum = 0;
                                   for ( int l = 0; l < L; l++ )
                                   {    if ( stack.Base( j, i+l ) 
                                             != hb.EdgeObject(e)[i+l+r] )
                                        {    qsum += stack.Qual( j, i + l );    
                                                  }    }
                                   if ( r == 0 ) baseline = qsum;
                                   else alt = Min( alt, qsum );    }
                              if ( alt >= min_ratio * baseline 
                                   && alt - baseline >= min_add ) 
                              {    anchors.push_back(i);
                                   break;    }    }
                         i += ( xpass == 1 ? +1 : -1 );    }    }
               if ( anchors.nonempty( ) )
               {    int low = Min(anchors), high = Max(anchors) + L;
                    cov.push( low, high );    }
               Sort(cov);
               vec<ho_interval> cov2;
               for ( int i = 0; i < cov.isize( ); i++ )
               {    int low = cov[i].Start( ), high = cov[i].Stop( ), j;
                    for ( j = i + 1; j < cov.isize( ); j++ )
                    {    if ( high - cov[j].Start( ) >= L ) 
                              high = Max( high, cov[j].Stop( ) );
                         else break;    }
                    cov2.push( low, high );
                    i = j - 1;    }
               cov = cov2;    }

          // Compute uncovered.

          vec<ho_interval> uncov;
          Uncovered( ne, cov, uncov );
          for ( int i = 0; i < uncov.isize( ); i++ )
          for ( int u = uncov[i].Start( ); u < uncov[i].Stop( ); u++ )
               flagged[e].Set( u, True );
          if ( VERBOSITY >= 1 )
          {    for ( int i = 0; i < uncov.isize( ); i++ )
               {    int start = uncov[i].Start( ), stop = uncov[i].Stop( );
                    out << "unanchored: " << start << "-" << stop << endl;    }    }
     
          // Report suspicious quality score sums.
     
          for ( int i = 0; i < ne; i++ )
          {    if ( flagged[e][i] ) continue;

               // Check for zero-coverage range.

               int cov = 0;
               for ( int j = 0; j < stack.Rows( ); j++ )
               {    if ( !stack.Def( j, i ) ) continue;
                    cov += stack.Qual(j,i);    }
               if ( cov == 0 )
               {    int k;
                    for ( k = i + 1; k < ne; k++ )
                    {    int cov = 0;
                         for ( int j = 0; j < stack.Rows( ); j++ )
                         {    if ( !stack.Def( j, k ) ) continue;
                              cov += stack.Qual(j,k);    }
                         if ( cov != 0 ) break;    }
                    if ( VERBOSITY >= 1 ) 
                         out << "zero coverage: " << i << "-" << k << endl;
                    uncov.push( i, k );
                    i = k - 1;
                    continue;    }

               // Check qsums.

               vec<int> qsum( 4, 0 ), ids( 4, vec<int>::IDENTITY );
               for ( int j = 0; j < stack.Rows( ); j++ )
               {    if ( !stack.Def( j, i ) ) continue;
                    qsum[ stack.Base(j,i) ] += stack.Qual(j,i);    }
               ReverseSortSync( qsum, ids );
               Bool ok = ( ids[0] == hb.EdgeObject(e)[i] );
               if ( !ok ) uncov.push( i, i+1 );
               if ( VERBOSITY <= 1 && ok && qsum[0] >= min_qsum 
                    && qsum[0] >= min_qsum_ratio * qsum[1] ) 
               {    continue;    }
               if ( VERBOSITY >= 1 )
               {    out << i << " -->";
                    for ( int j = 0; j < 4; j++ )
                    {    out << " " << as_base( ids[j] ) << "[" << qsum[j] << "]";
                         if ( ids[j] == hb.EdgeObject(e)[i] ) out << "*";    }
                    out << "\n";    }    }

          // Report results.

          vec<ho_interval> uncov2;
          ExtractGivenCoverage( ne, 1, uncov, uncov2 );
          for ( int i = 0; i < uncov2.isize( ); i++ )
          {    int start = uncov2[i].Start( ), stop = uncov2[i].Stop( );
               int j;
               for ( j = i + i; j < uncov2.isize( ); j++ )
               {    if ( uncov2[j].Start( ) - stop > gap_close ) break;
                    stop = Max( stop, uncov2[j].Stop( ) );    }
               if ( stop > start )
                    out << "flagged: " << start << "-" << stop << endl;
               for ( int u = start; u < stop; u++ )
                    flagged[e].Set( u, True );
               i = j - 1;    }

          // Save.

          #pragma omp critical
          {    reports.push( e, out.str( ) );    }    }

     // Generate output.
     
     Sort(reports);
     for ( int i = 0; i < reports.isize( ); i++ )
          cout << reports[i].second;
     flagged.WriteAll( fin_dir + "/a.flagged" );
     cout << "\n" << TimeSince(clock) << " spent in main loop" << endl;
     Scram(0);    }
