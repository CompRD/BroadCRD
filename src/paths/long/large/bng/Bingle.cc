///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Bingle.  Experimental code to align a DISCOVAR assembly to a BNG map assembly.
// This code has only been tested on human data and has some hardwired
// assumptions that will need to be relaxed to get it to work more generally.
//
// For reporting purposes we assume for now that a reference sequence is available.
//
// Code to be sensibly renamed at some point.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "FastIfstream.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "TokenizeString.h"
#include "graph/Digraph.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/long/large/Lines.h"
#include "paths/long/large/bng/BngAlign.h"
#include "paths/long/large/bng/Texps.h"

int main( int argc, char** argv )
{    RunTime( );

     BeginCommandArguments;

     /*
     CommandArgument_String_OrDefault_Doc(DISCO_DIR, 
          "/wga/scr4/jaffe/GapToy/50871", "path for DISCOVAR assembly directory");
     CommandArgument_String_OrDefault_Doc(BNG_ASSEMBLY,
          "/wga/scr4/vendor/bng/Aseembly/human_NA12878_5v1b_cmb/"
          "contigs/exp_refineFinal1", "path for BNG assembly");
     CommandArgument_String_OrDefault_Doc(BNG_CONTIG_HEAD, "exp_refineFinal1_contig",
          "head for contig files within BNG_ASSEMBLY" );
     */
     CommandArgument_String_OrDefault_Doc(SAMPLE, "", "NA12878 or tumor");
     CommandArgument_String_OrDefault_Doc(OUT_HEAD, "Bingle.51400",
          "head for output files");

     EndCommandArguments;

     // Configure hardcoded samples.

     String cmap, DISCO_DIR;
     if ( SAMPLE == "NA12878" )
     {    DISCO_DIR = "/wga/scr4/jaffe/GapToy/51400.newchem";
          cmap = "/wga/scr4/vendor/bng/NA12878/EXP_REFINEFINAL1.cmap";    }
     else if ( SAMPLE == "tumor" )
     {    DISCO_DIR = "/wga/scr4/jaffe/GapToy/50921.HCC1143+BL";
          cmap = "/wga/scr4/vendor/bng/HCC1143/EXP_REFINEFINAL1_HCC1143.cmap";    }
     else
     {    cout << "Unknown sample." << endl;
          exit(1);    }

     // Heuristics.  There are a bunch more embedded in the code.

     const int min_line = 10000;
     const int min_line2 = 25000;
     const int min_line3 = 20000;
     const int max_ignored_indel = 120;
     const int min_seed = 6500;
     const double max_diff = 0.025;
     const double max_d = 0.1;
     const int flank = 5;
     const int local_flank = 5;
     const int min_overlap = 10000;

     // Define enzyme.

     String cut = "GCTCTTC", rcut;
     StringReverseComplement( cut, rcut );

     // Load assembly.

     cout << Date( ) << ": loading assembly" << endl;
     String dir = DISCO_DIR + "/a.final";
     HyperBasevectorX hb;
     BinaryReader::readFile( dir + "/a.hbx", &hb );
     vec<int> inv;
     BinaryReader::readFile( dir + "/a.inv", &inv );
     vec<vec<vec<vec<int>>>> lines;
     BinaryReader::readFile( dir + "/a.lines", &lines );
     vec<int> lens;
     GetLineLengths( hb, lines, lens );
     const int K = hb.K( );

     // Load map of assembly to reference sequence.  Note that this would need
     // to be changed for a genome that does not have a reference sequence.  This
     // is also a very ugly way to load the map information.  

     vec<String> maplines;
     vec<int> maplines_id;
     fast_ifstream mapin( dir + "/a.lines.map" );
     String line;
     while(1)
     {    getline( mapin, line );
          if ( mapin.fail( ) ) break;
          if ( line.Contains( "line[" ) )
                maplines_id.push_back( line.Between( "line[", "]" ).Int( ) );
          else maplines_id.push_back(-1);
          maplines.push_back(line);    }

     // Go through the DISCOVAR lines, generating a data structure texps, that
     // we will use for mapping the lines to the BNG maps.
     //
     // This code also prints out a human-readable representation of the 
     // translation of each line to "map space".
     //
     // Example:
     //
     // 40630:
     // {{11426,2888,7665,382,540,2717,1527,375,2627,723}}
     // 
     // This is the simplest case: the line is unambiguously expanded as a list
     // of restrictions-site-free intervals.  Note that the first and last invervals
     // are incomplete and thus special.
     //
     // (Identical documentation in two other places.)

     vec<vec<vec<int>>> texps;
     Bool verbose = True;
     ostringstream out;
     Texps( hb, lines, cut, min_line, max_ignored_indel, texps, verbose, out );
     cout << out.str( );

     // Set up all_data, which is texps, organized differently.

     cout << Date( ) << ": setting up all_data" << endl;
     vec< quad<double,int,int,int> > all_data;
     int M;
     for ( M = 0; M < lines.isize( ); M++ )
          if ( lens[M] < min_line2 ) break;
     for ( int li = 0; li < M; li++ )
     {    for ( int i = 0; i < texps[li].isize( ); i++ )
          {    for ( int j = 0; j < texps[li][i].isize( ); j++ )
               {    if ( texps[li][i].size( ) >= 3 )
                    {    int x = texps[li][i][j];
                         all_data.push( x, li, i, j );    }    }    }    }
     ParallelSort(all_data);

     // Find the BNG assembly contigs.

     cout << Date( ) << ": loading BNG data" << endl;
     vec<int> mids;
     vec< vec<double> > POS;
     {    vec<double> pos;
          fast_ifstream xin(cmap);
          int last_mid = -1;
          int mid = -1;
          while(1)
          {    getline( xin, line );
               if ( xin.fail( ) ) break;
               if ( line.Contains( "#", 0 ) ) continue;
               vec<String> x;
               Tokenize( line, x );
               mid = x[0].Int( );
               if ( mid != last_mid && last_mid >= 0 )
               {    POS.push_back(pos);
                    pos.clear( );
                    mids.push_back(mid);     }
               pos.push_back( x[5].Double( ) );
               last_mid = mid;    }
          mids.push_back(mid);
          POS.push_back(pos);    }

     /*
     {    vec<String> all_files = AllFiles(BNG_ASSEMBLY);
          for ( int mi = 0; mi < all_files.isize( ); mi++ )
          {    String fn = all_files[mi];
               if ( !fn.Contains( ".cmap", -1 ) ) continue;
               if ( !fn.Contains( BNG_CONTIG_HEAD, 0 ) ) continue;
               if ( !fn.Between( BNG_CONTIG_HEAD, ".cmap" ).IsInt( ) ) continue;
               int mid = fn.Between( BNG_CONTIG_HEAD, ".cmap" ).Int( );
               mids.push_back(mid);    }    }
     */
     SortSync( mids, POS );
     int nmids = mids.size( );

     // Map DISCOVAR lines to BNG assembly.  This is where virtually all the run
     // time is spent.

     vec< vec< triple< pair<int,int>, int, double > > > all_aligns(nmids);
     vec<vec<int>> all_alen(nmids);

     for ( int mi = 0; mi < mids.isize( ); mi++ )
     {    int mid = mids[mi];

          // Save separations as a vector of integers b.

          // String amap = BNG_ASSEMBLY + "/" + BNG_CONTIG_HEAD + ToString(mid) 
          //      + ".cmap";
          cout << "\nPROCESSING MAP " << mid << endl << endl;
          // fast_ifstream in(amap);
          vec<double> pos = POS[mi];
          vec<int> b;
          for ( int j = 1; j < pos.isize( ); j++ )
               b.push_back( int( round( pos[j] - pos[j-1] ) ) );
          cout << "mean sep = " << Mean(b) << endl;
     
          // Find pairs (r,a) where r is a read loc, a is an assembly loc, 
          // r >= 20000, and |r-a|/a <= 2%.  (DESCRIPTION NO LONGER ACCURATE.)

          cout << Date( ) << ": start new approach" << endl;
          vec< triple< pair<int,int>, int, double > > aligns;

          // Set up data.

          vec< quad<int,int,int,int> > data;
          for ( int u = 0; u < b.isize( ); u++ )
          {    if ( b[u] < min_seed ) continue;
               double low = b[u] * ( 1.0 - max_diff );
               double high = b[u] * ( 1.0 + max_diff );
               int64_t L = LowerBound1( all_data, low );
               int64_t H = UpperBound1( all_data, high );
               for ( int64_t m = L; m < H; m++ )
               {    data.push( all_data[m].second, u, all_data[m].third,
                         all_data[m].fourth );    }    }
          
          // Start alignment loop.

          cout << Date( ) << ": start parallel loop" << endl;
          vec<Bool> ok( data.size( ), True );
          for ( int pass = 1; pass <= 2; pass++ )
          {
               #pragma omp parallel for schedule(dynamic, 1)
               for ( int di = 0; di < data.isize( ); di++ )
               {    if ( pass == 2 && !ok[di] ) continue;
     
                    int li = data[di].first, u = data[di].second;
                    int i = data[di].third, j = data[di].fourth;
     
                    int a_m = j;
                    vec<int> A = texps[li][i];
                    if ( pass == 1 )
                    {    if ( a_m < flank || A.isize( ) - a_m < flank ) continue;
                         vec<int> B;
                         for ( int l = a_m - flank; l < a_m + flank; l++ )
                              B.push_back( A[l] );
                         A = B;
                         a_m = flank;    }
                    // Why are we ADDING to left?  Doesn't make sense.
                    int left = 3*a_m/2 + local_flank;
                    int right = 3 * ( A.isize( ) - a_m ) / 2 + local_flank;
                    int r1 = Max( 0, u - left ), r2 = Min( b.isize( ), u + right );
                    vec<int> R;
                    for ( int v = r1; v < r2; v++ )
                         R.push_back( b[v] );
                    int r_m = u - r1;
                    ostringstream out;
                    int rstart, rstop;
                    int ndirect;
                    vec<int> p;
                    int rsum, asum;
                    double d = BngAlign( R, A, r_m, a_m, out, p, rstart, rstop, 
                         ndirect, rsum, asum, 1 );

                    if ( pass == 2 && ( rstart == 0 || rstop == R.isize( ) - 1 ) )
                         continue;

                    rstart += r1, rstop += r1;
                    int alen = 0;
                    for ( int z = rstart; z <= rstop; z++ )
                         alen += b[z];
                    if ( pass == 1 )
                    {    if ( d > max_d ) ok[di] = False;    }
                    else if ( d < max_d && alen >= min_line3 
                         && rstop - rstart + 1 >= 3 )
                    {    
                         #pragma omp critical
                         {    aligns.push( make_pair(rstart,rstop), 
                                   li, d );    }    }    }    }

          // Filter alignments.

          PRINT( aligns.size( ) );
          cout << Date( ) << ": sorting" << endl;
          UniqueSort(aligns);
          vec<int> alen( aligns.size( ), 0 );
          for ( int i = 0; i < aligns.isize( ); i++ )
          {    int rstart = aligns[i].first.first, rstop = aligns[i].first.second;
               for ( int j = rstart; j <= rstop; j++ )
                    alen[i] += b[j];    }
     
          // First round of deletions.
     
          vec<Bool> to_delete1( aligns.size( ), False );
          for ( int i1 = 0; i1 < aligns.isize( ); i1++ )
          for ( int i2 = 0; i2 < aligns.isize( ); i2++ )
          {    if ( i1 == i2 ) continue;
               Bool same_line = aligns[i1].second == aligns[i2].second;
               if ( same_line && aligns[i1].third < aligns[i2].third )
               {    to_delete1[i2] = True;    }    }
          EraseIf( aligns, to_delete1 );
          EraseIf( alen, to_delete1 );
     
          // Second round of deletions.
     
          vec<Bool> to_delete( aligns.size( ), False );
          for ( int i1 = 0; i1 < aligns.isize( ); i1++ )
          for ( int i2 = 0; i2 < aligns.isize( ); i2++ )
          {    if ( i1 == i2 ) continue;
               Bool subsumes = aligns[i1].first.first <= aligns[i2].first.first
                    && aligns[i1].first.second >= aligns[i2].first.second;
               Bool same_line = aligns[i1].second == aligns[i2].second;
     
               int overlap = 0;
               int low = Max( aligns[i1].first.first, aligns[i2].first.first );
               int high = Min( aligns[i1].first.second, aligns[i2].first.second );
               for ( int j = low; j <= high; j++ )
                    overlap += b[j];

               int l1 = aligns[i1].second, l2 = aligns[i2].second;
               double d1 = aligns[i1].third, d2 = aligns[i2].third;
               double score1 
                    = d1 / ( aligns[i1].first.second - aligns[i1].first.first );
               double score2 
                    = d2 / ( aligns[i2].first.second - aligns[i2].first.first );

               if ( l1 != l2 && low <= high && score1 < score2 )
               {    Bool PRINT_DELS = False;
                    if (PRINT_DELS)
                    {    cout << l1 << " beats " << l2 << ", d1 = " << d1 
                              << ", d2 = " << d2 << ", score1 = " << score1 
                              << ", score2 = " << score2 << endl;     }
                    to_delete[i2] = True;    }

               if ( !same_line && overlap >= min_overlap && alen[i1] > alen[i2]
                    && aligns[i1].third < aligns[i2].third )
               {    to_delete[i2] = True;    }
               if ( aligns[i1].third * 2 <= aligns[i2].third
                    && aligns[i2].third >= 0.05 && subsumes )
               {    to_delete[i2] = True;    }
               if ( aligns[i1].third < aligns[i2].third && same_line
                    && aligns[i1].first == aligns[i2].first )
               {    to_delete[i2] = True;    }
               if ( alen[i1] >= 5 * alen[i2] && subsumes
                    && aligns[i1].third * 0.5 <= aligns[i2].third )
               {    to_delete[i2] = True;    }    }
     
          EraseIf( aligns, to_delete );
          EraseIf( alen, to_delete );
          all_aligns[mi] = aligns;
          all_alen[mi] = alen;

          PRINT( aligns.size( ) );
          cout << Date( ) << ": done" << endl;    }

     // Get line inversion.  In rare cases a line does not have an "inverse line",
     // in which case the value of linv is -1.

     vec<int> linv( lines.size( ), -1 );
     {    vec<int> ids( lines.size( ), vec<int>::IDENTITY );
          vec< pair<int,int> > line_ends( lines.size( ) );
          #pragma omp parallel for
          for ( int i = 0; i < lines.isize( ); i++ )
          {    line_ends[i] = make_pair(
                    lines[i].front( )[0][0], lines[i].back( )[0][0] );    }
          ParallelSortSync( line_ends, ids );
          #pragma omp parallel for
          for ( int i = 0; i < line_ends.isize( ); i++ )
          {    int e1 = line_ends[i].first, e2 = line_ends[i].second;
               int p = BinPosition( line_ends, make_pair( inv[e2], inv[e1] ) );
               if ( p >= 0 ) linv[ ids[i] ] = ids[p];    }    }

     // Remove aligns where we have multiple locations, one much better.
     // First: flag inferior placements.

     const double flag_thresh = 1.8;
     cout << "\nANALYSIS OF DUPLICATES" << endl;
     vec< vec< triple<double,int,int> > > places( lines.size( ) );
     vec< vec< Bool > > flag(nmids);
     for ( int i = 0; i < all_aligns.isize( ); i++ )
          flag[i].resize( all_aligns[i].size( ), False );
     for ( int i = 0; i < all_aligns.isize( ); i++ )
     for ( int j = 0; j < all_aligns[i].isize( ); j++ )
     {    double d = all_aligns[i][j].third;
          places[ all_aligns[i][j].second ].push( d, i, j );    }
     for ( int i = 0; i < lines.isize( ); i++ )
     {    if ( linv[i] >= 0 && linv[i] < i ) continue;
          if ( linv[i] >= 0 ) places[i].append( places[ linv[i] ] );
          if ( places[i].isize( ) <= 1 ) continue;
          Sort( places[i] );
          cout << "\nline " << i; 
          if ( linv[i] >= 0 ) cout << "/" << linv[i];
          cout << " placed at:\n";
          for ( int j = 0; j < places[i].isize( ); j++ )
          {    double d = places[i][j].first;
               int rid = places[i][j].second, pos = places[i][j].third;
               cout << "[" << j+1 << "] " << rid << "." << pos << ", d = " 
                    << setprecision(2) << 100*d << "%";
               if ( d > flag_thresh * places[i][0].first ) 
               {    cout << " *****";
                    flag[rid][pos] = True;     }
               cout << endl;    }    }

     // Output flagged map.

     for ( int fpass = 1; fpass <= 2; fpass++ )
     {    Ofstream( out, OUT_HEAD + ".out" + ToString(fpass) );
          for ( int mi = 0; mi < nmids; mi++ )
          {    int mid = mids[mi];
               out << "\n" << ( fpass == 1 ? "FULL" : "CLEAN" );
               out << " MAP " << mid << "\n" << endl;
               vec< triple< pair<int,int>, int, double > > aligns = all_aligns[mi];
               vec<int> alen = all_alen[mi];
               if ( fpass == 2 )
               {    EraseIf( aligns, flag[mi] );
                    EraseIf( alen, flag[mi] );    }
               int last_p = -1, last_rp = -1;
               int unmapped = 0;
               for ( int i = 0; i < aligns.isize( ); i++ )
               {    int li = aligns[i].second;
                    int rstart = aligns[i].first.first; 
                    int rstop = aligns[i].first.second;
                    double d = aligns[i].third;
                    int p = Position( maplines_id, li ), rp = -1;
                    if ( linv[li] >= 0 ) rp = Position( maplines_id, linv[li] );
                    const int last_see = 10;
                    if ( p >= 0 )
                    {    if ( last_p >= 0 && p - last_p <= last_see )
                         {    for ( int j = last_p + 1; j < p; j++ )
                              {    out << "+" << maplines[j] 
                                        << " [UNMAPPED]" << endl;
                                   unmapped += lens[ maplines_id[j] ];    }    }
                         out << "+" << maplines[p];
                         out << " [r" << rstart << "-" << rstop << ", d = " 
                              << setprecision(2) << 100*d << "%" << "]";
                         last_p = p;    }
                    else if ( rp >= 0 )
                    {    if ( last_rp >= 0 && last_rp >= rp 
                              && last_rp - rp <= last_see )
                         {    for ( int j = last_rp - 1; j > rp; j-- )
                              // for ( int j = rp - 1; j >= last_rp + 1; j-- )
                              {    out << "-" << maplines[j] << " [UNMAPPED]" 
                                        << endl;
                                   unmapped += lens[ maplines_id[j] ];    }    }
                         out << "-" << maplines[rp];
                         out << " [r" << rstart << "-" << rstop << ", d = " 
                              << setprecision(2) << 100*d << "%" << "]";
                         last_rp = rp;    }
                    else
                    {    out << "line[" << li << "] aligned to r" << rstart << "-" 
                              << rstop << " with d = " << setprecision(2) << 100*d 
                              << "%" << ", len = " << all_alen[mi][i];    }
                    if ( fpass == 1 && flag[mi][i] ) out << "  *****";
                    out << endl;    }
               out << "\n";
               PRINT_TO( out, unmapped );    }    }

     Scram(0);     }

     /*

     // Commented out code to load the raw reads instead of the assembled maps.
     // Hardcoded location.

     vec<vec<int>> bng;
     String bng_data = "/wga/scr4/vendor/bng/BNXfile/all1.bnx";
     fast_ifstream bin(bng_data);
     vec<String> bs;
     while(1)
     {    getline( bin, line );
          if ( bin.fail( ) ) break;
          if ( line.Contains( "#", 0 ) ) continue;
          getline( bin, line );
          bs.push_back(line);    
          if ( bs.size( ) > 30000 ) break; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               }
     bng.resize( bs.size( ) );
     cout << Date( ) << ": parsing" << endl;
     #pragma omp parallel for
     for ( int i = 0; i < bng.isize( ); i++ )
     {    vec<String> s;
          Tokenize( bs[i], '\t', s );
          bng[i].resize( s.size( ) - 2 );
          for ( int j = 2; j < s.isize( ); j++ )
               bng[i][j-2] = int( round( s[j].Double( ) - s[j-1].Double( ) ) );    }
     */
