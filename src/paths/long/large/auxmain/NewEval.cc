///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "FastIfstream.h"
#include "FetchReads.h"
#include "MainTools.h"
#include "PrintAlignment.h"
#include "TokenizeString.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "paths/HyperBasevector.h"
#include "paths/long/ReadPath.h"
#include "paths/long/fosmid/Fosmids.h"
#include "paths/long/large/GapToyTools.h"
#include "paths/long/large/Lines.h"
#include "paths/long/large/MakeGaps.h"
#include "random/Shuffle.h"

#include "random/Random.h"

void TestEvent( const HyperBasevector& hb, const vecbitvector& flagged,
     const vec<int>& inv, const vec<vec<vec<vec<int>>>>& lines, const vec<int>& tol,
     const vec<int>& lens, const vec<vec<covcount>>& covs,
     const ReadPathVec& paths, const VecULongVec& paths_index,
     const vecbasevector& genome, vecString& Gnames, const int g, 
     const vecbasevector& bases, const VecPQVec& quals, const int gstop1, 
     const int gstart2, const int estop1, const int estart2, const int e, 
     Bool& categorized, int& cats, int& subs, int& indel1, vec<int>& indels, 
     vec<int>& subclust, const Bool TEST_SUPPORT )
{

     const Bool hide_unanchored = True;

     const int flank_add = 30;
     int flank = max( { gstop1-gstart2, estop1-estart2, 0 } )/2 + flank_add;
     int gstop1x = gstop1 - flank;
     int gstart2x = gstart2 + flank;
     int estop1x = estop1 - flank;
     int estart2x = estart2 + flank;
     const basevector &E = hb.EdgeObject(e), &G = genome[g];
     if ( gstop1x < gstart2x && estop1x < estart2x && gstop1x >= 0 
          && gstart2x <= G.isize( )
          && estop1x >= 0 && estart2x <= E.isize( ) )
     {    basevector b1( E, estop1x, estart2x - estop1x );
          basevector b2( genome[g], gstop1x, gstart2x - gstop1x );
          alignment al;
          SmithWatAffine( b1, b2, al );
          align a = al;
          if ( a.pos1( ) != 0 || a.Pos1( ) != b1.isize( )
               || a.pos2( ) != 0 || a.Pos2( ) != b2.isize( ) )
          {    return;    }
          cout << "\nEVENT from " << estop1x << " to " << estart2x << " on edge "
               << e << ", cov = " << covs[0][e].Cov( ) << endl;

          int li = tol[e];
          Bool reg = False;
          if ( lens[li] >= 2000 )
          {    for ( int j = 0; j < lines[li].isize( ); j += 2 )
                    if ( lines[li][j][0][0] == e ) reg = True;    }
          int ne = hb.Bases(e);
          if ( ne >= 2000 ) reg = True; // very flaky!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if ( !reg ) cout << "This event is not at a regular site." << endl;

          Bool flag = False;
          if ( flagged.size( ) > 0 )
          {    for ( int j = estop1x; j < estart2x; j++ )
                    if ( flagged[e][j] ) flag = True;    }
          if (flag) cout << "flagged" << endl;
          PrintVisualAlignment( False, cout, b1, b2, a );
          cats++;
          categorized = True;    
          vec<int> mgg = a.MutationsGap1Gap2( b1, b2 );
          if ( mgg[0] == 1 && Sum(mgg) == 1 )
          {    subs++;
               cats--;    }
          else if ( Sum(mgg) == 1 )
          {    indel1++;
               cats--;    }    
          else if ( mgg[0] == 0 && a.Nblocks( ) == 2 )
          {    indels.push_back( Abs( a.Gaps(1) ) );
               cats--;    }
          else if ( mgg[1] == 0 && mgg[2] == 0 )
          {    subclust.push_back( mgg[0] );    
               cats--;    }    

          // Heuristics.

          const int min_diff = 30;
          const int max_qsum = 50;
          const double min_ratio = 5.0;
          const int min_support_ratio = 5;

          // Test support.

          if ( !TEST_SUPPORT ) return;
          int ass_support = 0, ref_support = 0;

          for ( int pass = 1; pass <= 2; pass++ )
          {    int ex = ( pass == 1 ? e : inv[e] );
               for ( int j = 0; j < (int) paths_index[ex].size( ); j++ )
               {    int64_t id = paths_index[ex][j];
                    const ReadPath& p = paths[id];

                    // Find estart, the starting position of the read on the
                    // edge ( e or inv[e] ), and estop, the inferred stop position
                    // of the read on that same edge.

                    int estart = p.getOffset( );
                    for ( int l = 0; l < (int) p.size( ); l++ )
                    {    if ( p[l] == ex ) break;
                         estart -= hb.Kmers( p[l] );    }
                    int estop = estart + bases[id].isize( );

                    // In the reverse complement case ( inv[e] ), flip the start
                    // and stop positions.

                    if ( pass == 2 )
                    {    int restart = ne - estop, restop = ne - estart;
                         estart = restart, estop = restop;    }

                    // Check for out of bounds.  Possibly too stringent.

                    if ( estart > estop1x || estop < estart2x ) continue;

                    // Get the read, and in the reverse complement case, flip it.

                    basevector b = bases[id];
                    qualvector q;
                    quals[id].unpack(&q);
                    if ( pass == 2 )
                    {    b.ReverseComplement( );
                         q.ReverseMe( );    }

                    // Test anchoring.

                    vec<int> anchors;
                    const int L = 20;
                    const double min_ratio = 3.0;
                    const int min_add = 8;
                    const int radius = 5;
                    // effectively executing two passes:
                    // for ( int i = 0; i <= ne - L; i++ )
                    // for ( int i = ne - L; i >= 0; i-- )
                    for ( int xpass = 1; xpass <= 2; xpass++ )
                    {    int start = Max( 0, estart ), stop = Min( ne, estop );
                         int i = ( xpass == 1 ? start : stop - L );
                         while(1)
                         {    if ( i > stop - L || i < start ) break;
                              int baseline = 0, alt = 1000000000;
                              for ( int r = -radius; r <= +radius; r++ )
                              {    if ( i + r < start || i + r + L > stop ) continue;
                                   int qsum = 0;
                                   for ( int l = 0; l < L; l++ )
                                   {    if ( b[i+l-start] 
                                             != hb.EdgeObject(e)[i+l+r] )
                                        {    qsum += q[i+l-start];    }    }
                                   if ( r == 0 ) baseline = qsum;
                                   else alt = Min( alt, qsum );    }
                              if ( alt >= min_ratio * baseline 
                                   && alt - baseline >= min_add ) 
                              {    anchors.push_back(i);
                                   break;    }
                              i += ( xpass == 1 ? +1 : -1 );    }    }
                    Bool anchored = False;
                    int low = -1, high = -1;
                    if ( anchors.nonempty( ) )
                    {    low = Min(anchors), high = Max(anchors) + L;
                         if ( low + L <= estop1 && high - L >= estart2 )
                              anchored = True;    }

                    // Compute quality score sums.  We have to work from the
                    // beginning of the read, which requires separate calculations
                    // for the two passes.
     
                    int ass_qsum = 0, ref_qsum = 0;
                    int nb = b.isize( );
                    int r;
                    if ( pass == 1 )
                    {    r = estop1x - estart;
                         for ( int l = r; l < nb; l++ )
                         {    int epos = estop1x + l - r, gpos = gstop1x + l - r;
                              if ( epos >= ne || gpos >= G.isize( ) ) break;
                              if ( b[l] != E[epos] ) ass_qsum += q[l];
                              if ( b[l] != G[gpos] ) ref_qsum += q[l];    }    }
                    else
                    {    r = estop - estart2x;
                         for ( int l = nb - r - 1; l >= 0; l-- )
                         {    int epos = estart2x - 1 - ((nb-r-1)-l);
                              int gpos = gstart2x - 1 - ((nb-r-1)-l);
                              if ( epos < 0 || gpos < 0 ) break;
                              if ( b[l] != E[epos] ) ass_qsum += q[l];
                              if ( b[l] != G[gpos] ) ref_qsum += q[l];    }    }
     
                    // Test alignment.
     
                    if ( Abs( ass_qsum - ref_qsum ) < min_diff ) continue;
                    if ( ass_qsum > max_qsum && ref_qsum > max_qsum ) continue;
                    /*
                    if ( min_ratio * ass_qsum >= ref_qsum
                         && min_ratio* ref_qsum >= ass_qsum )
                    {    continue;    }
                    */

                    // Hide unanchored.

                    if ( ass_qsum < ref_qsum && !anchored && hide_unanchored )
                         continue;

                    // Record support.

                    if ( ass_qsum < ref_qsum ) 
                    {    if (anchored) ass_support++;    }
                    else ref_support++;

                    // Print support.
     
                    cout << "\npossible support from read " << id
                         << " " << ( pass == 1 ? "fw" : "rc" )
                         << ", starting at " << estart << endl;
                    PRINT2( ass_qsum, ref_qsum );
                    if ( anchors.nonempty( ) ) PRINT2( low, high );
                    if ( !anchored ) cout << "not anchored" << endl;
                    if ( ass_qsum == 0 ) continue;
     
                    // Print alignments.
     
                    avector<int> gaps(1), lengths(1);
                    qualvector qx;
                    int ng = G.isize( );
                    cout << "\nalignment to assembly:\n";
                    if ( pass == 1 )
                    {    basevector rx( b, r, nb - r );
                         qx.SetToSubOf( q, r, nb - r );
                         basevector a( E, estop1x, ne - estop1x );
                         gaps(0) = 0, lengths(0) = rx.size( );
                         align al( 0, 0, gaps, lengths );
                         if ( rx.size( ) > a.size( ) )
                         {    cout << "can't display" << endl;
                              continue;    }
                         PrintVisualAlignment( True, cout, rx, a, al, qx );    
                         cout << "\nalignment to reference:\n";
                         basevector gx( G, gstop1x, ng - gstop1x );
                         align alx( 0, 0, gaps, lengths );
                         if ( rx.size( ) > gx.size( ) )
                         {    cout << "can't display" << endl;
                              continue;    }
                         PrintVisualAlignment( True, cout, rx, gx, alx, qx );    }
                    else
                    {    basevector rx( b, 0, nb-r );
                         qx.SetToSubOf( q, 0, nb-r );
                         if ( nb-r > estart2x )
                         {    cout << "can't display" << endl;
                              continue;    }
                         basevector a( E, estart2x - (nb-r), nb-r );
                         gaps(0) = 0, lengths(0) = rx.size( );
                         align al( 0, 0, gaps, lengths );
                         if ( rx.size( ) > a.size( ) )
                         {    cout << "can't display" << endl;
                              continue;    }
                         PrintVisualAlignment( True, cout, rx, a, al, qx );    
                         cout << "\nalignment to reference:\n";
                         basevector gx( G, gstart2x - (nb-r), nb-r );
                         align alx( 0, 0, gaps, lengths );
                         if ( rx.size( ) > gx.size( ) )
                         {    cout << "can't display" << endl;
                              continue;    }
                         PrintVisualAlignment( 
                              True, cout, rx, gx, alx, qx );    }    }    }

          // Decide if we accept assembly variant.

          Bool accept = ( ass_support >= min_support_ratio 
               * Max( double(ref_support), 1.0 ) );

          if ( !reg ) accept = False;

          cout << "\nsupport: assembly[" << ass_support << "], reference["
               << ref_support << "]\n\n";    

          if ( ref_support > ass_support )
          {    cout << "reference support exceeds assembly support for event "
                    << "from " << estop1x << " to " << estart2x << " on edge " << e 
                    << ( !flag ? ", NOT FLAGGED" : "" ) << "\n" << endl;    }

          if (accept) 
          {    cout << "accepting assembly variant " << "from " << estop1x << " to "
                    << estart2x << " on edge " << e << ( flag ? ", FLAGGED" : "" ) 
                    << ", coverage = " << covs[0][e].Cov( ) << "\n" << endl;

               // Define variant for VCF.

               int A_start = 0, A_stop = b1.size( );
               int G_start = gstop1x, G_stop = gstart2x;
               while( A_start < A_stop && G_start < G_stop 
                    && b1[A_start] == G[G_start] )
               {    A_start++;
                    G_start++;    }
               while( A_start < A_stop && G_start < G_stop 
                    && b1[A_stop-1] == G[G_stop-1] )
               {    A_stop--;
                    G_stop--;    }
               if ( G_start == G_stop || A_start == A_stop )
               {    A_start--;
                    G_start--;    }
               ForceAssertGe( A_start, 0 );
               ForceAssertGe( G_start, 0 );
               basevector v1( G, G_start, G_stop - G_start );
               basevector v2( b1, A_start, A_stop - A_start );
               cout << "at " << Gnames[g] << "." << G_start+1 << ", replace "
                    << v1.ToString( ) << " by " << v2.ToString( ) 
                    << "\n" << endl;    
               cout << "VCF: " << Gnames[g] << "\t" << G_start+1 << "\t"
                    << "." << "\t" << v1.ToString( ) << "\t" << v2.ToString( )
                    << "\t" << 100 << "\t" << "PASS" << "\t" << "." << "\n" 
                    << endl;    }    }    }

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault_Doc(DIR, ".",
          "looks for DIR/a.{fastb,inv,paths,paths.inv}");
     CommandArgument_Bool_OrDefault(DETAILS, False);
     CommandArgument_Bool_OrDefault_Doc(F, False,
          "use Fosmid reference sequences rather than genome.fastb");
     CommandArgument_Int_OrDefault_Doc(SHOW_COUNT, 0,
          "note edges have less reads supporting them than this number");
     CommandArgument_Bool_OrDefault_Doc(FIND_PATHS, True,
          "try to find paths between perfect stretches");
     CommandArgument_Bool_OrDefault_Doc(TEST_SUPPORT, False,
          "assay read support for differences, and generate VCF");
     EndCommandArguments;

     // Heuristics.

     const int max_paths = 10;
     const int max_iterations = 1000;
     const int max_gap = 300;

     // Track stats.

     int cats = 0, uncats = 0, extreme_gaps = 0, uncircles = 0, missing = 0;
     int subs = 0, indel1 = 0;
     vec<int> indels, subclust;

     // Output VCF header lines.

     cout << "VCF: ##fileformat=VCF4.1" << endl;
     cout << "VCF: " << "#CHROM" << "\t" << "POS" << "\t" << "ID" << "\t" << "REF"
          << "\t" << "ALT" << "\t" << "QUAL" << "\t" << "FILTER"
          << "\t" << "INFO" << endl;

     // Load genome.

     vecbasevector genome;
     vecString Gnames; 
     if (F)
     {    vec<int> fosmids = AllFosmids( );
          genome.reserve( fosmids.size( ) );
          Gnames.reserve( fosmids.size( ) );
          for ( int i = 0; i < fosmids.isize( ); i++ )
          {    vecbasevector g;
               int f = fosmids[i];
               FetchReads( g, 0, "/wga/dev/references/Homo_sapiens/"
                    "NA12878_Fosmid_Pool.regions.fin/fos." + ToString(f) 
                    + ".fasta" );
               if ( g.size( ) == 1 )
               {    genome.push_back( g[0] ); 
                    Gnames.push_back( "F" + ToString(f) );    }    }    }
     else
     {    genome.ReadAll( DIR + "/../genome.fastb" );
          if ( IsRegularFile( DIR + "/../genome.names" ) )
          {    fast_ifstream in( DIR + "/../genome.names" );
               String line;
               while(1)
               {    getline( in, line );
                    if ( in.fail( ) ) break;
                    Gnames.push_back(line);    }    }
          else
          {    for ( int g = 0; g < (int) genome.size( ); g++ )
                    Gnames.push_back( ToString(g) );    }    }

     // Load data.

     HyperBasevector hb;
     BinaryReader::readFile( DIR + "/a.hbv", &hb );
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     vec<int> inv;
     BinaryReader::readFile( DIR + "/a.inv", &inv );
     ReadPathVec paths;
     VecULongVec paths_index;
     vecbasevector bases;
     VecPQVec quals;
     vecbitvector flagged;
     if (TEST_SUPPORT)
     {    paths.ReadAll( DIR + "/a.paths" );
          invert( paths, paths_index, hb.EdgeObjectCount( ) );
          bases.ReadAll( DIR + "/../data/frag_reads_orig.fastb" );
          quals.ReadAll( DIR + "/../data/frag_reads_orig.qualp" );
          if ( IsRegularFile( DIR + "/a.flagged" ) )
               flagged.ReadAll( DIR + "/a.flagged" );    }
     vec<vec<covcount>> covs;
     BinaryReader::readFile( DIR + "/a.covs", &covs );
     vec<vec<vec<vec<int>>>> lines;
     BinaryReader::readFile( DIR + "/a.lines", &lines );
     vec<int> tol, lens;
     GetTol( hb, lines, tol );
     GetLineLengths( hb, lines, lens );

     // Load edge counts.

     vec<int> xcount;
     if ( SHOW_COUNT > 0 )
     {    xcount.resize( hb.E( ), 0 );
          Ifstream( in, DIR + "/a.counts" );
          String line;
          getline( in, line );
	  ForceAssert( !in.fail( ) );
          vec<String> x;
          for ( int e = 0; e < hb.E( ); e++ )
 	  {    getline( in, line );
	       Tokenize( line, x );
               for ( int j = 1; j < x.isize( ); j++ )
                    xcount[e] += x[j].Int( );    }    }

     // Align to genome.

     vec< triple< pair<int,int>, pair<int,int>, int > > perfs;
     vec<perf_place> places;
     AlignToGenomePerf( hb, genome, perfs, places );

     // Go through each genome contig.

     vec<vec<perf_place>> pl( genome.size( ) );
     for ( int i = 0; i < places.isize( ); i++ )
          pl[ places[i].g ].push_back( places[i] );
     for ( int g = 0; g < (int) genome.size( ); g++ )
     {    const vec<perf_place>& P = pl[g];

          // Form perfect placements output.

          vec<vec<String>> rows;
          vec<String> row;
          row.push_back( "g.#", "gstart", "gstop", "len", "estart", "e" );
          rows.push_back(row);
          ostringstream out;
          int count = 0;
          for ( int i = 0; i < P.isize( ); i++ )
          {    vec<String> row;
               String name;
               if ( Gnames.size( ) == 0 ) name = ToString( P[i].g );
               else name = Gnames[ P[i].g ];
               row.push_back( name + "." + ToString(++count) );
               row.push_back( ToString( P[i].gstart ) );
               row.push_back( ToString( P[i].gstart + P[i].len ) );
               row.push_back( ToString( P[i].len ) );
               row.push_back( ToString( P[i].estart ) );
               ostringstream o;
               o << printSeq( P[i].e );
               row.push_back( ToString( o.str( ) ) );
               rows.push_back(row);    }
          PrintTabular( out, rows, 2 );
          vec<String> xlines;
          Tokenize( out.str( ), '\n', xlines );
     
          // Print, along with categorizations.

          cout << "\n==================================================="
               << "=============================\n";
          cout << "\n" << xlines[0] << endl;
          if ( xlines.size( ) == 1 )
          {    missing++;
               cout << "\n[MISSING]\n\n";
               continue;    }
          if ( P[0].Gstart( ) > 0 ) 
          {    cout << "\n[GAP AT BEGINNING]\n\n";
               extreme_gaps++;    }
          cout << "\n" << xlines[1] << endl;
          for ( int i = 1; i < P.isize( ); i++ )
          {    int gstop1 = P[i-1].Gstop( ), gstart2 = P[i].Gstart( );
               int estop1 = P[i-1].Estop(hb), estart2 = P[i].Estart( );
               int e = P[i].E( ).front( );
               Bool categorized = False;
               const int flank_add = 30;
               if ( Abs( gstart2 - gstop1 ) <= max_gap && P[i-1].E( ).back( ) == e )
               {    TestEvent( hb, flagged, inv, lines, tol, lens, covs, paths, 
                         paths_index, genome, Gnames, g, bases, quals, gstop1, 
                         gstart2, estop1, estart2, e, categorized, cats, subs, 
                         indel1, indels, subclust, TEST_SUPPORT );    }

               // Try to classify gap by finding paths between.

               int flank = flank_add;
               int e1 = P[i-1].E( ).back( ), e2 = P[i].E( ).front( );
               int eleft = hb.Bases(e1) - estop1, eright = estart2;
               int K = hb.K( );
               if ( FIND_PATHS && !categorized && eleft + flank <= K-1 
                    && eright + flank <= K-1 )
               {    int v = to_right[e1], w = to_left[e2];
                    vec<vec<int>> paths;
                    Bool ok = hb.EdgePaths( to_left, to_right, v, w, 
                         paths, -1, max_paths, max_iterations );
                    if (ok)
                    {    for ( int j = 0; j < paths.isize( ); j++ )
                         for ( int i = 0; i < paths[j].isize( ); i++ )
                              if ( hb.Bases( paths[j][i] ) == 0 ) ok = False;    }
                    if (ok)
                    {    for ( int j = 0; j < paths.isize( ); j++ )
                         {    if ( j == 0 ) cout << endl;
                              cout << "path " << j+1 << " of " << paths.size( )
                                   << endl;
                              int gleft = gstop1 - flank, gright = gstart2 + flank;

                              // Working around a bug....
                              if ( !( gleft < gright ) ) continue;

                              basevector b1 = hb.Cat( paths[j] );
                              int trim1 = K-1-eleft-flank, trim2 = K-1-eright-flank;
                              if ( trim1 + trim2 > b1.isize( ) )
                              {    ok = False;
                                   continue;    }
                              b1.SetToSubOf(b1, trim1, b1.isize( ) - trim1 - trim2);
                              basevector b2( genome[g], gleft, gright - gleft );
                              alignment al;
                              SmithWatAffine( b1, b2, al );
                              align a = al;
                              if ( !( a.pos1( ) == 0 && a.Pos1( ) == b1.isize( )
                                   && a.pos2( ) == 0 && a.Pos2( ) == b2.isize( ) ) )
                              {    ok = False;
                                   continue;    }
                              PrintVisualAlignment( 
                                   False, cout, b1, b2, a );    }    }    }

               if ( !categorized ) 
               {    cout << "\n[GAP: " << gstop1 << "-" << gstart2 << "]\n\n";
                    uncats++;    }
               cout << xlines[i+1] << endl;    

               if ( SHOW_COUNT > 0 )
               {    for ( int j = 0; j < P[i].E( ).isize( ); j++ )
                    {    int e = P[i].E( )[j];
                         if ( xcount[e] < SHOW_COUNT )
                         {    cout << "\nNote: edge " << e << " has only "
                                   << xcount[e] << " reads supporting it." 
                                   << endl;    }    }    }

               if ( i == P.isize( ) - 1 )
               {    if ( P[i].Gstop( ) < genome[g].isize( ) )
                    {    cout << "\n[GAP AT END]\n";
                         extreme_gaps++;    }
                    else if ( !F && P[0].Gstart( ) == 0 )
                    {    Bool circled = False;
                         int e1 = P.back( ).E( ).back( );
                         int e2 = P.front( ).E( ).front( );

                         if ( e1 == e2 )
                         {    if ( P.back( ).Estop(hb) == P.front( ).Estart( ) )
                                   circled = True;    }
                         else if ( to_right[e1] == to_left[e2] )
                         {    if ( hb.EdgeLengthBases(e1) - P.back( ).Estop(hb) 
                                   + P.front( ).Estart( ) == hb.K( ) - 1 )
                              {    circled = True;    }    }

                         if ( !circled )
                         {    cout << "\n[BAD CIRCLE]\n";
                              uncircles++;    }    }    }    }    }

     // Print perfects.

     if (DETAILS)
     {    cout << "\nDETAILS\n";
          vec<vec<String>> rows;
          vec<String> row = { "g.gstart", "e.estart", "len", "glen", "elen" };
          rows.push_back(row);
          for ( int i = 0; i < perfs.isize( ); i++ )
          {    vec<String> row;
               int g = perfs[i].first.first;
               int e = perfs[i].second.first;
               row.push_back( 
                    ToString(g) + "." + ToString( perfs[i].first.second ) );
               row.push_back( 
                    ToString(e) + "." + ToString( perfs[i].second.second ) );
               row.push_back( ToString( perfs[i].third ) );
               row.push_back( ToString( genome[g].size( ) ) );
               row.push_back( ToString( hb.EdgeLengthBases(e) ) );
               rows.push_back(row);    }
          PrintTabular( cout, rows, 2 );    }

     // Summarize.

     cout << "\nSUMMARY\n";
     if ( subs > 0 ) cout << "subs = " << subs << endl;
     if ( subclust.nonempty( ) )
     {    Sort(subclust);
          cout << "substitution clusters of size " 
               << printSeq(subclust) << endl;    }
     if ( indel1 > 0 ) cout << "single-base indels = " << indel1 << endl;
     if ( indels.nonempty( ) )
     {    Sort(indels);
          cout << "indels of size " << printSeq(indels) << endl;    }
     if ( cats > 0 ) cout << "other aligned = " << cats << endl;
     if ( uncats > 0 ) cout << "gaps = " << uncats << endl;
     if ( extreme_gaps > 0 ) cout << "extreme gaps = " << extreme_gaps << endl;
     if ( uncircles > 0 ) cout << "bad circles = " << uncircles << endl;
     if ( missing > 0 ) cout << "completely missing = " << missing << endl;

     // Compute N50 perfect stretch.

     cout << "\nN50 perfect stretch = " << N50PerfectStretch( hb, genome, F ) << endl;

     // Done.
     
     Scram(0);    }
