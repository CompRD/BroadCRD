///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// IsolatedErrorView.  For K = 20, in a given read set, find 2K+1-mers, of the
// form K-1-K, and for each K-K pair, find the middle bases that occur and their 
// quality scores, ignoring bases of quality <= 2.  Report fishy cases, defined as 
// follows for a given K-K pair:
// (a) ignore middle base with the most occurrences;
// (b) ignore middle bases whose quality score sum is less than 30 or exceeds 200;
// (c) ignore middle bases supported by less than three reads.
//
// The fishiness rate is useful because it could be used to calibrate error
// correction algorithms.  For reads aligned to part of the human genome, the 
// observed rate of fishy cases is about one per kb (example chr1:0-80Mb --> 70257).
// However this rate is likely exaggerated by misalignment, and thus one might
// anticipate that the genome-wide rate is substantially lower.  But this is not
// in fact the case (see below).  The most likely explanation is that highly
// repetitive sequences (and artifacts associated with them) dominate in the
// genome-wide case.
//
// Note that read ids are currently stored in an int, which would overflow if there
// were more than ~2 billion reads.
//
// SAMPLE RESULTS
//
// IsolatedErrorView TMP=tmp.xxx SAMPLE=human X=1:0-80M
//
// found 70257 loci
// loading used 1860 seconds
// setting up used 403 seconds
// indexing used 351 seconds
// sorting used 102 seconds
// looking used 96 seconds
// time used = 47.1 minutes, peak mem used = 10.3 GB
//
// IsolatedErrorView SAMPLE=human
//
// found 6261312 loci
// loading used 2409 seconds
// setting up used 16016 seconds
// indexing used 12591 seconds
// sorting used 6329 seconds
// looking used 3432 seconds
// time used = 11.3 hours, peak mem used = 449.9 GB

#define _GLIBCXX_PARALLEL

#include <omp.h>

#include "Basevector.h"
#include "MainTools.h"
#include "PackAlign.h"
#include "ParallelVecUtilities.h"
#include "Qualvector.h"
#include "kmers/KmerRecord.h"
#include "lookup/LookAlign.h"
#include "math/HoInterval.h"
#include "paths/long/CreateGenome.h"
#include "paths/long/DataSpec.h"
#include "paths/long/LoadAndCorrect.h"
#include "paths/long/Logging.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/local/Setup.h"
#include "system/SortInPlace.h"

inline void GetPassIds( const basevector& u, const int L, const int j, 
     const int npasses, int& pass1, int& pass2 )
{    if ( npasses == 4 )
     {    pass1 = u[j];
          pass2 = 3 - u[j+L-1];    }
     else if ( npasses == 16 )
     {    pass1 = 4 * u[j] + u[j+1];
          pass2 = 4 * ( 3 - u[j+L-1] ) + ( 3 - u[j+L-2] );    }
     else if ( npasses == 64 )
     {    pass1 = 16 * u[j] + 4 * u[j+1] + u[j+2];
          pass2 = 16 * ( 3 - u[j+L-1] ) + 4 * ( 3 - u[j+L-2] )
               + ( 3 - u[j+L-3] );    }
     else if ( npasses == 256 )
     {    pass1 = 64 * u[j] + 16 * u[j+1] + 4 * u[j+2] + u[j+3];
          pass2 = 64 * ( 3 - u[j+L-1] ) + 16 * ( 3 - u[j+L-2] )
               + 4 * ( 3 - u[j+L-3] ) + ( 3 - u[j+L-4] );    }
     else
     {    FatalErr("npasses illegal");    }    }

typedef triple< triple< kmer<20>, kmer<20>, char >, int, int > sandwich_pie_20;
typedef byte_pac< 12, 8 > squished_sandwich_pie_20;

int loc1 = 0, loc2 = 5, loc3 = 10, loc4 = 11, loc5 = 12, loc6 = 16;

void Compress( const sandwich_pie_20& x, squished_sandwich_pie_20& y )
{    memcpy( &y + loc1, &x.first.first, 5 );
     memcpy( (char*) &y + loc2, &x.first.second, 5 );
     memcpy( (char*) &y + loc3, &x.first.third, 1 );
     char z(0);
     memcpy( (char*) &y + loc4, &z, 1 );
     memcpy( (char*) &y + loc5, &x.second, 4 );
     memcpy( (char*) &y + loc6, &x.third, 4 );    }

void Uncompress( const squished_sandwich_pie_20& y, sandwich_pie_20& x )
{    memcpy( &x.first.first, (char*) &y + loc1, 5 );
     memcpy( &x.first.second, (char*) &y + loc2, 5 );
     memcpy( &x.first.third, (char*) &y + loc3, 1 );
     memcpy( &x.second, (char*) &y + loc5, 4 );
     memcpy( &x.third, (char*) &y + loc6, 4 );    }

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(SAMPLE, "presently human or rhody");
     CommandArgument_String_OrDefault_Doc(DATASET, "", "as in LongProto");
     CommandArgument_String_OrDefault_Doc(X, "", "region definition, as in LongProto");
     CommandArgument_String_OrDefault_Doc(TMP, "", "scratch directory");
     CommandArgument_Bool_OrDefault(VERBOSE, False);
     EndCommandArguments;

     // Initial set up.

     double clock = WallClockTime( );
     double load_clock = 0, setup_clock = 0, index_clock = 0, sort_clock = 0;
     double look_clock = 0;
     if ( !( SAMPLE == "human" && X == "" ) ) ForceAssert( TMP != "" );

     // Fetch reads.

     cout << Date( ) << ": loading data" << endl;
     load_clock -= WallClockTime( );
     vecbasevector bases;
     vecqualvector quals;
     if ( X == "" )
     {    String dir = "/wga/scr4/human_data/NA12878_H01UJ";
          bases.ReadAll( dir + "/frag_reads_orig.fastb" );
          quals.ReadAll( dir + "/frag_reads_orig.qualb" );    }
     else
     {    Mkdir777(TMP);
          String picard = "/wga/scr4/picard";
          vec<picard_spec> specs;
          String region;
          vec<String> chrlist;
          double SELECT_FRAC = 1.0;
          double COVERAGE_CAP = -1;
          vec<String> regions;
          ParseRegions( X, regions );
          long_data_spec spec( "" );
          long_logging logc( "" );
          SetupIlluminaData( spec, SAMPLE, DATASET, regions, COVERAGE_CAP, False, 
               specs, chrlist, region, SELECT_FRAC, logc );
          LoadIlluminaReads( SAMPLE, specs, region, TMP, logc, False, False );
          bases.ReadAll( "tmp.xxx/frag_reads_orig.fastb" );
          quals.ReadAll( "tmp.xxx/frag_reads_orig.qualb" );    }
     load_clock += WallClockTime( );
     cout << Date( ) << ": loading done, peak mem used = " 
          << setiosflags(ios::fixed) 
          << setprecision(1) << PeakMemUsageBytes( ) / 1000000000.0 
          << resetiosflags(ios::fixed) << " GB" << endl;

     // Define heuristics.

     const int flank = 20;
     const int min_qual = 3;
     const int min_qual_total = 30;
     const int min_group = 3;
     const int max_qual_total = 200;

     // Event counter.

     int xcount = 0;

     // Go through passes.

     const int npasses = 256;
     for ( int pass = 0; pass < npasses; pass++ )
     {
          // Index kmers in the reads.

          setup_clock -= WallClockTime( );
          const int K = flank;
          //typedef triple< kmer<K>, kmer<K>, char > sandwich; // (left,right,middle)
          const int L = 2*K + 1;
          vec<squished_sandwich_pie_20> kmers_plus;
          vec<int> counts( bases.size( ), 0 );
          #pragma omp parallel for
          for ( size_t i = 0; i < bases.size( ); i++ )
          {    const basevector& u = bases[i];
               int count = 0;
               for ( int j = 0; j <= u.isize( ) - L; j++ )
               {    int pass1 = 0, pass2 = 0;
                    GetPassIds( u, L, j, npasses, pass1, pass2 );
                    if ( pass1 < pass || pass2 < pass ) continue;
                    if ( pass1 > pass && pass2 > pass ) continue;
                    counts[i]++;    }    }
          vec<int64_t> starts( bases.size( ) + 1 );
          starts[0] = 0;
          for ( size_t i = 0; i < bases.size( ); i++ )
               starts[i+1] = starts[i] + counts[i];
          kmers_plus.resize( starts.back( ) );
          setup_clock += WallClockTime( ), index_clock -= WallClockTime( );
          #pragma omp parallel for
          for ( size_t i = 0; i < bases.size( ); i++ )
          {    const basevector& u = bases[i];
               int count = 0;
               kmer<K> x1, x1rc, x2, x2rc;
               sandwich_pie_20 sp;
               basevector b1, b2, b1rc, b2rc;
               for ( int j = 0; j <= u.isize( ) - L; j++ )
               {    int pass1 = 0, pass2 = 0;
                    GetPassIds( u, L, j, npasses, pass1, pass2 );
                    if ( pass1 < pass || pass2 < pass ) continue;
                    if ( pass1 > pass && pass2 > pass ) continue;
                    Bool type1 = ( pass1 == pass ), type2 = ( pass2 == pass );
                    int64_t r = starts[i] + count++;
                    x1.SetToSubOf( u, j );
                    x1rc = x1;
                    x1rc.ReverseComplement( );   
                    x2.SetToSubOf( u, j + K + 1 );
                    x2rc = x2;
                    x2rc.ReverseComplement( ); 

                    Bool fw;
                    if ( type1 && !type2 ) fw = True;
                    else if ( !type1 && type2 ) fw = False;
                    else 
                    {    x1.GetBasevector(b1), x2.GetBasevector(b2);
                         x1rc.GetBasevector(b1rc), x2rc.GetBasevector(b2rc);
                         fw = ( make_pair(b1,b2) < make_pair(b2rc,b1rc) );    }

                    char m = u[j+K];
                    char mrc = 3 - m;
                    sp.first = ( fw 
                         ? make_triple(x1,x2,m) : make_triple(x2rc,x1rc,mrc) );
                    sp.second = i; 
                    sp.third = ( fw ? j : -j-1 );
                    Compress( sp, kmers_plus[r] );    }    }
          index_clock += WallClockTime( );
          sort_clock -= WallClockTime( );
          // sortInPlaceParallel( kmers_plus.begin( ), kmers_plus.end( ) );
          ParallelSort(kmers_plus);
          sort_clock += WallClockTime( );

          // Define batches for alignment.

          look_clock -= WallClockTime( );
          int nproc = omp_get_max_threads( );
          vec<int64_t> bstart(nproc+1);
          for ( int i = 0; i <= nproc; i++ )
               bstart[i] = ( (int64_t) kmers_plus.size( ) * i ) / nproc;
          for ( int i = 1; i < nproc; i++ )
          {    int64_t& s = bstart[i];
               while( s > bstart[i-1] 
                    && memcmp( &kmers_plus[s], &kmers_plus[s-1], 10 ) == 0 )
               {    s--;    }    }

          // Look for events.

          #pragma omp parallel for
          for ( int bi = 0; bi < nproc; bi++ )
          for ( int64_t i = bstart[bi]; i < bstart[bi+1]; i++ )
          {    int64_t j;
               for ( j = i + 1; j < bstart[bi+1]; j++ )
               {    if ( memcmp( &kmers_plus[j], &kmers_plus[i], 10 ) != 0 ) 
                         break;    }
               vec<int> count(4,0), qsum(4,0);
               vec< vec<int> > qs(4), idsx(4);
               vec< vec<Bool> > fw(4);
               sandwich_pie_20 x;
               for ( int64_t k = i; k < j; k++ )
               {    Uncompress( kmers_plus[k], x );
                    count[ x.first.third ]++;
     
                    int id = x.second, pos = x.third;
     
                    if ( pos < 0 ) pos = -pos-1;
                    int q = quals[id][pos+K];
     
                    if ( q >= min_qual ) 
                    {    qsum[ x.first.third ] += q;
                         qs[ x.first.third ].push_back(q);    
                         fw[ x.first.third ].push_back( x.third >= 0 );    
                         idsx[ x.first.third ].push_back(id);    }    }

               vec<int> count_sort(count), qsum_sort(qsum);
               vec<int> ids( 4, vec<int>::IDENTITY );
               ReverseSortSync( count_sort, qsum_sort, ids, qs, fw );

               for ( int l = 1; l < 4; l++ )
               {    if ( count_sort[l] >= min_group && qsum_sort[l] >= min_qual_total
                         && qsum_sort[l] <= max_qual_total )
                    {    
                         #pragma omp critical
                         {    ++xcount;    }
                         if (VERBOSE)
                         {    ReverseSort( qs[l] );
                              int fwx = 0, rcx = 0;
                              for ( int m = 0; m < fw[l].isize( ); m++ )
                              {    if ( fw[l][m] ) fwx++;
                                   else rcx++;    }
                              basevector b1, b2;
                              sandwich_pie_20 x;
                              Uncompress( kmers_plus[i], x );
                              x.first.first.GetBasevector(b1);
                              x.first.second.GetBasevector(b2);
                              #pragma omp critical
                              cout << "\n[" << xcount << "] " << as_base( ids[0] ) 
                                   << " --> " << as_base( ids[l] ) << " :: "
                                   << printSeq(count) << " :: " << printSeq(qsum) 
                                   << " :: {" << printSeq(qs[l]) << "}" 
                                   << " :: f/r = " << fwx << "/" << rcx << "\n"
                                   << String( ToString(xcount).size( ) + 3, ' ' )
                                   << "ids = " << idsx[ ids[l] ][0] << ", ..."
                                   << " :: seq = " << b1.ToString( ) 
                                   << char( tolower( as_base( ids[l] ) ) ) 
                                   << b2.ToString( ) << endl;    }    }    }
               i = j - 1;    }
          look_clock += WallClockTime( );    }
     if (VERBOSE) cout << endl;
     cout << Date( ) << ": found " << xcount << " loci" << endl;
     cout << Date( ) << ": loading used " << int(round(load_clock)) << " seconds" 
          << endl;
     cout << Date( ) << ": setting up used " << int(round(setup_clock)) 
          << " seconds" << endl;
     cout << Date( ) << ": indexing used " << int(round(index_clock)) << " seconds" 
          << endl;
     cout << Date( ) << ": sorting used " << int(round(sort_clock)) << " seconds" 
          << endl;
     cout << Date( ) << ": looking used " << int(round(look_clock)) << " seconds" 
          << endl;
     cout << Date( ) << ": done, time used = " << TimeSince(clock)
          << ", peak mem used = " << setiosflags(ios::fixed) << setprecision(1)
          << PeakMemUsageBytes( ) / 1000000000.0 << resetiosflags(ios::fixed)
          << " GB" << endl;
     Scram(); // ******************************************************************

     // Old alignment-based approach, possibly defunct (rest of code).

     SystemSucceed( "QueryLookupTable K=12 MM=12 MC=0.15 "
          "SEQS=tmp.xxx/frag_reads_orig.fastb QUALS=tmp.xxx/frag_reads_orig.qualb "
          "L=genome.lookup SMITH_WAT=True PARSEABLE=True > aligns" );

     // Load genome.
     
     vecbasevector genome( "genome.fastb" );

     // Fetch aligns.

     vec<look_align> aligns;
     LoadLookAligns( "aligns", aligns );

     // Find bad sites.

     vec< triple< pair<int,char>, int, look_align > > bads;
     for ( int i = 0; i < aligns.isize( ); i++ )
     {    const look_align& la = aligns[i];
          int id = la.query_id;
          basevector b = bases[id];
          qualvector q = quals[id];
          if ( la.rc1 ) 
          {    b.ReverseComplement( );
               q.ReverseMe( );    }
          const align& a = la.a;
          vec<ho_interval> p1, p2;
          a.PerfectIntervals1( b, genome[0], p1 );
          a.PerfectIntervals2( b, genome[0], p2 );
          for ( int j = 0; j < p1.isize( ) - 1; j++ )
          {    if ( p1[j].Length( ) < flank ) continue;
               if ( p1[j+1].Length( ) < flank ) continue;
               if ( p1[j+1].Start( ) - p1[j].Stop( ) != 1 ) continue;
               if ( p2[j+1].Start( ) - p2[j].Stop( ) != 1 ) continue;

               int pos1 = p1[j].Stop( ), pos2 = p2[j].Stop( );
               int qual = q[pos1];
               if ( qual < min_qual ) continue;
               char ref_base = as_base( genome[0][pos2] );
               char read_base = as_base( b[pos1] );

               bads.push( make_pair(pos2,read_base), pos1, la );    }    }

     Sort(bads);
     int count = 0;
     for ( int i = 0; i < bads.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < bads.isize( ); j++ )
               if ( bads[j].first != bads[i].first ) break;
          int qual_total = 0;
          for ( int k = i; k < j; k++ )
          {    int pos1 = bads[k].second;
               const look_align& la = bads[k].third;
               int id = la.query_id;
               qualvector q = quals[id];
               if ( la.rc1 ) q.ReverseMe( );
               qual_total += q[pos1];    }
          if ( j - i >= min_group && qual_total >= min_qual_total )
          {    int pos2 = bads[i].first.first;
               cout << "\n[" << ++count << "] pos2 = " << pos2 << "\n";
               for ( int k = i; k < j; k++ )
               {    int pos1 = bads[k].second;
                    const look_align& la = bads[k].third;
                    int id = la.query_id;
                    basevector b = bases[id];
                    qualvector q = quals[id];
                    if ( la.rc1 ) 
                    {    b.ReverseComplement( );
                         q.ReverseMe( );    }
                    int qual = q[pos1];

                    char read_base = as_base( b[pos1] );
                    char ref_base = as_base( genome[0][pos2] );

                    if ( la.rc1 ) pos1 = b.isize( ) - pos1 - 1;
                    String dir = ( la.rc1 ? "rc" : "fw" );

                    PRINT6( id, pos1, dir, qual, ref_base, read_base );    }    }
          i = j - 1;    }    }

/*
The MUC1 case that we would like to recover:
CCGGGCTCCACCGCCCCCCCgGCCCACGGTGTCACCTCGGC
ALIGN 0 1007 - 16 37     28
ALIGN 0 1179 + 52 73      5
ALIGN 0 2940 - 116 137    7
ALIGN 0 2980 - 101 122    7
Change g to a: 102 hits.
*/

/*
X=11:12M-12.1M

[1] G --> C :: 0,3,73,0 :: 0,69,2442,0 :: {35,34} :: f/r = 2/0
    ids = 5919, ... :: seq = TTCCAATCAATAGAAAAAGAcGGAATCCTCCCTAACTCATT
    SearchHuman BAMS=free: >= 72.

[2] G --> A :: 3,0,44,0 :: 99,0,1360,0 :: {39,37,23} :: f/r = 0/3
    ids = 13090, ... :: seq = GGCTTTGGTATCAGAATGATaCTGGCCTCATAAAATGAGTT
    SearchHuman BAMS=free: >= 72.

[3] C --> A :: 4,37,0,0 :: 103,1195,0,0 :: {39,35,29} :: f/r = 0/3
    ids = 645, ... :: seq = ATGTATTTCAAAATAATAAGaGCTATCTATGACAAACCCAC
    SearchHuman BAMS=free: >= 144.

[4] C --> A :: 3,86,0,0 :: 66,2606,0,0 :: {33,33} :: f/r = 2/0
    ids = 665, ... :: seq = AAGCCAGGCAGAGACACAACaAAAAAAGAGAATTTTAGACC

[5] A --> G :: 11,0,9,0 :: 231,0,188,0 :: {30,30,28,28,26,26,20} :: f/r = 7/0
    ids = 67, ... :: seq = TCCAGCCTGGGCGACAGAGCgAGACTCCTTCTCAAAAAAAA

[6] A --> T :: 50,0,0,3 :: 1615,0,0,71 :: {34,30,7} :: f/r = 1/2
    ids = 14890, ... :: seq = TACAATGAGAACTACAAAACtCTGCTGAAAGAAATCATAGA

[7] T --> G :: 0,0,3,46 :: 0,0,30,1353 :: {16,7,7} :: f/r = 3/0
    ids = 1850, ... :: seq = GCGCAGGTCCCTGGGAAGAGgTCTAGAAAGATTCATATTGC
    SearchHuman BAMS=free: 3.

[8] C --> T :: 1,73,0,3 :: 34,2380,0,95 :: {35,35,25} :: f/r = 3/0
    ids = 3899, ... :: seq = GGATTACATTTATTGATTTGtGTATATTGAACCAGCCTTGC

[9] T --> G :: 0,1,3,36 :: 0,0,43,1080 :: {29,7,7} :: f/r = 3/0
    ids = 354, ... :: seq = ACGTAGCTCTCGTGCCTTGGgTTTCAGCTCCATCAGGTCCT

[10] G --> T :: 0,0,38,3 :: 0,0,1273,64 :: {39,15,10} :: f/r = 3/0
     ids = 5215, ... :: seq = CACAGAGATGTAATGCACATtACAACTAACTAATTTGACAT

[11] G --> A :: 3,0,64,0 :: 96,0,2078,0 :: {34,33,29} :: f/r = 1/2
     ids = 11274, ... :: seq = TCAATGTACAAAAATCACAAaCATTCTTATACACCAATAAC
     SearchHuman BAMS=free: >= 36.

[12] G --> T :: 0,0,40,3 :: 0,0,1090,51 :: {36,15} :: f/r = 1/1
     ids = 20442, ... :: seq = TAGAACCCAGGTAGTCGGGGtATGTAGAGTGCACGCATATA

[13] C --> A :: 3,49,0,0 :: 34,1709,0,0 :: {27,7} :: f/r = 0/2
     ids = 8423, ... :: seq = TCATCTATTCAATCAGACTCaTTTGAATAGGCATTTGGCAG

[14] T --> G :: 0,0,4,13 :: 0,0,30,425 :: {12,8,5,5} :: f/r = 4/0
     ids = 3522, ... :: seq = GGCACATTGCAGGCGGGGGGgAATACACTCTGCCCCAGGAG

[15] T --> G :: 0,1,3,41 :: 0,5,30,1079 :: {20,5,5} :: f/r = 3/0
     ids = 679, ... :: seq = TGAGTACTATAGCTGGAGGGgATAGATGGAAATAGTACTCA
*/
