///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Try to capture the sequence within a gap.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "Intvector.h"
#include "MainTools.h"
#include "paths/HyperBasevector.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/RemodelGapTools.h"
#include "paths/Unipath.h"
#include "paths/long/LongReadsToPaths.h"
#include "paths/long/ReadPath.h"
#include "paths/long/large/tenx/TenxDirs.h"

int main( )
{    RunTime( );
     
     int N = 1;

     // Control.

     const Bool tenx = False;

     // Define examples.

     vec<int64_t> pids;

     // #1
     // const int e1 = 9776648, e2 = 9727806, K = 80; // captured, closed

     // #2
     // const int e1 = 9695423, e2 = 9734276, K = 48; // captured, closed

     // #3
     const int e1 = 7988288, e2 = 9920679, K = 24;    // uncaptured, closed
     pids = {9207508,113885979,185087549,
          270372233,386904880,392793776,415390484,102007768,345235013,39958119};
     vec<int64_t> pids2( {227426911,66147254,80335964,109897383,125065499,153005786,
          195371131,213760442,230717026,263422260,276157194,285771365,293842023,
          302257825,308746894,382978287,408651912,438629761,159595663,222759187} );
     // pids.append(pids2);

/*

seems like read 227771958 should be aligned to edge 7988288:
(first 80-mer multiply placed but first 120-mer is placed essentially uniquely)
518832,578078,7988288

0fw vs 0, 3 mismatches/0 indels (of 251), from 0-186 to 3292-3478 (of 3478)

33333333333333333333333333332333333333333333333333333333332333333333333333333333
44444888888888888888688888882788888888888888888888888888882778788888887888888778

                                                                                
TAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACAAGGTCAGGAGATCGAGACCATCCTGGCTAACACAGTGAAA
TAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACAAGGTCAGGAGATCGAGACCATCCTGGCTAACACAGTGAAA


31233333333333331333333333333333333123333333332333333333333133233213239233333333
8130777767788846167584788887888777509747784888978887777478710754751074 378461658

0123456789012345678901234567890123456789
                                                                  *             
CCCCGTCTCTACTAAAAATACAAAAAATTAGCTGGGCGAGGTGGCGGGCGCCTGTAGTCCCAGCTAATCGGGAGGCTGAG
CCCCGTCTCTACTAAAAATACAAAAAATTAGCTGGGCGAGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAG


32223339212323321211113123
7663446 357414406586574456

                    *  *  
GCAGGAGAATGGCGTGAACCTCGTGG
GCAGGAGAATGGCGTGAACCCCGGGG

matches ref 409, missing bases are 3965429 -- 3965803

GCGGAGCCTGCAGTGAGCCGAGATCGCGCCACTGCACTCCGGCCTGGGCGACAGCGAGACTCCGTCTCAAAAAAAAAAAA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTATTGACATGTACATTATCTATTCTGTAATGAGGCCATCC
CTCCTAGTTTCCATTGCAGAGGATTGGATCTGGAAATTGTGTTACTAAGAAAAATGCAGGAGAAGGTTTGAGTGTCCCTA
TTCCCATATGTAGTGGATTGTCATGCAACATACCTCTCAACCTCTTCCAGTGCATTTCTCCTGTACTGCAAAAGATGTGC
AACTGAAAAGAACATTTCCTGGATAGCATTTGATGTCATTTAGATTTAGCCAATC


line 3 in a.200/26695289

[6]   ..26695289                     [pid=+9207508,+113885979,+185087549,-270372233,+386904880,+392793776]
[1]   26495132,26535728..26695289    [pid=+227426911]
[1]   26535728,2670447..26695289     [pid=-415390484]
[1]   26670675..26695289             [pid=-102007768]
[1]   26695289..                     [pid=+345235013]
[17]  26695289..17580202             [pid=-66147254,-80335964,-109897383,+125065499,+153005786,-195371131,-213760442,-230717026,-263422260,+276157194,-285771365,+293842023,-302257825,-308746894,+382978287,+408651912,+438629761]
[1]   26695289..29512571             [pid=+39958119]
[1]   26859125..26695289             [pid=+159595663]
[1]   26859125,26495132..26695289    [pid=+222759187]


QueryLookupTable K=12 MM=12 MC=0.05 SEQS=../data/frag_reads_orig.fastb L=/wga/scr4/bigrefs/grch38/genome.lookup VISUAL=True AI=True TARGETS_TO_PROCESS=409 SMITH_WAT=True VISUAL_ABBR=False MF=5K SEQS_TO_PROCESS=227771958

227771958fw vs 409, 14 mismatches/0 indels (of 251), from 0-251 to 3965243-3965494 (of 4604811)

                                                                                
TAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACAAGGTCAGGAGATCGAGACCATCCTGGCTAACACAGTGAAA
TAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACAAGGTCAGGAGATCGAGACCATCCTGGCTAACACAGTGAAA


                                                                  *             
CCCCGTCTCTACTAAAAATACAAAAAATTAGCTGGGCGAGGTGGCGGGCGCCTGTAGTCCCAGCTAATCGGGAGGCTGAG
CCCCGTCTCTACTAAAAATACAAAAAATTAGCTGGGCGAGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAG


                    *  *                 * *    *  **    *       *      *    *  
GCAGGAGAATGGCGTGAACCTCGTGGGCGGAGCCTGCAGTGTGTCGAGTTCTGGCCAGTGCACTCTGGCCTGTGCGAGAG
GCAGGAGAATGGCGTGAACCCCGGGGGCGGAGCCTGCAGTGAGCCGAGATCGCGCCACTGCACTCCGGCCTGGGCGACAG


    **     
CGAGCATCCGT
CGAGACTCCGT

227771959rc vs 409, 0 mismatches/0 indels (of 251), from 0-251 to 3965518-3965769

243-387 read match
338-526 read match 540744467
429-803 gap
518-769 perfect read
529-780 perfect read 540744466

*/


     // Constants.

     const int read_len = 250; // *************************************************
     const int end_back = 400;

     // Hardcoded directories.

     String dir, odir, tdir;
     SetTenxDirs( N, dir, odir, tdir );

     // Load assembly.

     cout << Date( ) << ": loading assembly" << endl;
     HyperBasevectorX hb;
     BinaryReader::readFile( dir + "/a.hbx", &hb );
     vec<int> inv;
     BinaryReader::readFile( dir + "/a.inv", &inv );

     // Load alignments and barcodes.

     vec< pair<int,int> > places;
     vec<uint32_t> bcs;
     if (tenx)
     {    cout << Date( ) << ": loading alignments" << endl;
          BinaryReader::readFile( odir + "/10X.aligns", &places );
          cout << Date( ) << ": loading barcodes" << endl;
          BinaryReader::readFile( tdir + "/10X.bc", &bcs );    }

     // Identify the barcodes.

     vec<uint32_t> bx;
     cout << Date( ) << ": identifying barcodes" << endl;
     for ( int64_t id = 0; id < places.jsize( ); id++ )
     {    int e = places[id].first, epos = places[id].second;
          if ( e == e1 || e == inv[e1] || e == e2 || e == inv[e2] )
          {    bx.push_back( bcs[id] );    }    }
     UniqueSort(bx);
     cout << "using " << bx.size( ) << " barcodes" << endl;

     // Identify the reads.

     cout << Date( ) << ": identifying reads" << endl;
     vec<vec<int64_t>> hits( bx.size( ) );
     #pragma omp parallel for
     for ( int64_t id = 0; id < places.jsize( ); id++ )
     {    int p = BinPosition( bx, bcs[id] );
          if ( p < 0 ) continue;
          #pragma omp critical
          {    hits[p].push_back(id);    }    }
     for ( int i = 0; i < bx.isize( ); i++ )
          Sort( hits[i] );

     // Load the reads.

     cout << Date( ) << ": loading tenx reads" << endl;
     vec<int64_t> all;
     for ( int i = 0; i < bx.isize( ); i++ )
     for ( int j = 0; j < hits[i].isize( ); j++ )
     {    int64_t pid = hits[i][j]/2;
          all.push_back( 2*pid, 2*pid+1 );    }
     UniqueSort(all);
     cout << "total reads in barcode set = " << ToStringAddCommas( all.size( ) )
          << endl;
     vecbasevector tenx_bases;
     tenx_bases.Read( tdir + "/10X.fastb", all );

     // Find the ids of the reads that are incident upon e1 and e2.

     cout << Date( ) << ": finding read ids" << endl;
     vec<int> s = {e1, inv[e1], e2, inv[e2]};
     VecULongVec x;
     x.Read( dir + "/a.paths.inv", s );
     vec<int64_t> xall;
     for ( int i = 0; i < (int) x.size( ); i++ )
     for ( int j = 0; j < (int) x[i].size( ); j++ )
     {    int64_t pid = x[i][j]/2;
          xall.push_back( 2*pid, 2*pid+1 );    }

     // Load the paths for these reads.

     cout << Date( ) << ": loading " << xall.size( ) << " paths" << endl;
     ReadPathVec paths;
     paths.Read( dir + "/a.paths", xall );

     // Find the reads that might be in the gap, and their orientations.

     cout << Date( ) << ": looking for gappers" << endl;
     vec< pair<int64_t,Bool> > gappers;
     for ( int i = 0; i < (int) paths.size( ); i++ )
     {    const ReadPath& p = paths[i];
          if ( p.size( ) == 0 ) continue;
          int64_t id1 = xall[i];
          int64_t id2 = ( id1 % 2 == 0 ? id1 + 1 : id1 - 1 );
          int start = p.getOffset( );
          for ( int j = 0; j < (int) p.size( ); j++ )
          {    int e = p[j];
               if ( e == e1 )
               {    int stop = start + read_len;
                    if ( stop > hb.Bases(e) )
                    {    gappers.push( id1, True );
                         gappers.push( id2, False );    }
                    else if ( stop > hb.Bases(e) - end_back )
                    {    gappers.push( id2, False );    }    }
               if ( e == inv[e1] )
               {    if ( start < 0 )
                    {    gappers.push( id1, False );    }    }
               if ( e == e2 )
               {    if ( start < 0 )
                         gappers.push( id1, True );    }
               if ( e == inv[e2] )
               {    int stop = start + read_len;
                    if ( stop > hb.Bases(e) )
                    {    gappers.push( id1, False );
                         gappers.push( id2, True );    }
                    else if ( stop > hb.Bases(e) - end_back )
                    {    gappers.push( id2, True );    }    }
               start -= hb.Kmers(e);    }     }
     UniqueSort(gappers);
     for ( int i = 0; i < pids.isize( ); i++ )
     {    int64_t pid = pids[i];
          gappers.push( 2*pid, True );
          gappers.push( 2*pid, False );
          gappers.push( 2*pid+1, True );
          gappers.push( 2*pid+1, False );    }
     UniqueSort(gappers);

     // Load and orient the gappers.

     cout << Date( ) << ": loading " << gappers.size( ) << " reads" << endl;
     vecbasevector bases;
     vec<int64_t> gids( gappers.size( ) );
     for ( int i = 0; i < (int) gappers.size( ); i++ )
          gids[i] = gappers[i].first;
     bases.Read( dir + "/../data/frag_reads_orig.fastb", gids );
     for ( int i = 0; i < (int) gappers.size( ); i++ )
          if ( !gappers[i].second ) bases[i].ReverseComplement( );
     
     // Form the graph from the two edges and the reads.

     cout << Date( ) << ": forming graph" << endl;
     HyperBasevector hbg;
     /*
     {    vecbasevector stuff;
          stuff.push_back( hb.EdgeObject(e1) );
          stuff.push_back( hb.EdgeObject(e2) );
          stuff.Append(bases);
          stuff.Append(tenx_bases);
          HyperKmerPath h;
          vecKmerPath paths, paths_rc;
          unsigned const COVERAGE = 50u;
          LongReadsToPaths( stuff, K, COVERAGE, 0, False, 
               &hbg, &h, &paths, &paths_rc );    }
     */

     {    vecbasevector stuff;
          stuff.push_back( hb.EdgeObject(e1) );
          stuff.push_back( hb.EdgeObject(e2) );
          stuff.Append(bases);
          stuff.Append(tenx_bases);
          vecKmerPath paths, paths_rc, unipaths;
          // vec<big_tagged_rpint> pathsdb, unipathsdb;
          vec<tagged_rpint> pathsdb, unipathsdb;
          ReadsToPathsCoreY( stuff, K, paths );
          CreateDatabase( paths, paths_rc, pathsdb );
          Unipath( paths, paths_rc, pathsdb, unipaths, unipathsdb );
          digraph A;
          BuildUnipathAdjacencyGraph( paths, paths_rc, pathsdb, unipaths,
               unipathsdb, A );
          HyperKmerPath h;
          BuildUnipathAdjacencyHyperKmerPath( K, A, unipaths, h );
          KmerBaseBroker kbb( K, paths, paths_rc, pathsdb, stuff );
          hbg = HyperBasevector( h, kbb );    }

     // Identify edges sharing kmers with e1 and e2.

     cout << Date( ) << ": mapping back" << endl;
     vec<Bool> m1( hbg.E( ), False ); 
     vec<Bool> m2( hbg.E( ), False );
     {    vecbasevector woof;
          woof.push_back( hb.EdgeObject(e1) );
          woof.push_back( hb.EdgeObject(e2) );
          for ( int e = 0; e < hbg.E( ); e++ )
               woof.push_back( hbg.EdgeObject(e) );
          vec< triple<kmer<K>,int,int> > kmers_plus;
          MakeKmerLookup0( woof, kmers_plus );
          for ( int64_t i = 0; i < kmers_plus.jsize( ); i++ )
          {    int64_t j;
               for ( j = i + 1; j < kmers_plus.jsize( ); j++ )
                    if ( kmers_plus[j].first != kmers_plus[i].first ) break;
               for ( int64_t k1 = i; k1 < j; k1++ )
               for ( int64_t k2 = i; k2 < j; k2++ )
               {    int id1 = kmers_plus[k1].second;
                    if ( id1 > 1 ) continue;
                    int id2 = kmers_plus[k2].second - 2;
                    if ( id2 < 0 ) continue;
                    if ( id1 == 0 ) m1[id2] = True;
                    else m2[id2] = True;    }
               i = j - 1;    }    }

     // Dump fasta file.

     Ofstream( fout, "/wga/dev/jaffe/BroadCRD/yyy.fasta" );
     for ( int e = 0; e < hbg.E( ); e++ )
          hbg.EdgeObject(e).Print( fout, e );

     // Display the graph.

     cout << Date( ) << ": displaying graph" << endl;
     vec<String> edge_color( hbg.E( ) );
     for ( int e = 0; e < hbg.E( ); e++ )
     {    if ( m1[e] && m2[e] ) edge_color[e] = "purple";
          else if ( m1[e] ) edge_color[e] = "red";
          else if ( m2[e] ) edge_color[e] = "blue";    }
     vec<String> edge_names( hbg.E( ) );
     for ( int e = 0; e < hbg.E( ); e++ )
          edge_names[e] = ToString(e);
     vec<double> lengths( hbg.E( ) );
     for ( int i = 0; i < hbg.E( ); i++ )
          lengths[i] = hbg.Kmers(i);
     Ofstream( dout, "/wga/dev/jaffe/BroadCRD/yyy.dot" );
     const Bool LABEL_CONTIGS = False;
     const Bool DOT_LABEL_VERTICES = False;
     hbg.PrettyDOT( dout, lengths, HyperBasevector::edge_label_info(
          HyperBasevector::edge_label_info::DIRECT, &edge_names ),
          LABEL_CONTIGS, DOT_LABEL_VERTICES, NULL, NULL, NULL, NULL,
          NULL, &edge_color );    }
