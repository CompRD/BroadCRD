///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// ScaffoldUnipaths.  This is just a hack for testing scaffolding.
//


// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
// MakeDepend: dependency QueryLookupTable

#include <omp.h>

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "Superb.h"
#include "efasta/EfastaTools.h"
#include "feudal/BinaryStream.h"
#include "kmers/KmerRecord.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "paths/AlignReadsOnUnibases.h"
#include "paths/AssemblyCleanupTools.h"
#include "paths/BigMapTools.h"
#include "paths/GetNexts.h"
#include "paths/LongReadTools.h"
#include "paths/PdfEntry.h"
#include "paths/ProcessGap.h"
#include "paths/UnibaseUtils.h"
#include "paths/UnipathScaffold.h"
#include "Vec.h"
#include "graph/FindCells.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_Bool_OrDefault(USE_ALTERNATIVE, True);
     CommandArgument_String_OrDefault(SUBDIR, "test");
     //CommandArgument_Int_OrDefault_Doc(K1, 96, "little K");
     CommandArgument_Int_OrDefault_Doc(K2, 96, "big K");
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
       "Number of threads to use (use all available processors if set to 0)");
     //CommandArgument_String_Doc(HEAD1, "head for K1 unibases");
     CommandArgument_String_Doc(HEAD2, "head for K2 unibases");
     CommandArgument_String_OrDefault(JUMPS_IN, "jump_reads_ec");
     CommandArgument_Bool_OrDefault(WRITE, True);
     CommandArgument_Int_OrDefault(FORCE_FIRST, -1);
     CommandArgument_Bool_OrDefault(RAISE_VERBOSE, False);
     CommandArgument_Bool_OrDefault(DELETE_VERBOSE, False);
     CommandArgument_Bool_OrDefault(FORCE_JALIGNS, False);
     // algorithmic options

     CommandArgument_Bool_OrDefault(FILTER_BY_JUMPS, True);
     CommandArgument_Int_OrDefault_Doc( MIN_KMERS, 100,
	"Min number of k-mers in the unipath used in scaffolding." );

     // logging options

     CommandArgument_Int_OrDefault(VERBOSITY, 0);
     CommandArgument_Bool_OrDefault(VALIDATE, False);
     CommandArgument_Bool_OrDefault(SHOW_TRANSLATION, False);
     //CommandArgument_Bool_OrDefault(SHOW_LINKS, False);
     CommandArgument_Bool_OrDefault(SHOW_LINKS2, False);
     //CommandArgument_Bool_OrDefault(SHOW_INITIAL_LINKS, False);
     CommandArgument_Bool_OrDefault(SHOW_INITIAL_LINKS2, False);
     //CommandArgument_Bool_OrDefault(PRINT_GRAPH1, False);
     CommandArgument_Bool_OrDefault(PRINT_GRAPH2, False);
     CommandArgument_Bool_OrDefault(SCAFFOLDING_VERBOSE, False);

     EndCommandArguments;

     // Thread control, etc.

     double clock = WallClockTime( );
     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads( NUM_THREADS );

     // Heuristics.

     //const int min_kmers1 = 800;
     const int min_kmers2 = MIN_KMERS;

     
     // Define directories.
     //String KS1 = ToString(K1), 
     String KS2 = ToString(K2);
     

     String data_dir = PRE + "/" + DATA;
     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;

     String uscaffolds_out_head = sub_dir + "/" + HEAD2 + ".unibases.k" + KS2;

     // Load unibases and generate ancillary data.

     cout << Date( ) << ": loading unibases" << endl;

     //vecbasevector unibases1( run_dir + "/" + HEAD1 + ".unibases.k" + KS1 );
     vecbasevector unibases2( run_dir + "/" + HEAD2 + ".unibases.k" + KS2 );
     //int nuni1 = unibases1.size( ), 
     int nuni2 = unibases2.size( );
     //vec<int> to_rc1, 
     vec<int> to_rc2;
     //UnibaseInvolution( unibases1, to_rc1 );
     UnibaseInvolution( unibases2, to_rc2 );

     // Define alignments between unibases2.  There are two sets.  The first set
     // set (nexts2) consists of K2-1 base alignments.  The second set (nexts2x) 
     // consists of alignments of size between K1-1 and K2-1.

//      vec< vec<int> > nexts2;
//      GetNexts( K2, unibases2, nexts2 );
//      const int K0 = 80;
//      ForceAssertLe( K0, K1 );
//      vec<int64_t> starts;
//      starts.push_back(0);
//      for ( size_t i = 0; i < unibases2.size( ); i++ )
//      {    const basevector& u = unibases2[i];
//           starts.push_back( starts.back( ) + u.isize( ) - K0 + 1 );    }
//      vec< triple<kmer<K0>,int,int> > kmers_plus( starts.back( ) );
//      #pragma omp parallel for schedule( dynamic, 1 )
//      for ( size_t i = 0; i < unibases2.size( ); i++ )
//      {    const basevector& u = unibases2[i];
//           for ( int j = 0; j <= u.isize( ) - K0; j++ )
//           {    int64_t r = starts[i] + j;
//                kmers_plus[r].first.SetToSubOf( u, j );
//                kmers_plus[r].second = i, kmers_plus[r].third = j;    }    }
//      ParallelSort(kmers_plus);
//      vec< kmer<K0> > kmers( kmers_plus.size( ) );
//      for ( size_t i = 0; i < kmers.size( ); i++ )
//           kmers[i] = kmers_plus[i].first;
//      cout << Date( ) << ": mapping from 2 to 2" << endl;
//      vec< vec< pair<int,int> > > nexts2x( unibases2.size( ) );
//      #pragma omp parallel for
//      for ( int u1 = 0; u1 < (int) unibases2.size( ); u1++ )
//      {    kmer<K0> x;
//           x.SetToSubOf( unibases2[u1], unibases2[u1].isize( ) - K0 );
//           int64_t low = LowerBound( kmers, x ), high = UpperBound( kmers, x );
//           for ( int64_t z = low; z < high; z++ )
//           {    int u2 = kmers_plus[z].second, pos2 = kmers_plus[z].third;
//                if ( u2 == u1 && pos2 == unibases2[u1].isize( ) - K0 ) continue;
//                int overlap = pos2 + K0;
//                if ( overlap < K1 || overlap > unibases2[u1].isize( ) ) continue;
//                Bool match = True;
//                for ( int j = 0; j < overlap - K0; j++ )
//                {    if ( unibases2[u2][j] 
//                          != unibases2[u1][ unibases2[u1].isize( ) - overlap + j ] )
//                     {    match = False;
//                          break;    }    }
//                if (match) nexts2x[u1].push( u2, overlap );    }    }
//      vec< vec< pair<int,int> > > nexts2y(nexts2x);
//      for ( int i = 0; i < nexts2y.isize( ); i++ )
//      {    vec<Bool> to_remove( nexts2y[i].size( ), False );
//           for ( int j = 0; j < nexts2y[i].isize( ); j++ )
//           {    if ( nexts2y[i][j].second < K2 - 1 )
//                     to_remove[j] = True;    }
//           EraseIf( nexts2y[i], to_remove );    }

//      // Locate unibases on genome.

//      const int LG = 12;
//      vec< vec< pair<int,int> > > Glocs;
//      vecbasevector genome;
//      vec< vec<placementy> > Ulocs1, Ulocs2;
//      if (VALIDATE)
//      {    cout << Date( ) << ": locating on genome" << endl;
//           Glocs.resize( IPow( 4, LG ) );
//           genome.ReadAll( data_dir + "/genome.fastb" );
//           for ( size_t i = 0; i < genome.size( ); i++ )
//           {    for ( int j = 0; j <= genome[i].isize( ) - LG; j++ )
//                {    int n = KmerId( genome[i], LG, j );
//                     Glocs[n].push( i, j );    }    }
//           Ulocs1.resize(nuni1);
//           for ( int u = 0; u < nuni1; u++ )
//                Ulocs1[u] = FindGenomicPlacementsY( unibases1[u], LG, genome, Glocs );
//           Ulocs2.resize(nuni2);

//           // Find perfect placements of unibases2.  Then use QueryLookupTable
//           // to find imperfect placements.

//           int count = 0;
//           {    Ofstream( out, run_dir + "/tmp/BigMap.ids" );
//                for ( int u = 0; u < nuni2; u++ )
//                {    Ulocs2[u] = FindGenomicPlacementsY( 
//                          unibases2[u], LG, genome, Glocs );    
//                     if ( Ulocs2[u].empty( ) ) 
//                     {    count++;
//                          out << u << "\n";    }    }    }
//           if ( count > 0 )
//           {    fast_pipe_ifstream qlt( "QueryLookupTable K=12 MM=12 MC=0.15 SEQS="
//                     + run_dir + "/" + HEAD2 + ".unibases.k" + KS2 + " L=" + data_dir 
//                     + "/genome.lookup SEQS_IS_FASTB=True PARSEABLE=True "
//                     + "SEQS_TO_PROCESS=@" + run_dir + "/tmp/BigMap.ids" );
//                String line;
//                look_align la;
//                while(1)
//                {    getline( qlt, line );
//                     if ( qlt.fail( ) ) break;
//                     if ( line.Contains( "QUERY", 0 ) )
//                     {    la.ReadParseable(line);
//                          Ulocs2[la.query_id].push( la.target_id, la.pos2( ), 
//                               la.Pos2( ), la.Fw1( ) );    }    }    }
//           Remove( run_dir + "/tmp/BigMap.ids" );    }

     // Get copy number.

     cout << Date( ) << ": getting copy number" << endl;
     VecPdfEntryVec CN2( ( run_dir + "/" + HEAD2 + ".unipaths.predicted_count.k"
          + KS2 ).c_str( ) );
     vec<int> predicted_CN2( nuni2, -1 );
     for ( int i = 0; i < nuni2; i++ )
          GetMostLikelyValue( predicted_CN2[i], CN2[i] );
     vec<double> CN_raw;
     BinaryReader::readFile( run_dir + "/" + HEAD2 + ".unipaths.cn_raw.k" + KS2,
                             &CN_raw );

     // Load jumps and align to unibases2.

     vec< triple<int64_t,int,int> > jaligns;
     vec<basevector> jbases_sorted;
     vec<int64_t> jbases_sorted_id;
     PairsManager jpairs;
     if (FILTER_BY_JUMPS){    
       cout << Date( ) << ": loading jump data" << endl;
       vecbasevector jbases( run_dir + "/" + JUMPS_IN + ".fastb" );
       jpairs.Read( run_dir + "/" + JUMPS_IN + ".pairs" );
       jpairs.makeCache( );
       String jaligns_file = run_dir + "/" + JUMPS_IN + ".on." + HEAD2 + ".unibases.k" + KS2 + ".aligns";
       String jbases_sorted_file = run_dir + "/" + JUMPS_IN + ".sorted.fastb";
       String jbases_sorted_id_file = run_dir + "/" + JUMPS_IN + ".sorted_id";
       if ( FORCE_JALIGNS || ! IsRegularFile( jaligns_file ) ){
	 cout << Date() << ": aligning jumps on unibases" << endl;
	 AlignReadsOnUnibases( jbases, unibases2, jaligns, jbases_sorted, 
			       jbases_sorted_id, cout );

	 vecbvec vjbases_sorted;
	 for ( size_t jbi = 0; jbi < jbases_sorted.size(); jbi++ )
	   vjbases_sorted.push_back_reserve( jbases_sorted[jbi] );
	 vjbases_sorted.WriteAll( jbases_sorted_file );
	 BinaryWriter::writeFile( jaligns_file, jaligns );
	 BinaryWriter::writeFile( jbases_sorted_id_file, jbases_sorted_id );
       }else{
	 cout << Date() << ": reading jumps on unibases data" << endl;
	 vecbvec vjbases_sorted( jbases_sorted_file );
	 jbases_sorted.reserve( vjbases_sorted.size() );
	 for ( size_t jbi = 0; jbi < vjbases_sorted.size(); jbi++ )
	   jbases_sorted.push_back( vjbases_sorted[jbi] );
	 BinaryReader::readFile( jaligns_file, &jaligns );
	 BinaryReader::readFile( jbases_sorted_id_file, &jbases_sorted_id );
       }
	  
     }

     // Map large unibases1 to unibases2.

//      vec< vec< pair<int,int> > > to_big;
//      if ( K2 == 640 ) ToBig<640>( unibases1, unibases2, to_big );
//      else
//      {    cout << "Not implemented for K2 = " << K2 << endl;
//           exit(1);    }
//      if (SHOW_TRANSLATION)
//      {    cout << "\nMap from unibases1 to unibases2\n\n";
//           for ( int u1 = 0; u1 < nuni1; u1++ )
//           {    int nkmers1 = unibases1[u1].isize( ) - K1 + 1;
//                if ( nkmers1 < min_kmers1 ) continue;
//                if ( to_big[u1].size( ) > 1 ) cout << u1 << " MULTIMAPS!\n";
//                for ( int j = 0; j < to_big[u1].isize( ); j++ )
//                {    int u2 = to_big[u1][j].first, pos2 = to_big[u1][j].second;
//                     int cn = predicted_CN1[u1];
//                     prob_t maxp = 0;
//                     for ( uint z = 0; z < CN1[u1].size( ); z++ )
//                          if ( CN1[u1][z].first == cn ) maxp = CN1[u1][z].second;
//                     String prob = ToString( 100.0 * maxp, 1 );
//                     cout << u1 << "[l=" << nkmers1 << ",rc=" << to_rc1[u1] << ",cn="
//                          << cn << "(" << prob << "%)" << ",cn_raw=" 
//                          << setiosflags(ios::fixed) << setprecision(2) << CN_raw[u1]
//                          << resetiosflags(ios::fixed) << "] --> " << u2 << "[l=" 
//                          << unibases2[u2].isize( ) - K2 + 1 << ",rc=" << to_rc2[u2] 
//                          << "] at " << pos2 << "\n";    }    }
//           cout << "\n";    }

     // Load predicted gaps and use them to make a digraph.

     digraphE<linklet> G;
     {    vec< vec<int> > from(nuni2), to(nuni2);
          vec< vec<int> > from_edge_obj(nuni2), to_edge_obj(nuni2);
          vec<linklet> edges;
          String line;
          fast_ifstream in( 
               run_dir + "/" + HEAD2 + ".unibases.k" + KS2 + ".predicted_gaps.txt" );
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               if ( line.Contains( "#", 0 ) ) continue;
               int u1, u2, sep, dev, nlinks;
               istringstream iline( line.c_str( ) );
               iline >> u1 >> u2 >> sep >> dev >> nlinks;
               int nkmers1 = unibases2[u1].isize( ) - K2 + 1;
               int nkmers2 = unibases2[u2].isize( ) - K2 + 1;
               if ( nkmers1 < min_kmers2 || nkmers2 < min_kmers2 ) continue;
               if ( u2 == to_rc2[u1] ) continue; 
               from[u1].push_back(u2), to[u2].push_back(u1);
               from_edge_obj[u1].push_back( edges.size( ) );
               to_edge_obj[u2].push_back( edges.size( ) );
               edges.push( sep, dev, nlinks, 0 );    }
          for ( int u = 0; u < nuni2; u++ )
          {    SortSync( from[u], from_edge_obj[u] );
               SortSync( to[u], to_edge_obj[u] );    }
          G.Initialize( from, to, edges, to_edge_obj, from_edge_obj );    
	  PRINT( G.NComponents() );
     }

     // Identify unipaths that appear not to have copy number one.

     vec<Bool> not_one;
     const int max_link_ratio = 5;
     const int dev_mult = 3;
     NotOne( K2, unibases2, to_rc2, G, min_kmers2, max_link_ratio, dev_mult, 
          RAISE_VERBOSE, not_one );

     // Additional test for not copy number one.

     for ( int u1 = 0; u1 < nuni2; u1++ )
     {    for ( int j = 0; j < G.From(u1).isize( ); j++ )
          {    int u2 = G.From(u1)[j];
               const linklet& l = G.EdgeObjectByIndexFrom( u1, j );
               int nkmers1 = unibases2[u1].isize( ) - K2 + 1;
               int nkmers2 = unibases2[u2].isize( ) - K2 + 1;
               double cn1 = CN_raw[u1], cn2= CN_raw[u2];
               double min_cn_mult = 1.5;
               double min_size_mult = 4.0;
               const int min_links = 10;
               if ( l.nlinks >= min_links && cn1 >= min_cn_mult * cn2
                    && double(nkmers1) <= min_size_mult * double(nkmers2) )
               {    if (RAISE_VERBOSE)
                    {    cout << "\n";
                         PRINT7( u1, u2, cn1, cn2, nkmers1, nkmers2, l.nlinks );
                         cout << u1 << " not of copy number one\n";    }
                    not_one[u1] = True;    }    }    }

     // for ( int u = 0; u < nuni2; u++ )
//      {    int nkmers1 = unibases1[u].isize( ) - K2 + 1;
//           if ( nkmers1 < min_kmers1 ) continue;
//           if ( !to_big[u].solo( ) ) 
//           {    not_one[u] = True;
//                if ( RAISE_VERBOSE ) 
//                     cout << u << " not of copy number one?\n";    }    }

     // Print initial links.

//      if (SHOW_INITIAL_LINKS)
//      {    cout << "\nInitial unibase1 links\n\n";
//           for ( int u1 = 0; u1 < nuni1; u1++ )
//           for ( int j = 0; j < G.From(u1).isize( ); j++ )
//           {    int u2 = G.From(u1)[j];
//                int nkmers1 = unibases1[u1].isize( ) - K1 + 1;
//                int nkmers2 = unibases1[u2].isize( ) - K1 + 1;
//                const linklet& l = G.EdgeObjectByIndexFrom( u1, j );
//                cout << "u1 = " << u1 << "[l=" << nkmers1
//                     /* << ",cn=" << predicted_CN1[u1] */
//                     << "], u2 = " << u2
//                     << "[l=" << nkmers2
//                     /* << ",cn=" << predicted_CN1[u2] */
//                     << "], ";
//                int sep = l.sep, dev = l.dev, nlinks = l.nlinks;
//                cout << "sep = " << sep << " +/- " << dev << ", nlinks = "
//                     << nlinks << "\n";    }
//           cout << "\n";    }

     // Derive initial links2.

//      {
//      vec< triple<int,int,linklet> > links2;
//      const int max_dist_mult = 15;
//      const int max_self_dist_mult = 3;
//      for ( int u1 = 0; u1 < nuni1; u1++ )
//      {    if ( !to_big[u1].solo( ) ) continue;
//           int v1 = to_big[u1][0].first, pos1 = to_big[u1][0].second;
//           for ( int j = 0; j < G.From(u1).isize( ); j++ )
//           {    int u2 = G.From(u1)[j];
//                linklet l = G.EdgeObjectByIndexFrom( u1, j );
//                if ( !to_big[u2].solo( ) ) continue;
//                int v2 = to_big[u2][0].first, pos2 = to_big[u2][0].second;
//                l.sep -= ( unibases2[v1].isize( ) - pos1 
//                     - unibases1[u1].isize( )) + pos2;
//                if ( l.sep < - K2 - 1 - max_dist_mult * l.dev ) continue;
//                if ( v1 == v2 && l.sep < - K2 - 1 - max_self_dist_mult * l.dev ) 
//                     continue;
//                links2.push( v1, v2, l );    }    }
//      Sort(links2);
//      digraphE<linklet> G2;
//      vec< vec<int> > from(nuni2), to(nuni2);
//      vec< vec<int> > from_edge_obj(nuni2), to_edge_obj(nuni2);
//      vec<linklet> edges;
//      for ( int i = 0; i < links2.isize( ); i++ )
//      {    
//           // In this version, for better or for worse, if there are multiple 
//           // link groups, we choose the biggest one, rather than combine them.

//           int j;
//           for ( j = i + 1; j < links2.isize( ); j++ )
//           {    if ( links2[j].first != links2[i].first ) break;
//                if ( links2[j].second != links2[i].second ) break;    }
//           int u1 = links2[i].first, u2 = links2[i].second;
//           const linklet& l = links2[i].third;
//           int sep = l.sep, dev = l.dev, nlinks = l.nlinks;
//           from[u1].push_back(u2), to[u2].push_back(u1);
//           from_edge_obj[u1].push_back( edges.size( ) );
//           to_edge_obj[u2].push_back( edges.size( ) );
//           edges.push( sep, dev, nlinks, 0 );
//           i = j - 1;    }
//      for ( int u = 0; u < nuni2; u++ )
//      {    SortSync( from[u], from_edge_obj[u] );
//           SortSync( to[u], to_edge_obj[u] );    }
//      G2.Initialize( from, to, edges, to_edge_obj, from_edge_obj );
//      if (SHOW_INITIAL_LINKS2)
//      {    cout << "\nInitial unibase2 links\n\n";
//           for ( int u1 = 0; u1 < nuni2; u1++ )
//           for ( int j = 0; j < G2.From(u1).isize( ); j++ )
//           {    int u2 = G2.From(u1)[j];
//                int nkmers1 = unibases2[u1].isize( ) - K2 + 1;
//                int nkmers2 = unibases2[u2].isize( ) - K2 + 1;
//                const linklet& l = G2.EdgeObjectByIndexFrom( u1, j );
//                cout << "u1 = " << u1 << "[l=" << nkmers1 << "], u2 = " << u2
//                     << "[l=" << nkmers2 << "], ";
//                int sep = l.sep, dev = l.dev, nlinks = l.nlinks;
//                cout << "sep = " << sep << " +/- " << dev << ", nlinks = " 
//                     << nlinks << "\n";    }
//           cout << "\n";    }
//      }

     // Remove weak links.

     RemoveWeakLinks( unibases2, to_rc2, not_one, G, max_link_ratio, 
          dev_mult, DELETE_VERBOSE );
     PRINT( G.NComponents() );
     // Print links.
     
//      if (SHOW_LINKS)
//      {    cout << "\nUnibase1 links\n\n";
//           for ( int u1 = 0; u1 < nuni1; u1++ )
//           for ( int j = 0; j < G.From(u1).isize( ); j++ )
//           {    int u2 = G.From(u1)[j];
//                int nkmers1 = unibases1[u1].isize( ) - K1 + 1;
//                int nkmers2 = unibases1[u2].isize( ) - K1 + 1;
//                const linklet& l = G.EdgeObjectByIndexFrom( u1, j );
//                cout << "u1 = " << u1 << "[l=" << nkmers1 << ",cn=" 
//                     << predicted_CN1[u1] << "], u2 = " << u2
//                     << "[l=" << nkmers2 << ",cn=" << predicted_CN1[u2] << "], ";
//                int sep = l.sep, dev = l.dev, nlinks = l.nlinks;
//                cout << "sep = " << sep << " +/- " << dev << ", nlinks = " 
//                     << nlinks << "\n";    }
//           cout << "\n";    }

     // Infer links between unibases2.

//      vec< triple<int,int,linklet> > links2;
//      const int max_dist_mult = 15;
//      const int max_self_dist_mult = 3;
//      for ( int u1 = 0; u1 < nuni1; u1++ )
//      {    if ( !to_big[u1].solo( ) ) continue;
//           int v1 = to_big[u1][0].first, pos1 = to_big[u1][0].second;
//           for ( int j = 0; j < G.From(u1).isize( ); j++ )
//           {    int u2 = G.From(u1)[j];
//                linklet l = G.EdgeObjectByIndexFrom( u1, j );
//                if ( !to_big[u2].solo( ) ) continue;
//                int v2 = to_big[u2][0].first, pos2 = to_big[u2][0].second;
//                l.sep -= ( unibases2[v1].isize( ) - pos1 
//                     - unibases1[u1].isize( )) + pos2;
//                if ( l.sep < - K2 - 1 - max_dist_mult * l.dev ) continue;
//                if ( v1 == v2 && l.sep < - K2 - 1 - max_self_dist_mult * l.dev ) 
//                     continue;
//                links2.push( v1, v2, l );    }    }
//      Sort(links2);
//      digraphE<linklet> G2;
//      vec< vec<int> > from(nuni2), to(nuni2);
//      vec< vec<int> > from_edge_obj(nuni2), to_edge_obj(nuni2);
//      vec<linklet> edges;
//      for ( int i = 0; i < links2.isize( ); i++ )
//      {    
//           // In this version, for better or for worse, if there are multiple 
//           // link groups, we choose the biggest one, rather than combine them.

//           int j;
//           for ( j = i + 1; j < links2.isize( ); j++ )
//           {    if ( links2[j].first != links2[i].first ) break;
//                if ( links2[j].second != links2[i].second ) break;    }
//           int u1 = links2[i].first, u2 = links2[i].second;
//           const linklet& l = links2[i].third;
//           int sep = l.sep, dev = l.dev, nlinks = l.nlinks;
//           from[u1].push_back(u2), to[u2].push_back(u1);
//           from_edge_obj[u1].push_back( edges.size( ) );
//           to_edge_obj[u2].push_back( edges.size( ) );
//           edges.push( sep, dev, nlinks, 0 );
//           i = j - 1;    }
//      for ( int u = 0; u < nuni2; u++ )
//      {    SortSync( from[u], from_edge_obj[u] );
//           SortSync( to[u], to_edge_obj[u] );    }
//      G2.Initialize( from, to, edges, to_edge_obj, from_edge_obj );

     // Test weak links.

//      const int max_to_call_weak = 4;
//      vec<int> to_removex;
//      for ( int u1 = 0; u1 < nuni2; u1++ )
//      {    int max_nlinks = 0;
//           for ( int j = 0; j < G2.From(u1).isize( ); j++ )
//           {    const linklet& l = G2.EdgeObjectByIndexFrom( u1, j );
//                max_nlinks = Max( max_nlinks, l.nlinks );    }
//           if ( max_nlinks > max_to_call_weak ) continue;
//           vec<Bool> keep( G2.From(u1).size( ), False );
//           for ( int j = 0; j < G2.From(u1).isize( ); j++ )
//           {    int u2 = G2.From(u1)[j];
//                const linklet& l = G2.EdgeObjectByIndexFrom( u1, j );
//                vec< vec< pair<int,int> > > walks1;
//                int bad;
//                GetWalks( u1, u2, l.sep, l.dev, unibases2, nexts2x, walks1, bad );
//                if ( walks1.nonempty( ) ) keep[j] = True;    }
//           if ( Sum(keep) == 0 ) continue;
//           for ( int j = 0; j < G2.From(u1).isize( ); j++ )
//           {    if ( keep[j] ) continue;
//                to_removex.push_back( 
//                     G2.EdgeObjectIndexByIndexFrom( u1, j ) );    }    }
//      G2.DeleteEdges(to_removex);

     // Print links.

//      if (SHOW_LINKS2)
//      {    cout << "\nUnibase2 links\n\n";
//           for ( int u1 = 0; u1 < nuni2; u1++ )
//           for ( int j = 0; j < G2.From(u1).isize( ); j++ )
//           {    int u2 = G2.From(u1)[j];
//                int nkmers1 = unibases2[u1].isize( ) - K2 + 1;
//                int nkmers2 = unibases2[u2].isize( ) - K2 + 1;
//                const linklet& l = G2.EdgeObjectByIndexFrom( u1, j );
//                cout << "u1 = " << u1 << "[l=" << nkmers1 << "], u2 = " << u2
//                     << "[l=" << nkmers2 << "], ";
//                int sep = l.sep, dev = l.dev, nlinks = l.nlinks;
//                cout << "sep = " << sep << " +/- " << dev << ", nlinks = " 
//                     << nlinks << "\n";    }
//           cout << "\n";    }

     // Build unipath scaffolds, allowing for the circular case.  We only use
     // unipaths that are predicted to have copy number one.

     cout << Date( ) << ": building unipath scaffolds" << endl;
     vec<superb> uscaffolds2;
     vec<Bool> circled;
     vec<Bool> used2( nuni2, False );
     //for ( int u1 = 0; u1 < nuni2; u1++ )
       //      {    if ( to_big[u1].solo( ) )
       //           {    int u2 = to_big[u1][0].first;
       //                used2[u2] = False;    }    }
       
     if ( USE_ALTERNATIVE ){
       UnipathScaffoldAlt( G, unibases2, K2, to_rc2, used2, min_kmers2,
			   FORCE_FIRST, SCAFFOLDING_VERBOSE, uscaffolds2, circled );
     }
     else
       UnipathScaffold( G, unibases2, K2, to_rc2, used2, min_kmers2,
			FORCE_FIRST, SCAFFOLDING_VERBOSE, uscaffolds2, circled );
     
     VecEFasta eunibases2;
     for ( size_t i = 0; i < unibases2.size() ; i++ )
       eunibases2.push_back( efasta( fastavector(unibases2[i] ) ) );
     

     Assembly A( uscaffolds2, eunibases2 );
     PRINT( A.scaffoldsNtigs() );
     A.remove_unused_contigs();   
     A.check_integrity();
     PRINT( A.scaffoldsNtigs() );
     //A.remove_small_contigs( 1000, 1 );
     //A.remove_small_scaffolds( 1 );
     PRINT( A.scaffoldsNtigs() );

     A.WriteAll( uscaffolds_out_head );
     
     

//      {
//        cout << Date( ) << ": writing unipath scaffolds to " << uscaffolds_out_head << ".superb"  << endl;
//        WriteSuperbs( uscaffolds_out_head + ".superb" , uscaffolds2 );
       
//        size_t count = 0;
//        Ofstream( ofa, uscaffolds_out_head + ".contigs.fasta" );
//        for ( size_t i = 0; i < unibases2.size() ; i++)
// 	 unibases2[i].PrintCol(ofa,ToString(count++),80);
//      }
     

//      // Set up global output data structures.

//      VecEFasta tigse;
//      vec<superb> scaffolds;

//      // Simultaneously patch gaps in and evaluate scaffolds.

//      String single_line = "-------------------------------------------------------"
//           "-----------------------------\n";
//      String double_line = "======================================================="
//           "=============================\n";
//      cout << Date( ) << ": patching and evaluating unipath scaffolds" << endl;
//      vec< vec<int> > extras( uscaffolds2.size( ) );
//      vec<String> reports( uscaffolds2.size( ) );
//      for ( int si = 0; si < uscaffolds2.isize( ); si++ )
//      {    const superb& s = uscaffolds2[si];
//           ostringstream out;

//           // Set up data structures to store the new scaffold.

//           VecEFasta SEQ;
//           vec<int> SEP, DEV;

//           // Go through the contigs.

//           placementy lastp;
//           int lastover = K2 - 1;
//           Bool have_lastp = False;
//           int gapcount = 0;
//           for ( int j = 0; j < s.Ntigs( ); j++ )
//           {    int u = s.Tig(j);

//                // Find walks from last.

//                Bool at_circ = ( circled[si] && j == s.Ntigs( ) - 1 );
//                if ( j > 0 )
//                {    int sep = s.Gap(j-1), dev = s.Dev(j-1), u1 = s.Tig(j-1), u2 = u;
//                     efasta epatch;
//                     Bool found_patch = False;
//                     vec<int> ex;
//                     ProcessGap( s, u1, u2, sep, dev, unibases2, nexts2x, nexts2y, 
//                          K2, jbases_sorted, jbases_sorted_id, jpairs, jaligns, out, 
//                          VERBOSITY, VALIDATE, LG, genome, Glocs, have_lastp, lastp, 
//                          gapcount, G2, found_patch, epatch, lastover, ex );
//                     extras[si].append(ex);
//                     if (found_patch)
//                     {    if ( epatch.size( ) > 0 )
//                          {    SEQ.back( ).resize( SEQ.back( ).isize( ) 
//                                    - lastover );    }
//                          SEQ.back( ) += epatch;
//                          if (at_circ)
//                          {    SEQ[0] = SEQ[0].substr( 
//                                    lastover, SEQ[0].isize( ) - lastover );    }    }
//                     else 
//                     {    if ( VERBOSITY >= 1 ) out << "GAP!\n";
//                          if ( !at_circ )
//                          {    SEP.push_back(sep), DEV.push_back(dev);    } 
//                          out << "\n";    }    }

//                // Save current.

//                if ( !at_circ )
//                {    if ( SEQ.size( ) > SEP.size( ) )
//                     {    basevector b( unibases2[u], lastover, 
//                               unibases2[u].isize( ) - lastover );
//                          SEQ.back( ) += b.ToString( );    }
//                     else SEQ.push_back( unibases2[u].ToString( ) );    }

//                // Print placement.

//                vec<placementy> goods;
//                if ( VALIDATE && VERBOSITY >= 1 )
//                {    goods = Ulocs2[u];
//                     Bool have_good = False;
//                     if (have_lastp)
//                     {    for (int l = 0; l < goods.isize( ); l++)
//                          {    if ( Follow( lastp, goods[l], lastover + 1 ) )
//                               {    have_good = True;
//                                    goods[0] = goods[l];
//                                    goods.resize(1);
//                                    break;    }    }    }
//                     if ( have_lastp && !have_good ) out << single_line;    }
//                if ( VERBOSITY >= 1 )
//                     out << u << "[" << unibases2[u].isize( ) - K2 + 1 << "]";
//                if ( VALIDATE && VERBOSITY >= 1 )
//                {    out << " (";
//                     for ( int l = 0; l < goods.isize( ); l++ )
//                          Print( out, goods[l], unibases2[u].size( ) );
//                     if ( goods.solo( ) )
//                     {    lastp = goods[0];
//                          have_lastp = True;    }
//                     else have_lastp = False;
//                     out << " )";    }
//                if ( VERBOSITY >= 1 ) out << "\n";    }

//           // Save efasta scaffold.

//           superb snew;
//           snew.SetNtigs( SEQ.size( ) );
//           for ( int j = 0; j < snew.Ntigs( ); j++ )
//           {    snew.SetTig( j, tigse.size( ) + j );
//                snew.SetLen( j, SEQ[j].Length1( ) );
//                if ( j < snew.Ngaps( ) )
//                {    snew.SetGap( j, SEP[j] );
//                     snew.SetDev( j, DEV[j] );    }    }
//           scaffolds.push_back(snew);
//           tigse.append(SEQ); 
//           int n = extras[si].size( );
//           for ( int j = 0; j < n; j++ )
//                extras[si].push_back( to_rc2[ extras[si][j] ] );
//           UniqueSort( extras[si] );
//           out << "\nlength excluding gaps = " 
//                << ToStringAddCommas( snew.ReducedLength( ) ) << endl;
//           reports[si] = out.str( );    }

//      // Add unipaths2 that have not been used.

//      used2.resize(0);
//      used2.resize( nuni2, False );
//      for ( int i = 0; i < scaffolds.isize( ); i++ )
//      {    for ( int j = 0; j < uscaffolds2[i].Ntigs( ); j++ )
//                used2[ uscaffolds2[i].Tig(j) ] = True;
//           for ( int j = 0; j < extras[i].isize( ); j++ )
//                used2[ extras[i][j] ] = True;    }
//      for ( int u = 0; u < nuni2; u++ )
//      {    if ( used2[u] || used2[ to_rc2[u] ] ) continue;
//           int nkmers = unibases2[u].isize( ) - K2 + 1;
//           if ( nkmers < 1000 ) continue;
//           if ( to_rc2[u] < u ) continue;
//           superb s;
//           s.SetNtigs(1);
//           s.SetTig( 0, u );
//           s.SetLen( 0, unibases2[u].size( ) );
//           uscaffolds2.push_back(s);
//           superb sx;
//           sx.SetNtigs(1);
//           sx.SetTig( 0, tigse.size( ) );
//           sx.SetLen( 0, unibases2[u].size( ) );
//           scaffolds.push_back(sx);
//           tigse.push_back( efasta( unibases2[u].ToString( ) ) );
//           circled.push_back(False);
//           reports.push_back( ToString(u) + "[" + ToString(nkmers) + "]\n" );
//           extras.push_back( vec<int>( ) );    }

//      // Remove scaffolds whose base unipaths are subsumed in other scaffolds.
//      // Currently implemented as O(n^2), which is nuts.

//      vec<Bool> to_remove( scaffolds.size( ), False );
//      for ( int i1 = 0; i1 < scaffolds.isize( ); i1++ )
//      {    if ( to_remove[i1] ) continue;
//           for ( int i2 = 0; i2 < scaffolds.isize( ); i2++ )
//           {    if ( i2 == i1 ) continue;
//                {    vec<int> u2;
//                     for ( int j = 0; j < uscaffolds2[i2].Ntigs( ); j++ )
//                          u2.push_back( uscaffolds2[i2].Tig(j) );
//                     Sort(u2);
//                     if ( BinSubset( u2, extras[i1] ) ) 
//                          to_remove[i2] = True;    }    }    }
//      EraseIf( scaffolds, to_remove );
//      EraseIf( circled, to_remove );
//      EraseIf( reports, to_remove );

//      // Dump logging.

//      if ( VERBOSITY >= 1 ) 
//      {    for ( int i = 0; i < reports.isize( ); i++ )
//           {    cout << "\n" << double_line << "\nscaffold [" << i << "]";
//                if ( circled[i] ) cout << " - circular";
//                cout << "\n" << reports[i];    }    }

//      // Print unibase graphs.

//      if (PRINT_GRAPH1)
//      {    vec< vec<int> > nexts1;
//           GetNexts( K1, unibases1, nexts1 );
//           cout << "\nunibases1 graph:\n\n";
//           for ( int u = 0; u < nuni1; u++ )
//           {    cout << u << "[l=" << unibases1[u].isize( ) - K1 + 1
//                     << ",rc=" << to_rc1[u] << "] -->";
//                for ( int j = 0; j < nexts1[u].isize( ); j++ )
//                     cout << " " << nexts1[u][j];
//                cout << "\n";    }
//           cout << "\n";    }
//      if (PRINT_GRAPH2)
//      {    cout << "\nunibases2 graph:\n\n";
//           for ( int u = 0; u < nuni2; u++ )
//           {    cout << u << "[l=" << unibases2[u].isize( ) - K2 + 1
//                     << ",rc=" << to_rc2[u] << "] -->";
//                for ( int j = 0; j < nexts2[u].isize( ); j++ )
//                     cout << " " << nexts2[u][j];
//                cout << "\n";    }
//           cout << "\n";    }

//      // Write assembly.

//      if (WRITE)
//      {    cout << "\n" << Date( ) << ": writing assembly" << endl;
//           Assembly A( scaffolds, tigse );
//           A.remove_unused_contigs( );
//           A.WriteAll( sub_dir + "/" + SCAFFOLDS_OUT );    }

//      // Print summary statistics.

//      int ntigs = 0;
//      int64_t total_bases = 0;
//      for ( int j = 0; j < scaffolds.isize( ); j++ )
//      {    const superb& s = scaffolds[j];
//           ntigs += s.Ntigs( );
//           for ( int r = 0; r < s.Ntigs( ); r++ )
//                total_bases += tigse[ s.Tig(r) ].Length1( );    }
//      cout << "\nSummary:\ntotal scaffolds = " << scaffolds.size( )
//           << ", total contigs = " << ntigs 
//           << ", total bases = " << ToStringAddCommas(total_bases) << endl;
//      cout << "time used = " << TimeSince(clock) << endl;    }
}
