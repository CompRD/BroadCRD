///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/**
   Program: ComputeReadPairSeps

   Determines average separation and standard deviation of paired read libraries.
   This is done by placing read k-mer paths on a set of unipaths. The process has two 
   stages: i) try to place each pair of reads on the same unipath. If enough pairs
   (SAMPLE_SIZE) were placed this way compute separation stats. Else ii) try to 
   compute distances between different unipaths on which reads of each pair were 
   placed using unipath adjacency graph. This step is not currently computationally 
   efficient. In most cases it is not going to be needed.

   WARNING!!!!!!!!!! pathing of paired reads and unipaths must be consistent!!!!!!!!!!

   INPUT FILES:
     (unipaths reads).unipaths.k*
     (unipaths reads).unipathsdb.k*
     (unipaths reads).unipath_adjgraph.k*
     (paired reads).paths.k*
     (paired reads).pairs

   OUTPUT FILES:
     OUT_HEAD.pairs

   @file

*/


#include "MainTools.h"
#include "Map.h"
#include "Set.h"
#include "paths/KmerPath.h"
#include "graph/Digraph.h"
#include "PairsManager.h"
#include "graph/FindCells.h"
#include "paths/SeedNeighborhood.h"
#include "feudal/BinaryStream.h"

// auxiliary routines    -------------------------------
//
inline void compute_edges( const digraph& AG, digraphE<fsepdev>& LG, ulonglong nunipaths, const vecKmerPath& paths , longlong& nCells){


  const int max_cell_size = 10;
  vec< vec<int> > cells;
  cout << Date( ) << ": Begin FindCells" << endl;
  FindCells( AG, max_cell_size, cells );
  cout << Date( ) << ": End FindCells" << endl;
  PRINT( cells.size() );
  nCells = cells.size();
  
  {
    vec< vec<int> > d( nunipaths );
    vec< fsepdev > edges;
    LG.Initialize( d, d, edges, d, d);
  }
  
  for ( int i = 0; i < cells.isize( ); i++ ){    
    vec< vec<int> > allpaths;
    int v = cells[i].front( ), w = cells[i].back( );
    if ( LG.HasEdge( v, w ) ) continue;
    set< pair< vec<int>, set<int> > > partials;
    vec<int> just;
    set<int> justs;
    just.push_back(v), justs.insert(v);
    partials.insert( make_pair( just, justs ) );
    set< vec<int> > pathss;
    while( !partials.empty( ) ){    
      pair< vec<int>, set<int> > p = *partials.begin( );
      partials.erase( partials.begin( ) );
      int x = p.first.back( );
      Bool to_process = True;
      if ( x == w ) {    
	pathss.insert( p.first );
	to_process = False;    
      }
      if (to_process){    
	for ( int i = 0; i < AG.From(x).isize( ); i++ ){    
	  int y = AG.From(x)[i];
	  if ( Member( p.second, y ) ) continue;
	  pair< vec<int>, set<int> > q = p;
	  q.first.push_back(y);
	  q.second.insert(y);
	  partials.insert(q);    
	}    
      }    
    }
    allpaths.reserve( pathss.size( ) );
    for ( set< vec<int> >::iterator i = pathss.begin( ); 
	  i != pathss.end( ); ++i ){    
      allpaths.push_back(*i);    
    }
    
    if ( allpaths.empty( ) ) continue;
    vec<int> dist;
    for ( int i = 0; i < allpaths.isize( ); i++ ){    
      int d = 0;
      for ( int j = 0; j < allpaths[i].isize( ) - 1; j++ ){
	ulonglong pathi = allpaths[i][j];
	if ( pathi < paths.size() )
	  d += paths[ allpaths[i][j] ].KmerCount( );
      }
      dist.push_back(d);    
    }
    int m = Min(dist), M = Max(dist);
    Sort(dist);
    LG.AddEdge( v, w, fsepdev( (m+M)/2, (M-m)/2 ) );    
  }  
}

int main( int argc, char *argv[] ){
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_String_Doc(UREADS_HEAD, 
    "looks for RUN_DIR/UREADS_HEAD.unipath_adjgraph");
  CommandArgument_String_Doc(READS_HEAD, 
    "looks for RUN_DIR/READS_HEAD.{paths, paths_rc, pathsdb, pairs}");
  CommandArgument_Int(K);
  CommandArgument_Int_OrDefault_Doc(TARGET_SAMPLE_SIZE, 10000, "desired sample size");
  CommandArgument_Int_OrDefault_Doc(MIN_SAMPLE_SIZE, 500, "minimum sample size");
  CommandArgument_Int_OrDefault_Doc(ORIG_DEV_BRACKET, 5, 
    "Radius (in std dev) of data censoring around original (lab specified) mean paired read separation.");
  CommandArgument_Int_OrDefault_Doc(RANDOM_SEED, 133333,
    "Seed value for the random generator." );
  CommandArgument_String_OrDefault(OUT_HEAD, "");
  CommandArgument_Bool_OrDefault(VERBOSE, False);
  EndCommandArguments;
  

  ForceAssertGe( TARGET_SAMPLE_SIZE, MIN_SAMPLE_SIZE );
  
  String KS = ToString( K );
  
  String data_dir = PRE + "/" + DATA;
  String run_dir = data_dir + "/" + RUN;
  
  String full_ureads_head       = run_dir + "/" + UREADS_HEAD;
  String unipaths1_file         = full_ureads_head + ".unipaths.k" + KS;
  String unipaths1db_file       = full_ureads_head + ".unipathsdb.k" + KS;
  String unipath_adjgraph_file  = full_ureads_head + ".unipath_adjgraph.k" + KS;
  
  String full_reads_head        = run_dir + "/" + READS_HEAD;
  String paths_file             = full_reads_head + ".paths.k" + KS;
  String pairs_file             = full_reads_head + ".pairs";
  
  vecKmerPath unipaths1( unipaths1_file );
  BREAD2( unipaths1db_file, vec<tagged_rpint>, unipaths1db );
    
  digraph AG; 
  digraphE<fsepdev> LG;

  cout << Date() << ": reading paths" << endl;
  vecKmerPath paths( paths_file );


  cout << Date() << ": loading pairing information..." << endl;
  longlong n_reads_total = 
    MastervecFileObjectCount( full_reads_head + ".fastb" );
  PairsManager pairs( pairs_file );
  size_t npairs     = pairs.nPairs();
  size_t nreads     = pairs.nReads();
  size_t nLibraries = pairs.nLibraries();

  PRINT3( nLibraries, nreads, npairs );
  

  cout << Date() << ": randomizing data set" << endl;
  srand(RANDOM_SEED);
  vec<longlong> perm( npairs );
  {
    if ( 1 ){
      for ( size_t i = 0; i < npairs; i++)
	perm[ i ] = (double) rand() / RAND_MAX * npairs; 

    }else{
      vec<longlong> randoms( npairs );
      for ( size_t i = 0; i < npairs; i++){
	randoms[i] = rand();
	perm[i] = i;
      }
      cout << Date() << ": sorting" << endl; 
      SortSync(randoms, perm);
    }
  }



  cout << Date() << ": computing separations" << endl;
  vec< vec< float> >   pairSepData( nLibraries );
  vec< StdMap<String,longlong> > libStats( nLibraries );
  longlong nInfoCheck = npairs/20;
  PRINT( perm.size() );
  vec<char> algnd2diffUnis( perm.size(), 0 );
  for ( int pass = 1; pass <= 2; pass++ ){

    if ( pass == 2 ){
      
      uint nfound = 0;
      for ( size_t il = 0; il < nLibraries; il++ ){
	if ( pairSepData[ il ].isize() >= TARGET_SAMPLE_SIZE )
	  nfound++;
      }
      if ( nfound == nLibraries ){
	cout << "no need to consider unipath adjacency graph, there is enough data already." << endl;
	for ( size_t il = 0; il < nLibraries; il++ )
	  PRINT2( il, pairSepData[il].size() );
	break;
      }else{
	cout << "have to consider unipath adjacency graph, not enough data yet." << endl;
	for ( size_t il = 0; il < nLibraries; il++ )
	  PRINT2( il, pairSepData[il].size() );
      }

      cout << Date() << Date() << ": reading unipath adjacency graph" << endl;
      BinaryReader::readFile( unipath_adjgraph_file, &AG );
      
      // compute unipath graph edges for graph cells
      longlong nCells = 0;
      compute_edges( AG, LG, unipaths1.size(), paths, nCells );
      if ( nCells == 0 ){
	cout << "WARNING: no cells found in the adjacency graph -> could not compute graph based separation of read pairs." << endl;
	break;
      }
    }

    for ( ulonglong  mi= 0;  mi < perm.size(); mi++ ){
      if ( perm.size( ) >= 100 && mi % (perm.size()/100) == 0 ){ 
	Dot( cout, mi / (double)perm.size() * 100.0);
      }
      if ( mi == perm.size() -1 ) cout << endl;
      if ( pass == 2 && ! algnd2diffUnis[ mi ] )
	continue;
      int pi   = perm[ mi ];
      int lib  = pairs.libraryID( pi );
      
      if ( pairSepData[ lib ].isize() >= TARGET_SAMPLE_SIZE )
	continue;  // enough samples accumulated for this library

      longlong id1 = pairs.ID1(pi);
      vec< pair<int,int> > uo1; 
      KmerPath& p1 = paths[ id1 ];
      
      kmer_id_t k01 = p1.Segment(0).Start();
      vec<longlong> locs1;
      Contains( unipaths1db, k01, locs1 );
      ForceAssert( locs1.size() <= 1 );
      

      longlong id2 = pairs.ID2(pi);
      vec< pair<int,int> > uo2;
      KmerPath p2 = paths[ id2 ];
      p2.Reverse();
     
      kmer_id_t k02 = p2.Segment(0).Start();
      vec<longlong> locs2;
      Contains( unipaths1db, k02, locs2 );
      ForceAssert( locs2.size() <= 1 );
      
      if ( pass == 1 ) libStats[lib]["nPairsSmpld"]++;

      if ( locs1.size() != 1 || locs2.size() != 1 ) continue; // one of the reads not aligned

      if ( pass == 1 ) libStats[lib]["nPairsAlgnd"]++;

      Bool sameUni = 
	unipaths1db[locs1[0]].PathId() == unipaths1db[ locs2[0] ].PathId() ? True : False; 

      if ( pass == 1 ){
	if ( sameUni ){
	  libStats[lib]["nPairsSameUni"]++;
	  libStats[lib]["nPairsSameCntg"]++;
	}
	else libStats[lib]["nPairsDiffUni"]++;

	algnd2diffUnis[mi] = 1;
      }
      
      
      if ( VERBOSE && mi > 0 && 
	   (mi % nInfoCheck == 0 || mi == 10000 || mi == perm.size() -1) ) {
	for ( size_t library = 0; library < nLibraries; library++)
	  PRINT2( library, pairSepData[library].size() );
      }

      tagged_rpint& t1 = unipaths1db[ locs1[0] ];
      longlong uid1    = t1.PathId( );
      longlong offset1 = k01 - t1.Start();
      for ( int r = 0; r < t1.PathPos( ); r++ )
	offset1 += unipaths1[uid1].Segment(r).Length( );
      uo1.push_back( make_pair( uid1, int(offset1) ) );  

      tagged_rpint& t2 = unipaths1db[ locs2[0] ];
      longlong uid2    = t2.PathId( );
      longlong offset2 = k02 - t2.Start( );
      for ( int r = 0; r < t2.PathPos( ); r++ )
	offset2 += unipaths1[uid2].Segment(r).Length( );
      uo2.push_back( make_pair( uid2, int(offset2) ) );
          
      int origSep = pairs.getLibrarySep( lib );
      int origDev = pairs.getLibrarySD( lib );

      //cout << "computing distance" << endl;
      
      int dist = 0, insertSize = 0, sep = 0;
      Bool distComputed = False;
      if ( pass == 1 && uid1 == uid2 ){
	dist = uo2[ 0 ].second - uo1[ 0 ].second;

	if ( dist < 0 ) libStats[ lib ]["nOuties"]++;
	else{
	  libStats[ lib ]["nInnies"]++;
	  sep = dist - paths[ id1 ].TotalLength() - K;
	  if ( abs(sep - origSep) < ORIG_DEV_BRACKET * origDev )
	    pairSepData.at( lib ).push_back( sep );
	}
	
      }else if ( pass == 2 ){  
  
	vec< vec<int> > lpaths;
	if ( ! LG.HasEdge(uid1, uid2) )
	  continue;	   

	libStats[lib]["nPairsSameCntg"]++;
	LG.EdgePaths( uid1, uid2, lpaths );
	
	vec<int> lpathsLens( lpaths.isize(), 0 );
	for ( int p = 0; p < lpaths.isize(); p++ ){
	  for ( int k = 1; k < lpaths[p].isize() -1; k++ )
	    lpathsLens[p] += unipaths1[ lpaths[p][k] ].TotalLength();
	}
	Sort( lpathsLens);
	
	if ( lpathsLens.size() == 0 ) 
	  dist = 0;
	else if ( lpathsLens.size() == 1 )
	  dist = lpathsLens[ 0 ];
	else {
	  int n = lpathsLens.isize();
	  dist = ( n%2 == 1 ) ? lpathsLens[n/2] : 
	    ( lpathsLens[n/2] + lpathsLens[n/2 -1] )/ 2;
	}
	
	dist += unipaths1[uid1].TotalLength() 
	  - uo1[ 0 ].second + uo2[ 0 ].second;
	sep = dist - paths[ id1 ].TotalLength() - K;

	if ( dist < 0 ) libStats[ lib ]["nOuties"]++;
	else{
	  libStats[ lib ]["nInnies"]++;
	  sep = dist - paths[ id1 ].TotalLength() - K;
	  if ( abs(sep - origSep) < ORIG_DEV_BRACKET * origDev )
	    pairSepData.at( lib ).push_back( sep );
	}
	
      }      
    }
  }

  // compute library stats
  cout << Date() << "\nprocessing " << pairSepData.size( ) << " libraries" << endl;
  vec<String> lib2warning( nLibraries, "" );
  for ( int lib = 0; lib < pairSepData.isize(); lib++ ){
    
    int origSep = pairs.getLibrarySep( lib );
    int origDev = pairs.getLibrarySD( lib );

    libStats[lib]["libNo"]       = lib;
    libStats[lib]["origSep"]     = origSep;
    libStats[lib]["origDev"]     = origDev;
    libStats[lib]["fnlSmplSize"]  = pairSepData[lib].size();
    
    if ( (double) libStats[lib]["fnlSmplSize"] / (double) libStats[lib]["nPairsSameCntg"] < 0.5 )
      lib2warning[lib] += " *SMALL FRACTION OF CONTIG PAIRS IN SAMPLE*";

    if ( libStats[lib]["nInnies"] < libStats[lib]["nOuties"] )
      lib2warning[lib] += " *FEWER INNIES THAN OUTIES*";


    if ( pairSepData[lib].isize() < MIN_SAMPLE_SIZE ){
      cout << "WARNING: cannot compute separation stats, too few samples in library= " << lib << ".\nCopying input values.";
      libStats[lib]["newSep"]    = libStats[lib]["origSep"];
      libStats[lib]["newDev"]    = libStats[lib]["origDev"];
      libStats[lib]["newMedian"] = -99999;
      lib2warning[lib]          += " *SAMPLE TOO SMALL*";
    }else{
      
      if ( pairSepData[lib].isize() < TARGET_SAMPLE_SIZE )
	lib2warning[lib] += " *SAMPLE SIZE BELOW TARGET*";

      longlong n  = pairSepData[ lib ].size();
      ForceAssertGt( n, 0 );

      Sort(pairSepData[ lib ]);   
      
      NormalDistribution res = SafeMeanStdev( pairSepData[lib] );
      float sepAve = res.mu_;
      float sepDev = res.sigma_;
      pairs.changeLibrarySepSd( lib, sepAve, sepDev );

      float sepMedian = ( n%2 == 1 ) ? pairSepData[lib][n/2] : 
	( pairSepData[lib][n/2] + pairSepData[lib][n/2 -1] )/ 2;
      
      libStats[lib]["newSep"]    = sepAve;
      libStats[lib]["newDev"]    = sepDev;
      libStats[lib]["newMedian"] = sepMedian;

    }
  }
  
  
  vec<String> titles;
  titles.push_back("libNo", "origSep", "origDev", "newSep", "newDev", "newMedian" );
  titles.push_back("nPairsSmpld", "nPairsAlgnd", "nPairsSameUni", "nPairsDiffUni", "nPairsSameCntg" );
  titles.push_back("nInnies", "nOuties", "fnlSmplSize");


  cout << "                     --------------- STATS -------------------\n";
  titles.Print( cout ); cout << " warning"; cout << endl;
  for ( size_t lib = 0; lib < nLibraries; lib++ ){
    
    for (unsigned it = 0; it < titles.size(); it++ ){
      String value = 
	libStats[lib][ titles[it] ] != -99999 ? ToString( libStats[lib][ titles[it] ] ) : "NA";
      cout.width( titles[it].size() ); cout.fill(' '); 
      cout << value << " ";
    }
    cout << lib2warning[lib];
    cout << endl;
  }
  cout << "--------------------------------------------------------------\n"; cout << endl;


  
  if ( ! OUT_HEAD.empty() ){
    cout << Date() << ": writing pairs" << endl;
    //vec<read_pairing> new_pairs = pairs.convert_to_read_pairings();
    //WritePairs( run_dir, new_pairs, nreads, False, OUT_HEAD );
    pairs.Write( run_dir + "/" + OUT_HEAD + ".pairs" );
    cout << Date() << ": done!" << endl;
  }

}





   
