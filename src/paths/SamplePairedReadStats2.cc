///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/**
   Program: SamplePairedReadStats

   Determines separation stats, orientation, etc. of paired reads.
   This is done by aligning a random sample of pairs to a reference

   INPUT FILES:
     reference.fastb
     (paired reads).pairs

   OUTPUT FILES:
     READS_HEAD.outies (stats for incorrectly oriented pairs)
     if OUT_HEAD nonempty then following files are written
     OUT_HEAD.pairs
     if WRTE_STATS writes READS_HEAD.sample_stats containing detailed stats
     
   @file

*/


#include "MainTools.h"
#include "PairsManager.h"
#include "lookup/FirstLookupFinderECJ.h"
#include "lookup/LookAlign.h"
#include "lookup/LookupTabBuilder.h"
#include "random/Shuffle.h"
#include "paths/UnibaseUtils.h"
#include "VecUtilities.h"
#include "Map.h"
#include "paths/GetNexts.h"
#include "paths/HyperFastavector.h"
#include "paths/CLinFastavec.h"
#include "paths/Unipath.h"
#include "paths/ReadsToPathsCoreX.h"

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

int N50hbv( HyperBasevector& hbv ){
  vec<int> lens( hbv.Edges().size() );
  for ( size_t u = 0; u < hbv.Edges().size(); u++ )
    lens[u] = hbv.Edges()[u].isize();
  PRINT( Max(lens) );
  return N50(lens);
}

// auxiliary routines    -------------------------------
void pair_alignment_data( const PairsManager& pairs, const vecbasevector& reads, const VecQualNibbleVec &quals,
			  FirstLookupFinderECJ& lfinder, FirstLookupFilterECJ& lfilter,
			  const longlong pi, const int LookupK, const Bool use_quals, 
			  Bool& Read1Algnd, Bool& Read2Algnd, Bool& PairAlgnd, Bool& PairUnique, 
			  Bool& PairSameContig, Bool& PairLogical, Bool& PairInnie, Bool& PairOutie,
			  Bool FLIP, Bool TRIM, int64_t& sep );
double gm( vec<float> data, double& mean, double& sigma, double& stdErr, const Bool& verbose);
void actual_gm( const vec<float>& data, vec<double>& means, vec<double>& sigmas, vec<double>& weights, const Bool& verbose );
//----------------------------------------------------------

// Fraction of seps falling within N sigma of mean. --bruce
double FractionInDist(const vec<float> &data, const double mean, const double sigma, const double N=3)
{
  ForceAssertGt(sigma, 0.0);
  int count = 0;
  for (u_int i = 0; i < data.size(); ++i) {
    double x = (data[i]-mean) / sigma;
    if (abs(x) < N) ++count;
  }
  return (double)count / (double)data.size();
}

size_t filter_outliers( vec<float>& data, double zet_start, double zet_end ){
  
  ForceAssert( zet_start >= zet_end );
  Sort(data);
  {
    NormalDistribution distr0 = SafeMeanStdev( data );
    PRINT2( distr0.mu_, distr0.sigma_ ); 
  }
  PRINT3( data.front(), data.back(), data.size() );
  { // remove border values; sometimes there could be many data points of the same value 
    // at the edges, perhaps an alignment artifact
    cout << "Removing border values" << endl;
    size_t i2 = data.size() -1;
    float middleV = data.at( round( 0.5 * data.size() ) );
    if ( data.back() > middleV )
      while( i2 > 0 && data[i2] == data.back() ) i2--;
    size_t i1 = 0;
    if ( data.front() < middleV )
      while( i1 < data.size() -1 && data[i1] == data.front() ) i1++;
    
    ForceAssertGe( i2, i1 );
    PRINT3( i2, i1, data.size() );
    data.SetToSubOf( data, i1, i2 - i1 + 1 );

    {
      NormalDistribution distr1 = SafeMeanStdev( data );
      int median = Median( data );
      PRINT3( distr1.mu_, distr1.sigma_, median ); 
    }

  }


  // filter outliers
  size_t nremoved = 0;
  cout << "filtering outliers" << endl;
  PRINT3( data.front(), data.back(), data.size() );
  double mean, sigma;
  for ( int i = zet_start; i >= zet_end; i-- ){
    double zet_cut = (double)i;
    PRINT( zet_cut );
    Bool RemovedData = True;
    while( RemovedData == True ){
      RemovedData = False;
      vec<float> cdata = data;
      size_t k = 0;
      for ( k = 0; k < data.size() -1; k++ )
	if ( data[k+1] != data.front() || data[k+1] == data.back() )
	  break;
      cdata.Erase(0, k+1 );
      NormalDistribution distr = SafeMeanStdev( cdata );
      int median = Median( cdata );
      mean  = distr.mu_;
      sigma = distr.sigma_;
      double zetB = abs( (data.front() - median)/sigma );
      if ( zetB > zet_cut ){
	//PRINT( zetB );
	nremoved += k+1;
	RemovedData = True;
	data = cdata;
      }

      cdata = data;
      for ( k = data.size() -1; k > 0; k-- )
	if ( data[k-1] != data.back() || data[k-1] == data.front() )
	  break;
      cdata.Erase(k, data.size() );
      distr = SafeMeanStdev( cdata );
      mean  = distr.mu_;
      sigma = distr.sigma_;
      median - Median( cdata );
      double zetE = abs( (data.back() - median)/sigma );
      if ( zetE > zet_cut ){
	//PRINT( zetE );
	nremoved += data.size() - k -1;
	RemovedData = True;
	data = cdata;
      }

    }
    
    {
      PRINT2( data.front(), data.back() );
      NormalDistribution distr2 = SafeMeanStdev( data );
      PRINT2( distr2.mu_, distr2.sigma_ ); 
    }
  }
    
  PRINT4( data.front(), data.back(), data.size(), nremoved );

  return nremoved;
}


int main( int argc, char *argv[] ){
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_String_OrDefault_Doc(REF,"", 
    "looks for RUN_DIR/REF.fastb");
  CommandArgument_String_OrDefault_Doc(UNIBASES,"", 
    "looks for RUN_DIR/UNIBASES as a reference");
  CommandArgument_UnsignedInt_OrDefault_Doc(UNIBASES_K,0, 
    "needed by involution algorithm");
  CommandArgument_String_Doc(READS_HEAD, 
    "looks for RUN_DIR/READS_HEAD.{fastb, pairs/pairs}");
  CommandArgument_Bool_OrDefault_Doc( FLIP, True,
    "reverse-complement reads for alignment; use for jumps, not frags.");
  CommandArgument_String_OrDefault_Doc(OUT_HEAD, "", 
    "Name for output files. No output files written if not specified");
  CommandArgument_Bool_OrDefault_Doc(WRITE_STATS, True, 
    "Write sample stats into READS_HEAD.sample_stats file.");
  CommandArgument_Bool_OrDefault_Doc(WRITE_SEPS, True, 
    "Write sampled separations into READS_HEAD.seps.{in,out} files.");
  CommandArgument_Int_OrDefault_Doc(ORIG_DEV_BRACKET, 5, 
    "Radius (in std dev) of data censoring around original (lab specified) mean paired read separation.");

  CommandArgument_Bool_OrDefault_Doc( NEW_LOOKUP_TABLE, True,
    "build new lookup table even if a file already exists.");

  CommandArgument_Int_OrDefault_Doc(TRIM, 0, 
    "Trim that many bases from the beginning of a read (just before alignment).");
  CommandArgument_Int_OrDefault_Doc(MAX_TARGET_SAMPLE_SIZE, 300000, 
    "don't try to get bigger sample size");
  CommandArgument_Int_OrDefault_Doc(MIN_TARGET_SAMPLE_SIZE, 10000, 
    "try to get at least that sample size");
  CommandArgument_Int_OrDefault_Doc(MIN_SAMPLE_SIZE, 100, 
    "minimum sample size for computing stats (note: this could be samller than MIN_TARGET_SAMPLE_SIZE)");
  CommandArgument_Double_OrDefault_Doc(TARGET_UNIBASE_COVERAGE, .2,
    "unibases covering this fraction of total are enough");
  CommandArgument_Int_OrDefault_Doc(MAX_SEPARATION, 100000, 
    "maximum aligned pair separation considered");
  CommandArgument_UnsignedInt_OrDefault_Doc(MIN_UNIBASE_LENGTH, 0, 
    "only include unibases of this length or longer for alignments");
  CommandArgument_Bool_OrDefault_Doc(USE_QUALS, False, 
    "Use quality scores in alignments");
  CommandArgument_Bool_OrDefault_Doc(PERFECT, True, 
    "Require perfect alignments (we usually use corrected reads on a subset of a reference)");
   CommandArgument_Bool_OrDefault_Doc(FLATTEN, True, 
    "Flatten unibase graph to get longer unibases");
  CommandArgument_Int_OrDefault_Doc(RANDOM_SEED, 133333,
    "Seed value for the random generator." );
  CommandArgument_Bool_OrDefault(VERBOSE, False);
  CommandArgument_Bool_OrDefault_Doc(WRITE_NOTHING, False, 
    "Don't write any ouput files.");
  EndCommandArguments;
  
  ForceAssertGe( MAX_TARGET_SAMPLE_SIZE, MIN_TARGET_SAMPLE_SIZE );
  ForceAssertGe( MIN_TARGET_SAMPLE_SIZE, MIN_SAMPLE_SIZE );
  ForceAssertGt( MIN_TARGET_SAMPLE_SIZE, 0 );

  ForceAssert( !REF.empty() || !UNIBASES.empty() );
  if ( !UNIBASES.empty() )
    ForceAssertNe( UNIBASES_K, 0u );
  
  int K = UNIBASES_K;
  int LookupK = 12;

  String data_dir = PRE + "/" + DATA;
  String run_dir = data_dir + "/" + RUN;
  
  String full_ref_head          = run_dir + "/" + REF;
  String ref_fastb_file         = full_ref_head + ".fastb";
  
  String full_reads_head        = run_dir + "/" + READS_HEAD;
  String reads_fastb            = full_reads_head + ".fastb";
  String quals_file             = full_reads_head + ".qualb";
  String pairs_file             = full_reads_head + ".pairs";
  String pairto_file            = full_reads_head + ".pairto";
  
 //  String results_file           = OUT_HEAD.empty() ? full_reads_head + ".sample_stats" : run_dir + "/" + OUT_HEAD + ".sample_stats";
//   String seps_file              = OUT_HEAD.empty() ? full_reads_head + ".sample_seps" : run_dir + "/" + OUT_HEAD + ".sample_seps"; 
//   String outies_file            = OUT_HEAD.empty() ? full_reads_head + ".outies" : run_dir + "/" + OUT_HEAD + ".outies";
  String results_file           = full_reads_head + ".sample_stats";
  String seps_file              = full_reads_head + ".sample_seps"; 
  String outies_file            = full_reads_head + ".outies";


  // get reference
  vecbasevector ref;
  String out_lookup_file;
  if ( UNIBASES.empty() ){
    cout << "Using reference file" << endl;
    out_lookup_file = full_ref_head + ".lookuptab";
    if ( ! IsRegularFile( out_lookup_file ) || NEW_LOOKUP_TABLE ) {
      cout << Date( ) << ": creating and saving lookup table " << out_lookup_file << endl;
      LookupTabBuilder lookup( LookupK );
      lookup.addFastb( ref_fastb_file.c_str( ) );
      lookup.write( out_lookup_file.c_str( ) );
    }
    cout << Date() << ": reading reference" << endl;
    ref.ReadAll( ref_fastb_file );
  }else{
    cout << "Using unibases as reference" << endl;

    
    String uni_fastb_file = run_dir + "/" + UNIBASES;
    
    if ( FLATTEN ){

      vecbasevector unibases(uni_fastb_file);
      int nuni = unibases.size();

      size_t maxUniLen = 0; 
      vec<int> ulens( unibases.size() );
      for ( size_t u = 0; u < unibases.size(); u++ ){
	ulens[u] = unibases[u].isize();
	if ( unibases[u].size() > maxUniLen )
	  maxUniLen = unibases[u].size();
      }
      PRINT2( maxUniLen, N50(ulens) );
      
      cout << Date() << ": getting nexts for unibases" << endl;
      vec< vec<int> > nexts;
      GetNexts( K, unibases, nexts );
      ForceAssertEq( nexts.isize(), nuni);
      
      vec< vec<int> > from(nuni), to(nuni);
      for ( int i = 0; i < nexts.isize(); i++ ){
	for ( int j = 0; j < nexts[i].isize(); j++ ){
	  from[i].push_back(nexts[i][j]);
	  to[ nexts[i][j] ].push_back( i );
	}
      }
      
      for ( int i = 0; i < from.isize(); i++){
	Sort( from[i] );
	ForceAssertLe( from[i].size(), 4u );
	Sort( to[i] );
	ForceAssertLe( to[i].size(), 4u );
      }
      
      cout << Date() << ": building adjacency graph" << endl;
      digraph AG;
      AG.Initialize( from, to );
      
      cout << Date() << ": building hyperbasevector" << endl;
      HyperBasevector hbv;
      BuildUnibaseAdjacencyHyperBasevector( K-1, AG, unibases, hbv );
      AG.Clear();
      PRINT( N50hbv( hbv ) );
      PRINT2( hbv.N(), hbv.TotalEdgeLength() );
      hbv.RemoveSmallComponents( 1000 );
      PRINT2( hbv.N(), hbv.TotalEdgeLength() );
      PRINT( N50hbv( hbv ) );
      
      cout << Date() << ": popping bubbles" << endl;
      int npopped = 0;
      while( 1 ) {
	vec<int> bubbleLocs; bubbleLocs.reserve( hbv.N() );
	for ( int v = 0; v < hbv.N(); v++ ){
	  if ( hbv.From(v).size() == 2 && hbv.From(v)[0] == hbv.From(v)[1] ){
	    int w = hbv.From(v)[0];
	    if ( hbv.To(v).size() > 2 || hbv.From(w).size() > 2 )
	      continue;
	    int u0 = hbv.EdgeObjectIndexByIndexFrom( v, 0 );
	    int u1 = hbv.EdgeObjectIndexByIndexFrom( v, 1 );
	    if ( unibases[u0].size() == unibases[u1].size() )
	      bubbleLocs.push_back( v );
	  }
	}
	if ( bubbleLocs.size() == 0 )
	  break;
	npopped += bubbleLocs.size();
	hbv.PopBubbles( bubbleLocs );
      }
      PRINT( npopped );
      
      hbv.RemoveUnneededVertices();
      hbv.RemoveDeadEdgeObjects();
      hbv.RemoveEdgelessVertices();
      PRINT( N50hbv( hbv ) );
      {
	vec<int> to_delete;
	// remove incoming
	for ( int v = 0; v < hbv.N(); v++ ){
	  if ( hbv.To(v).size() == 0 && hbv.From(v).size() == 1 ){ 
	    size_t w = hbv.From(v)[0];
	    if ( hbv.To(w).size() > 1 ){
	      for ( size_t s = 0; s < hbv.To(w).size(); s++ ){
		if ( hbv.To(s).size() > 0 ){ 
		  for ( size_t i = 0; i < hbv.EdgesBetween( v, w ).size(); i++ ){
		    size_t uid = hbv.EdgesBetween(v,w)[i];
		    if ( hbv.Edges()[uid].size() < 200 )
		      to_delete.push_back( uid );
		  }
		  break;
		}
	      }
	    }
	  }
	}
	int nBadIncoming = to_delete.size();
	// remove outgoing 
	for ( int w = 0; w < hbv.N(); w++ ){
	  if ( hbv.From(w).size() == 0 && hbv.To(w).size() == 1 ){
	    size_t v = hbv.To(w)[0];
	    if ( hbv.From(v).size() > 1 ){
	      for ( size_t s = 0; s < hbv.From(v).size(); s++ ){
		if ( hbv.From(s).size() > 0 ){
		  for ( size_t i = 0; i < hbv.EdgesBetween( v, w ).size(); i++ ){
		    size_t uid = hbv.EdgesBetween(v,w)[i];
		    if ( hbv.Edges()[uid].size() < 200 )
		      to_delete.push_back( uid );
		  }
		  break;
		}
	      }
	    }
	  }
	}
	
	 cout << "Shaving function identified " << to_delete.size()
	      << " edges to delete (" << nBadIncoming
	      << " incoming and " << to_delete.size() -nBadIncoming
	      << " outgoing)" << endl;
	
	 hbv.DeleteEdges( to_delete );
	 hbv.RemoveUnneededVertices();
	 hbv.RemoveDeadEdgeObjects();
	 hbv.RemoveEdgelessVertices();
      }
      PRINT3( hbv.N(), hbv.TotalEdgeLength(), hbv.EdgeObjectCount() );
      PRINT( N50hbv( hbv ) );

      vecbasevector lunibases; 
      lunibases.reserve( hbv.Edges().size() );
      Bool REPATH = False;
      if ( REPATH ){
	{	
	  vecbasevector reads( hbv.Edges().size() );
	  for ( size_t i = 0; i < hbv.Edges().size(); i++ )
	    reads[i] = hbv.Edges()[i];
	  
	  cout << Date() << ": pathing" << endl;
	  vecKmerPath paths, paths_rc;
	  vec<tagged_rpint> pathsdb;
	  int n_threads = 12;
	  ReadsToPathsCoreY( reads, K, paths, paths_rc, pathsdb, "", n_threads);
	  
	  cout << Date() << ": unipathing" << endl;
	  vecKmerPath unipaths;
	  vec<tagged_rpint> unipathsdb;
	  Unipath( paths, paths_rc, pathsdb, unipaths, unipathsdb, True);
	  
	  cout << Date() << ": building KmerBaseBroker (needed for unibases)" << endl;
	  KmerBaseBroker kbb( K, paths, paths_rc, pathsdb, reads );
	  
	  
	  cout << Date() << ": building unibases" << endl;
	  lunibases.reserve(unipaths.size());
	  for ( size_t i = 0; i < unipaths.size(); i++ )
	    lunibases.push_back( kbb.Seq( unipaths[i] ) );
	}
      }else{
	lunibases.reserve( 2 * lunibases.size() );
	for ( size_t i = 0; i < hbv.Edges().size(); i++ ){
	  lunibases.push_back( hbv.Edges()[i] );
	  lunibases.push_back( hbv.Edges()[i] );
	  lunibases.back().ReverseComplement();
	}
	lunibases.UniqueSort();
      }

      vec<int> lulens(lunibases.size() );
      size_t maxLUniLen = 0;
      for ( size_t u = 0; u < lunibases.size(); u++ ){
	lulens[u] = lunibases[u].isize();
	if ( lunibases[u].size() > maxLUniLen )
	  maxLUniLen = lunibases[u].size();
      }
      PRINT2( maxLUniLen, N50(lulens) );
      
      uni_fastb_file = run_dir + "/" + UNIBASES + ".flattened";
      lunibases.WriteAll( uni_fastb_file );
    }


    String strand_uni_fastb_file = uni_fastb_file + ".oneway.tmp2";
    out_lookup_file = strand_uni_fastb_file + ".lookuptab";
    if ( ! IsRegularFile( uni_fastb_file ) || NEW_LOOKUP_TABLE ) {
      cout << Date() << ": loading unibases" << flush;
      vecbasevector unibases(uni_fastb_file);
      cout << ": " << unibases.size() << " loaded" << endl;      
      vec<int> toRc;
      cout << Date() << ": removing rc copies" << endl;
      UnibaseInvolution( unibases, toRc );
      
      cout << Date() << ": building unibase reference file" << endl;
      uint64_t ulensSum = 0; int ulongest = 0;
      vec<int> suids( unibases.size(), vec<int>::IDENTITY );
      vec<int> ulens( unibases.size() );
      vec<int> ustatus( unibases.size(), 1 );
      for ( size_t i = 0; i < unibases.size(); i++ ){
	ulens[i] = unibases[i].size();
	ulensSum += ulens[i];
	if ( ulongest < ulens[i] )
	  ulongest = ulens[i];
	if ( ustatus[i] == 1 )
	  ustatus[ toRc[i] ] = -1; //don't use rc, disregard palindromes
      }
      vec<int> sulens( ulens );
      
      ReverseSortSync( sulens, suids ); 
      
      uint64_t approxLength = (ulensSum - unibases.size() * UNIBASES_K) / 2;
      uint64_t targetLength = round(approxLength * TARGET_UNIBASE_COVERAGE);
      
      cout << Date() << ": total unibase length " << approxLength
	   << ", targetting " << (ToString(100*TARGET_UNIBASE_COVERAGE))
	   << "% (" << targetLength << ")" << endl;
      
      uint64_t targetLensSum = 0;
      vecbasevector strand_unibases;
      for ( size_t i = 0; i < suids.size(); i++ ){
	size_t uid = suids[i];
	if ( ustatus[uid] > 0 && unibases[uid].size() > (max(MIN_UNIBASE_LENGTH, 2 * UNIBASES_K))) {
	  strand_unibases.push_back( unibases[uid] );
	  targetLensSum += ulens[uid];
	  if ( targetLensSum >= targetLength)
	    break;
	}
      }
      cout << Date() << ": selected " << strand_unibases.size() << " unibases of size range "
	   << strand_unibases.back().size() << "-" << strand_unibases[0].size()
	   << " covering " << targetLensSum << " bases" << endl;
      
      strand_unibases.WriteAll( strand_uni_fastb_file );
    }
    if ( ! IsRegularFile( out_lookup_file ) || NEW_LOOKUP_TABLE ){
      cout << Date() << ": building lookup table" << endl;
      LookupTabBuilder lookup( LookupK );
      lookup.addFastb( strand_uni_fastb_file.c_str( ) );
      cout << Date() << ": writing lookup table" << endl;
      lookup.write( out_lookup_file.c_str( ) );
    }
    cout << Date() << ": reading reference=" << strand_uni_fastb_file << endl;
    ref.ReadAll( strand_uni_fastb_file );
  } // -----------------------------------------------------


  // find longest contig/unipath in the reference
  int maxRefTigLen = 0;
  for ( size_t u = 0; u < ref.size(); u++ )
    maxRefTigLen = ref[u].isize() > maxRefTigLen ? ref[u].isize() : maxRefTigLen;
  PRINT( maxRefTigLen );

  uint64_t refLength = 0;
  vec<int64_t> srefLens( ref.size(), 0 );
  for (size_t i = 0; i < ref.size(); i++){
    srefLens[i] = ref[i].size();
    refLength += srefLens[i];
  }
  ReverseSort( srefLens );
  
  vec<double>  lenDataWeights( vec<double>(srefLens.at(0) +1, 1.0 ) );
  for ( size_t sri = 0; sri < srefLens.size(); sri++ ){    
    for ( longlong is = 1; is <= srefLens.at(0); is++ ){
      longlong ldiff = srefLens[sri] - is + 1;
      if ( ldiff > 0 )
	lenDataWeights.at(is) += ldiff;
      else
	break;
    }
  }
  
  for ( longlong is = 1; is <= srefLens.at(0); is++ )
    lenDataWeights.at(is) /= (double)refLength;
  
  


  cout << Date( ) << ": loading lookup table " << out_lookup_file << endl;
  LookupTab lookup_tab( out_lookup_file.c_str( ) );

  size_t nreads = MastervecFileObjectCount( reads_fastb );
   
  cout << Date() << ": loading pairs" << endl;
  PairsManager pairs;
  if ( IsRegularFile( pairs_file ) )
    pairs.Read( pairs_file );
  else if ( IsRegularFile( pairto_file ) )
    pairs.ReadFromPairtoFile( pairto_file, nreads );
  else
    FatalErr(" ERROR: did not find pairs file ");

  uint64_t npairs   = pairs.nPairs();
  size_t nLibraries = pairs.nLibraries();
  cout << Date( ) << ": found " << npairs << " pairs in " << nLibraries << " libraries" << endl;
  vec<String> titles;
  titles.push_back("libNo", "readLen", "origSep", "origDev");
  titles.push_back("newInSep", "newInDev", "3sigma%", "inStdErr");
  titles.push_back("newOutSep", "newOutDev", "outStdErr");
  titles.push_back("nLibraryPairs", "nPairsSampled", "nRead1Algnd", "nRead2Algnd");
  titles.push_back("nPairsAlgnd", "nPairsUniq", "nPairsSameContig", "nPairsLogical", "nPairsInnies", "nPairsOuties");
  titles.push_back("nSampleInnies", "nSampleOuties");
  vec< StdMap<String,longlong> > libStats( nLibraries ); // for storing results
  for ( size_t lib = 0; lib < nLibraries; lib++ ){
    libStats[lib]["readLen"] = 0;
    for ( size_t ti = 0; ti < titles.size(); ti++ )
      libStats[ lib ][ titles[ti] ] = 0;
  }

  vec<String> libID2libName( nLibraries );
  vec<longlong> libID2trgtSmplSize( nLibraries, MAX_TARGET_SAMPLE_SIZE );
  for (size_t lib = 0; lib < nLibraries; lib++ ){
    longlong origSep = pairs.getLibrarySep( lib );
    longlong origDev = pairs.getLibrarySD( lib );
    if ( origSep + 5 * origDev > MAX_SEPARATION )
      MAX_SEPARATION = origSep + 5 * origDev;
    libStats[ lib ]["origSep"] = origSep;
    libStats[ lib ]["origDev"] = origDev;
    libID2libName[ lib ]       = pairs.getLibraryName( lib );
    libStats[ lib ]["libNo"]   = lib;
    // require sample size large enough to reduce standard error below 0.5 
    // NOTE1: this assumes that original estimate of std. dev. is not too bad
    longlong trgtSmplSize = (longlong) pow( (double)origDev/0.45, 2.0 );
    if ( trgtSmplSize < MIN_TARGET_SAMPLE_SIZE)
      libID2trgtSmplSize[lib] = MIN_TARGET_SAMPLE_SIZE;
    else if ( trgtSmplSize > MAX_TARGET_SAMPLE_SIZE )
      libID2trgtSmplSize[lib] = MAX_TARGET_SAMPLE_SIZE;

    //else
    //cout << "WARNING: library " << lib << " named " << libID2libName[lib] << 
    //" will likely have standard error estimation bigger than 1." << endl;
  }
  cout << "target sample sizes: " << endl;
  libID2trgtSmplSize.Print( cout ); cout << endl;
  

  cout << Date() << ": computing pairs in each library" << endl;
  vec< vec<size_t> > lib_pairs(nLibraries);
  for (size_t pi = 0; pi < npairs; pi++ ) {
    libStats[pairs.libraryID(pi)]["nLibraryPairs"]++;
    lib_pairs[pairs.libraryID(pi)].push_back(pi);
  }
  
  srand(RANDOM_SEED);
  for (size_t lib = 0; lib < nLibraries; lib++ ){
    random_shuffle(lib_pairs[lib].begin(), lib_pairs[lib].end());
  }

  cout << Date() << ": reading reads from " << reads_fastb << endl;
  vecbasevector reads( reads_fastb );
  
  VecQualNibbleVec quals;
  if ( USE_QUALS ){
    cout << Date() << ": reading quals from " << quals_file << endl;
    LoadQualNibbleVec(quals_file, &quals);
    ForceAssertEq(nreads, quals.size());
  }



  cout << Date() << " Computing separations for libraries:" << endl;
  vec< vec< float> >   pairSepDataInnies( nLibraries );
  vec< vec< float> >   pairSepDataOuties( nLibraries );
  vec< vec< float> >   readLengthData( nLibraries );
  
  vec< Bool > found_data( nLibraries, False );
  
#pragma omp parallel for
  for (size_t lib = 0; lib < nLibraries; lib++ ) {
    // set alignment filtering parameters
    if ( libStats[lib]["origSep"] >  maxRefTigLen )
      continue;
    FirstLookupFilterECJ lfilter;
    lfilter.orientation           = FirstLookupFilterECJ::ALL;
    lfilter.max_kmer_freq         = 10000;
    lfilter.max_extend            = lfilter.max_kmer_freq;
    lfilter.max_placements        = 2; // we only care if there are multiple plausible placements
    if ( PERFECT ){
      lfilter.score_delta           = 0;
      lfilter.score_max             = -1;
      lfilter.mismatch_threshhold   = 1000000;
      lfilter.mismatch_neighborhood = 0;
      lfilter.mismatch_backoff      = 0;
      lfilter.max_error_rate        = 0;
    }else{
      lfilter.score_delta           = 20;
      lfilter.score_max             = 100;
      lfilter.mismatch_threshhold   = 3;
      lfilter.mismatch_neighborhood = 8;
      lfilter.mismatch_backoff      = 3;
      lfilter.max_placements        = 2; // we only care if there are multiple plausible placements
    }
    FirstLookupFinderECJ lfinder(lfilter, lookup_tab, ref, UNIBASES_K);
    for ( size_t  mi= 0;  mi < lib_pairs[lib].size() && !found_data[lib]; mi++ ) {
        size_t pi = lib_pairs[lib][mi];
      ForceAssertEq(lib, (size_t)pairs.libraryID(pi));
    
      if ( found_data[ lib ] ) continue;

      if ( pairSepDataInnies[ lib ].isize() >= libID2trgtSmplSize[lib] ){
	double outiePercLoc = (double)libStats[lib]["nPairsOuties"] / 
	  ( double)(libStats[lib]["nPairsOuties"] + 
		    libStats[lib]["nPairsInnies"] ) * 100.0;
	
	found_data[ lib ] = True;
	if (VERBOSE) {
	  cout << "\nfound enough data for "; PRINT2( lib, libID2libName[lib] );
	  PRINT3(  pairSepDataInnies[ lib ].isize(), 
		   pairSepDataOuties[ lib ].isize(), outiePercLoc );
	} else cout << " " << lib << flush;
      }
      
      libStats[ lib ]["nPairsSampled"]++;
      
      if ( PERFECT )
	lfilter.min_size = Min( reads[ pairs.ID1(pi) ].size(), reads[ pairs.ID2(pi) ].size() );
      Bool Read1Algnd = False, Read2Algnd = False, PairAlgnd = False, PairUnique = False,
	PairSameContig = False, PairLogical = False, PairInnie = False, PairOutie = False;
      int64_t sep = MAX_SEPARATION + 1;  
      pair_alignment_data( pairs, reads, quals, lfinder, lfilter,
			   pi, LookupK, USE_QUALS,
			   Read1Algnd, Read2Algnd, PairAlgnd, PairUnique, 
			   PairSameContig, PairLogical, PairInnie, PairOutie,
			   FLIP, TRIM, sep );

      if ( Read1Algnd )     libStats[ lib ]["nRead1Algnd"]++;
      if ( Read2Algnd )     libStats[ lib ]["nRead2Algnd"]++;
      if ( PairAlgnd )      libStats[ lib ]["nPairsAlgnd"]++;
      if ( PairUnique )     libStats[ lib ]["nPairsUniq"]++;
      if ( PairSameContig ) libStats[ lib ]["nPairsSameContig"]++;
      if ( PairLogical )    libStats[ lib ]["nPairsLogical"]++;
      if ( PairInnie )      libStats[ lib ]["nPairsInnies"]++;
      if ( PairOutie )      libStats[ lib ]["nPairsOuties"]++;

      
      if ( abs(sep) <= MAX_SEPARATION ){
	if ( PairInnie )
	  pairSepDataInnies[lib].push_back( sep );
	else if ( PairOutie )
	  pairSepDataOuties[lib].push_back( sep );
	
	float rlen = 
	  (float)(reads[ pairs.ID1(pi) ].size() + reads[ pairs.ID2(pi) ].size())/2.0;
	readLengthData[lib].push_back( rlen );
      }

    }
  }

  cout << endl;

  // compute library stats
  cout << Date() << " computing stats from data" << endl;
  for ( size_t lib = 0; lib < nLibraries; lib++ ){
    if ( VERBOSE )
      PRINT2( lib, libID2libName[lib] );
    if (  libStats[ lib ]["nPairsInnies"] < MIN_SAMPLE_SIZE ){
      cout << "WARNING: library " << lib << " named " << libID2libName[lib] << 
	": unable to compute stats; sample too small (" << 
	pairSepDataInnies[lib].isize() << " < " << MIN_SAMPLE_SIZE << 
	"), preserving original sep = " << libStats[lib].at("origSep") << ", dev = " <<
	libStats[lib].at("origDev") << endl;
      
      // return original values if computation not possible
      libStats[ lib ]["newInSep"]  = libStats[ lib ]["origSep"]; 
      libStats[ lib ]["newInDev"]  = libStats[ lib ]["origDev"];
      libStats[ lib ]["newOutSep"] = 0; 
      libStats[ lib ]["newOutDev"] = 0;
      continue;
    }
    
    
    libStats[ lib ]["nSampleInnies"] = pairSepDataInnies[lib].size();
    libStats[ lib ]["nSampleOuties"] = pairSepDataOuties[lib].size();
    libStats[ lib ]["readLen"] = Mean( readLengthData[lib] );
    

    vec<float> dataInnies = pairSepDataInnies[lib];
    Sort( dataInnies );
    //size_t nremoved = filter_outliers( dataInnies, 100.0, 100.0 );
    
    vec<double> histoInnies( dataInnies.back() + 1, 0.0 );
    for ( size_t i = 0; i < dataInnies.size(); i++  )
      histoInnies[ abs(dataInnies[i]) ]++;
    
    cout << Date() << " adding innie data to compensate for contig sizes" << endl;
    PRINT( dataInnies.size() );
    for ( size_t lsep = 0; lsep < histoInnies.size(); lsep++ ){
      int lins   = lsep + 2 * libStats[lib]["readLen"]; 
      int factor = floor( histoInnies[lsep] * ( 1.0/lenDataWeights.at(lins) - 1.0 ) );
      for ( int i = 1; i <= factor; i++ )
	dataInnies.push_back( lsep );
    }
    PRINT( dataInnies.size() );

    cout << Date() << " computing stats for lib = " << lib << endl;
    double meanInnies        = (double) libStats[lib].at("origSep");
    double devInnies         = (double) libStats[lib].at("origDev");
    double stdErrInnies      = 0.0;
    double max_weight_innies = gm(dataInnies, meanInnies, devInnies, stdErrInnies, VERBOSE );

    double sigma3 = FractionInDist(pairSepDataInnies[lib], meanInnies, devInnies, 3);
    libStats[ lib ]["3sigma%"]   = round(sigma3*100.0);

    // sanity check samples against mu,sigma
    if ( sigma3 < .80)
      cout << "WARNING: library " << lib << " named " << libID2libName[lib] << 
	": only " << libStats[lib]["3sigma%"] << "% of pairs are within 3 sigma of mean. Possible problem with data." << endl;
    
    libStats[ lib ]["newInSep"]   = round(meanInnies);
    libStats[ lib ]["newInDev"]   = round(devInnies + stdErrInnies);
    libStats[ lib ]["inStdErr"]   = round(stdErrInnies);

    // update pairs library information
    pairs.changeLibrarySepSd( lib, round(meanInnies), round(devInnies) );
    
    if ( pairSepDataOuties[lib].size() >= (unsigned) MIN_SAMPLE_SIZE ) {
      Sort( pairSepDataOuties[lib] );

      vec<float> dataOuties = pairSepDataOuties[lib];
      Sort( dataOuties );
      //size_t nremoved = filter_outliers( dataOuties, 100.0, 100.0 );
      
      vec<double> histoOuties( dataOuties.back() + 1, 0.0 );
      for ( size_t i = 0; i < dataOuties.size(); i++  )
	histoOuties[ abs(dataOuties[i]) ]++;
      
      cout << Date() << " adding outie data to compensate for contig sizes" << endl;
      PRINT( dataOuties.size() );
      for ( size_t lsep = 0; lsep < histoOuties.size(); lsep++ ){
	int lins   = lsep + 2 * libStats[lib]["readLen"]; 
	int factor = floor( histoOuties[lsep] * ( 1.0/lenDataWeights.at(lins) - 1.0 ) );
	for ( int i = 1; i <= factor; i++ )
	  dataOuties.push_back( lsep );
      }
      PRINT( dataOuties.size() );


      double meanOuties        = (double) dataOuties.front() - 1.0;
      double devOuties         = 1.0;
      double stdErrOuties      = 0.0;
      double max_weight_outies = gm(dataOuties, meanOuties, devOuties, stdErrOuties, VERBOSE );
	
      double sigma3 = FractionInDist(pairSepDataOuties[lib], meanOuties, devOuties, 3);
      // sanity check samples against mu,sigma
      if ( sigma3 < .80)
	cout << "WARNING: library " << lib << " named " << libID2libName[lib] << 
	  ": only " << round(sigma3*100) << "% of outie pairs are within 3 sigma of mean. Possible problem with data." << endl;

      libStats[ lib ]["newOutSep"]    = round(meanOuties);
      libStats[ lib ]["newOutDev"]    = round(devOuties + stdErrOuties);
      libStats[ lib ]["outStdErr"]    = round(stdErrOuties);
    } else {
      libStats[ lib ]["outSkewness"]   = 0;
      libStats[ lib ]["outKurtosis"]   = 0;
    }

    
    if ( abs( libStats[lib].at("newInSep") - libStats[lib].at("origSep") ) > 
	 2.0 * libStats[lib].at("origDev") )
      cout << "WARNING: library " << lib << " named " << libID2libName[lib] << 
	": computed separation (" << libStats[lib].at("newInSep") << " +/- " << 
	libStats[lib].at("newInDev") << ") significantly different from the original one (" << libStats[lib].at("origSep") << " +/- 2 * " << 
	libStats[lib].at("origDev") << ")" << endl;

  }

 
  
  vec<double> refInWeights( nLibraries, 0.0 );
  vec<double> refOutWeights( nLibraries, 0.0 );
  for ( size_t lib = 0; lib < libStats.size(); lib++ ){
    for ( size_t sri = 0; sri < srefLens.size(); sri++ ){
      if ( srefLens[sri] >= libStats[lib]["newInSep"] + 2 * libStats[lib]["readLen"])
	refInWeights[lib] += 
	  srefLens[sri]  - libStats[lib]["newInSep"] - 2 * libStats[lib]["readLen"];
      if ( srefLens[sri] >= libStats[lib]["newOutSep"] + 2 * libStats[lib]["readLen"])
	refOutWeights[lib] += 
	  srefLens[sri] - libStats[lib]["newOutSep"] - 2 * libStats[lib]["readLen"];
    }
    refInWeights[lib] /= (double)refLength;
    refOutWeights[lib] /= (double)refLength;
  }
  if ( VERBOSE ){
    cout << "reference weights for innies:" << endl;
    refInWeights.Print(cout); cout << endl;
    cout << "reference weights for outies:" << endl;
    refOutWeights.Print(cout); cout << endl;
  }

  for ( size_t lib = 0; lib < libStats.size(); lib++ ){
    libStats[ lib ]["percInnies"]   = 0;
    libStats[ lib ]["percOuties"]   = 0;
    if ( refInWeights[lib] > 0.0 && refOutWeights[lib] > 0.0 ){
      double sumBoth                       = 
	(double)libStats[ lib ].at("nPairsInnies")/refInWeights[lib] + 
	(double)libStats[ lib ].at("nPairsOuties")/refOutWeights[lib];
      if ( sumBoth > 0.0 ){
	libStats[ lib ]["percInnies"]   =  
	  round( (double)libStats[ lib ].at("nPairsInnies")/refInWeights[lib]/sumBoth * 100.0 );
	libStats[ lib ]["percOuties"]   =  
	  round( (double)libStats[ lib ].at("nPairsOuties")/refOutWeights[lib]/sumBoth * 100.0 );
      }
    }

    if (libStats[ lib ]["percOuties"] < 1 ){
      // very small sample size, likely there are no real outies present
      libStats[ lib ]["newOutSep"] = libStats[ lib ]["newOutDev"] = 0;
    }

    longlong nPairsInniesWghtd = refInWeights[lib] > 0.0 ? libStats[ lib ].at("nPairsInnies")/refInWeights[lib] : libStats[ lib ].at("nPairsInnies");
    longlong nPairsOutiesWghtd = refOutWeights[lib] > 0.0 ? libStats[ lib ].at("nPairsOuties")/refOutWeights[lib] : libStats[ lib ].at("nPairsOuties");
    if ( nPairsInniesWghtd < nPairsOutiesWghtd )
      cout << "WARNING: library " << lib << " named " << libID2libName[lib] << 
	": found more outies (" << nPairsOutiesWghtd << ") than innies (" 
	   << nPairsInniesWghtd << "); possible problem with reads  orientation." 
	   << endl;
    
  }


  if ( ! WRITE_NOTHING ){
    // Write .outies output file
    ofstream dotout( outies_file.c_str() );
    for ( size_t lib = 0; lib < nLibraries; lib++ )
      dotout << libStats[ lib ].at("newOutSep") << " " << libStats[lib].at("newOutDev") << " " << libStats[ lib ].at("percOuties") << endl;
    
    
    // WRITING NEW PAIRS FILE
    if ( ! OUT_HEAD.empty() ){
      cout << Date() << ": writing pairs" << endl;
      pairs.Write( run_dir + "/" + OUT_HEAD + ".pairs" );
    } 
  }
  
  cout << Date() << ": preparing stats output" << endl;
  // compute the character length of longest library name for formatting purposes
  string longestLibName = "libraryName";
  size_t longestLibNameSize = longestLibName.size();
  for ( size_t lib = 0; lib < nLibraries; lib++ )
    if ( libID2libName[lib].size() > longestLibNameSize )
      longestLibNameSize = libID2libName[lib].size();

  cout << Date() << ": preparing short stats output" << endl;
  /// ----------------- BUILDING SHORT STATS OUTPUT ----------------------------------
  stringstream sout;
  vec<String> short_titles; 
  short_titles.push_back("libNo", "readLen", "origSep", "origDev");
  short_titles.push_back("newInSep", "newInDev");
  short_titles.push_back("3sigma%", "inStdErr");
  short_titles.push_back("newOutSep", "newOutDev", "outStdErr");
  short_titles.push_back("percInnies", "percOuties");
  sout << "                              ------ SHORT STATS ------\n";
  String libNameTag = "libraryName";
  sout.width( longestLibNameSize ); sout.fill(' '); 
  sout << libNameTag; sout << " ";
  short_titles.Print( sout ); sout << endl;
  for ( size_t lib = 0; lib < nLibraries; lib++ ){
    sout.width( libNameTag.size() ); sout.fill(' '); sout << libID2libName[ lib ]; sout << " ";
    for (unsigned it = 0; it < short_titles.size(); it++ ){
      String value = 
      //libStats[lib].at( short_titles.at(it) ) != -1 ? 
      //ToString( libStats[lib].at( short_titles.at(it) ) ) : "NA";
	ToString( libStats[lib].at( short_titles.at(it) ) );
      sout.width( short_titles[it].size() ); sout.fill(' '); 
      sout << value << " ";
    }
    sout << endl;
    
  }
  sout << "--------------------------------------------------------------\n\n\n";
  
  string shold;
  if (!VERBOSE) {
    shold = sout.str();
    cout << shold;
  }

  cout << Date() << ": preparing long stats output" << endl;
  /// ----------------- BUILDING LONG STATS OUTPUT----------------------------------  
  sout << "                         ------ MORE DETAILED STATS ------\n";
  sout.width( longestLibNameSize );  sout.fill(' '); sout << libNameTag; sout << " ";
  titles.Print( sout ); sout << endl;
  for ( size_t lib = 0; lib < nLibraries; lib++ ){
    sout.width( libNameTag.size() ); sout.fill(' '); sout << libID2libName[ lib ];  sout << " ";
    for (unsigned it = 0; it < titles.size(); it++ ){
      String value = 
	//libStats[lib].at( titles[it] ) != -1 ? ToString( libStats[lib][ titles[it] ] ) : "NA";
	ToString( libStats[lib][ titles[it] ] );
      sout.width( titles[it].size() ); sout.fill(' '); 
      sout << value << " ";
    }
    sout << endl;
  }
  sout << "--------------------------------------------------------------\n\n";

  // PRINTING AND WRITING STATS
  cout << Date() << ": writing stats" << endl;
  shold = sout.str();
  if (VERBOSE)
    cout << shold;
  if ( ! WRITE_NOTHING && WRITE_STATS ){ 
    ofstream rout( results_file.c_str() );
    rout << shold;
  }
  
  if ( ! WRITE_NOTHING && WRITE_SEPS) {
    PRINT4( full_reads_head, seps_file, READS_HEAD, OUT_HEAD );
    cout << Date() << ": writing seps" << endl;
    ofstream in_seps((seps_file+".in").c_str());
    ofstream out_seps((seps_file+".out").c_str());
    for ( size_t lib = 0; lib < nLibraries; lib++ ) {
      for (u_int i = 0; i < pairSepDataInnies[lib].size(); ++i)
	in_seps << libStats[lib]["libNo"] << " " << pairSepDataInnies[lib][i] << endl;
      for (u_int i = 0; i < pairSepDataOuties[lib].size(); ++i)
	out_seps << libStats[lib]["libNo"] << " " << pairSepDataOuties[lib][i] << endl;
    }
    in_seps.close();
    out_seps.close();
  }
      

  cout << Date() << ": Done!" << endl;

} // main() ends here.







//-------------------------------
void pair_alignment_data( const PairsManager& pairs, const vecbasevector& reads, const VecQualNibbleVec& quals,
			  FirstLookupFinderECJ& lfinder, FirstLookupFilterECJ& lfilter,
			  const longlong pi, const int LookupK, const Bool use_quals,
			  Bool& Read1Algnd, Bool& Read2Algnd, Bool& PairAlgnd, Bool& PairUnique, 
			  Bool& PairSameContig, Bool& PairLogical, Bool& PairInnie, Bool& PairOutie,
			  Bool FLIP, Bool TRIM, int64_t& sep ){


  ulonglong id1 = pairs.ID1(pi);
  bvec read1 = reads[ id1 ];
  if ( FLIP ) {
    read1.ReverseComplement();
  }
  if ( TRIM > 0 ) read1.SetToSubOf( read1, TRIM, -1 );
  
  list<first_look_align> hits1;
  lfilter.orientation = FirstLookupFilterECJ::ALL;
  //lfilter.debug_reads.resize(1);
  //lfilter.debug_reads[0] = id1;
  if ( use_quals ){
    QualNibbleVec qual1 = quals[id1];
    if ( FLIP ) qual1.ReverseMe();
    //if ( TRIM > 0 ) qual1.SetToSubOf( qual1, TRIM, -1 );
    lfinder.getAlignments( read1, qual1, id1, &hits1);
  }else{
    lfinder.getAlignments( read1, id1, &hits1);
  }
  
  
  ulonglong id2 = pairs.ID2(pi);
  bvec    read2 = reads[ id2 ];
  
  if ( FLIP ) {
    read2.ReverseComplement();
  }
  if ( TRIM > 0 )read2.SetToSubOf( read2, TRIM, -1 );
  
  list<first_look_align> hits2;
  lfilter.orientation = FirstLookupFilterECJ::ALL;
  //lfilter.debug_reads[0] = id1;
  if ( use_quals ){
    QualNibbleVec qual2 = quals[id2];
    if ( FLIP ) qual2.ReverseMe();
    //if ( TRIM > 0 ) qual2.SetToSubOf( qual2, TRIM, -1 );
    lfinder.getAlignments( read2, qual2, id2, &hits2);
  }else{
    lfinder.getAlignments( read2, id2, &hits2);
  }
  
  
  if ( hits1.size() > 0 ) Read1Algnd = True;
  if ( hits2.size() > 0 ) Read2Algnd = True;
  
  if ( hits1.size() > 0 && hits2.size() > 0 )
    PairAlgnd = True;
  
  if ( hits1.size() != 1 || hits2.size() != 1 )
    return;  // take only unique alignments
  
  PairUnique = True;
 
  first_look_align &fla1 = *hits1.begin();
  first_look_align &fla2 = *hits2.begin();
  
  if ( fla1.target_loc.getContig() != fla2.target_loc.getContig() ) 
    return;  // not the same contig
  
  PairSameContig = True;
  
  // orientations
  Bool fd1 = fla1.is_FW();
  Bool fd2 = fla2.is_FW();
  
  if ( (fd1 && fd2) || (!fd1 && !fd2) ) // require correct (logical) orientation
    return;
  PairLogical = True;
  
  // if we flipped for alinment, reverse sense of forwardness
  if (FLIP) {
    fd1 = !fd1;
    fd2 = !fd2;
  }

  // use convention [s1,e1] refer to forward read, [s2,e2] to reverse
  int64_t s1,e1,s2,e2;
  if (fd1) {
    s1 = fla1.get_start_on_target(read1.isize(),LookupK);
    e1 = s1 + read1.isize();
    e2 = fla2.get_start_on_target(read2.isize(),LookupK);
    s2 = e2 + read2.isize();
  } else {
    e2 = fla1.get_start_on_target(read1.isize(),LookupK);
    s2 = e2 + read1.isize();
    s1 = fla2.get_start_on_target(read2.isize(),LookupK);
    e1 = s1 + read2.isize();
  }

  // positive->innie, negative->outie
  int64_t dist = s2 - s1;

  // fragment length
  int64_t length;

  // This is a bit of a hack which assumes we will be handling jump
  // libraries with FLIP=True, but not regular frag libraries.  Jump
  // library "innies" appear as outies here, because they are flipped
  // earlier in the pipeline.  However, those fragments smaller than the
  // read length cause confusion. --bruce
  if (FLIP) {
    int64_t readlengths = read1.isize() + read2.isize();
    dist -= readlengths;
    length = (dist > 0) ? (dist+readlengths) : (e1 - e2);
  } else {
    length = (dist > 0) ? dist : (e1 - e2);
  }

#ifdef notdef
  // show bases of reads beyond end of fragment
  if (length < read1.isize()) {
    bvec b;
    b.SetToSubOf(read1,length,read1.size()-length);
    b.Print(cout);
    b.SetToSubOf(read2,length,read2.size()-length);
    b.Print(cout);
  }
#endif

  // we define separation as fragment length minus read portion,
  // regardless of orientation (per David). --bruce
  sep = length - (read1.isize() + read2.isize());

  if ( dist >= 0 )
    PairInnie = True;
  else
    PairOutie = True;
}




// ----------------------
double gm( vec<float> data, double& emean, double& esigma, double& stdErr, const Bool& verbose ){
  cout << "\n  --------- gm -------" << endl;
  Sort( data );
  size_t nremoved = filter_outliers( data, 10.0, 10.0 );
  size_t N = data.size();

  int n10 = round( data.size() / 10.0 );
  int n1 = round( data.size() * 0.20 );
  int n2 = round( data.size() * 0.70 );
  
  double interv = data.back() - data.front();
  int intervSt = data.front();
  for ( int i = n1; i <= n2; i++ )
    if ( data.at( i + n10 ) - data.at( i ) < interv ){
      interv = data[i + n10] - data[i];
      intervSt = data[i];
    }
  double sigmaI = interv / 2.0 / 0.126;
  PRINT6( n1, n2, n10, interv, intervSt, sigmaI);

  //for ( double i = 0.12; i <= 0.13; i+=0.001){
  //double erfv = erf( (double)i/sqrt(2) );
  // PRINT2( i, erfv );
  //}

  double sigma = 1.0;
  ForceAssertGt( sigma, 0 );
  size_t K = 20;
  vec<double> means(K,0.0), sigmas(K, sigma), weights(K, 1.0/(K+1.0) );
  
  means[0] = data.front(); means[K -1] = data.back();
  
  for ( size_t i = 1; i < K -1; i++ ){
    int j = i * (N - 1.0) / (K -1.0);
    means[i] = data[ j ];
  }

  means.push_back( Median( data ) );
  sigmas.push_back( sigmaI );
  weights.push_back( 1.0 / (K +1.0) );

  if (verbose){
    cout << "meansOrig:   "; means.Print(cout); cout << endl;
    cout << "sigmasOrig:  "; sigmas.Print(cout); cout << endl;
    cout << "weightsOrig: "; weights.Print(cout); cout << endl;
  }
  
  actual_gm( data, means, sigmas, weights, verbose );    
  int n_data = data.size();


  int d1 = floor(data.front()), d2 = ceil(data.back());
  map<int,long double> measurep;
  map<int,long double> errorsf;
  for ( size_t i = 0; i < means.size(); i++ ){
    long double lsigma = sigmas[i];
    long double lmean = means[i];
    for ( int l = d1; l <= d2; l++ ){
      measurep[l] 
	+= (double)data.size() * weights[i] / (sqrt(2 * M_PI ) * lsigma ) * exp( - pow((l - lmean),2.0)/(2 * pow(lsigma,2.0) ) );
    }
  }

  

  map<int,long double> measurec;
  int cd1 = round( d1 + 0.5 * sigmaI);
  int cd2 = round(d2 - 0.5 * sigmaI );
  for ( int l = cd1; l <= cd2; l++ ){
    double v = 0.0;
    for ( int j = round( l - 0.5 * sigmaI); j <= round(l + 0.5 * sigmaI); j++ )
      measurec[l] += measurep[j];
  } 
      

  

  long double max = 0.0;
  int v1 = d1, v2 = d2, v0 = v1 + (v2 -v1)/2;
  for ( int l = d1; l < d2; l++)
    if ( measurec[l] > max ){
      max = measurec[l];
      v0 = l;
    }
    
  double thresh = 0.001 * (double)data.size();
  for ( int l = v0; l >= cd1; l-- )
    if (measurec[l] < thresh ){
      v1 = l;
      break;
    }
  
  for ( int l = v0; l <= cd2; l++ )
    if (measurec[l] < thresh ){
      v2 = l;
      break;
    }

  PRINT3( v1, v0, v2 );

  
  int v1l = v1, v2l = v2;
  for ( int l = v0; l > cd1; l-- )
    if (measurec[l-1] > measurec[l] ){
      v1l = l;
      break;
    }
  
  for ( int l = v0; l < cd2; l++ )
    if (measurec[l+1] > measurec[l] ){
      v2l = l;
      break;
    }
  
  PRINT5( v0, v1, v2, v1l, v2l );
  PRINT5( data[v0], data[v1], data[v2], data[v1l], data[v2l] );

  vec<float> cdata;
  cdata.reserve( data.size() );
  
  for ( size_t i = 0; i < data.size(); i++ )
    if ( data[i] >= v1 && data[i] <= v2 )
      cdata.push_back( data[i] );
  
  data = cdata; 
  cdata.clear();
  
  
  NormalDistribution distr0 = SafeMeanStdev( data );
 
  emean  = distr0.mu_;
  esigma = distr0.sigma_;

  stdErr = esigma/sqrt( (double)data.size() );
  PRINT2( emean, esigma );
  cout << "  ----------------- gm done ---------------------\n\n" << endl;
  return Max( weights );
}

void actual_gm( const vec<float>& data, vec<double>& means, vec<double>& sigmas, vec<double>& weights, const Bool& verbose ){
  size_t N = data.size();
  size_t K = means.size();
  double ratio = 9999;
  size_t nsteps = 0;
  double oldloglike = 1.0, loglike = 1.0;
  while( nsteps < 2  || abs(1.0 - ratio) > 0.001 && nsteps < 100){
    nsteps++;
    oldloglike = loglike; 
    vec< vec<double> > ddat( N, vec<double> (K, 0.0 ) );
    for ( size_t k = 0; k < K; k++ )
      for ( size_t n = 0; n < N; n++ )
	if ( weights[k] > 0.0 )
	  ddat[n][k] = -0.5 * pow( (data[n] - means[k]) / sigmas[k], 2.0 ) +
	    log( weights[k] )  + log( 1.0 / sqrt(2.0 * M_PI ) ) 
	    - log( sigmas[k] );
    
    
    loglike = 0.0; 
    for ( size_t n = 0; n < N; n++ ){
      double max = ddat[n][0];
      for ( size_t k = 0; k < K; k++ )
	if ( weights[k] > 0.0 )
	  if ( ddat[n][k] > max ) max = ddat[n][k];

      double sum = 0.0;
      for ( size_t k = 0; k < K; k++ )
	if ( weights[k] > 0.0 )
	  sum += exp( ddat[n][k] - max );
      
      double tmp = max + log( sum );
      for ( size_t k = 0; k < K; k++ )
	if ( weights[k] > 0.0 )
	  ddat[n][k] = exp( ddat[n][k] - tmp );
      
      loglike += tmp;
    }
    
    for ( size_t k = 0; k < K; k++ ){
      if ( weights[k] > 0.0 ){
	double wgt = 0.0;
	for ( size_t n = 0; n < N; n++ )
	  wgt += ddat[n][k];
	
	weights[k] = wgt / N; 
	
	double sum = 0.0;
	for ( size_t n = 0; n < N; n++ )
	  sum += ddat[n][k] * data[n];
	
	means[k] = sum / wgt;
	double sum2 = 0.0;
	for ( size_t n = 0; n < N; n++ )
	  sum2 += ddat[n][k] * pow( data[n] - means[k], 2.0 );
	// since we are dealing with integer seps, it's not improbable to get a peak
	// with only one value (sigma = 0)...this causes all sorts of trouble, so
	// putting floor of 1.0 on sigma. --bruce
	sigmas[k] = max(sqrt( sum2 / wgt ), 1.0);
      }else{
	means[k] = data.front() - 100000;
	sigmas[k] = 1e-99;
      }
    }
    
    ratio = abs( loglike / oldloglike );
  }
  if (verbose){
    PRINT2(ratio, nsteps);
    cout << "means:   "; means.Print(cout); cout << endl;
    cout << "sigmas:  "; sigmas.Print(cout); cout << endl;
    cout << "weights: "; weights.Print(cout); cout << endl;
  }
 
  if (verbose) cout << "\n\n";
}
