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

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

// auxiliary routines    -------------------------------
void pair_alignment_data( const PairsManager& pairs, const vecbasevector& reads, const VecQualNibbleVec &quals,
			  FirstLookupFinderECJ& lfinder, FirstLookupFilterECJ& lfilter,
			  const longlong pi, const int LookupK,
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
      PRINT2( distr1.mu_, distr1.sigma_ ); 
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
      NormalDistribution distr = SafeMeanStdev( data );
      mean  = distr.mu_;
      sigma = distr.sigma_;
      double zetB = abs( (data.front() - mean)/sigma );
      double zetE = abs( (data.back() - mean)/sigma );
      size_t i_max_zet = 0;
      if ( zetE >= zetB && zetE > zet_cut ){
	RemovedData = True;
	data.Erase( data.size() -1, data.size() );
      }else if ( zetB > zet_cut ){
	RemovedData = True;
	data.Erase( 0, 1 );
      }

    }
    
    {
      NormalDistribution distr2 = SafeMeanStdev( data );
      PRINT2( distr2.mu_, distr2.sigma_ ); 
    }
  }
    
  PRINT3( data.front(), data.back(), data.size() );

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
  CommandArgument_Int_OrDefault_Doc(MAX_TARGET_SAMPLE_SIZE, 1000000, 
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
  CommandArgument_Int_OrDefault_Doc(RANDOM_SEED, 133333,
    "Seed value for the random generator." );
  CommandArgument_Bool_OrDefault(VERBOSE, False);
  EndCommandArguments;
  
  ForceAssertGe( MAX_TARGET_SAMPLE_SIZE, MIN_TARGET_SAMPLE_SIZE );
  ForceAssertGe( MIN_TARGET_SAMPLE_SIZE, MIN_SAMPLE_SIZE );
  ForceAssertGt( MIN_TARGET_SAMPLE_SIZE, 0 );

  ForceAssert( !REF.empty() || !UNIBASES.empty() );
  if ( !UNIBASES.empty() )
    ForceAssertNe( UNIBASES_K, 0u );

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

  String results_file           = full_reads_head + ".sample_stats";
  String seps_file              = full_reads_head + ".sample_seps";
  String outies_file            = full_reads_head + ".outies";

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
	  ustatus[ toRc[i] ] = -1; //don't use rc, we will not worry about palindromes
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
    cout << Date() << ": reading reference" << endl;
    ref.ReadAll( strand_uni_fastb_file );
  }
  

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
    else if ( trgtSmplSize < MAX_TARGET_SAMPLE_SIZE )
      libID2trgtSmplSize[lib] = trgtSmplSize;
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

  cout << Date() << ": reading quals from " << quals_file << endl;
  VecQualNibbleVec quals;
  LoadQualNibbleVec(quals_file, &quals);
  ForceAssertEq(nreads, quals.size());

  // set alignment filtering parameters
  FirstLookupFilterECJ lfilter;
  lfilter.orientation = FirstLookupFilterECJ::ALL;
  lfilter.min_size = 20;
  lfilter.max_kmer_freq = 10000;
  lfilter.max_extend = lfilter.max_kmer_freq;
  lfilter.score_delta = 20;
  lfilter.score_max = 100;
  lfilter.min_match = 20;
  lfilter.mismatch_threshhold = 3;
  lfilter.mismatch_neighborhood = 8;
  lfilter.mismatch_backoff = 3;
  lfilter.max_placements = 2; // we only care if there are multiple plausible placements
  FirstLookupFinderECJ lfinder(lfilter, lookup_tab, ref, UNIBASES_K);

  cout << Date() << " Computing separations for libraries:" << endl;
  vec< vec< float> >   pairSepDataInnies( nLibraries );
  vec< vec< float> >   pairSepDataOuties( nLibraries );
  vec< vec< float> >   readLengthData( nLibraries );
  
  vec< Bool > found_data( nLibraries, False );
  
#pragma omp parallel for
  for (size_t lib = 0; lib < nLibraries; lib++ ) {
    for ( size_t  mi= 0;  mi < lib_pairs[lib].size() && !found_data[lib]; mi++ ) {
        size_t pi = lib_pairs[lib][mi];
      ForceAssertEq(lib, (size_t)pairs.libraryID(pi));
    
      //if ( mi % (perm.size()/100) == 0 )
      //Dot( cout, mi / (double)perm.size() * 100.0);
      //if ( mi == perm.size() -1 ) cout << endl;
    
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
      
      //lfilter.min_size = Min( reads[ pairs.ID1(pi) ].size(), reads[ pairs.ID2(pi) ].size() );
      Bool Read1Algnd = False, Read2Algnd = False, PairAlgnd = False, PairUnique = False,
	PairSameContig = False, PairLogical = False, PairInnie = False, PairOutie = False;
      int64_t sep = MAX_SEPARATION + 1;  
      pair_alignment_data( pairs, reads, quals, lfinder, lfilter,
			   pi, LookupK,
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
    size_t nremoved = filter_outliers( dataInnies, 10.0, 8.0 );
    
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
      size_t nremoved = filter_outliers( dataOuties, 10.0, 8.0 );
      
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

    if ( libStats[lib].at("nPairsInnies") < libStats[lib].at("nPairsOuties") )
      cout << "WARNING: library " << lib << " named " << libID2libName[lib] << 
	": found more outies (" << libStats[lib].at("nPairsOuties") << ") than innies (" 
	   << libStats[lib].at("nPairsInnies") << "); possible problem with reads  orientation." 
	   << endl;
    
  }


  // Write .outies output file
  ofstream dotout( outies_file.c_str() );
  for ( size_t lib = 0; lib < nLibraries; lib++ )
    dotout << libStats[ lib ].at("newOutSep") << " " << libStats[lib].at("newOutDev") << " " << libStats[ lib ].at("percOuties") << endl;

  
  // WRITING NEW PAIRS FILE
  if ( ! OUT_HEAD.empty() ){
    cout << Date() << ": writing pairs" << endl;
    pairs.Write( run_dir + "/" + OUT_HEAD + ".pairs" );
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
  if ( WRITE_STATS ){ 
    ofstream rout( results_file.c_str() );
    rout << shold;
  }
  
  if (WRITE_SEPS) {
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
			  const longlong pi, const int LookupK,
			  Bool& Read1Algnd, Bool& Read2Algnd, Bool& PairAlgnd, Bool& PairUnique, 
			  Bool& PairSameContig, Bool& PairLogical, Bool& PairInnie, Bool& PairOutie,
			  Bool FLIP, Bool TRIM, int64_t& sep ){


  ulonglong id1 = pairs.ID1(pi);
  bvec read1 = reads[ id1 ];
  QualNibbleVec qual1 = quals[id1];
  if ( FLIP ) {
    read1.ReverseComplement();
    qual1.ReverseMe();
  }
  if ( TRIM > 0 ) read1.SetToSubOf( read1, TRIM, -1 );
  
  list<first_look_align> hits1;
  lfilter.orientation = FirstLookupFilterECJ::ALL;
  //lfilter.debug_reads.resize(1);
  //lfilter.debug_reads[0] = id1;
  lfinder.getAlignments( read1, qual1, id1, &hits1);
  
  ulonglong id2 = pairs.ID2(pi);
  bvec    read2 = reads[ id2 ];
  QualNibbleVec qual2 = quals[id2];
  if ( FLIP ) {
    read2.ReverseComplement();
    qual2.ReverseMe();
  }
  if ( TRIM > 0 )read2.SetToSubOf( read2, TRIM, -1 );
  
  list<first_look_align> hits2;
  lfilter.orientation = FirstLookupFilterECJ::ALL;
  //lfilter.debug_reads[0] = id1;
  lfinder.getAlignments( read2, qual2, id2, &hits2);
 
  
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
  Sort( data );
  size_t nremoved = filter_outliers( data, 10.0, 7.0 );
  size_t N = data.size();

  double sigma = 1.0;
  ForceAssertGt( sigma, 0 );
  size_t K = 3;
  vec<double> means(K,0.0), sigmas(K, sigma), weights(K, 1.0/K );
  
  means[0] = data.front(); means[K -1 ] = data.back();
  
  for ( size_t i = 1; i < K -1; i++ ){
    int j = i * (N - 1.0) / (K -1.0);
    means[i] = data[ j ];
  }

  if (verbose){
    cout << "meansOrig:   "; means.Print(cout); cout << endl;
    cout << "sigmasOrig:  "; sigmas.Print(cout); cout << endl;
    cout << "weightsOrig: "; weights.Print(cout); cout << endl;
  }
  
  Bool RemovedGaussian = True;
  while( K >= 1 && RemovedGaussian ){
    RemovedGaussian = False;
    actual_gm( data, means, sigmas, weights, verbose );    
    vec<int> indices( K );
    for ( size_t i = 0; i < K; i++ )
      indices[i] = i;
    vec<double> weights_sorted( weights );
    ReverseSortSync( weights_sorted, indices );
    for ( size_t i = 0; i < indices.size() -1; i++ ){
      int c = indices[i];
      int n = indices[i+1];
      if ( abs( means[c] - means[n] ) <= 3.0 * sigmas[c] || 
	   abs( means[c] - means[n] ) < 3.0 * sigmas[n] ){
	K--;
	RemovedGaussian = True;
	means.Erase(n,n+1); sigmas.Erase(n,n+1); weights.Erase(n,n+1); indices.Erase(i,i+1);
	break;
      }
    }

  }
  
  vec<int> indices( K );
  for ( size_t i = 0; i < K; i++ )
    indices[i] = i;
  vec<double> weights_sorted( weights );
  cout << "weights: "; weights.Print(cout); cout << endl;
  ReverseSortSync( weights_sorted, indices );
  double mean0 = means[indices[0]], sigma0 = sigmas[indices[0]];
  double maxweight = weights[ indices[0] ];
  if (maxweight > 0) {
    emean  = mean0;
    esigma = sigma0;
  }
  if ( verbose ){
    PRINT4( K, emean, esigma, maxweight );
    cout << endl;
  }

  stdErr = esigma/sqrt( (double)data.size() );
  PRINT(stdErr);
  return maxweight;
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
    if (verbose){
      PRINT2(ratio, nsteps);
      cout << "means:   "; means.Print(cout); cout << endl;
      cout << "sigmas:  "; sigmas.Print(cout); cout << endl;
      cout << "weights: "; weights.Print(cout); cout << endl;
    }
  }
  if (verbose) cout << "\n\n";
}
