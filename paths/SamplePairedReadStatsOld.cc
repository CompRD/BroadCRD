///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/**
   Program: ComputeReadPairSeps

   Determines separation stats, orientation, etc. of paired reads.
   This is done by aligning a random sample of pairs to a reference.

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
#include "lookup/FirstLookupFinder.h"
#include "lookup/LookAlign.h"
#include "lookup/LookupTabBuilder.h"
#include "paths/UnibaseUtils.h"
#include "VecUtilities.h"
#include "Map.h"


// auxiliary routines    -------------------------------
void pair_alignment_data( const PairsManager& pairs, const vecbasevector& reads, 
			  FirstLookupFinder& lfinder, FirstLookupFilter& lfilter,
			  const longlong pi, const int LookupK,
			  Bool& Read1Algnd, Bool& Read2Algnd, Bool& PairAlgnd, Bool& PairUnique, 
			  Bool& PairSameContig, Bool& PairLogical, Bool& PairInnie, Bool& PairOutie,
			  Bool FLIP, Bool TRIM, int64_t& sep );
double gm( vec<float> data, double& mean, double& sigma, const Bool& verbose);
void actual_gm( const vec<float>& data, vec<double>& means, vec<double>& sigmas, vec<double>& weights, const Bool& verbose );
//----------------------------------------------------------

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
  CommandArgument_Int_OrDefault_Doc(UNIBASES_K,-1, 
    "needed by involution algorithm");
  CommandArgument_String_Doc(READS_HEAD, 
    "looks for RUN_DIR/READS_HEAD.{fastb, pairs/pairto}");
  CommandArgument_Bool_OrDefault_Doc( FLIP, True,
    "reverse-complement reads for alignment (not for actual stats reporiting).");
  CommandArgument_String_OrDefault_Doc(OUT_HEAD, "", 
    "Name for output files. No output files written if not specified");
  CommandArgument_Bool_OrDefault_Doc(WRITE_STATS, False, 
    "Write sample stats into READS_HEAD.sample_stats file.");
  CommandArgument_Int_OrDefault_Doc(ORIG_DEV_BRACKET, 5, 
    "Radius (in std dev) of data censoring around original (lab specified) mean paired read separation.");
  CommandArgument_Int_OrDefault_Doc(MEAN_JUMP_CONSTRUCT_SIZE, 300, 
    "mean size of fragments constructud from jumping libraries");
  CommandArgument_Int_OrDefault_Doc(MEAN_JUMP_CONSTRUCT_DEV, 60, 
    "size std. deviation of fragments constructud from jumping libraries");

  CommandArgument_Int_OrDefault_Doc(TRIM, 0, 
    "Trim that many bases from the beginning of a read (just before alignment).");
  CommandArgument_Int_OrDefault_Doc(TARGET_SAMPLE_SIZE, 10000, 
    "sample size");
  CommandArgument_Int_OrDefault_Doc(MIN_SAMPLE_SIZE, 100, 
    "minimum sample size");
  CommandArgument_Int_OrDefault_Doc(MAX_SEPARATION, 100000, 
    "maximum aligned pair separation considered");
  CommandArgument_Int_OrDefault_Doc(RANDOM_SEED, 133333,
    "Seed value for the random generator." );
  CommandArgument_Bool_OrDefault(VERBOSE, False);
  EndCommandArguments;
  
  ForceAssertGe( TARGET_SAMPLE_SIZE, MIN_SAMPLE_SIZE );
  
  ForceAssert( !REF.empty() || !UNIBASES.empty() );
  if ( !UNIBASES.empty() )
    ForceAssertNe( UNIBASES_K, -1 );

  int LookupK = 12;

  String data_dir = PRE + "/" + DATA;
  String run_dir = data_dir + "/" + RUN;
  
  String full_ref_head          = run_dir + "/" + REF;
  String ref_fastb_file         = full_ref_head + ".fastb";
  
  String full_reads_head        = run_dir + "/" + READS_HEAD;
  String reads_fastb            = full_reads_head + ".fastb";
  String pairs_file             = full_reads_head + ".pairs";
  String pairto_file            = full_reads_head + ".pairto";

  String results_file           = full_reads_head + ".sample_stats";
  String outies_file            = full_reads_head + ".outies";

  vecbasevector ref;
  String out_lookup_file; 
  if ( UNIBASES.empty() ){
    cout << "Using reference file" << endl;
    out_lookup_file = full_ref_head + ".lookuptab";
    if ( ! IsRegularFile( out_lookup_file ) ) {
      cout << Date( ) << ": creating and saving lookup table " << out_lookup_file << endl;
      LookupTabBuilder lookup( LookupK );
      lookup.addFastb( ref_fastb_file.c_str( ) );
      lookup.write( out_lookup_file.c_str( ) );
    }
    ref.ReadAll( ref_fastb_file );
  }else{
    cout << "Using unibases as reference" << endl;
    String uni_fastb_file = run_dir + "/" + UNIBASES;
    String strand_uni_fastb_file = uni_fastb_file + ".oneway.tmp";
    out_lookup_file = strand_uni_fastb_file + ".lookuptab";
    if ( ! IsRegularFile( out_lookup_file ) ) {
      vecbasevector unibases(uni_fastb_file);
      vec<int> toRc;
      UnibaseInvolution( unibases, toRc );
      
      vec<int> ustatus( unibases.size(), 1 );
      for ( size_t i = 0; i < unibases.size(); i++ ){
	if ( ustatus[i] == 1 )
	  ustatus[ toRc[i] ] = -1; //don't use rc, we will not worry about palindromes
      }
      vecbasevector strand_unibases;
      for ( size_t uid = 0; uid < unibases.size(); uid++ ){
	if ( ustatus[uid] > 0 )
	  strand_unibases.push_back( unibases[uid] );
      }    
      strand_unibases.WriteAll( strand_uni_fastb_file );
      LookupTabBuilder lookup( LookupK );
      lookup.addFastb( strand_uni_fastb_file.c_str( ) );
      lookup.write( out_lookup_file.c_str( ) );
    }
    ref.ReadAll( strand_uni_fastb_file );
  }
  
  cout << Date( ) << ": loading lookup table " << out_lookup_file << endl;
  LookupTab lookup_tab( out_lookup_file.c_str( ) );


  cout << Date() << ": reading reads from " << reads_fastb << endl;
  vecbasevector reads( reads_fastb );
  longlong nreads = reads.size();
  
  cout << Date() << " loading pairs" << endl;
  PairsManager pairs;
  if ( IsRegularFile( pairs_file ) )
    pairs.Read( pairs_file );
  else if ( IsRegularFile( pairto_file ) )
    pairs.ReadFromPairtoFile( pairto_file, nreads );
  else
    FatalErr(" ERROR: did not find pairs file ");

  size_t npairs     = pairs.nPairs();
  size_t nLibraries = pairs.nLibraries();
  PRINT2( npairs, nLibraries );

  vec<String> titles;
  titles.push_back("libNo", "origSep", "origDev");
  titles.push_back("newInSep", "newInDev");
  titles.push_back("newOutSep", "newOutDev");
  titles.push_back("nLibraryPairs", "nPairsSampled", "nRead1Algnd", "nRead2Algnd");
  titles.push_back("nPairsAlgnd", "nPairsUniq", "nPairsSameContig", "nPairsLogical", "nPairsInnies", "nPairsOuties");
  titles.push_back("nSampleInnies", "nSampleOuties");
  vec< StdMap<String,longlong> > libStats( nLibraries ); // for storing results
  for ( size_t lib = 0; lib < nLibraries; lib++ )
    for ( size_t ti = 0; ti < titles.size(); ti++ )
      libStats[ lib ][ titles[ti] ] = 0;


  vec<String> libID2libName( nLibraries );
  for (size_t lib = 0; lib < nLibraries; lib++ ){
    libStats[ lib ]["origSep"] = pairs.getLibrarySep( lib );
    libStats[ lib ]["origDev"] = pairs.getLibrarySD( lib );
    libID2libName[ lib ]       = pairs.getLibraryName( lib );
    libStats[ lib ]["libNo"]   = lib;
  }
  for (ulonglong pi = 0; pi < npairs; pi++ )
    libStats[ pairs.libraryID( pi ) ]["nLibraryPairs"]++;
 

  cout << Date() << ": randomizing data set" << endl;
  srand(RANDOM_SEED);
  vec<longlong> perm( npairs );
  for ( size_t i = 0; i < npairs; i++)
    perm[ i ] = (double) rand() / RAND_MAX * npairs;  // allow repetition, faster
  

  // set alignment filtering parameters
  FirstLookupFilter lfilter;
  lfilter.min_size         = 20;
  lfilter.max_kmer_freq    = 2000;
  lfilter.max_extend       = 2000;
  lfilter.max_error_rate   = 0.3;
  FirstLookupFinder lfinder( lfilter, lookup_tab, ref );
 

  cout << Date() << " Computing separations" << endl;
  vec< vec< float> >   pairSepDataInnies( nLibraries );
  vec< vec< float> >   pairSepDataOuties( nLibraries );
  
  for ( ulonglong  mi= 0;  mi < perm.size(); mi++ ){
    if ( mi % (perm.size()/100) == 0 )
      Dot( cout, mi / (double)perm.size() * 100.0);
    if ( mi == perm.size() -1 ) cout << endl;
    
    longlong pi  = perm[ mi ];
    int lib      = pairs.libraryID( pi );
    
    if ( pairSepDataInnies[ lib ].isize() >= TARGET_SAMPLE_SIZE && 
	 pairSepDataOuties[ lib ].isize() >= MIN_SAMPLE_SIZE )
      continue;    
   
    
    libStats[ lib ]["nPairsSampled"]++;
    
    int rlen1   = reads[ pairs.ID1(pi) ].size();
    int rlen2   = reads[ pairs.ID2(pi) ].size();
    int rlenMin = rlen1 < rlen2 ? rlen1 : rlen2;
    lfilter.min_size         = (int)round( 0.5 * rlenMin );
    Bool Read1Algnd = False, Read2Algnd = False, PairAlgnd = False, PairUnique = False,
      PairSameContig = False, PairLogical = False, PairInnie = False, PairOutie = False;
    int64_t sep = MAX_SEPARATION + 1;  
    pair_alignment_data( pairs, reads, lfinder, lfilter,
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

    if ( PairInnie && abs(sep) <= MAX_SEPARATION  )
      pairSepDataInnies[lib].push_back( sep );

    if ( PairOutie && abs(sep) <= MAX_SEPARATION )
      pairSepDataOuties[lib].push_back( sep );
    
  }


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
    
    if ( libStats[lib].at("nPairsInnies") < libStats[lib].at("nPairsOuties") )
      cout << "WARNING: library " << lib << " named " << libID2libName[lib] << 
	": found more outies (" << libStats[lib].at("nPairsOuties") << ") than innies (" 
	   << libStats[lib].at("nPairsInnies") << "); possible problem with reads  orientation." 
	   << endl;
    
    libStats[ lib ]["nSampleInnies"] = pairSepDataInnies[lib].size();

    double meanInnies        = (double) libStats[lib].at("origSep");
    double devInnies         = (double) libStats[lib].at("origDev");
    double max_weight_innies = gm(pairSepDataInnies[lib], meanInnies, devInnies, VERBOSE );
    

    if ( max_weight_innies < 0.5 )
      cout << "WARNING: library " << lib << " named " << libID2libName[lib] << 
	": main gaussian has small weight = " << max_weight_innies << 
	   ". Possible problem with data." << endl;

    libStats[ lib ]["newInSep"]   = round(meanInnies);
    libStats[ lib ]["newInDev"]   = round(devInnies);
 
    
    // update pairs library information
    pairs.changeLibrarySepSd( lib, round(meanInnies), round(devInnies) );
    
    if ( pairSepDataOuties[lib].size() >= (unsigned) MIN_SAMPLE_SIZE ){
      Sort( pairSepDataOuties[lib] );
      double meanOuties        = (double) pairSepDataOuties[lib].front() - 1.0;
      double devOuties         = 1.0;
      double max_weight_outies = gm(pairSepDataOuties[lib], meanOuties, devOuties, VERBOSE );
      
      libStats[ lib ]["newOutSep"]    = round(meanOuties);
      libStats[ lib ]["newOutDev"]        = round(devOuties);
    }
    
    
    if ( abs( libStats[lib].at("newInSep") - libStats[lib].at("origSep") ) > 
	 2.0 * libStats[lib].at("origDev") )
      cout << "WARNING: library " << lib << " named " << libID2libName[lib] << 
	": computed separation (" << libStats[lib].at("newInSep") << " +/- " << 
	libStats[lib].at("newInDev") << ") significantly different from the original one (" << libStats[lib].at("origSep") << " +/- 2 * " << 
	libStats[lib].at("origDev") << ")" << endl;
    
  }
  
  
  
  for ( size_t lib = 0; lib < libStats.size(); lib++ ){
    double sumBoth                       = 
      (double)libStats[ lib ].at("nPairsInnies") + 
      (double)libStats[ lib ].at("nPairsOuties");
    if ( sumBoth > 0.0 ){
      libStats[ lib ]["percInnies"]   =  
	round( (double)libStats[ lib ].at("nPairsInnies")/sumBoth * 100.0 );
      libStats[ lib ]["percOuties"]   =  
	round( (double)libStats[ lib ].at("nPairsOuties")/sumBoth * 100.0 );
    }else{
      libStats[ lib ]["percInnies"]   = 0;
      libStats[ lib ]["percOuties"]   = 0;
    }

    if (libStats[ lib ]["percOuties"] < 1 ){
      // very small sample size, likely there are no real outies present
      libStats[ lib ]["newOutSep"] = libStats[ lib ]["newOutDev"] = 0;
    }
    
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

  cout << Date() << ": preparing short stats output";
  /// ----------------- BUILDING SHORT STATS OUTPUT ----------------------------------
  stringstream sout;
  vec<String> short_titles; 
  short_titles.push_back("libNo", "origSep", "origDev");
  short_titles.push_back("newInSep", "newInDev", "newOutSep", "newOutDev");
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
	libStats[lib].at( short_titles.at(it) ) != -1 ? 
	ToString( libStats[lib].at( short_titles.at(it) ) ) : "NA";
      sout.width( short_titles[it].size() ); sout.fill(' '); 
      sout << value << " ";
    }
    sout << endl;
    
  }
  sout << "--------------------------------------------------------------\n\n\n";
  
  
  cout << Date() << ": preparing long stats output";
  /// ----------------- BUILDING LONG STATS OUTPUT----------------------------------  
  sout << "                         ------ MORE DETAILED STATS ------\n";
  sout.width( longestLibNameSize );  sout.fill(' '); sout << libNameTag; sout << " ";
  titles.Print( sout ); sout << endl;
  for ( size_t lib = 0; lib < nLibraries; lib++ ){
    sout.width( libNameTag.size() ); sout.fill(' '); sout << libID2libName[ lib ];  sout << " ";
    for (unsigned it = 0; it < titles.size(); it++ ){
      String value = 
	libStats[lib].at( titles[it] ) != -1 ? ToString( libStats[lib][ titles[it] ] ) : "NA";
      sout.width( titles[it].size() ); sout.fill(' '); 
      sout << value << " ";
    }
    sout << endl;
  }
  sout << "--------------------------------------------------------------\n\n";


  // PRINTING AND WRITING STATS
  cout << Date() << ": writing stats" << endl;
  string shold = sout.str();;
  cout << shold;
  if ( WRITE_STATS ){ 
    ofstream rout( results_file.c_str() );
    rout << shold;
  }
  

  cout << Date() << ": Done!" << endl;

} // main() ends here.







//-------------------------------
void pair_alignment_data( const PairsManager& pairs, const vecbasevector& reads, 
			  FirstLookupFinder& lfinder, FirstLookupFilter& lfilter,
			  const longlong pi, const int LookupK,
			  Bool& Read1Algnd, Bool& Read2Algnd, Bool& PairAlgnd, Bool& PairUnique, 
			  Bool& PairSameContig, Bool& PairLogical, Bool& PairInnie, Bool& PairOutie,
			  Bool FLIP, Bool TRIM, int64_t& sep ){

  ulonglong id1 = pairs.ID1(pi);
  bvec read1 = reads[ id1 ];
  if ( FLIP ) read1.ReverseComplement();
  if ( TRIM > 0 ) read1.SetToSubOf( read1, TRIM, -1 );
  
  list<first_look_align> hits1;
  lfilter.orientation = FirstLookupFilter::ALL;
  lfinder.getAlignments( read1, id1, hits1);
  
  ulonglong id2 = pairs.ID2(pi);
  bvec    read2 = reads[ id2 ];
  if ( FLIP ) read2.ReverseComplement();
  if ( TRIM > 0 )read2.SetToSubOf( read2, TRIM, -1 );
  
  list<first_look_align> hits2;
  lfilter.orientation = FirstLookupFilter::ALL;
  lfinder.getAlignments( read2, id2, hits2);
 
  
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
  
  Bool fd1 = True ? fla1.n_mismatches >= 0 : False;
  Bool fd2 = True ? fla2.n_mismatches >= 0 : False;
  
  // compute offsets
  longlong off1 = fd1 ? fla1.target_loc.getOffset() : Max( (int)(fla1.target_loc.getOffset() + LookupK - read1.size() ), 0 );
  
  longlong off2 = fd2 ? fla2.target_loc.getOffset() : Max( (int)(fla2.target_loc.getOffset() + LookupK - read2.size() ), 0 );
  
  
  if ( FLIP ){ // reads were flipped, reverse orientation information
    fd1 = ! fd1;
    fd2 = ! fd2;
  }
  
  if ( (fd1 && fd2) || (!fd1 && !fd2) ) // require correct (logical) orientation
    return;
  
  PairLogical = True;
  
  int dist = 0;
  if ( fd1 && ! fd2 ){
    dist = off2 - off1;
    if ( dist >= 0 ) sep = dist - read1.size();
    else             sep = dist + read2.size();
  }
  else if ( ! fd1 && fd2 ){
    dist = off1 - off2;
    if ( dist >= 0 ) sep = dist - read2.size();
    else             sep = dist + read1.size();
  }
  else{
    cout << "ERROR: unexpected placement of reads" << endl;
    ForceAssert( 0 == 1 );
  }
  
  if ( dist >= 0 )
    PairInnie = True;
  else
    PairOutie = True;
  
}
// ----------------------
double gm( vec<float> data, double& emean, double& esigma, const Bool& verbose ){
  Sort( data );
  size_t N = data.size();
  int i1 = round( N * 1.0/100.0 ), i2 = round( N * ( 1.0 - 1.0/100.0 ) );
  data.SetToSubOf( data, i1, i2 - i1 + 1 );

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
  while( K > 1 && RemovedGaussian ){
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
  ReverseSortSync( weights_sorted, indices );  
  emean  = means[ indices[0] ];
  esigma = sigmas[ indices[0] ];
  double maxweight = weights[ indices[0] ];
  if ( verbose ){
    PRINT4( K, emean, esigma, maxweight );
    cout << endl;
  }

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
	sigmas[k] = sqrt( sum2 / wgt );
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
  cout << "\n\n";
}
