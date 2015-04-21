///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/** 
   Program: GcBiasCurveLG

   Estimate and store the <GC bias> curve for a given K value.  This shows,
   for each possible GC content, the positive or negative bias for reading
   genome regions with that gc content: do such genome regions get relatively
   more reads than average, or relatively fewer reads than average?

   This is a modified version of the original script. It does not assume that 
   the bias computation is possible for each GC content. 

   The GC bias curve computed by this program can be used in the following ways:

      - the <SeqToPaths> program can use it to generate simulated reads with GC bias

      - the <UnibaseCopyNumber> program can use it to calculate unbibase copy nos

   INPUT FILES:
     reads.fastb

   OUTPUT FILES:
     reads.freq_table.k*
     reads.gc_bias.k*

   @file
*/

#include "MainTools.h"
#include "kmers/KmerShape.h"
#include "PairsManager.h"
#include "kmer_freq/KmerFrequencyTable.h"
#include "kmer_freq/WriteKmerFrequencies.h"




/**
   Function: smooth_histogram
   This is a local function that smooths the histogram between start and end positions
   by computing an average value in a section specified by radiusLeft and radiusRight
   (to allow nonsymmetric computations)

*/
void smooth_histogram( const vec<longlong>& histo, vec<double>& smoothHisto, 
		       int radiusLeft, int radiusRight, longlong start, longlong end ){
  ForceAssertGe(radiusLeft, 0);
  ForceAssertGe(radiusRight, 0);
  ForceAssertGe(start, 0);
  ForceAssertLe(end, histo.isize() -1);
  if ( histo.isize() != smoothHisto.isize() ){
    //cout << "WARNING: size of smoothed histogram not equal to the size of an original. Resizing\n";
    smoothHisto.resize( histo.isize() );
  }
  for (longlong i = start; i <= end; i++ )
    smoothHisto[i] = 0;
  for (longlong i = 0; i < start; i++)
    smoothHisto[i] = histo[i];
  for (longlong i = end + 1; i < histo.isize(); i++)
    smoothHisto[i] = histo[i];
  for (longlong i = start; i <= end; i++ ){
    longlong leftRange = (i - radiusLeft >= start ) ? i - radiusLeft : start;
    longlong rightRange = ( i + radiusRight <= end ) ? i + radiusRight : end;
    longlong fullRange = rightRange - leftRange + 1;
    longlong rangeSum = 0;
    for ( longlong j = leftRange; j <= rightRange; j++ )
      rangeSum += histo[j];
    double smoothVal = (double) rangeSum / (double) fullRange;
    smoothHisto[i] = smoothVal;
  }
  return;
} 


int main( int argc, char** argv ) {
  
  RunTime();
  
  BeginCommandArguments;
  CommandDoc("Compute GC Bias curve from reads.");

  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_KShape(K);
  
  CommandArgument_String_OrDefault(READS_IN,"reads");
  CommandArgument_Int_OrDefault(BEGIN, 0 );
  CommandArgument_Int_OrDefault(END, -1 );
  CommandArgument_UnsignedInt_OrDefault(TRUNC, 0 );
  CommandArgument_Double_OrDefault_Doc(TAIL_FRACTION, 0.40,"fraction of data in the histogram tail" );
  CommandArgument_String_OrDefault(BIAS_CURVE_OUT,"");
  
  CommandArgument_Bool_OrDefault(VERBOSE, False);
  
  EndCommandArguments;

  if ( BIAS_CURVE_OUT == "" )
    BIAS_CURVE_OUT = READS_IN + ".gc_bias";

  String dataDir = PRE + "/" + DATA;
  String runDir = PRE + "/" + DATA + "/" + RUN;
  String freqTableFile = runDir + "/" + READS_IN + ".freq_table.k" + ToString(K);

  if ( VERBOSE )
    cout << Date() << " reading reads " << endl;
  vecbasevector reads;
  if ( BEGIN == 0 && END == -1 )
    reads.ReadAll( runDir + "/" + READS_IN + ".fastb");
  else {
    ForceAssertLe( BEGIN, END );
    reads.SparseReadRange( runDir + "/" + READS_IN + ".fastb", BEGIN, END );
  }
 
  
  if ( TRUNC > 0 ){
    if ( VERBOSE ) 
      cout << Date() << " truncating reads" << endl;
    for ( size_t r = 0; r < reads.size(); r++ )
      if ( reads[r].size() > TRUNC )
        reads[r].resize( TRUNC );
  }
   
  // Compute kmer frequencies if they were not computed before
#define CASE(K) WriteKmerFrequencies<K>( reads, freqTableFile, true )
  if ( ! IsRegularFile( freqTableFile ) ){
    if ( VERBOSE )
      cout << "computing kmer frequency table" << endl;
    DISPATCH_ON_KSHAPE( K, CASE );
  }
  
  int Ksize = GetKmerSize( K );  // kmer size
    
  if ( VERBOSE )
    cout << Date() << " reading the frequency table file" << endl; 
  KmerFrequencyTable freqTable( K, freqTableFile );

  // gcKmerHist cotains numbers of kmers having the same frequency of occurance
  vec< vec<longlong> > gcKmerHist;
  freqTable.GetGCHistograms( gcKmerHist );

  // clean up the histogram by removing low count and orphaned frequency bins
  // Sometimes there are single kmer types occuring very many times (perhaps
  // an experimental artifact.
  // Lots of heruristics.
  int lowKmerCount = 2;  /// HEURISTIC!!!!!! this should be justified or derived somehow
  for ( unsigned gc = 0; gc < gcKmerHist.size(); gc++){
    longlong lastFreq = gcKmerHist[gc].isize() -1;
    
    // zoroeing low count bins
    for ( longlong freq = 1; freq <= lastFreq; freq++ ){
      if ( gcKmerHist[gc][freq] <= lowKmerCount ) 
	gcKmerHist[gc][freq] = 0;
    }
    
    // zeroeing isolated bins
    for ( longlong freq = 1; freq < lastFreq; freq++){
      if ( gcKmerHist[gc][ freq ] > 0 && gcKmerHist[gc][ freq -1 ] == 0 &&
	   gcKmerHist[gc][ freq +1 ] == 0 ){
	gcKmerHist[gc][ freq ] = 0;
      }
    } 
    if ( lastFreq > 2 )
      if ( gcKmerHist[gc][ lastFreq ] > 0 && gcKmerHist[gc][ lastFreq -1 ] == 0 )
	gcKmerHist[gc][ lastFreq ] = 0;
  }

  // clean up by getting rid of zeroed high frequency bin
  longlong maxFreq = 0;
  for ( int gc = 0; gc < gcKmerHist.isize(); gc++){
    for ( longlong freq = gcKmerHist[gc].isize() -1; freq >= 0; freq-- ){
      if ( gcKmerHist[gc][freq] > 0 ){
	gcKmerHist[gc].resize( freq + 1 );
	if ( freq > maxFreq ) maxFreq = freq;
	break;
      }
      if ( freq == 0 )
	gcKmerHist[gc].resize( 0 ); // no data !!!!!!!!!!
    }
  }
   

  
  // gcHist contains a histogram of numbers of kmer occurences (instances)
  vec< vec<longlong> > gcHist( gcKmerHist.isize() );
  vec<longlong> Hist( maxFreq + 1 );
  // gcCount contains numbers of all kmer occurences (instances) having a given gc content
  vec<longlong> gcCount( gcHist.isize(), 0 );
  vec<longlong> gcKmerCount( gcKmerHist.isize(), 0 );
  longlong allCount = 0;
  for ( int gc = 0; gc < gcKmerHist.isize(); gc++){
    gcHist[gc].resize( gcKmerHist[ gc ].isize() );
    for (unsigned i = 0; i < gcKmerHist[ gc ].size(); i++){
      gcHist[gc][i]    = gcKmerHist[ gc ][ i ] * i;
      gcCount[gc]     += gcHist[ gc ][ i ];
      gcKmerCount[gc] += gcKmerHist[ gc ][ i ]; 
      Hist[ i ] += 
	gcHist[gc][i];
      allCount += gcHist[gc][i];
    }
  }
   
  vec<double> bias( Ksize + 1, 1.0 ); // initialize it to 1.0 - if values cannot be computed than 
                                      // flat curve is going to be written out
  vec<double> tailFreq( Ksize +1, 0.0 );
  
  // NOTE: we use a frequency above which a specified fraction of the data (TAIL_FRACTION) is 
  // located in the histogram. This is not a precise approach, since for example such fraction
  // can be located at different distances from average in different GC bins. It can also be 
  // influance by varying numbers of higher copy-number kmers. However, the computations
  // following this module in the assembly pipeline are not very sensitive to small changes in 
  // the bias curve. It is more important to get the picture approximately right.
  // Computing locations of frequency peaks is tricky because of the nature of the data.

  // first find (global) frequency of the histogram tail
  longlong globalTailFreq = 0;
  longlong sumLoc = 0;
  for ( longlong i = Hist.isize() -1; i > 0; i--){
    sumLoc += Hist[ i ];
    if ( (double)sumLoc/allCount >= TAIL_FRACTION ){
      globalTailFreq = i;
      break;
    }
  }
  if ( VERBOSE ) PRINT( globalTailFreq );
  
  
  // fore each GC value find a frequency with most kmers in the histogram tail and perform sanity checks
  double MAX_REASONABLE_FREQ       = 1.5 * globalTailFreq; // HERUSTIC!!!!!!!!!!
  unsigned MINIMUM_KMERS           = 500;      // HEURISTIC!!!!!!!!
  int gcMaxTail                    = -1;
  longlong gcMaxTailCount         = 0;
  vec<longlong> gcTailCount( gcHist.isize() ); // approximate number of kmer occurences with copy number 1
  for ( int gc = 0; gc <= Ksize; ++gc ) {
    longlong sumLoc = 0;
    if ( gcKmerCount[ gc ] < MINIMUM_KMERS ) continue;
    for ( ulonglong i = gcHist[gc].isize() -1; i > 0; i--){
      sumLoc += gcHist[gc][ i ];
      if ( (double)sumLoc / gcCount[gc] >= TAIL_FRACTION && i < MAX_REASONABLE_FREQ ){
	tailFreq[ gc ] = i;
	gcTailCount[ gc ] = sumLoc;
	if ( sumLoc > gcMaxTailCount ){
	  gcMaxTail = gc;
	  gcMaxTailCount = sumLoc;
	}
	break;
      }
    }
  }
  
  
  
  int upperCurveLimit = gcMaxTail;
  int lowerCurveLimit = gcMaxTail;
  for ( int gc = gcMaxTail +1; gc <= Ksize; gc++){
    if( tailFreq[gc] <= 0 ){ 
      upperCurveLimit = gc -1;
      break;
    }else{ upperCurveLimit = gc; }
  }
  for ( int gc = gcMaxTail -1; gc >= 0; --gc ){
    if( tailFreq[gc] <= 0 ){
      lowerCurveLimit = gc +1;
      break;
    }else{ lowerCurveLimit = gc; }
  }
  
  if ( VERBOSE ) cout << "lowerCurveLimit = " << lowerCurveLimit << endl;
  if ( VERBOSE ) cout << "upperCurveLimit = " << upperCurveLimit << endl;
  
  
  // smooth the frequencies
  vec<longlong> tailFreqTmp( tailFreq.isize(), 0);
  for (int i=0; i < tailFreq.isize(); i++) 
    tailFreqTmp[i] = (longlong) tailFreq[i];
  smooth_histogram( tailFreqTmp, tailFreq, 1, 1, lowerCurveLimit, upperCurveLimit );
  
  // extend over gc's with undefined maximum values
  for (int gc = 0; gc < lowerCurveLimit; gc++ )
    tailFreq[ gc ] = tailFreq[ lowerCurveLimit ];
  for (int gc = upperCurveLimit + 1; gc < tailFreq.isize(); gc++ )
    tailFreq[ gc ] = tailFreq[ upperCurveLimit ];
  
 
  longlong totalSamples = 0;
  double sumMaximums    = 0;
  for ( int gc = 0; gc <= Ksize; ++gc ) {
    totalSamples   += gcCount[gc];
    sumMaximums    += 1.0 / tailFreq[gc] * (double)gcCount[gc];       
  }
  double aveMaximum = sumMaximums / (double) totalSamples;
 
  for ( int gc = 0; gc <= Ksize; ++gc ) 
    bias[gc] = (double) tailFreq[gc] * aveMaximum;  


  if ( VERBOSE ){
    cout << "BIAS\n";
    for ( int gc = 0; gc <= Ksize; ++gc )
      PRINT3( gc, tailFreq[gc], bias[gc] );
  }
  
   
  // Save the file.
  String out = runDir + "/" + BIAS_CURVE_OUT + ".k" + ToString(K);
  ofstream outstrm( out.c_str() );
  outstrm << bias;
  outstrm.close();
  
  cout << Date( ) << ": Done with GcBiasCurveLG!" << endl;
  return 0;
}
  
       
       
       
  

