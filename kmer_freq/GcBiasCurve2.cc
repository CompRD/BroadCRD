/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/** 
   Program: GcBiasCurve

   Estimate and store the <GC bias> curve for a given K value.  This shows,
   for each possible GC content, the positive or negative bias for reading
   genome regions with that gc content: do such genome regions get relatively
   more reads than average, or relatively fewer reads than average?

   This is a modified version of the original script that assumed that does not 
   assume that bias computation is possible for each GC content. It also requires
   additional input of a genome size.

   The GC bias curve computed by this program can be used in the following ways:

      - the <SeqToPaths> program can use it to generate simulated reads with GC bias

      - the <UnibaseCopyNumber> program can use it to calculate unbibase copy nos

   INPUT FILES:
     reads.fastb

   OUTPUT FILES:
     reads.fastb.unpaired.freq.k*
     reads.gc_bias.k*

   @file
*/

#include "MainTools.h"
#include "kmers/KmerShape.h"
#include "ReadPairing.h"
#include "kmer_freq/KmerFrequencyTable.h"
#include "kmer_freq/WriteKmerFrequencies.h"
#include <map>

/**
   Function: smooth_histogram
   This is a local function that smooths the histogram between start and end positions
   by computing an average value in a section specified by radiiusLeft and radiusRight
   (to allow nonsymmetric computations)

*/
void smooth_histogram( const vec<longlong>& histo, vec<double>& smoothHisto, 
		       int radiusLeft, int radiusRight, int start, int end ){
  ForceAssertGe(radiusLeft, 0);
  ForceAssertGe(radiusRight, 0);
  ForceAssertGe(start, 0);
  ForceAssertLe(end, histo.isize() -1);
  if ( histo.isize() != smoothHisto.isize() ){
    //cout << "WARNING: size of smoothed histogram not equal to the size of an original. Resizing\n";
    smoothHisto.resize( histo.isize() );
  }
  for (int i = start; i <= end; i++ )
    smoothHisto[i] = 0;
  for (int i = 0; i < start; i++)
    smoothHisto[i] = histo[i];
  for (int i = end + 1; i < histo.isize(); i++)
    smoothHisto[i] = histo[i];
  for (int i = start; i <= end; i++ ){
    int leftRange = (i - radiusLeft >= start ) ? i - radiusLeft : start;
    int rightRange = ( i + radiusRight <= end ) ? i + radiusRight : end;
    int fullRange = rightRange - leftRange + 1;
    longlong rangeSum = 0;
    for ( int j = leftRange; j <= rightRange; j++ )
      rangeSum += histo[j];
    double smoothVal = (double) rangeSum / (double) fullRange;
    smoothHisto[i] = smoothVal;
  }
  return;
} 


int main( int argc, char** argv ) {
  
  RunTime();
  
  BeginCommandArguments;
  
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_KShape(K);
  
  CommandArgument_String_OrDefault(READS_IN,"reads.fastb");
  CommandArgument_Int_OrDefault(BEGIN, 0 );
  CommandArgument_Int_OrDefault(END, -1 );
  CommandArgument_UnsignedInt_OrDefault(TRUNC, 0 );
  CommandArgument_String_OrDefault(BIAS_CURVE_OUT,"reads.gc_bias");

  CommandArgument_Bool_OrDefault(EXCLUDE_PAIRED, True);
  
  EndCommandArguments;

  String dataDir = PRE + "/" + DATA;
  String runDir = PRE + "/" + DATA + "/" + RUN;
  String freqTableFile = runDir + "/" + READS_IN + ".unpaired.freq.k" + ToString(K);



  
  vecbasevector reads;
  if ( BEGIN == 0 && END == -1 )
    reads.ReadAll( runDir + "/" + READS_IN );
  else {
    ForceAssertLe( BEGIN, END );
    reads.SparseReadRange( runDir + "/" + READS_IN, BEGIN, END );
  }
 
  
  // Get genome size.  Note that this is cheating. Approximate genome size should be used.
  longlong genomeSize = StringOfFile( runDir + "/genome.size", 1 ).Int( );

  /**
  String QualsFile = dataDir + "/" + READS_IN + ".qualb";
  int nqualReads = MastervecFileObjectCount(QualsFile);
  if ( nqualReads != reads.isize() ){
    cout << "ERROR: number of reads with qualities not the same as the number of reads ";
    cout << nqualReads << " != " << reads.isize() << " Stopping!" << endl;
    exit(1);
  }
  /// this is a, currently not used, part that attempts to compute the GC-bias from 
  quality scores. This could account for the read-out bias but perhaps not for the PCR-bias
  related to read library preparation with PCR.
  
  cout << Date() << " computing average misscall probabilities" << endl;
  map< int,map<char,double> > base2missProb;
  map< int,map<char,longlong> >   base2count;
  int batchsize = 100000;
  {    
    vecqualvector quals;                                                
    for ( int i = 0; i < reads.isize(); i += batchsize ){
      int start = i, stop = Min( reads.isize(), i + batchsize );            
      quals.ReadRange( QualsFile, start, stop );
      for ( int j = 0; j < quals.isize( ); j++ ){
	if( reads[ i + j ].isize() != quals[j].isize() ){
	  cout << "ERROR: read and qual vector sizes not equal. Stopping" << endl;
	  exit(1);
	}
	String sread = reads[i+j].ToString();
	if( sread.isize() != quals[j].isize() ){
	  cout << "ERROR: sread and qual vector sizes not equal. Stopping" << endl;
	  exit(1);
	} 
	int GC=0;
	for ( int pos = 0; pos < reads[i+j].isize(); pos++)
	  if( reads[i+j].at(pos) == 'G' || reads[i+j].at(pos) == 'C' )
	    GC++;
	for ( int pos = 0; pos < quals[j].isize(); pos++){
	  //if ( quals[j][pos] < 30) continue;
	  char base   = reads[i+j].at(pos);
	  double prob = pow( 10, -1 * (double)quals[j][pos]/10.0 );
	  base2missProb[GC][base] += prob;
	  base2count[GC][base]++;
	}
      }
      quals.clear( );    
      }    
  }
  
  cout << "base misscall probabilities: " << endl;
  for ( map< int, map<char,double> >::iterator i1 = base2missProb.begin(); i1 != base2missProb.end(); i1++){
    int gc = i1->first;
    for ( map<char,double>::iterator i2 = base2missProb[gc].begin(); i2 != base2missProb[gc].end(); i2++){
      i2->second = i2->second / (double)base2count[gc][ i2->first];
      cout << i1->first << "\t" << i2->first << "\t" << i2->second << endl;
    }
  }
  cout << Date() << " misscall probabilites computed" << endl;
  */
  
  if ( EXCLUDE_PAIRED ) {
    vec<read_pairing> pairs;
    ReadPairsFile( runDir + "/reads.pairto", pairs );
    for ( unsigned int p = 0; p < pairs.size(); ++p )
      if ( ! pairs[p].Dead() ) {
        reads[ pairs[p].id1 ].clear().shrink_to_fit();
        reads[ pairs[p].id2 ].clear().shrink_to_fit();
      }
  }

  if ( TRUNC > 0 )
    for ( size_t r = 0; r < reads.size(); r++ )
      if ( reads[r].size() > TRUNC )
        reads[r].resize( TRUNC );

  size_t sizeReadSum = 0;
  for ( size_t i = 0; i < reads.size(); i++ ) 
    sizeReadSum += reads[i].size();
  
  
  size_t nreads = reads.size();
  cout << "number of reads = " << nreads << endl;
  cout << "number of bases in reads = " << sizeReadSum << endl;

  size_t aveReadLen = sizeReadSum/nreads;
  cout << "average size of a read = " << aveReadLen << endl;

  int expectedCoverage = (int)( nreads * aveReadLen / (double)(genomeSize * 2) );
  cout << "expected coverage per position = " << expectedCoverage << endl;

  int coverageDev = (int) ( sqrt( expectedCoverage ) );
  cout << "standard deviation of coverage (Poisson) = " << coverageDev << endl;


#define CASE(K) WriteKmerFrequencies<K>( reads, freqTableFile )
  DISPATCH_ON_KSHAPE( K, CASE );
    
  KmerFrequencyTable freqTable( K, freqTableFile );

  // Get the peaks of the kmer frequency curve for each GC content.
  // If the peak is not calculated (usually due to insufficient sample
  // size), set the peak to 1.
  /// those peaks are no longer uses in calcullation. Kept here temporarily.
  vec<int> peaks( GetKmerSize( K )+1);
  for ( int gc = 0; gc <= GetKmerSize( K ); ++gc ) {
    peaks[gc] = freqTable.GetPeak( gc );
    cout << "before max(): gc = " << gc << " peak = " << peaks[gc] << endl;
    peaks[gc] = max( 1, peaks[gc] );
    cout << "after: gc = " << gc << " peak = " << peaks[gc] << endl;
  }

  // Calculate the number of *instances* of kmers with each GC
  // content.  These become the weights for a weighted average of the
  // peaks of the frequency curves.

  
  vec< double > averages( GetKmerSize( K ) +1,0.0);

  /// gcKmerHist contains information about numbers of different kmers with a 
  /// given frequency of occurance
  vec< vec<longlong> > gcKmerHist;
  freqTable.GetGCHistograms( gcKmerHist );

  /// gcHist contains information about numbers of all kmers in the genome with
  /// a given frequency of occurance
  vec< vec<longlong> > gcHist( gcKmerHist.isize() );
  vec<long> gcCount( gcHist.isize() );
  for ( int gc = 0; gc < gcHist.isize(); gc++){
    gcHist[gc].resize( gcKmerHist[gc].isize() );
    for (int i = 0; i < gcHist[gc].isize(); i++){
      gcHist[gc][i] = gcKmerHist[gc][i] * i;
      gcCount[gc]  += gcHist[gc][i];
    }
  }



  /// smooth the histogram
  vec< vec<double> > gcHistSmoothLeft( gcHist.isize() );
  vec< vec<double> > gcHistSmoothRight( gcHist.isize() );
  vec< vec<double> > gcHistSmoothCenter( gcHist.isize() );
  for ( int gc = 0; gc <= GetKmerSize( K ); ++gc ) {
    longlong gcTotal = 0;
    /// do not consider rather unreasonable repeat frequencies (we are looking for a copy number one)
    int maxFreq2consid =
      (gcHist[gc].isize() -1 < expectedCoverage *2) ? gcHist[gc].isize() -1 : expectedCoverage *2; 
  
    /// unique kmers are not included in the histogram (because of the way WriteKmerFrequencies was
    /// called. Therefore, we do not want to include those initial zeros in smoothing
    int minPresentFreq = 0;
    for ( int freq = 1; freq <= maxFreq2consid; ++freq ){
      if ( gcHist[gc][freq] > 0 ){
	minPresentFreq = freq;
	break;
      }
    } 
    if ( minPresentFreq == 0){
      cout << "ERROR: minPresentFreq is 0 for gc = " << gc << ", stopping at line " << __LINE__ << endl; 
      exit(1);
    }
    int smoothingRadiusLeft  = 2;
    int smoothingRadiusRight = 0;
    smooth_histogram( gcHist[gc], gcHistSmoothLeft[gc], smoothingRadiusLeft, smoothingRadiusRight, 
		      minPresentFreq, gcHist[gc].isize() -1 );
    smoothingRadiusLeft  = 0;
    smoothingRadiusRight = 2;
    smooth_histogram( gcHist[gc], gcHistSmoothRight[gc], smoothingRadiusLeft, smoothingRadiusRight, 
		      minPresentFreq, gcHist[gc].isize() -1 );
    smoothingRadiusLeft  = 1;
    smoothingRadiusRight = 1;
    smooth_histogram( gcHist[gc], gcHistSmoothCenter[gc], smoothingRadiusLeft, smoothingRadiusRight, 
		      minPresentFreq, gcHist[gc].isize() -1 );
  }  

 

  
  // find a frequency with most kmers and perform sanity checks
  vec< double > cn1freq( GetKmerSize( K ) +1, 0.0 );
  for ( int gc = 0; gc <= GetKmerSize( K ); ++gc ) {
    longlong gcTotal = 0;
    int maxFreq2consid =
      (gcHist[gc].isize() -1 < expectedCoverage *2) ? gcHist[gc].isize() -1 : expectedCoverage *2; 

    int minPresentFreq = 0;
    for ( int freq = 1; freq <= maxFreq2consid; ++freq ){
      if ( gcHist[gc][freq] > 0 ){
	minPresentFreq = freq;
	break;
      }
    } 
    if ( minPresentFreq == 0){
      cout << "ERROR: minPresentFreq is 0 for gc = " << gc << ", stopping at " << __LINE__ << endl; 
      exit(1);
    }
    //cout << "starting frequency for gc = " << gc << " is " << minPresentFreq << endl;
    double maxVal = gcHistSmoothRight[gc][minPresentFreq];
    int maxFreq = minPresentFreq;
    for ( int freq = minPresentFreq +1; freq <= maxFreq2consid; ++freq ) {
      if ( gcHistSmoothRight[gc][freq] > maxVal ){
	maxFreq = freq;
	maxVal  = gcHistSmoothCenter[gc][freq];
      }
    }
    
    if ( maxFreq == minPresentFreq ){
      /// copy number one maximum is too sparsely populated, so that the lowest frequency bin
      /// contains most kmers
      cn1freq[gc]  = -1;
    }else if ( maxVal <= maxFreq2consid ){
      /// single kmer with limiting frequency would have more instances than all copy number one kmers
      cn1freq[gc]  = -1;
    }else{
      cn1freq[gc]  = maxFreq;
    }
  }

  /*  diagnostic output
  cout << "cn1freq original: " << endl;
  for ( int i = 0; i < cn1freq.isize(); i++ ){
    cout << i << "\t" << cn1freq[i] << endl;
  }
  */

  /// examine and smooth the maximums curve
  ///    first find GC bin with maximum number of kmers
  vec<double> newbias( GetKmerSize( K ) + 1, 0.0 );
  vec<double> bias( GetKmerSize( K ) + 1 );  /// old way of calculating bias
  long maxCountForGc = -1;
  int gcWithMaxCount = -1;
  for ( int gc = 0; gc <= GetKmerSize( K ); ++gc ) {
    if ( gcCount[gc] > maxCountForGc){
      gcWithMaxCount = gc;
      maxCountForGc  = gcCount[gc];
    }
  }
  if ( gcWithMaxCount == -1 ){
    cout << "ERROR: did not find any kmers in histogram, stopping" << endl;
    exit(1);
  }
  cout << "GC value with maximum number of kmers is " << gcWithMaxCount << endl;
  cout << " and the estimated copy number one for this GC is " << cn1freq[gcWithMaxCount] << endl;


  Bool FlatCurve = False;
  if ( cn1freq[gcWithMaxCount] < 0 ){
    cout << "WARNING: could not make any GC-bias estimates. Assuming that the reason for this is "<< endl;
    cout << "too few samples. Will output a flat curve (no GC bias)\n";
    for ( int i = 0; i < newbias.isize(); i++) newbias[i] = 1.0;
  }else{

    int upperCurveLimit = gcWithMaxCount;
    int lowerCurveLimit = gcWithMaxCount;
    for (int gc = gcWithMaxCount +1; gc <= GetKmerSize( K ); gc++){
      if( cn1freq[gc] < 0 ){ 
	upperCurveLimit = gc -1;
	break;
      }else{ upperCurveLimit = gc; }
    }
    for (int gc = gcWithMaxCount -1; gc >= 0; gc--){
      if( cn1freq[gc] < 0 ){
	lowerCurveLimit = gc +1;
	break;
      }else{ lowerCurveLimit = gc; }
    }
    
    cout << "lowerCurveLimit = " << lowerCurveLimit << endl;
    cout << "upperCurveLimit = " << upperCurveLimit << endl;
    
    for (int gc = 0; gc < lowerCurveLimit; gc++ )
      cn1freq[gc] = cn1freq[lowerCurveLimit];
    for (int gc = upperCurveLimit + 1; gc < cn1freq.isize(); gc++ )
      cn1freq[gc] = cn1freq[upperCurveLimit];
    
    
    
    vec<longlong> cn1freqTmp( cn1freq.isize(), 0);
    for (int i=0; i<cn1freq.isize(); i++) cn1freqTmp[i] = (longlong) cn1freq[i];
    //cout << "Before smoothing\n";
    //for (int i=0; i<cn1freq.isize(); i++) cout << i << "\t" << cn1freq[i] << "\n";
    //cout << endl;
    smooth_histogram( cn1freqTmp, cn1freq, 1, 1, lowerCurveLimit, upperCurveLimit );
    //cout << "\nAfter smoothing\n";
    //for (int i=0; i<cn1freq.isize(); i++) cout << i << "\t" << cn1freq[i] << "\n";
    //cout << endl;
    
    longlong totalSamples = 0;
    longlong sumPeaks     = 0;
    double sumMaximums  = 0;
    for ( int gc = 0; gc <= GetKmerSize( K ); ++gc ) {
      longlong totalSamplesForGc = 0;
      double aveFreqForGc = 0.0;
      int maxFreq2consid =(gcHist[gc].isize() < expectedCoverage * 2) ? gcHist[gc].size() : expectedCoverage * 2; 
      for ( int freq = 0; freq <= maxFreq2consid; ++freq ) {
	aveFreqForGc += gcHist[gc][freq];
      }
      for ( int freq = 0; freq < gcHist[gc].isize(); ++freq ) {
	totalSamplesForGc += gcHist[gc][freq];
      }
      totalSamples   += totalSamplesForGc;
      sumPeaks       += (longlong) ( 1.0 / (double)peaks[gc] * (double)totalSamplesForGc );
      sumMaximums    += 1.0 / cn1freq[gc] * (double)totalSamplesForGc;
    }
    double avePeak = (double)sumPeaks / (double)totalSamples;
    double aveMaximum = sumMaximums / (double) totalSamples;
    
    /* diagnostic output
    for ( int gc = 0; gc <= GetKmerSize( K ); ++gc ) {
      double frac = (double)gcCount[gc]/totalSamples;
      cout << "gc = " << gc << " fraction " << frac << endl;
    }
    */
    for ( int gc = 0; gc <= GetKmerSize( K ); ++gc ){ 
      newbias[gc] = (double) cn1freq[gc] * aveMaximum;
      bias[gc] = (double) peaks[gc] * avePeak;
    }
  }
  
  cout << "NEW BIAS\n";
  for ( int gc = 0; gc <= GetKmerSize( K ); ++gc )
    PRINT3( gc, cn1freq[gc], newbias[gc] );
  
  // The bias is calculated as the ratio of the peak for that GC to
  // the weighted average peak across all GCs.
  cout << "\nORIGINAL BIAS (for comparison)\n";
  for ( int gc = 0; gc <= GetKmerSize( K ); ++gc )
    PRINT3( gc, peaks[gc], bias[gc] );
  
  // Save the file.
  String out = runDir + "/" + BIAS_CURVE_OUT + ".k" + ToString(K);
  ofstream outstrm( out.c_str() );
  outstrm << newbias;
  outstrm.close();
  
  return 0;
}
  
       
       
       
  

