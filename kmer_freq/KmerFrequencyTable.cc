/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "kmer_freq/KmerFrequencyTable.h"

#include "math/Functions.h"
#include "kmers/KmerRecord.h"
#include "Vec.h"

/**
   \file
   \copydoc KmerFrequencyTable
   \ingroup grp_edits
 */

KmerFrequencyTable::KmerFrequencyTable( const KmerShapeId& KSHAPE, const String& filename )
  : m_kmerMap( KSHAPE, filename ),
    m_pCache( new cache )
{}

KmerFrequencyTable::~KmerFrequencyTable() {
  delete m_pCache;
}

void 
KmerFrequencyTable::GetGCBands( const vec<longlong>& totalByGC,
                                const int min_samples,
                                const int maxSampleRatio,
                                const bool includeTails,
                                vec< pair<int,int> >& bands ) const {
  vec<longlong>::const_iterator maxIter = 
    max_element( totalByGC.begin(), totalByGC.end() );
  int maxIdx = distance( totalByGC.begin(), maxIter );
  longlong maxValue = *maxIter;
  longlong minTotal = max<longlong>( maxValue/maxSampleRatio, min_samples );
  
  bands.push_back( make_pair( maxIdx, maxIdx ) );
  
  int gcLow = maxIdx+1, gcHigh = maxIdx+1;
  while ( gcHigh < totalByGC.isize() ) {
    int bandTotal = 0;
    for ( int gc = gcLow; gc <= gcHigh; ++gc )
      bandTotal += totalByGC[gc];
    if ( bandTotal < minTotal ) {
      if ( gcHigh < totalByGC.isize()-1 ) {
        gcHigh++;
        continue;
      } 
      else {
        if ( includeTails )
          bands.back().second = gcHigh;
        break;
      }
    }
    bands.push_back( make_pair( gcLow, gcHigh ) );
    gcLow = ++gcHigh;
  }
  
  gcLow = maxIdx-1, gcHigh = maxIdx-1;
  while ( gcLow >= 0 ) {
    int bandTotal = 0;
    for ( int gc = gcLow; gc <= gcHigh; ++gc )
      bandTotal += totalByGC[gc];
    if ( bandTotal < minTotal ) {
      if ( gcLow > 0 ) {
        gcLow--;
        continue;
      }
      else {
        if ( includeTails )
          bands.back().first = gcLow;
        break;
      }
    }
    bands.push_back( make_pair( gcLow, gcHigh ) );
    gcHigh = --gcLow;
  }
}

// Returns the first frequency with the highest count.
int
KmerFrequencyTable::GetPeak( const int gc_content ) const {
  if ( gc_content < 0 ) {
    if ( m_pCache->peak == -1 ) {
      vec<longlong> hist;
      this->GetHistogram( hist );
      m_pCache->peak = 1;
      for ( int freq = 2; freq < hist.isize(); ++freq )
        if ( freq*hist[freq] > (m_pCache->peak)*hist[m_pCache->peak] )
          m_pCache->peak = freq;
    }
    return m_pCache->peak;
  } else {
    const int K = this->GetKmerSize();
    ForceAssertLe( gc_content, K );
    if ( m_pCache->peakGC.empty() ) {
      m_pCache->peakGC.resize(K+1,0);
      
      vec< vec<longlong> > histByGC;
      this->GetGCHistograms( histByGC );

      vec<longlong> totalByGC( histByGC.size(), 0 );
      for ( unsigned int i = 0; i < histByGC.size(); ++i ) 
        for ( unsigned int j = 0; j < histByGC[i].size(); ++j )
          totalByGC[i] += histByGC[i][j];

      vec< pair<int,int> > gcBands;
      const int minSamples = 0;
      const int maxSampleRatio = 10;
      const bool includeTails = false;
      this->GetGCBands( totalByGC, minSamples, maxSampleRatio, includeTails, gcBands );

      const int maxFreq = m_kmerMap.GetMaxValue();
      for ( unsigned int b = 0; b < gcBands.size(); ++b ) {
        int gcLow = gcBands[b].first, gcHigh = gcBands[b].second;
        int peak = 0;
        longlong peakValue = 0;
        for ( int gc = gcLow; gc <= gcHigh; ++gc )
          peakValue += histByGC[gc][peak];
        for ( int freq = 1; freq < maxFreq+1; ++freq ) {
          longlong freqValue = 0; 
          for ( int gc = gcLow; gc <= gcHigh; ++gc )
            freqValue += histByGC[gc][freq];
          freqValue *= freq;
          if ( freqValue > peakValue ) {
            peak = freq;
            peakValue = freqValue;
          }
        }
        for ( int gc = gcLow; gc <= gcHigh; ++gc )
          m_pCache->peakGC[gc] = peak;
      }          
    }
    return m_pCache->peakGC[gc_content];
  }
}

int
KmerFrequencyTable::GetFLM( const vec<longlong>& hist, 
                            const double threshold, 
                            const int minSamples ) const {
  longlong total = 0, totalInstances = 0;
  for ( int freq = 1; freq < hist.isize(); ++freq ) {
    total += hist[freq];
    totalInstances += freq*hist[freq];
  }
  if ( total < minSamples )
    return -1;
  
  int median = 0;
  longlong medianTotal = 0;
  longlong medianTarget = totalInstances/2;
  for ( ; median < hist.isize(); ++median ) {
    medianTotal += median*hist[median];
    if ( medianTotal > medianTarget )
      break;
  }
  if ( median == hist.isize() )
    return -1;
  
  int lowerBound = 0;
  while ( hist[lowerBound] == 0 && lowerBound <= median )
    ++lowerBound;
  int localMin = median;
  int localMax = median;
  for ( int i = median-1; i >= lowerBound; --i ) {
    if ( hist[i] <= hist[localMin] &&
         (double)hist[i]*threshold <= (double)hist[localMax] )
      localMin = i;
    else if ( hist[i] > hist[localMax] )
      localMax = i;
  }
  
  if ( localMin == median )
    return -1;
  
  return localMin;
}

int
KmerFrequencyTable::GetFirstLocalMin( const double threshold, 
                                      const int gc_content, 
                                      const int min_samples ) const {
  if ( gc_content < 0 ) {
    if ( m_pCache->localMin == -2 ) {
      vec<longlong> hist;
      this->GetHistogram( hist );
      m_pCache->localMin = this->GetFLM( hist, threshold, min_samples );
    }
    return m_pCache->localMin;
  } else {
    const int K = this->GetKmerSize();
    ForceAssertLe( gc_content, K );
    if ( m_pCache->localMinGC.empty() ) {
      m_pCache->localMinGC.resize(K+1,-1);
      
      vec< vec<longlong> > histByGC;
      this->GetGCHistograms( histByGC );
      
      vec<longlong> totalByGC( histByGC.size(), 0 );
      for ( unsigned int gc = 0; gc < histByGC.size(); ++gc ) 
        for ( unsigned int freq = 0; freq < histByGC[gc].size(); ++freq )
          totalByGC[gc] += histByGC[gc][freq];
      
      vec< pair<int,int> > gcBands;
      const int maxSampleRatio = K;
      const bool includeTails = true;
      this->GetGCBands( totalByGC, min_samples, maxSampleRatio, includeTails, gcBands );
      sort( gcBands.begin(), gcBands.end() );
      
      for ( unsigned int b = 0; b < gcBands.size(); ++b ) {
        int gcLow = gcBands[b].first, gcHigh = gcBands[b].second;
        PRINT2( gcLow, gcHigh );
        
        vec<longlong> bandHist( histByGC[gcLow].size(), 0 );
        for ( int gc = gcLow; gc <= gcHigh; ++gc )
          for ( int freq = 0; freq < histByGC[gc].isize(); ++freq )
            bandHist[freq] += histByGC[gc][freq];
        
        int localMin = this->GetFLM( bandHist, threshold, min_samples );
        
        PRINT( localMin );
        for ( int gc = gcLow; gc <= gcHigh; ++gc )
          m_pCache->localMinGC[gc] = localMin;
      }
    }
    return m_pCache->localMinGC[gc_content];
  }
}
