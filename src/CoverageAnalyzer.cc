// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#include <fstream>
#include <functional>
#include <map>

#include "CoverageAnalyzer.h"
#include "Qualvector.h"
#include "SeqInterval.h"
#include "STLExtensions.h"

CoverageAnalyzer::CoverageAnalyzer(ostream* pOStream) 
  : m_pSeqLengths( 0 ),
    m_pOStream(pOStream)
{
}

CoverageAnalyzer::CoverageAnalyzer(const vec<seq_interval>& vecIntervals,
                                   const vec<int>* pVecSeqLengths,
                                   ostream* pOStream)
  : m_pSeqLengths( 0 ),
    m_pOStream(pOStream)
{
  CreateCoverages(vecIntervals, pVecSeqLengths);
}

void
CoverageAnalyzer::CreateCoverages(const vec<seq_interval>& vecIntervals,
                                  const vec<int>* pVecSeqLengths)
{
  m_covChanges.clear();

  if ( m_pSeqLengths )
    delete m_pSeqLengths;

  if ( pVecSeqLengths )
    m_pSeqLengths = new vec<int>( *pVecSeqLengths );
  else
    m_pSeqLengths = 0;

  for ( vec<seq_interval>::const_iterator iter = vecIntervals.begin();
	iter != vecIntervals.end(); ++iter )
    { 
      int idSeq  = iter->SeqId();
      int nStart = iter->Begin();
      int nStop  = iter->End();

      if ( nStart > nStop )
        {
	  if ( m_pOStream )
	    *m_pOStream << "seq_interval "
			<< iter->IntervalId()
			<< " is of negative length "
			<< "and will not be counted" << endl;
	  continue;
        }

      if (m_pSeqLengths)
        {
	  int seqLength = (*m_pSeqLengths)[idSeq];

	  if ( (nStart >= seqLength || nStop < 0) )
            {
	      if ( m_pOStream )
		*m_pOStream << "seq_interval "
			    << iter->IntervalId()
			    << " is out side of coverage region "
			    << "and will not be counted" << endl;
	      continue;
            }
                
	  if (nStart < 0)
            {
	      if (m_pOStream)
		*m_pOStream << "seq_interval "
			    << iter->IntervalId()
			    << " begins at " << nStart
			    << " (and will be reset to 0)." << endl;
	      nStart = 0;
            }

	  if (nStop > seqLength)
            {
	      if (m_pOStream)
		*m_pOStream << "seq_interval "
			    << iter->IntervalId()
			    << " ends "
			    << (nStop - seqLength)
			    << " base(s) over the end of sequence "
			    << "(and will be reset to"
			    << " length of sequence)." << endl;
	      nStop = seqLength;
            }
        }
        
      m_covChanges.push_back( CoverageChange( iter->SeqId(), nStart, +1 ) );
      m_covChanges.push_back( CoverageChange( iter->SeqId(), nStop,  -1 ) );
    }

  sort( m_covChanges.begin(), m_covChanges.end() );

  for ( unsigned int i = 1; i < m_covChanges.size(); ++i )
    {
      if ( m_covChanges[i].seqId == m_covChanges[i-1].seqId &&
	   m_covChanges[i].coord == m_covChanges[i-1].coord )
        {
	  m_covChanges[i].amount += m_covChanges[i-1].amount;
	  m_covChanges[i-1].amount = 0;
        }
    }

  m_covChanges.erase( remove_if( m_covChanges.begin(),
				 m_covChanges.end(),
				 mem_fun_ref( &CoverageChange::IsEmpty ) ),
		      m_covChanges.end() );
}


void
CoverageAnalyzer::GetCoveragesBetween(unsigned int minValue,
                                      unsigned int maxValue,
                                      vec<seq_interval>& vecIntervals) const
{
  vecIntervals.clear();

  if ( m_pSeqLengths )
    {
      vec<CoverageChange>::const_iterator change = m_covChanges.begin();
      
      for ( int seqId = 0; seqId < (int)m_pSeqLengths->size(); ++seqId )
        {
	  int seqLength = (*m_pSeqLengths)[seqId];
	  
	  // If we have gone through an entire contig without any coverage changes,
	  // assume the coverage is 0 throughout
	  if ( change == m_covChanges.end() || seqId < change->seqId )
            {
	      if ( minValue <= 0 )
		vecIntervals.push_back( seq_interval( -1, seqId, 0, seqLength ) );
                
	      continue;
            }
            
	  int begin = 0;
	  unsigned int cov = 0;
	  bool covIsInWindow = ( cov >= minValue && cov <= maxValue );

	  for ( ; change != m_covChanges.end() && change->seqId == seqId; ++change )
            {
	      unsigned int newCov = cov + change->amount;
	      bool newCovIsInWindow = ( newCov >= minValue && newCov <= maxValue );
                
	      if ( covIsInWindow && ! newCovIsInWindow )
                {
		  int end = change->coord;
		  if ( end - begin > 0 )
		    vecIntervals.push_back( seq_interval( -1, seqId, begin, end ) );
                }
	      else if ( ! covIsInWindow && newCovIsInWindow )
                {
		  begin = change->coord;
                }

	      cov = newCov;
	      covIsInWindow = newCovIsInWindow;
            }
            
	  if ( covIsInWindow && begin != seqLength )
	    vecIntervals.push_back( seq_interval( -1, seqId, begin, seqLength ) );
        }
    }

  else
    {
      int lastSeqId = -1;
      int begin = 0;
      unsigned int cov = 0;
      bool covIsInWindow = false;

      vec<CoverageChange>::const_iterator change = m_covChanges.begin();
      for ( ; change != m_covChanges.end(); ++change )
        {
	  if ( change->seqId == lastSeqId )
            {
	      unsigned int newCov = cov + change->amount;
	      int newCovIsInWindow = ( newCov >= minValue && newCov <= maxValue );
                
	      if ( covIsInWindow && ! newCovIsInWindow )
                {
		  int end = change->coord;
		  vecIntervals.push_back( seq_interval( -1, change->seqId, begin, end ) );
                }
	      else if ( ! covIsInWindow && newCovIsInWindow )
                {
		  begin = change->coord;
                }

	      cov = newCov;
	      covIsInWindow = newCovIsInWindow;
            }
	  else
            {
	      lastSeqId = change->seqId;
	      begin = change->coord;
	      cov = change->amount;
	      covIsInWindow = ( cov >= minValue && cov <= maxValue );
            }
        }
    }
}

void
CoverageAnalyzer::GetAllCoverages(vec<seq_interval>& vecIntervals) const
{
  vecIntervals.clear();

  if ( m_pSeqLengths )
    {
      vec<CoverageChange>::const_iterator change = m_covChanges.begin();
        
      for ( int seqId = 0; seqId < (int)m_pSeqLengths->size(); ++seqId )
        {
	  int seqLength = (*m_pSeqLengths)[seqId];

	  if ( change == m_covChanges.end() || seqId < change->seqId )
            {
	      vecIntervals.push_back( seq_interval( 0, seqId, 0, seqLength ) );
                
	      continue;
            }
            
	  int begin = 0;
	  int cov = 0;

	  for ( ; change != m_covChanges.end() && change->seqId == seqId; ++change )
            {
	      int end = change->coord;
                    
	      if ( end - begin > 0 )
		vecIntervals.push_back( seq_interval( cov, seqId, begin, end ) );

	      cov += change->amount;
	      begin = change->coord;
            }
            
	  if ( begin != seqLength )
	    vecIntervals.push_back( seq_interval( cov, seqId, begin, seqLength ) );
        }
    }

  else
    {
      int lastSeqId = -1;
      int begin = 0;
      int cov = 0;

      vec<CoverageChange>::const_iterator change = m_covChanges.begin();
      for ( ; change != m_covChanges.end(); ++change )
        {
	  if ( change->seqId == lastSeqId )
            {
	      int end = change->coord;
	      vecIntervals.push_back( seq_interval( cov, change->seqId, begin, end ) );

	      begin = change->coord;
	      cov += change->amount;
            }
	  else
            {
	      lastSeqId = change->seqId;
	      begin = change->coord;
	      cov = change->amount;
            }
        }
    }
}

void CoverageAnalyzer::GetAllCoverageInWindows(vec<seq_interval>& intervals,
                                               const unsigned int windowSize,
                                               const unsigned int scale ) const
{
  const double fScale = scale;

  GetAllCoverages(intervals);

  if ( m_pSeqLengths ) {
    vec<seq_interval>::iterator iInterval = intervals.begin();
    
    vec<seq_interval> windows;
    for ( int seq = 0; seq < (int) m_pSeqLengths->size(); ++seq ) {
      
      if ( iInterval != intervals.end() && iInterval->SeqId() > seq )
        continue;

      int seqLen = (*m_pSeqLengths)[seq];

      for ( int begin = 0; begin < seqLen; begin += windowSize ) {
        int end = min<int>( begin + windowSize, seqLen );
        seq_interval window( 0, seq, begin, end );
        int sumOfValues = 0;
        while ( iInterval != intervals.end() &&
                iInterval->SeqId() == seq &&
                iInterval->Begin() < end ) {
          int partInWindow = iInterval->HasAmountOfOverlapWith( window );
          sumOfValues += partInWindow * iInterval->IntervalId();
          if ( iInterval->End() <= end )
            ++iInterval;
          else
            break;
        }
        window.SetIntervalId( (int) round( fScale * (double)sumOfValues / (double)(end-begin) ) );
        windows.push_back( window );
      }
    }
    
    intervals.swap( windows );
  }

  else {
    cout << "CoverageAnalyzer::GetAllCoverageInWindows() is not defined when " << endl;
    cout << "sequence lengths have not been provided." << endl;
    
    TracebackThisProcess();
  }
}

void CoverageAnalyzer::CountAllCoverages(vec<longlong> & counts) const {
  vec<seq_interval> intervals;
  GetAllCoverages(intervals);
  int maxcov=0;
  for (int i=0; i != intervals.isize(); ++i) {
    maxcov = max(maxcov, intervals[i].IntervalId());
  }
  counts.assign(maxcov+1, 0);
  for (int i=0; i != intervals.isize(); ++i) {
    counts[intervals[i].IntervalId()] += intervals[i].Length();
  }
}
  



void CoverageAnalyzer::PrintCoverageStats(ostream & os,
					  bool printZeroCoverageContigs) {
  vec<seq_interval> intervals;
  GetAllCoverages(intervals);
  vec<seq_interval>::iterator inter;

  if (!m_pSeqLengths) {
     os << "Error in CoverageAnalyzer::PrintCoverageStats!  Null m_pSeqLengths.\n";
     return;
  }
  longlong numContigs = m_pSeqLengths->size();

  // Below is a somewhat slow way of computing Q3/Q1.  It is like computing the
  // mean.  One could compute the mean of numbers by the usual sorting and 
  // picking the middle number. This is what is being done below.  map<int,int>
  // automatically sorts in increasing size of the first pair of numbers. 
  // This works quickly since the number of different coverages is quite small.
  // However, a faster way to compute the mean of numbers is to pick a random
  // pivot number, partition the set of numbers, then one knows the median lies
  // in the larger partition.  One can continue recursively.  The average
  // running time is O(n), although the worst case is O(n^2).

  // One could also compute the mean below, which is simply Q2.
  // Below map has pairs (# of coverage of a base,
  //                       # times that coverage comes up)[contig]
  //   For example statCount[3][2] = 5 means that in contig 3 (contigs come 1st)
  //                                   a coverage of 2 came up 5 times.
  // We used a map since it's an associative container that sorts automatically.
  map<longlong,longlong> *statCount = new map<longlong,longlong>[numContigs];

  longlong tempSeqId, tempIntervalId, intervalSize;
  for (inter= intervals.begin(); inter != intervals.end(); ++inter) {
      tempSeqId = inter->SeqId();
      tempIntervalId = inter->IntervalId();

      // We are assuming that within an interval, the IntervalId() is constant.
      // This assumption is not used in PrintTextCoverage method below since
      // there they want one line of output per base.
      intervalSize = abs((inter->End()) - (inter->Begin()));
      if (tempIntervalId) {
          if (statCount[tempSeqId].find(tempIntervalId) !=
	      statCount[tempSeqId].end())
	         statCount[tempSeqId][tempIntervalId] += intervalSize;
          else
	         statCount[tempSeqId][tempIntervalId] = intervalSize;
      }
  }

  map<longlong,longlong>::iterator miter;

  os << "\n\n";
  os << "Contig     Mean     StdDev     Median     Q3/Q1\n";
  os << "------     ----     ------     ------     -----\n";
  for (longlong contig = 0; contig < numContigs; contig++) {
     longlong sumCover = 0, sumNums = 0;
     longlong Q1=0, Q2=0, Q3=0;
     bool noCompQ1 = false, noCompQ2 = false, noCompQ3 = false;
     char Q3divQ1[65];
     double Mean, StdDev = 0.0;
     longlong quartileIndex, tempNum;

     for (miter = statCount[contig].begin();
          miter != statCount[contig].end();
	  ++miter) {
        sumNums += miter->second;
        sumCover += (miter->first)*(miter->second);
     }
     int totalLength = (*m_pSeqLengths)[contig];
     Mean = ((double)sumCover)/totalLength;

     int tempVal;
     for (miter = statCount[contig].begin();
          miter != statCount[contig].end();
	  ++miter) {
        // This adds the non-zero values
        // tempVal is the non-zero value (# times covered)
        // miter->second is the number of times this value comes up
        tempVal = miter->first;
        StdDev += (miter->second)*(Mean - tempVal)*(Mean - tempVal);
     }
     // This adds the zero values, i.e. tempVal = 0;
     StdDev += (totalLength - sumNums)*(Mean*Mean);
     StdDev = StdDev/totalLength;
     StdDev = sqrt(StdDev);

     // This is the length of the contig number inter->SeqId()
     // Note that totalLength = sumNums + <# bases w/0 coverage>
     quartileIndex = totalLength/4; 

     // Below initially holds the number of zeros
     tempNum = totalLength - sumNums;

     if (tempNum > quartileIndex)
             noCompQ1 = true;
     if (tempNum > 2*quartileIndex)  // Will be true if above is true.
             noCompQ2 = true;
     if (tempNum > 3*quartileIndex)  // Will be true if above is true.
             noCompQ3 = true;

     miter = statCount[contig].begin();

     while (tempNum < quartileIndex && miter != statCount[contig].end() && !noCompQ1) {
        tempNum += miter->second;
        ++miter;
     }
     if (!noCompQ3) {
       --miter;   // back up to the previous one output then move forward
       Q1 = miter->first;
       ++miter;
     }
     while (tempNum < 2*quartileIndex && miter != statCount[contig].end() && !noCompQ2) {
        tempNum += miter->second;
	++miter;
     }
     if (!noCompQ2) {
       --miter;   // back up to the previous one output then move forward
       Q2 = miter->first;
       ++miter;
     }
     while (tempNum < 3*quartileIndex && miter != statCount[contig].end() && !noCompQ1) {
        tempNum += miter->second;
	++miter;
     }
     if (!noCompQ1) {
       --miter;
       Q3 = miter->first;
     }


     if (Q3 != 0 && Q1 == 0)
             sprintf(Q3divQ1, "infinity");
     else if (Q3 == 0 && Q1 == 0)
             sprintf(Q3divQ1, "0/0");
     else 
             sprintf(Q3divQ1, "%f", ((double)Q3)/Q1);


       os << "  " << contig << "       " << Mean << "   " << StdDev << "      " << Q2 << "       " << Q3divQ1 <<  "\n";
  }

  os << "\n\n\n";
}



void CoverageAnalyzer::PrintTextCoverage(ostream & os,
					 bool printLastBaseInContig,
					 bool printZeros) const {
  if (!m_pSeqLengths) printLastBaseInContig = false; //don't know length!
  vec<seq_interval> intervals;
  GetAllCoverages(intervals);

  vec<seq_interval>::iterator inter;
  for (inter= intervals.begin(); inter != intervals.end(); ++inter) {
    if (printZeros || inter->IntervalId()) { //coverage is not 0
      for (int i= inter->Begin(); i != inter->End(); ++i) {
	os << inter->SeqId() << "\t" << i << "\t" 
	   << inter->IntervalId() << "\n";
      }
    } 
    else if (!printZeros && printLastBaseInContig 
	     && inter->End() == (*m_pSeqLengths)[inter->SeqId()]) {
	os << inter->SeqId() << "\t" << inter->End() - 1 << "\t" 
	   << inter->IntervalId() << "\n";
    }      
  }
}

