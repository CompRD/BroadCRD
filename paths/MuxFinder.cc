// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#include "paths/MuxFinder.h"


void MuxFinder::StudyReads( ) {
  m_equalLengthGaplessReads = true;

  size_t pathId=0;
  // Find the first nonempty Fw read:
  while( pathId < m_pathsFw.size() && m_pathsFw[pathId].NSegments()==0 )
    ++pathId;
  if( pathId == m_pathsFw.size() )
    FatalErr("MuxFinder: no nonempty reads!  Surely not what you want.");

  int pathMinLength = m_pathsFw[pathId].MinLength();
  int pathMaxLength = m_pathsFw[pathId].MaxLength();

  if ( pathMinLength != pathMaxLength )
    m_equalLengthGaplessReads = false;

  // Check remaining Fw reads:
  for ( ++pathId; pathId < m_pathsFw.size() && m_equalLengthGaplessReads; ++pathId )
    if ( m_pathsFw[pathId].NSegments() != 0 &&
	 ( m_pathsFw[pathId].MinLength() != pathMinLength ||
	   m_pathsFw[pathId].MaxLength() != pathMaxLength ) )
      m_equalLengthGaplessReads = false;
  // Check all Rc reads:
  for ( pathId=0 ; pathId < m_pathsRc.size() && m_equalLengthGaplessReads; ++pathId )
    if ( m_pathsRc[pathId].NSegments() != 0 &&
	 ( m_pathsRc[pathId].MinLength() != pathMinLength ||
	   m_pathsRc[pathId].MaxLength() != pathMaxLength ) )
      m_equalLengthGaplessReads = false;
  
  if ( m_equalLengthGaplessReads && m_verbose )
  {
    cout << "All reads are of equal length and have no gaps. "
         << " No transitive dominance check needed." << endl;
  }
  
}



// Establishes which extension dominates the other, if any.
void MuxFinder::EstablishDominance( LeftExtension& extA, 
				    LeftExtension& extB )
{
  KmerPathLoc leftScanA = extA.GetLoc(), rightScanA = extA.GetLoc();
  KmerPathLoc leftScanB = extB.GetLoc(), rightScanB = extB.GetLoc();

  // If they disagree at all, neither can dominate the other.


  bool leftScanGood = IsPerfectMatchToLeft( leftScanA, leftScanB );

  if ( ! leftScanGood )
    return;

  bool rightScanGood = IsPerfectMatchToRight( rightScanA, rightScanB );

  if ( ! rightScanGood )
    return;

  // If their overlap is less than the minimum, neither dominates.
  if ( KmersInInterval( leftScanA, rightScanA ) < m_minOverlap )
    return;

  // If extA ended before extB to the left, it dominates extB.
  // 
  // A       --------------...
  // B    -----------------...
  if ( ! leftScanB.atBegin() )
    extB.SetIsDominated( true );

  // Otherwise, if extB path ended before extA to the left, it
  // dominates extA.
  // 
  // A    -----------------...
  // B       --------------...
  else if ( ! leftScanA.atBegin() )
    extA.SetIsDominated( true );

  // If extA and extB end at the same time to the left, whichever one
  // ends first to the right is submissive.
  //
  // dom   --------------------
  // sub   -----------------
  
  // (We know that leftScanA == leftScanA.GetPath().Begin() and
  // leftScanB == leftScanB.GetPath().Begin() after previous ifs.)

  else if ( ! rightScanB.atEnd() )
    extA.SetIsDominated( true );

  else if ( ! rightScanA.atEnd() )
    extB.SetIsDominated( true );

  // If they're coterminal on both ends, the one with the lowest
  // pathId is considered dominant.
  else 
    if ( extA.GetPathId() < extB.GetPathId() )
      extB.SetIsDominated( true );
    else if ( extB.GetPathId() < extA.GetPathId() )
      extA.SetIsDominated( true );

}



// Establishes which of many extensions are not dominated by any other.
int MuxFinder::EstablishDominance( vec<LeftExtension>& extensions,
				   longlong& dominanceCounter1,
				   TaskTimer& dominanceTimer1, 
				   longlong& dominanceCounter2,
				   TaskTimer& dominanceTimer2 )
{
  dominanceTimer1.Start();
  
  // First we try to assert the dominance of each extension over
  // each of the extensions of the same number of segments or more.
  for ( unsigned int domIdx = 0; domIdx < extensions.size(); ++domIdx )
  {
    LeftExtension& dominant = extensions[domIdx];
    
    // If this extension is already dominated by some other
    // extension, then any extension that it might dominate has
    // probably been marked as dominated by that other extension.
    // We'll check for exceptions to this in a separate pass.
    if ( dominant.IsDominated() )
      continue;
    
    for ( unsigned int subIdx = domIdx+1; subIdx < extensions.size(); ++subIdx )
    {
      LeftExtension& submissive = extensions[subIdx];
      
      // If this has already been dominated, don't check it again.
      if ( submissive.IsDominated() )
        continue;

      EstablishDominance( dominant, submissive );
      
      // If we determine that the would-be dominator is
      // actually dominated, stop checking.
      if ( dominant.IsDominated() )
        break;
    }
  }
  
  dominanceTimer1.Stop();
  
  int numDominated1 = count_if( extensions.begin(), extensions.end(),
                                mem_fun_ref( &LeftExtension::IsDominated ) );
  
  dominanceCounter1 += numDominated1;

  // If all reads are of equal length and have no gaps, then dominance
  // is transitive, i.e. if A dominates B, and B dominates C, then A
  // dominates C, so we won't encounter the situation described below.
  if ( m_equalLengthGaplessReads )
    return ( extensions.size() - numDominated1 );

  dominanceTimer2.Start();
  
  // Now we check for once-or-more-removed dominance.  This can occur
  // like so: consider read A with three extensions B, C, and D.
  //
  // A          ------------------
  // B        -----------------------~~~~
  // C      -------------------
  // D    -------------------------+++++
  //
  // Clearly, B dominates C and C dominates D, but since B and D don't
  // match to the right, B does not dominate D.  In the first pass, C
  // will be flagged as dominated, and so we won't check if C
  // dominates D.  Thus, we have to do a second pass where we check
  // each remaining undominated extension against all the dominated
  // extensions of equal or lesser number of segments.
  //
  // We don't need to check against non-dominated extensions, since we
  // know that in the loop above each surviving extension was checked
  // against all of the surviving extensions of equal or lesser number
  // of segments.
  for ( unsigned int currIdx = 0; currIdx < extensions.size(); ++currIdx )
  {
    LeftExtension& current = extensions[currIdx];
    
    // Only check surviving extensions.
    if ( current.IsDominated() )
      continue;
    
    for ( unsigned int domIdx = 0; domIdx < extensions.size(); ++domIdx )
    {
      // It can't dominate itself.
      if ( domIdx == currIdx )
        continue;
      
      LeftExtension& dominant = extensions[domIdx];
      
      // If this one has more segments than the current one (and
      // therefore all the ones after it do), it can't dominate.
      if ( dominant.GetSegment() > current.GetSegment() )
        break;
      
      // If this extension survived the first round, then since
      // the current one also did, we know that we've already
      // tried to establish the dominance of this one over the
      // current one and failed, so don't try again.
      if ( ! dominant.IsDominated() )
        continue;
      
      EstablishDominance( dominant, current );
      
      // If we found a dominator, move on.
      if ( current.IsDominated() )
        break;
    }
  }
  
  dominanceTimer2.Stop();

  int numDominated2 = count_if( extensions.begin(), extensions.end(),
                                mem_fun_ref( &LeftExtension::IsDominated ) );
    
  dominanceCounter2 += numDominated2 - numDominated1;

  int numFoundMuxes = extensions.size() - numDominated2;

  return numFoundMuxes;
}


longlong MuxFinder::FindMuxes( MuxGraph& theMuxGraph, const int partition )
{
  return this->FindMuxes( theMuxGraph, 1000000, 0, m_pathsDb.size() - 1, partition );
}


longlong MuxFinder::FindMuxes( MuxGraph& theMuxGraph,
			       const longlong chunkSize,
			       const longlong start,
			       const longlong stop,
                               const int partition )
{
  ForceAssertLe( start, stop );

  longlong numExtensionsThisChunk = 0;
  longlong numPaths = 0;
  longlong numMuxes = 0;

  histogram<int> numMuxesHisto;
  numMuxesHisto.HideUnderflow();
  numMuxesHisto.AddLinearProgressionOfBins(0,1,10); // 0,1,...,9
  numMuxesHisto.AddLinearProgressionOfBins(10,10,9); // 10,20,...,90
  numMuxesHisto.AddLogProgressionOfBins(100,10000,1,10,vec<int>(1,1));//10^{3-5}

  const longlong numPasses = ( m_pathsDb.size() + chunkSize - 1 ) / chunkSize;

  if ( m_verbose )
    cout << Date() << ": finding muxes for " << 2 * m_pathsFw.size() << " paths" 
	 << " in " << numPasses << " passes." << endl;

  TaskTimer rangeTimer;

  TaskTimer coterminalFullTimer, 
    coterminalFindTimer, coterminalSortTimer, 
    coterminalDominanceTimer1, coterminalDominanceTimer2;

  TaskTimer properFullTimer, 
    properFindTimer, properSortTimer, 
    properDominanceTimer1, properDominanceTimer2;

  longlong coterminalDominanceCounter1 = 0, coterminalDominanceCounter2 = 0;
  longlong properDominanceCounter1 = 0, properDominanceCounter2 = 0;

  vec<LeftExtension> extensions;

  longlong cachedStart = -1;
  longlong cachedUpper = 0;
  longlong cachedLower = 0;

  longlong chunkNum = 1;
  longlong chunkStart = start;
  longlong chunkStop  = min( start+chunkSize-1, stop );

  //for ( longlong current = 0; current < m_pathsDb.size(); ++current )
  for ( longlong current = start; current <= stop; ++current )
  {
    if ( current == chunkStop )
    {
      if ( m_verbose )
	Dot( cout, chunkNum-1 );

      if ( mp_log )
      {
        *mp_log << Date() << ": pass " << chunkNum++ << endl;
        *mp_log << "  Entries " << chunkStart << " to " << chunkStop << ": " 
             << numPaths << " paths, "
             << numMuxes << " minimal extensions so far" << endl;
        *mp_log << "  Extensions for reads in this chunk: " << numExtensionsThisChunk << endl;
        *mp_log << "  Finding ranges:           " << rangeTimer << endl;
        *mp_log << "  Checking coterminal exts: " << coterminalFullTimer << endl;
        *mp_log << "    Finding extensions:       " << coterminalFindTimer << endl;
        *mp_log << "    Sorting extensions:       " << coterminalSortTimer << endl;
        *mp_log << "    Dominance (first pass):   " << coterminalDominanceTimer1 << endl;
        *mp_log << "      Removed " << coterminalDominanceCounter1 << " extensions." << endl;
        *mp_log << "    Dominance (second pass):  " << coterminalDominanceTimer2 << endl;
        *mp_log << "      Removed " << coterminalDominanceCounter2 << " extensions." << endl;
        *mp_log << "  Checking proper exts:     " << properFullTimer << endl;
        *mp_log << "    Finding extensions:       " << properFindTimer << endl;
        *mp_log << "    Sorting extensions:       " << properSortTimer << endl;
        *mp_log << "    Dominance (first pass):   " << properDominanceTimer1 << endl;
        *mp_log << "      Removed " << properDominanceCounter1 << " extensions." << endl;
        *mp_log << "    Dominance (second pass):  " << properDominanceTimer2 << endl;
        *mp_log << "      Removed " << properDominanceCounter2 << " extensions." << endl;
      }

      numExtensionsThisChunk = 0;
      coterminalDominanceCounter1 = 0;
      coterminalDominanceCounter2 = 0;
      properDominanceCounter1 = 0;
      properDominanceCounter2 = 0;
      
      rangeTimer.Reset();

      coterminalFullTimer.Reset();
      coterminalFindTimer.Reset();
      coterminalSortTimer.Reset();
      coterminalDominanceTimer1.Reset();
      coterminalDominanceTimer2.Reset();

      properFullTimer.Reset();
      properFindTimer.Reset();
      properSortTimer.Reset();
      properDominanceTimer1.Reset();
      properDominanceTimer2.Reset();

      chunkStart = chunkStop+1;
      chunkStop = min( chunkStart+chunkSize-1, stop );
    }

    // Look for entries in the paths database that correspond to
    // initial segments of paths.
    if ( m_pathsDb[current].PathPos() != 0 )
      continue;

    ++numPaths;
             
    OrientedKmerPathId thisOkpi( m_pathsDb[current].PathId() );
    const KmerPath* p_thisPath = thisOkpi.GetPathPtr( m_pathsFw, m_pathsRc );
    const int thisPathNumSegs = p_thisPath->NSegments();
    
    // If this read is too short, skip it (if it is a walking read).
    if ( p_thisPath->KmerCount() < m_minNumKmers && thisOkpi.GetId() >= partition )
      continue;

    //cout << current << ": " << *p_thisPath << " " << flush;
    
    // Find the range of entries that could intersect this segment.
    
    rangeTimer.Start();
    
    // The lower bound is easy.
    longlong from = current - m_pathsDb[current].Lookback();
    
    if ( cachedStart != m_pathsDb[current].Start() )
    {
      cachedStart = m_pathsDb[current].Start();
      cachedLower = distance( m_pathsDb.Begin(),
                              lower_bound( m_pathsDb.Begin(), m_pathsDb.End(), 
                                           m_pathsDb[current] ) );
      cachedUpper = distance( m_pathsDb.Begin(),
                              upper_bound( m_pathsDb.Begin(), m_pathsDb.End(),
                                           m_pathsDb[current] ) );
      
      // PRINT3( cachedStart, cachedLower, cachedUpper );
    }

    longlong to = cachedUpper;
   
    rangeTimer.Stop();
    
    //cout << "searching " << ( to - from ) << endl;
    
    extensions.clear();
    extensions.reserve( to - from );

    // We do two passes: In the first pass, we look only for
    // coterminal extensions.  If we find one, we're done.  In the
    // second pass, assuming we haven't found a coterminal extension,
    // we look at all the others to find the minimal one.

    // Walk through these entries, looking for entries that overlap the current one.
    for ( int pass = 0; pass < 2; ++pass )
    {
      bool coterminalOnly = ( pass == 0 );

      TaskTimer& fullTimer = ( coterminalOnly ? coterminalFullTimer : properFullTimer );
      TaskTimer& findTimer = ( coterminalOnly ? coterminalFindTimer : properFindTimer );
      TaskTimer& sortTimer = ( coterminalOnly ? coterminalSortTimer : properSortTimer );
      TaskTimer& dominanceTimer1
        = ( coterminalOnly ? coterminalDominanceTimer1 : properDominanceTimer1 );
      TaskTimer& dominanceTimer2
        = ( coterminalOnly ? coterminalDominanceTimer2 : properDominanceTimer2 );
      longlong& dominanceCounter1 
        = ( coterminalOnly ? coterminalDominanceCounter1 : properDominanceCounter1 );
      longlong& dominanceCounter2 
        = ( coterminalOnly ? coterminalDominanceCounter2 : properDominanceCounter2 );

      fullTimer.Start();

      extensions.clear();

      int passSpecificFrom, passSpecificTo;

      if ( coterminalOnly )
      {
        passSpecificFrom = cachedLower;
        passSpecificTo   = to;
      }
      else
      {
        passSpecificFrom = from;
        passSpecificTo   = to;
      }
 
      //PRINT3( pass, passSpecificFrom, passSpecificTo );

      findTimer.Start();

      for ( longlong other = passSpecificFrom; other < passSpecificTo; ++other )
      {
        ForceAssertLe( m_pathsDb[other].Start(), m_pathsDb[current].Start() );

	// Skip segments which do not even intersect the current segment
	if( m_pathsDb[other].Stop() < m_pathsDb[current].Start() )
	  continue;

        if ( coterminalOnly )
        {
          ForceAssertEq( m_pathsDb[other].Start(), m_pathsDb[current].Start() );
          if ( m_pathsDb[other].PathPos() > 0 )
            continue;
        }
        else
        {
          ForceAssertLe( m_pathsDb[other].Start(), m_pathsDb[current].Start() );
          if ( other >= cachedLower && m_pathsDb[other].PathPos() == 0 )
            continue;
        }

        // Eliminate segments that extend the current segment to the
        // right when the current path has more than one segment
        if ( thisPathNumSegs > 1 && m_pathsDb[other].Stop() > m_pathsDb[current].Stop() )
          continue;
      
        OrientedKmerPathId otherOkpi( m_pathsDb[other].PathId() );
        
        // Do not count alignments of a read to itself.
        if ( otherOkpi.GetId() == thisOkpi.GetId() )
          continue;
        
        // Now we make sure they agree perfectly until one path or the other ends.
        int thisSeg = 0;
        
        const KmerPath* p_otherPath = otherOkpi.GetPathPtr( m_pathsFw, m_pathsRc );
        int otherSeg = m_pathsDb[other].PathPos();

        // Eliminate too-short paths.
        if ( p_otherPath->KmerCount() < m_minNumKmers )
          continue;
        
        // Eliminate paths that obviously don't match to the right.
        if ( otherSeg < p_otherPath->NSegments()-1 &&
             m_pathsDb[other].Stop() < m_pathsDb[current].Stop() )
          continue;
        
        KmerPathLoc firstMatchingKmer( p_otherPath, otherSeg );
        firstMatchingKmer.SetKmer( m_pathsDb[current].Start() );
        
        KmerPathLoc thisRightScan( p_thisPath->Begin() );
        KmerPathLoc otherRightScan( firstMatchingKmer );
        
        //PRINT( *p_thisPath );
        //PRINT( *p_otherPath );

        // We have an extension if the two paths match perfectly to
        // the right AND (in pass 1) it is coterminal with the current
        // path on the left and coterminal with or shorter than the
        // current path on the right or (in pass 2) the other path
        // extends the current path to the left.
        if ( IsPerfectMatchToRight( thisRightScan, otherRightScan ) )
        {
          if ( coterminalOnly )
          {
            if ( ! otherRightScan.atEnd() )
              continue;
            if ( thisRightScan.atEnd() && otherOkpi < thisOkpi )
              continue;
          }

          // If the overlap is insufficient, don't consider this an
          // extension, unless the reads are partitioned and this is
          // an end read, in which case we'll consider anything that
          // entirely subsumes the read to be an extension.
          if ( KmersInInterval( firstMatchingKmer, otherRightScan ) < m_minOverlap &&
               ! ( thisOkpi.GetId() < partition && thisRightScan.atEnd() ) )
            continue;

          extensions.push_back( LeftExtension( otherOkpi, firstMatchingKmer ) );

          continue;
          
          // CHECK: The extended read should dominate the extending read.
          
          LeftExtension trivialExt( thisOkpi, p_thisPath->Begin() );
          LeftExtension actualExt( otherOkpi, firstMatchingKmer );
          
          EstablishDominance( trivialExt, actualExt );
          ForceAssert( ! trivialExt.IsDominated() );
          ForceAssert( actualExt.IsDominated() );
        }
      }
    
      findTimer.Stop();
      //cout << "found " << extensions.size() << " extensions" << endl;
      
      if ( extensions.empty() )
      {
        fullTimer.Stop();

        continue;
      }
    
      numExtensionsThisChunk += extensions.size();
    
      sortTimer.Start();
      
      // Sort the extensions by the number of segments involved.
      sort( extensions.begin(), extensions.end() );
      
      sortTimer.Stop();

      /*
      cout << "Checking dominance among " << extensions.size() << " extensions"
      //   << " of " << *p_thisPath 
           << endl;
      */
      
      int numFoundMuxes = EstablishDominance( extensions,
                                              dominanceCounter1, dominanceTimer1,
                                              dominanceCounter2, dominanceTimer2 );
    
      //cout << "found " << numFoundMuxes << " minimal extensions" << endl;
    
      ForceAssert( numFoundMuxes > 0 );
  
      numMuxes += numFoundMuxes;
      numMuxesHisto.AddDatum(numFoundMuxes);
      
      
      vec<Mux> newMuxes;
      newMuxes.reserve( numFoundMuxes );
      for ( unsigned int extIdx = 0; extIdx < extensions.size(); ++extIdx )
        if ( ! extensions[extIdx].IsDominated() )
          newMuxes.push_back( extensions[extIdx].ToMux() );
      
      theMuxGraph.SetMuxesOf( thisOkpi, newMuxes );
      
      if ( m_printHiMuxCount && numFoundMuxes >= 100 )
      {
        int highestSeg = 0;
        for ( unsigned int extIdx = 0; extIdx < extensions.size(); ++extIdx )
          if ( ! extensions[extIdx].IsDominated() )
            highestSeg = max( highestSeg, extensions[extIdx].GetSegment() );
        
        cout << "base ";
        for ( int i = 0; i < highestSeg; ++i )
          cout << KmerPathInterval::Blank();
        for ( int i = 0; i < p_thisPath->NSegments(); ++i )
          cout << p_thisPath->Segment(i);
        cout << endl;
        
        for ( unsigned int extIdx = 0; extIdx < extensions.size(); ++extIdx )
          if ( ! extensions[extIdx].IsDominated() )
          {
            cout << "mux  ";
            for ( int i = highestSeg - extensions[extIdx].GetSegment(); i > 0; --i )
              cout << KmerPathInterval::Blank();
            const KmerPath* p_otherPath = extensions[extIdx].GetLoc().GetPathPtr();
            for ( int i = 0; i < p_otherPath->NSegments(); ++i )
              cout << p_otherPath->Segment(i);
            cout << endl;
          }
        
        cout << endl;
      }


      fullTimer.Stop();

      // If we've gotten here, then we've found the one minimal
      // coterminal extension, so we can stop.
      if ( coterminalOnly )
        break;
    }

    if ( extensions.empty() )
      numMuxesHisto.AddDatum( 0 );
  }
  
  if ( m_verbose ) {
    cout << endl;
    cout << "Histogram of # muxes per insert:\n" << numMuxesHisto << endl;
  }
  
  return numMuxes;
}

