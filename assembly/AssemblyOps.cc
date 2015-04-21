// Copyright (c) 2004 The Broad Institute at MIT and Harvard

#include "assembly/AssemblyOps.h"
#include "CoverageAnalyzer.h"
#include "Equiv.h"

#include <map>
//------------------------------------------------------------

ContigSplit::ContigSplit( Assembly *p_assembly,
                          const Contig &theContig,
                          const vec< vec<ReadLocation> > &readLocSets )
  : AssemblyOp( p_assembly ),
    m_contig( theContig ),
    m_locSets( readLocSets )
{
  // Validate the inputs.

  if ( ! m_contig.IsValid() )
  {
    cout << "Invalid Contig given to ContigSplit." << endl;
    TracebackThisProcess();
  }
  
  bool badContigFound = false;

  vec< vec<ReadLocation> >::iterator locSetIter;
  vec<ReadLocation>::iterator locIter;
  for ( locSetIter = m_locSets.begin();
        locSetIter != m_locSets.end();
        ++locSetIter )
    for ( locIter = locSetIter->begin();
          locIter != locSetIter->end();
          ++locIter )
      if ( locIter->GetContig() != m_contig )
      {
        if ( ! badContigFound )
        {
          cout << "Location from bad contig in ContigSplit:" << endl;
          cout << "Expected: " << m_contig << endl;
        }
        cout << "Found: " << *locIter << endl;
      }
  
  if ( badContigFound )
    TracebackThisProcess();
}

void ContigSplit::Execute( ostream *pLog )
{
  vec<ContigLocation> contigLocs;
  m_contig.GetSelfLocations( contigLocs );

  basevector contigBases;
  qualvector contigQuals;

  if ( pLog )
    *pLog << "Splitting " << m_contig << " into: ";

  vec< vec<ReadLocation> >::iterator locSetIter;

  for ( locSetIter = m_locSets.begin();
        locSetIter != m_locSets.end();
        ++locSetIter )
  {
    // If the location set is empty, create an invalid resulting contig.
    if ( locSetIter->empty() )
    {
      m_resultingContigs.push_back( Contig() );
      continue;
    }

    // Find the extent of the read locations on the contig.
    int minReadBegin = locSetIter->front().Begin();
    int maxReadEnd   = locSetIter->front().End();

    vec<ReadLocation>::iterator locIter;
    for ( locIter = locSetIter->begin();
          locIter != locSetIter->end();
          ++locIter )
    {
      minReadBegin = min( minReadBegin, locIter->Begin() );
      maxReadEnd   = max( maxReadEnd,   locIter->End() );
    }

    minReadBegin = max( minReadBegin, 0 );
    maxReadEnd   = min( maxReadEnd,   m_contig.GetLength() );
    if ( minReadBegin >= maxReadEnd )
    {
      m_resultingContigs.push_back( Contig() );
      continue;
    }

    // Copy the bases and quals from the original contig.
    contigBases.SetToSubOf( m_contig.GetBases(), minReadBegin, maxReadEnd - minReadBegin );
    contigQuals.SetToSubOf( m_contig.GetQuals(), minReadBegin, maxReadEnd - minReadBegin );

    // Create the new contig and save it for later retrieval.
    Contig resultingContig = mp_assembly->NewContig( contigBases, contigQuals );
    m_resultingContigs.push_back( resultingContig );

    if ( pLog )
      *pLog << resultingContig
            << " (" << locSetIter->size() << " reads,"
            << " " << resultingContig.GetLength() << " bases)  ";

    // Copy the read locations to the new contig, offset by the correct amount.
    for ( locIter = locSetIter->begin();
          locIter != locSetIter->end();
          ++locIter )
    {
      ReadLocation newReadLoc( resultingContig,
                               locIter->GetRead(),
                               locIter->GetInterval().CopyShiftedBy( -minReadBegin),
                               locIter->GetOrientation() );
      mp_assembly->AddReadLocation( newReadLoc );
    }
 
    // Place the new contig in the old contig's supers, offset by the correct amount.
      vec<ContigLocation>::iterator contigLocIter;
      for ( contigLocIter = contigLocs.begin();
          contigLocIter != contigLocs.end();
          ++contigLocIter )
    {
      int resultingContigBegin = contigLocIter->GetInterval().Begin() + minReadBegin;
      int resultingContigEnd   = resultingContigBegin + resultingContig.GetLength();

      ContigLocation newContigLoc( contigLocIter->GetSuper(),
                                   resultingContig,
                                   Interval( resultingContigBegin,
                                             resultingContigEnd ),
                                   contigLocIter->GetOrientation() );
      mp_assembly->AddContigLocation( newContigLoc );
    }
  }

  if ( pLog )
    *pLog << endl;

  // Remove the contig's read locations and sequence from the
  // assembly.  This will signal the Assembly to purge it before it
  // saves itself to disk.
  mp_assembly->RemoveContig( m_contig );
}


//------------------------------------------------------------

SuperSplit::SuperSplit( Assembly *p_assembly,
                        const Super &theSuper,
                        const vec< vec<ContigLocation> > &contigLocSets )
  : AssemblyOp( p_assembly ),
    m_super( theSuper ),
    m_locSets( contigLocSets )
{
  // Validate the inputs.

  if ( ! m_super.IsValid() )
  {
    cout << "Invalid Super given to SuperSplit." << endl;
    TracebackThisProcess();
  }
  
  bool badSuperFound = false;

  vec< vec<ContigLocation> >::iterator locSetIter;
  vec<ContigLocation>::iterator locIter;
  for ( locSetIter = m_locSets.begin();
        locSetIter != m_locSets.end();
        ++locSetIter )
    for ( locIter = locSetIter->begin();
          locIter != locSetIter->end();
          ++locIter )
      if ( locIter->GetSuper() != m_super )
      {
        if ( ! badSuperFound )
        {
          cout << "Location from bad super in SuperSplit:" << endl;
          cout << "Expected: " << m_super << endl;
        }
        cout << "Found: " << *locIter << endl;
      }
  
  if ( badSuperFound )
    TracebackThisProcess();
}

void SuperSplit::Execute( ostream *pLog )
{
  // Remove the supercontig's contig locations.  The assembly will
  // purge the empty supercontig before saving itself to disk.
  mp_assembly->Clear( m_super );

  if ( pLog )
    *pLog << "Cutting " << m_super << " into: ";

  vec< vec<ContigLocation> >::iterator locSetIter;
  vec<ContigLocation>::iterator locIter;

  for ( locSetIter = m_locSets.begin();
        locSetIter != m_locSets.end();
        ++locSetIter )
  {
    // Create an invalid super for empty sub-collections.
    if ( locSetIter->empty() )
    {
      m_resultingSupers.push_back( Super() );
      continue;
    }

    // Find the start of the leftmost contig.
    int minContigBegin = locSetIter->front().Begin();

    for ( locIter = locSetIter->begin();
          locIter != locSetIter->end();
          ++locIter )
      minContigBegin = min( minContigBegin, locIter->Begin() );

    // Create a new supercontig and copy it for later retrieval.
    Super resultingSuper = mp_assembly->NewSuper();
    m_resultingSupers.push_back( resultingSuper );

    if ( pLog )
      *pLog << resultingSuper << " (" << locSetIter->size() << " contigs)  ";

    // Copy the sub-collection of contig locations to the new super,
    // offset such that the super is normalized (i.e., the leftmost
    // contig starts at zero).
    for ( locIter = locSetIter->begin();
          locIter != locSetIter->end();
          ++locIter )
    {
      ContigLocation newContigLoc( resultingSuper,
                                   locIter->GetContig(),
                                   locIter->GetInterval().CopyShiftedBy( -minContigBegin ),
                                   locIter->GetOrientation() );
      mp_assembly->AddContigLocation( newContigLoc );
    }
  }

  if ( pLog )
    *pLog << endl;
}


//------------------------------------------------------------

ContigSlice::ContigSlice( Assembly *p_assembly, 
                          const Contig &theContig,
                          const int excisionStart,
                          const int excisionEnd )
  : ContigSplit( p_assembly, theContig, vec< vec<ReadLocation> >( 2 ) )
{
  ForceAssertLe( excisionStart, excisionEnd );

  vec<ReadLocation> readLocs;
  m_contig.GetReadLocations( readLocs );

  // Calculate the set of read locations to be put in the "left" and
  // "right" sets.
  vec<ReadLocation>::iterator readLocIter;
  for ( readLocIter = readLocs.begin();
        readLocIter != readLocs.end();
        ++readLocIter )
  {
    if ( readLocIter->End() < excisionStart ||
         readLocIter->Begin() <= excisionStart &&
         readLocIter->IsReverse() )
      m_locSets[0].push_back( *readLocIter );
    
    else if ( readLocIter->Begin() > excisionEnd ||
              readLocIter->End() >= excisionEnd &&
              readLocIter->IsForward() )
      m_locSets[1].push_back( *readLocIter );
  }
}


//------------------------------------------------------------

ContigTrim::ContigTrim( Assembly *p_assembly, 
                        const Contig &theContig,
                        const int leftTrim,
                        const int rightTrim )
  : ContigSplit( p_assembly, theContig, vec< vec<ReadLocation> >( 3 ) )
{
  vec<ReadLocation> readLocs;
  m_contig.GetReadLocations( readLocs );

  int firstBaseSaved = leftTrim;
  int lastBaseSaved = m_contig.GetLength() - rightTrim;

  // Calculate the set of read locations to be put in the "to be
  // saved" set.
  vec<ReadLocation>::iterator readLocIter;
  for ( readLocIter = readLocs.begin();
        readLocIter != readLocs.end();
        ++readLocIter )
    if ( readLocIter->Begin() < firstBaseSaved )
      m_locSets[1].push_back( *readLocIter );

    else if ( readLocIter->End() > lastBaseSaved )
      m_locSets[2].push_back( *readLocIter );

    else
      m_locSets[0].push_back( *readLocIter );
}


//------------------------------------------------------------

ContigSanitize::ContigSanitize( Assembly *p_assembly, 
				const Contig &theContig,
                                const int endForgiveness )
  : ContigSplit( p_assembly, theContig, vec< vec<ReadLocation> >( 0 ) )
{
  vec<ReadLocation> readLocs;
  m_contig.GetReadLocations( readLocs );

  vector<seq_interval> readSeqIntervals;
  readSeqIntervals.reserve( readLocs.size() );

  // figure out the read coverage of the contig
  int readIndex = 0;
  vec<ReadLocation>::iterator readLocIter;
  for ( readLocIter = readLocs.begin();
        readLocIter != readLocs.end();
        ++readLocIter )
  {
    seq_interval readInterval( readIndex,
			       0,
			       readLocIter->Begin(),
			       readLocIter->End() );
    readSeqIntervals.push_back( readInterval );
    readIndex++;
  }

  vec<int> contigLength( 1, theContig.GetLength() );
  CoverageAnalyzer theAnalyzer( readSeqIntervals, &contigLength );

  vec<seq_interval> regionsCovered;
  theAnalyzer.GetCoveragesAtLeast(1, regionsCovered);

  m_locSets.resize( regionsCovered.size() );

  // loop through the covered regions
  int regionIndex = 0;
  vec<seq_interval>::iterator intervalIter;
  for ( intervalIter = regionsCovered.begin();
        intervalIter != regionsCovered.end();
        ++intervalIter )
  {
    int regionBegin = intervalIter->Begin();
    int regionEnd   = intervalIter->End();

    if ( regionBegin <= endForgiveness )
      regionBegin = 0;
    if ( m_contig.GetLength() - regionEnd <= endForgiveness )
      regionEnd = m_contig.GetLength();
    
    // extract the ReadLocations in the region
    for ( readLocIter = readLocs.begin();
          readLocIter != readLocs.end();
          ++readLocIter )
    {
      // We need to be a little forgiving in order to capture reads
      // that hang off the end of the contig.
      if ( readLocIter->End() > regionBegin &&
           readLocIter->Begin() < regionEnd )
      {
        m_locSets[regionIndex].push_back( *readLocIter );
      }
    }

    regionIndex++;
  }

  // If there is one contiguous coverage region that covers the
  // entirety of the contig and all reads fall into it, then there's
  // nothing to do.
  if ( regionsCovered.size() == 1 &&
       regionsCovered.front().Begin() <= endForgiveness &&
       regionsCovered.front().End()   >= m_contig.GetLength() - endForgiveness &&
       m_locSets.front().isize() == theContig.GetNumReadLocations() )
      m_needNewContigs = false;
  else
    m_needNewContigs = true;
}


void
ContigSanitize::Execute( ostream *pLog )
{
  if ( m_needNewContigs )
    this->ContigSplit::Execute( pLog );
}


//------------------------------------------------------------

SuperSlice::SuperSlice( Assembly *p_assembly,
                        const Super &theSuper,
                        const int excisionStart,
                        const int excisionEnd )
  : SuperSplit( p_assembly, theSuper, vec< vec<ContigLocation> >( 2 ) )
{
  ForceAssertLe( excisionStart, excisionEnd );

  vec<ContigLocation> allContigLocs;
  m_super.GetContigLocations( allContigLocs );
  
  // Calculate the set of contig locations to be put in the "left" and
  // "right" sets.
  for ( vec<ContigLocation>::iterator contigLocIter = allContigLocs.begin();
        contigLocIter != allContigLocs.end();
        ++contigLocIter )
  {
    if ( contigLocIter->End() <= excisionStart )
      //        start     end
      //  ------  |        |
      //   tig
      m_locSets[0].push_back( *contigLocIter );
    
    else if ( contigLocIter->Begin() >= excisionEnd )
      //        start     end
      //          |        |    -----
      m_locSets[1].push_back( *contigLocIter );
    
    else if ( contigLocIter->Begin() > excisionStart &&
              contigLocIter->End()   < excisionEnd )
      //        start     end
      //          | -----  |
      m_contigsToRemove.push_back( contigLocIter->GetContig() );
    
    else
      // one of:
      //        start     end
      //        --|---     |
      //          |     ---|--
      //        --|--------|--
    {
      // Create a ContigSlice object to remove the portion of the contig in
      // the excision region.
      m_contigSlices.push_back( ContigSlice( p_assembly,
                                             contigLocIter->GetContig(),
                                             excisionStart - contigLocIter->Begin(),
                                             excisionEnd   - contigLocIter->Begin() ) );
    }
  }
}

void SuperSlice::Execute( ostream *pLog )
{
  // Execute the ContigSlices calculated above.

  for ( vec<ContigSlice>::iterator excisionIter = m_contigSlices.begin();
        excisionIter != m_contigSlices.end();
        ++excisionIter )
  {
    excisionIter->Execute( pLog );

    vec<Contig> newContigs;
    excisionIter->GetResultingContigs( newContigs );

    // Find the contig location of each new valid contig that places
    // it in this super and copy it to the appropriate set.
    vec<ContigLocation> newContigLocs;
    for ( unsigned int side = 0; side < 2; ++side )
      if ( newContigs[side].IsValid() )
      {
        newContigs[side].GetSelfLocations( newContigLocs );
        for ( vec<ContigLocation>::iterator newContigLocIter = newContigLocs.begin();
              newContigLocIter != newContigLocs.end();
              ++newContigLocIter )
          if ( newContigLocIter->GetSuper() == m_super )
            m_locSets[side].push_back( *newContigLocIter );
      }
  }

  for ( vec<Contig>::iterator contigIter = m_contigsToRemove.begin();
        contigIter != m_contigsToRemove.end();
        ++contigIter )
    mp_assembly->RemoveContig( *contigIter );

  // Now that everything is neatly in whole contigs, perform the split.
  this->SuperSplit::Execute( pLog );
}


//------------------------------------------------------------

SuperTrim::SuperTrim( Assembly *p_assembly,
                      const Super &theSuper,
                      const int leftTrim,
                      const int rightTrim )
  : SuperSplit( p_assembly, theSuper, vec< vec<ContigLocation> >( 3 ) )
{
  vec<ContigLocation> allContigLocs;
  m_super.GetContigLocations( allContigLocs );

  int firstBaseSaved = leftTrim;
  int lastBaseSaved = m_super.GetLength() - rightTrim;
  
  // Calculate the set of contig locations to be put in the "to be
  // saved" set.
  for ( vec<ContigLocation>::iterator contigLocIter = allContigLocs.begin();
        contigLocIter != allContigLocs.end();
        ++contigLocIter )
  {
    if ( contigLocIter->End() <= firstBaseSaved )
      //          L        R
      //  -----   |        |
      //   tig
      m_locSets[1].push_back( *contigLocIter );
    
    else if ( contigLocIter->Begin() >= lastBaseSaved )
      //          L        R
      //          |        |    -----
      m_locSets[2].push_back( *contigLocIter );
    
    else if ( contigLocIter->Begin() > firstBaseSaved &&
              contigLocIter->End()   < lastBaseSaved )
      //          L        R
      //          | -----  |
      m_locSets[0].push_back( *contigLocIter );
    
    else
      // one of:
      //          L        R
      //        --|---     |
      //          |     ---|--
      //        --|--------|--
    {
      // Create a ContigTrim object to remove the portion of the contig in
      // the excision region.
      m_contigTrims.push_back( ContigTrim( p_assembly,
                                           contigLocIter->GetContig(),
                                           leftTrim  - contigLocIter->Begin(),
                                           rightTrim - contigLocIter->Begin() ) );
    }
  }
}

void SuperTrim::Execute( ostream *pLog )
{
  // Execute the ContigTrims calculated above.

  for ( vec<ContigTrim>::iterator trimIter = m_contigTrims.begin();
        trimIter != m_contigTrims.end();
        ++trimIter )
  {
    trimIter->Execute( pLog );

    vec<Contig> newContigs;
    trimIter->GetResultingContigs( newContigs );

    // Find the contig location of each new valid contig that places
    // it in this super and copy it to the appropriate set.
    vec<ContigLocation> newContigLocs;
    for ( int part = 0; part < 3; ++part )
      if ( newContigs[part].IsValid() )
      {
        newContigs[part].GetSelfLocations( newContigLocs );
        for ( vec<ContigLocation>::iterator newContigLocIter = newContigLocs.begin();
              newContigLocIter != newContigLocs.end();
              ++newContigLocIter )
          if ( newContigLocIter->GetSuper() == m_super )
            m_locSets[part].push_back( *newContigLocIter );
      }
  }

  // Now that everything is neatly in whole contigs, perform the split.
  this->SuperSplit::Execute( pLog );
}


//------------------------------------------------------------

SuperSanitize::SuperSanitize( Assembly *p_assembly, 
			      const Super &theSuper,
                              const int minLinks,
			      const double maxStretch,	
			      const int minShortLinks,
			      const int lenShortLinks )
  : SuperSplit( p_assembly, theSuper, vec< vec<ContigLocation> >( 0 ) )
{
  vec<ContigLocation> contigLocs;
  m_super.GetContigLocations( contigLocs );

  equiv_rel connected( contigLocs.size() );

  map<ContigLocation,int> locToPos;
  for ( int i = 0; i < contigLocs.isize(); ++i )
    locToPos.insert( make_pair( contigLocs[i], i ) );

  for ( int i = 0; i < contigLocs.isize(); ++i ) {
    Contig theContig = contigLocs[i].GetContig();
    vec<Link> links;
    theContig.GetExternalLinks( links );
    
    vec<int> linkLongCount( contigLocs.size(), 0 );
    vec<int> linkShortCount( contigLocs.size( ), 0 );
    vec<Link>::iterator linkIter;
    for ( linkIter = links.begin();
          linkIter != links.end();
          ++linkIter )
    {
      // omit links between distinct supers
      if ( linkIter->IsBetweenSupers() )
	continue;

      // omit illogical links
      if ( ! linkIter->IsLogical() )
	continue;

      if ( maxStretch > 0 &&
           ! linkIter->IsValid( maxStretch ) )
        continue;

      pair<ContigLocation,ContigLocation> linkContigLocs
        = linkIter->GetContigLocations();

      ContigLocation otherLoc = linkContigLocs.first;
      if ( otherLoc == contigLocs[i] )
        otherLoc = linkContigLocs.second;

      int otherPos = locToPos[ otherLoc ];
      int lenLink = linkIter->GetExpectedSep( );

      if ( lenLink < lenShortLinks ) linkShortCount[otherPos]++;
      else linkLongCount[otherPos]++;
    }

    for ( int j = 0; j < linkShortCount.isize(); ++j ) {
      int totCount = linkLongCount[j] + linkShortCount[j];
      
      // minShortLinks is not defined.
      if ( minShortLinks == 0 ) {
	if ( totCount >= minLinks ) connected.Join( i, j );
	continue;
      }

      // minShortLinks is defined, but there are some long links.
      if ( linkLongCount[j] > 0 ) {
	if ( totCount >= minLinks ) connected.Join( i, j );
	continue;
      }

      // minShortLinks is defined, and there are only short links.
      if ( linkShortCount[j] >= minShortLinks )
	connected.Join( i, j );
    }

  }

  // create new supers only when there is more than one connected
  // component in the super
  if ( connected.OrbitCount() > 1 ) {
    m_locSets.resize( connected.OrbitCount() );

    vec<int> reps;
    connected.OrbitRepsAlt( reps );

    for ( int i = 0; i < contigLocs.isize(); ++i ) {
      int component = BinPosition( reps, connected.ClassId( i ) );
      m_locSets[component].push_back( contigLocs[i] );
    }
  }
}


void
SuperSanitize::Execute( ostream *pLog )
{
  if ( ! m_locSets.empty() )
    this->SuperSplit::Execute( pLog );
}


//------------------------------------------------------------


SuperPop::SuperPop( Assembly *p_assembly,
		    const Super &theSuper,
		    const int minContigSize ) 
  : SuperSplit( p_assembly, theSuper, vec< vec<ContigLocation> >( 3 ) )
{
  vec<ContigLocation> allContigLocs;
  m_super.GetContigLocations( allContigLocs );

  int numContigLocs = allContigLocs.size();
  if ( numContigLocs > 3 ) {

  // Calculate the set of contig locations to be put in the "to be
  // saved" set.
    for ( int i=0; i<numContigLocs; ++i ) {
      if ( i==0 && allContigLocs[i].Length() < m_minContigSize )
	m_locSets[0].push_back( allContigLocs[i] );
      else if ( i==numContigLocs-1 && allContigLocs[i].Length() >= m_minContigSize ) 
	m_locSets[2].push_back( allContigLocs[i] );
      else
	m_locSets[1].push_back( allContigLocs[i] );
    }
  }
}

void SuperPop::Execute(  ostream *pLog )
{
  this->SuperSplit::Execute( pLog );
}


//------------------------------------------------------------

SuperAppend::SuperAppend(  Assembly *p_assembly,
			   const Super &super1,
			   const Super &super2,
			   const int gap, 
			   const Orientation theOrientation ) 
  : AssemblyOp(p_assembly), m_superToAppendTo(super1), m_superToAppend(super2), m_gap(gap), m_orientation(theOrientation)
{}

void SuperAppend::Execute( ostream *pLog )
{

  if ( m_orientation == orient_RC ) {
     mp_assembly->Reverse(m_superToAppend);
  }

  int superLen = m_superToAppendTo.GetLength();
  int startOnNewSuper = superLen + m_gap;

  vec<ContigLocation> appendLocs;
  m_superToAppend.GetContigLocations(appendLocs);
  sort( appendLocs.begin(), appendLocs.end() );

  int numLocs = appendLocs.size();
  for ( int i=0; i<numLocs; ++i ) {
    
    int beginOnNewSuper =  startOnNewSuper + appendLocs[i].Begin();
    Contig oldContig = appendLocs[i].GetContig();
    
    ContigLocation newLoc( m_superToAppendTo, oldContig, beginOnNewSuper, appendLocs[i].GetOrientation() );

    mp_assembly->AddContigLocation( newLoc );

  }

  mp_assembly->Clear(m_superToAppend);
}

    
//------------------------------------------------------------

SuperSort::SuperSort(  Assembly *p_assembly,
			 const Super &super)
  : AssemblyOp(p_assembly), m_superToSort(super)
{}

void SuperSort::Execute( ostream *pLog )
{

  vec<ContigLocation> origLocs;
  m_superToSort.GetContigLocations(origLocs);
  sort( origLocs.begin(), origLocs.end() );

  int numLocs = origLocs.size();
  mp_assembly->Clear(m_superToSort);

  for ( int i=0; i<numLocs; ++i ) {
    
    Contig oldContig = origLocs[i].GetContig();
    
    ContigLocation newLoc( m_superToSort, oldContig, origLocs[i].Begin(), origLocs[i].GetOrientation() );

    mp_assembly->AddContigLocation( newLoc );

  }

}

    
//------------------------------------------------------------

SuperInsert::SuperInsert( Assembly *p_assembly,
	       const Super &superToInsertInto,
	       const vec<Super> &supersToInsert,
	       const int insertionPoint,
	       const vec<int> gaps,
	       const vec<Orientation> theOrientations )
  : AssemblyOp(p_assembly), 
    m_superToInsertInto(superToInsertInto), 
    m_supersToInsert(supersToInsert), 
    m_insertionPoint(insertionPoint),
    m_gaps(gaps), 
    m_orientations(theOrientations) 
{}

void
SuperInsert::Execute( ostream *pLog  )
{

  ForceAssertEq(m_supersToInsert.size(), m_gaps.size());
  ForceAssertEq(m_supersToInsert.size(), m_orientations.size());
  
  // first we slice at the insertion point
  SuperSlice slicey(mp_assembly, m_superToInsertInto, m_insertionPoint, m_insertionPoint);
  slicey.Execute( pLog );

  // then we append 
  vec<Super> newSupers;
  slicey.GetResultingSupers(newSupers);

  int numToInsert( m_supersToInsert.size() );
  for ( int i=0; i<numToInsert; ++i ) {

    if ( pLog )
      *pLog << "Appending " << m_supersToInsert[i].GetId() << " to " << m_superToInsertInto.GetId() << " at " << m_insertionPoint <<" with gap " << m_gaps[i] << " and orientation " << m_orientations[i] << endl;

    SuperAppend appendx(mp_assembly, newSupers[0], m_supersToInsert[i], m_gaps[i], m_orientations[i]);
    appendx.Execute( pLog );

  }
  
  // finally we append the bit previously sliced off (gapping is weak!)
  SuperAppend appendx(mp_assembly, newSupers[0], newSupers[1], 100, orient_FW);
  appendx.Execute( pLog );
  

}

 
 
