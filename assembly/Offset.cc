// Copyright (c) 2004 The Broad Institute at MIT and Harvard


#include "assembly/Offset.h"

#include <functional>
#include <math.h>

// SuperOffset

template<>
Offset<Super>::Offset( const Link &theLink )
  : mDev( theLink.GetPair().GetExpectedStdDev() ),
    mPair( theLink.GetPair() )
{
  pair<ReadLocationInSuper,ReadLocationInSuper> pairReadLocs 
    = theLink.GetReadLocationsInSupers();
  
  mAnchored = pairReadLocs.first.GetSuper();
  mFloating = pairReadLocs.second.GetSuper();
  
  // We'll anchor the larger supercontig (since that will generally
  // reduce the expense should the floater need to be flipped).  If
  // they are of equal length, order them by id (so there is always a
  // consistent ordering for two given supers).
  if ( mAnchored.GetLength() < mFloating.GetLength() ||
       mAnchored.GetLength() == mFloating.GetLength() &&
       mAnchored.GetId() < mFloating.GetId() )
  {
    swap( pairReadLocs.first, pairReadLocs.second );
    swap( mAnchored, mFloating );
  }
  
  ReadLocationInSuper &anchoredReadLoc = pairReadLocs.first;

  mIntervalOnAnchored = Interval(anchoredReadLoc.Begin(), anchoredReadLoc.End());

  Read read2 = pairReadLocs.second.GetRead();
  
  int impliedReadBegin, impliedReadEnd;
  Orientation impliedReadOrient;
  if ( pairReadLocs.first.IsForward() )
  {
    impliedReadBegin  = anchoredReadLoc.End() + mPair.GetExpectedSep();
    impliedReadEnd    = impliedReadBegin + read2.GetLength();
    impliedReadOrient = orient_RC;
  }
  else
  {
    impliedReadEnd    = anchoredReadLoc.Begin() - mPair.GetExpectedSep();
    impliedReadBegin  = impliedReadEnd - read2.GetLength();
    impliedReadOrient = orient_FW;
  }
  
  ReadLocationInSuper &actualReadLoc = pairReadLocs.second;

  mIntervalOnFloating = Interval(actualReadLoc.Begin(), actualReadLoc.End());

  if ( impliedReadOrient == actualReadLoc.GetOrientation() )
  {
    mOffset = impliedReadBegin - actualReadLoc.Begin();
    mOrient = orient_FW;
  }
  else
  {
    mOffset = impliedReadBegin - ( mFloating.GetLength() - actualReadLoc.End() );
    mOrient = orient_RC;
  }

  mInsertSize = theLink.GetExpectedInsertSize();
}


// ContigOffset

template<>
Offset<Contig>::Offset( const Link &theLink )
  : mDev( theLink.GetPair().GetExpectedStdDev() ),
    mPair( theLink.GetPair() )
{
  pair<ReadLocation,ReadLocation> pairReadLocs 
    = theLink.GetReadLocations();
  
  mAnchored = pairReadLocs.first.GetContig();
  mFloating = pairReadLocs.second.GetContig();
  
  // We'll anchor the larger contig (since that will generally reduce
  // the expense should the floater need to be flipped).  If they are
  // of equal length, order them by id (so there is always a
  // consistent ordering for two given Contigs).
  if ( mAnchored.GetLength() < mFloating.GetLength() ||
       mAnchored.GetLength() == mFloating.GetLength() &&
       mAnchored.GetId() < mFloating.GetId() )
  {
    swap( pairReadLocs.first, pairReadLocs.second );
    swap( mAnchored, mFloating );
  }
  
  ReadLocation &anchoredReadLoc = pairReadLocs.first;

  mIntervalOnAnchored = Interval(anchoredReadLoc.Begin(), anchoredReadLoc.End());

  Read read2 = pairReadLocs.second.GetRead();
  
  int impliedReadBegin, impliedReadEnd;
  Orientation impliedReadOrient;
  if ( pairReadLocs.first.IsForward() )
  {
    impliedReadBegin  = anchoredReadLoc.End() + mPair.GetExpectedSep();
    impliedReadEnd    = impliedReadBegin + read2.GetLength();
    impliedReadOrient = orient_RC;
  }
  else
  {
    impliedReadEnd    = anchoredReadLoc.Begin() - mPair.GetExpectedSep();
    impliedReadBegin  = impliedReadEnd - read2.GetLength();
    impliedReadOrient = orient_FW;
  }
  
  ReadLocation &actualReadLoc = pairReadLocs.second;

  mIntervalOnFloating = Interval(actualReadLoc.Begin(), actualReadLoc.End());
  
  if ( impliedReadOrient == actualReadLoc.GetOrientation() )
  {
    mOffset = impliedReadBegin - actualReadLoc.Begin();
    mOrient = orient_FW;
  }
  else
  {
    mOffset = impliedReadBegin - ( mFloating.GetLength() - actualReadLoc.End() );
    mOrient = orient_RC;
  }

  mInsertSize = theLink.GetExpectedInsertSize();
}


template <class Entity>
WeightedOffset<Entity>::WeightedOffset( typename set< Offset<Entity> >::iterator begin,
                                        typename set< Offset<Entity> >::iterator end )
  : mWeight( 0 ), mStdev(-1)
{
  int numLinks = distance( begin, end );

  if ( numLinks == 0 )
    return;

  mPairs.resize( numLinks );
  transform( begin, end,
             mPairs.begin(),
             mem_fun_ref( &Offset<Entity>::GetPair ) );

  mAnchored = begin->GetAnchored();
  mFloating = begin->GetFloating();

  Spread( begin, end );

  GetWeightedMeanLoc( begin, end );
}


template <class Entity>
void WeightedOffset<Entity>::Spread( typename set< Offset<Entity> >::iterator begin,
				     typename set< Offset<Entity> >::iterator end )
{

  Interval onAnchored( begin->GetIntervalOnAnchored() );
  Interval onFloating( begin->GetIntervalOnFloating() );

  typename set< Offset<Entity> >::iterator offsetIter = begin;
  for ( ; offsetIter != end; ++offsetIter )
  {
    onAnchored.SetBegin( std::min( (offsetIter->GetIntervalOnAnchored()).Begin(), onAnchored.Begin() ) );
    onFloating.SetBegin( std::min( (offsetIter->GetIntervalOnFloating()).Begin(), onFloating.Begin() ) );

    onAnchored.SetEnd( std::max( (offsetIter->GetIntervalOnAnchored()).End(), onAnchored.End() ) );
    onFloating.SetEnd( std::max( (offsetIter->GetIntervalOnFloating()).End(), onFloating.End() ) );
  }
  
  mSpread = std::min( onAnchored.Length(), onFloating.Length() );

}

template <class Entity>
void
WeightedOffset<Entity>::GetWeightedMeanLoc( typename set< Offset<Entity> >::iterator begin,
                                            typename set< Offset<Entity> >::iterator end )
{
  mOrient = begin->GetOrientation();

  if ( distance( begin, end ) == 1 )
  {
    mOffset = begin->GetOffsetAmount();
    mStdev = begin->GetDev();
  }
  else
  {
    // The weighted mean offset, G, will be:
    //
    // sum( w[i] * offset[i] ) / sum( w[i] )
    //
    // where w[i] = 1 / ( stdev[i] * stdev[i] )
    
    float      sumOfWeights = 0.0;
    float      sumOfSquaredWeights = 0.0;
    float      weightedMeanOffset = 0.0;
    
    for ( typename set< Offset<Entity> >::iterator locIter = begin;
          locIter != end; ++locIter )
    {
      float dev = locIter->GetDev();
      float sigma = dev * dev;
      float weight = 1.0 / sigma;

      sumOfWeights += weight;
      sumOfSquaredWeights += weight * weight;
      weightedMeanOffset += weight * locIter->GetOffsetAmount();
    }
    weightedMeanOffset /= sumOfWeights;
    
    // The estimated variance of the offsets, V, will be:
    //
    // sum( weight[i] * ( offset[i] - G )^2 ) * sum( weight[i] )
    // ---------------------------------------------------------
    //        sum( weight[i]^2 ) - sum( weight[i] )^2
    //
    // and the weighted mean of the variances of the samples, MV, will be:
    //
    // 1 / ( sum( w[i] ) / N ) / sqrt(N)
    //
    // We will use the maximum of V and MV as the variance.
    
    float variance = 0.0;
    for ( typename set< Offset<Entity> >::iterator locIter = begin;
          locIter != end; ++locIter )
    {
      float difference = static_cast<float>( locIter->GetOffsetAmount() - weightedMeanOffset );
      
      float dev = locIter->GetDev();
      float sigma = dev * dev;
      float weight = 1.0 / sigma;
      
      variance += weight * difference * difference;
    } 
    
    variance *= sumOfWeights;
    variance /= ( sumOfWeights * sumOfWeights - sumOfSquaredWeights );
    
    float weightedVariance = 1.0 / ( sumOfWeights / static_cast<float>( mPairs.size() ) );
    weightedVariance /= sqrtf( mPairs.size() );

    mOffset = static_cast<int>( rintf( weightedMeanOffset ) );

    float chosenDev = sqrtf( max( variance, weightedVariance ) );
    mStdev = static_cast<int>( rintf( chosenDev ) );
  }

  mWeight = static_cast<float>( mPairs.size() );
  mWeight *= mWeight;
  mWeight /= static_cast<float>( mStdev );

  // Reward diversity of insert sizes.
  int minInsertSize = begin->GetInsertSize();
  int maxInsertSize = begin->GetInsertSize();
  for ( ; begin != end; ++begin )
  {
    minInsertSize = min( minInsertSize, begin->GetInsertSize() );
    maxInsertSize = max( maxInsertSize, begin->GetInsertSize() );
  }
  float reward = log10f( maxInsertSize ) - log10f( minInsertSize ) + 1.0;

  mWeight *= reward;
}


template <class Entity>
void
WeightedOffsetBuilder<Entity>::AddWeightedOffsetFromOffsets( vec< WeightedOffset<Entity> > &weightedOffsets,
                                                             set< Offset<Entity> > &currentOffsets )
{
  if ( mpFilter == 0 || (int) currentOffsets.size() >= mpFilter->GetMinLinks() )
  {
    WeightedOffset<Entity> weightedOffset( currentOffsets.begin(),
                                           currentOffsets.end() );
    if ( mpFilter == 0 || (*mpFilter)( weightedOffset ) )
    {
      weightedOffsets.push_back( weightedOffset );
      if ( mpLog )
      {
        copy( currentOffsets.begin(), currentOffsets.end(),
              ostream_iterator< Offset<Entity> >( *mpLog, "\n" ) );
        *mpLog << weightedOffset << "\n";
      }
    }
  }
}

template <class Entity>
void
WeightedOffsetBuilder<Entity>::BuildFromLinks( vec<Link>::iterator begin, vec<Link>::iterator end,
                                               vec< WeightedOffset<Entity> > &weightedOffsets )
{
  // First we make Offsets out of all the links and sort them.
  vec< Offset<Entity> > offsets;
  for ( ; begin != end; ++begin )
    offsets.push_back( Offset<Entity>( *begin ) );

  sort( offsets.begin(), offsets.end() );
  
  // We now find the maximal sets of "agreeing" offsets, where all
  // pairs of offsets in a set have some minimum probability of being
  // from the same distribution.  In practice, you can think of it as
  // transforming the offsets into intervals [offset-N*dev,
  // offset+N*dev) and finding the highest piles of these intervals.

  // To do this, we'll go through the offsets one at a time, keeping a
  // set of the currently overlapping offsets.  Whenever we encounter
  // an offset that doesn't overlap with something in the current set,
  // we make a weighted offset out of the current set then remove the
  // ones that don't overlap the new one, then add the new one.  The
  // key invariant is that at the end of each pass, all the offsets in
  // the set overlap.
  
  set< Offset<Entity> > currentOffsets;

  typename vec< Offset<Entity> >::iterator newOffsetIter = offsets.begin();
  for ( ; newOffsetIter != offsets.end(); ++newOffsetIter )
  {
    // Calculate the beginning of the interval we want to add.
    float newIntervalBegin = 
      newOffsetIter->GetOffsetAmount() - mNumDevs * newOffsetIter->GetDev();

    // If any offset in the current set doesn't overlap with the new
    // one, make a WeightedOffset out of the current set.  If it
    // passes the filter, add it to the output vector.
    typename set< Offset<Entity> >::iterator currOffsetIter = currentOffsets.begin();
    for ( ; currOffsetIter != currentOffsets.end(); ++currOffsetIter )
    {
      // Calculate the end of the next interval in the set.
      float currIntervalEnd = 
        currOffsetIter->GetOffsetAmount() + mNumDevs * currOffsetIter->GetDev();

      if ( currIntervalEnd <= newIntervalBegin )
      {
        AddWeightedOffsetFromOffsets( weightedOffsets, currentOffsets );
        break;
      }
    }

    // Now remove any offsets in the current set that don't overlap
    // with the new one.
    currOffsetIter = currentOffsets.begin();
    while ( currOffsetIter != currentOffsets.end() )
      if ( currOffsetIter->GetOffsetAmount() + mNumDevs * currOffsetIter->GetDev() <= 
           newIntervalBegin )
      {
        typename set< Offset<Entity> >::iterator deadIter = currOffsetIter;
        ++currOffsetIter;
        currentOffsets.erase( deadIter );
      }
      else
        ++currOffsetIter;

    // Insert the new one.
    currentOffsets.insert( *newOffsetIter );
  }

  // Get the last pile.
  AddWeightedOffsetFromOffsets( weightedOffsets, currentOffsets );
}

template <class Entity>
bool 
WeightedOffsetBuilder<Entity>::IsGoodContig( const Contig &theContig )
{
  if ( mMinGoodInternalLinks == 0 )
    return true;

  if ( mBadContigs.count( theContig ) )
    return false;

  if ( mGoodContigs.count( theContig ) )
    return true;

  vec<Link> contigLinks;
  theContig.GetLinks( contigLinks );

  int goodInternalLinkCount = 0;

  for ( vec<Link>::iterator linkIter = contigLinks.begin();
        linkIter != contigLinks.end() && 
          goodInternalLinkCount < mMinGoodInternalLinks;
        ++linkIter )
    if ( ! linkIter->IsBetweenContigs() &&
         linkIter->IsValid( 1.0 ) )
      ++goodInternalLinkCount;

  if ( goodInternalLinkCount < mMinGoodInternalLinks )
  {
    mBadContigs.insert( theContig );
    return false;
  }

  mGoodContigs.insert( theContig );
  return true;
}



// Various functors used to sort links in WeightedOffsetBuilder<Super>::Build().

struct order_links_by_supers_and_reads : public binary_function<Link,Link,bool>
{
    bool operator() ( const Link &lhs, const Link &rhs ) const
    {
        pair<ReadLocationInSuper,ReadLocationInSuper> lhsRLs, rhsRLs;
        lhsRLs = lhs.GetReadLocationsInSupers();
        rhsRLs = rhs.GetReadLocationsInSupers();

        if ( lhsRLs.first.GetSuper() > lhsRLs.second.GetSuper() ) 
          swap( lhsRLs.first, lhsRLs.second );
        if ( rhsRLs.first.GetSuper() > rhsRLs.second.GetSuper() ) 
          swap( rhsRLs.first, rhsRLs.second );

        if ( lhsRLs.first.GetSuper() < rhsRLs.first.GetSuper() ) return true;
        if ( lhsRLs.first.GetSuper() > rhsRLs.first.GetSuper() ) return false;
 
        if ( lhsRLs.second.GetSuper() < rhsRLs.second.GetSuper() ) return true;
        if ( lhsRLs.second.GetSuper() > rhsRLs.second.GetSuper() ) return false;

        if ( lhsRLs.first.GetRead() < rhsRLs.first.GetRead() ) return true;
        if ( lhsRLs.first.GetRead() > rhsRLs.first.GetRead() ) return false;
 
        if ( lhsRLs.second.GetRead() < rhsRLs.second.GetRead() ) return true;
        return false;
    }
};

struct order_links_by_supers_and_orient : public binary_function<Link,Link,bool>
{
    bool operator() ( const Link &lhs, const Link &rhs ) const
    {
        pair<ReadLocationInSuper,ReadLocationInSuper> lhsRLs, rhsRLs;
        lhsRLs = lhs.GetReadLocationsInSupers();
        rhsRLs = rhs.GetReadLocationsInSupers();

        if ( lhsRLs.first.GetSuper() > lhsRLs.second.GetSuper() ) 
          swap( lhsRLs.first, lhsRLs.second );
        if ( rhsRLs.first.GetSuper() > rhsRLs.second.GetSuper() ) 
          swap( rhsRLs.first, rhsRLs.second );

        if ( lhsRLs.first.GetSuper() < rhsRLs.first.GetSuper() ) return true;
        if ( lhsRLs.first.GetSuper() > rhsRLs.first.GetSuper() ) return false;
 
        if ( lhsRLs.second.GetSuper() < rhsRLs.second.GetSuper() ) return true;
        if ( lhsRLs.second.GetSuper() > rhsRLs.second.GetSuper() ) return false;

        if ( lhsRLs.first.GetOrientation() + lhsRLs.second.GetOrientation() <
             rhsRLs.first.GetOrientation() + rhsRLs.second.GetOrientation() ) return true;
        return false;
    }
};

// WeightedOffsetBuilder<Super> instantiation.

template
void 
WeightedOffsetBuilder<Super>::AddWeightedOffsetFromOffsets( vec< WeightedOffset<Super> > &weightedOffsets,
                                                            set< Offset<Super> > &theSet );
template
void 
WeightedOffsetBuilder<Super>::BuildFromLinks( vec<Link>::iterator begin, vec<Link>::iterator end,
                                              vec< WeightedOffset<Super> > &weightedOffsets );
template
bool
WeightedOffsetBuilder<Super>::IsGoodContig( const Contig &theContig );

template<>
void WeightedOffsetBuilder<Super>::Build( const set<Contig> &contigsToUse,
                                          vec< WeightedOffset<Super> > &offsets,
                                          bool expand )
{
  if ( mpLog )
    *mpLog << "Building offsets:" << "\n";
  //   {
  //     *pLog << "Looking for links " << ( expand ? "to" : "within" )
  //           << " the following " << contigsToUse.size() << " contigs:" << "\n";
  //     copy( contigsToUse.begin(), contigsToUse.end(),
  //           ostream_iterator<Contig>( *pLog, "\n" ) );
  //   }

  mInterentityLinks.clear();

  mGetLinksTimer.Start();
  for ( set<Contig>::iterator contigIter = contigsToUse.begin();
        contigIter != contigsToUse.end();
        ++contigIter )
  {
    contigIter->GetExternalLinks( mContigLinks );

    for ( vec<Link>::iterator linkIter = mContigLinks.begin();
          linkIter != mContigLinks.end();
          ++linkIter )
      if ( mUseInternalLinks || 
           linkIter->IsBetweenSupers() )
      {
        pair<Contig,Contig> pairContigs = linkIter->GetContigs();

        if ( ! this->IsGoodContig( pairContigs.first ) ||
             ! this->IsGoodContig( pairContigs.second ) )
          continue;

        if ( ! mUseRepetitiveReads &&
             ( Read( linkIter->GetRead1() ).IsRepetitive() ||
               Read( linkIter->GetRead2() ).IsRepetitive() ) )
          continue;

        if ( ! expand && 
             ( contigsToUse.count( pairContigs.first ) &&
               contigsToUse.count( pairContigs.second ) ) ||
             expand &&
             ( contigsToUse.count( pairContigs.first ) ||
               contigsToUse.count( pairContigs.second ) ) )
          mInterentityLinks.push_back( *linkIter );
      }
  }
  mGetLinksTimer.Stop();
    
  mSortLinksTimer1.Start();
  sort( mInterentityLinks.begin(), mInterentityLinks.end(),
        order_links_by_supers_and_reads() );
  mSortLinksTimer1.Stop();

  mEraseLinksTimer.Start();
  mInterentityLinks.erase( unique( mInterentityLinks.begin(),
                                  mInterentityLinks.end() ),
                          mInterentityLinks.end() );
  mEraseLinksTimer.Stop();

  mSortLinksTimer2.Start();
  sort( mInterentityLinks.begin(), mInterentityLinks.end(),
        order_links_by_supers_and_orient() );
  mSortLinksTimer2.Stop();
  
  vec<WeightedSuperOffset> newOffsets;

  mBuildTimer.Start();
  pair<vec<Link>::iterator,vec<Link>::iterator> range;
  range.first = range.second = mInterentityLinks.begin();
  for ( ; range.first != mInterentityLinks.end(); range.first = range.second )
  {
    range = equal_range( range.first, mInterentityLinks.end(),
                         *(range.first),
                         order_links_by_supers_and_orient() );

    this->BuildFromLinks( range.first, range.second, newOffsets );
  }
  mBuildTimer.Stop();

  if ( offsets.empty() )
    offsets.reserve( newOffsets.size() );

  mPushTimer.Start();
  for ( vec<WeightedSuperOffset>::iterator newOffsetIter = newOffsets.begin();
        newOffsetIter != newOffsets.end(); ++newOffsetIter )
  {
    offsets.push_back( *newOffsetIter );

    push_heap( offsets.begin(), offsets.end() );
  }
  mPushTimer.Stop();

  if ( mpLog && mGetLinksTimer.GetTasksTimed() == 100 )
  {
    PRINT_TO( *mpLog, mGetLinksTimer );
    PRINT_TO( *mpLog, mSortLinksTimer1 );
    PRINT_TO( *mpLog, mEraseLinksTimer );
    PRINT_TO( *mpLog, mSortLinksTimer2 );
    PRINT_TO( *mpLog, mBuildTimer );
    PRINT_TO( *mpLog, mPushTimer );

    mGetLinksTimer.Reset();
    mSortLinksTimer1.Reset();
    mEraseLinksTimer.Reset();
    mSortLinksTimer2.Reset();
    mBuildTimer.Reset();
    mPushTimer.Reset();
  }
}



// Various functors used to sort links in WeightedOffsetBuilder<Contig>::Build().

struct order_links_by_contigs_and_reads : public binary_function<Link,Link,bool>
{
    bool operator() ( const Link &lhs, const Link &rhs ) const
    {
        pair<ReadLocation,ReadLocation> lhsRLs, rhsRLs;
        lhsRLs = lhs.GetReadLocations();
        rhsRLs = rhs.GetReadLocations();

        // The ReadLocations in a Link are always stored such that the
        // first is less than the second.  ReadLocations are ordered
        // by Contig id.

        if ( lhsRLs.first.GetContig() < rhsRLs.first.GetContig() ) return true;
        if ( lhsRLs.first.GetContig() > rhsRLs.first.GetContig() ) return false;

        if ( lhsRLs.second.GetContig() < rhsRLs.second.GetContig() ) return true;
        if ( lhsRLs.second.GetContig() > rhsRLs.second.GetContig() ) return false;

        if ( lhsRLs.first.GetRead() < rhsRLs.first.GetRead() ) return true;
        if ( lhsRLs.first.GetRead() > rhsRLs.first.GetRead() ) return false;

        if ( lhsRLs.second.GetRead() < rhsRLs.second.GetRead() ) return true;
        return false;
    }
};

struct order_links_by_contigs_and_orient : public binary_function<Link,Link,bool>
{
    bool operator() ( const Link &lhs, const Link &rhs ) const
    {
        pair<ReadLocation,ReadLocation> lhsRLs, rhsRLs;
        lhsRLs = lhs.GetReadLocations();
        rhsRLs = rhs.GetReadLocations();

        // The ReadLocations in a Link are always stored such that the
        // first is less than the second.  ReadLocations are ordered
        // by Contig id.

        if ( lhsRLs.first.GetContig() < rhsRLs.first.GetContig() ) return true;
        if ( lhsRLs.first.GetContig() > rhsRLs.first.GetContig() ) return false;

        if ( lhsRLs.second.GetContig() < rhsRLs.second.GetContig() ) return true;
        if ( lhsRLs.second.GetContig() > rhsRLs.second.GetContig() ) return false;

        if ( lhsRLs.first.GetOrientation() + lhsRLs.second.GetOrientation() <
             rhsRLs.first.GetOrientation() + rhsRLs.second.GetOrientation() ) return true;
        return false;
    }
};


// WeightedOffsetBuilder<Contig> instantiation.

template
void 
WeightedOffsetBuilder<Contig>::AddWeightedOffsetFromOffsets( vec< WeightedOffset<Contig> > &weightedOffsets,
                                                             set< Offset<Contig> > &theSet );
template
void 
WeightedOffsetBuilder<Contig>::BuildFromLinks( vec<Link>::iterator begin, vec<Link>::iterator end,
                                              vec< WeightedOffset<Contig> > &weightedOffsets );
template
bool
WeightedOffsetBuilder<Contig>::IsGoodContig( const Contig &theContig );

template<>
void WeightedOffsetBuilder<Contig>::Build( const set<Contig> &contigsToUse,
                                           vec< WeightedOffset<Contig> > &offsets,
                                           bool expand )
{
  if ( mpLog )
    *mpLog << "Building offsets:" << "\n";
  //   {
  //     *pLog << "Looking for links " << ( expand ? "to" : "within" )
  //           << " the following " << contigsToUse.size() << " contigs:" << "\n";
  //     copy( contigsToUse.begin(), contigsToUse.end(),
  //           ostream_iterator<Contig>( *pLog, "\n" ) );
  //   }

  mInterentityLinks.clear();

  mGetLinksTimer.Start();
  for ( set<Contig>::iterator contigIter = contigsToUse.begin();
        contigIter != contigsToUse.end();
        ++contigIter )
  {
    if ( mUseInternalLinks )
      contigIter->GetLinks( mContigLinks );
    else
      contigIter->GetExternalLinks( mContigLinks );

    for ( vec<Link>::iterator linkIter = mContigLinks.begin();
          linkIter != mContigLinks.end();
          ++linkIter )
    {
      pair<Contig,Contig> pairContigs = linkIter->GetContigs();
      
      if ( ! this->IsGoodContig( pairContigs.first ) ||
           ! this->IsGoodContig( pairContigs.second ) )
        continue;
      
      if ( ! mUseRepetitiveReads &&
           ( Read( linkIter->GetRead1() ).IsRepetitive() ||
             Read( linkIter->GetRead2() ).IsRepetitive() ) )
        continue;
      
      if ( ! expand && 
           ( contigsToUse.count( pairContigs.first ) &&
             contigsToUse.count( pairContigs.second ) ) ||
           expand &&
           ( contigsToUse.count( pairContigs.first ) ||
             contigsToUse.count( pairContigs.second ) ) )
        mInterentityLinks.push_back( *linkIter );
    }
  }
  mGetLinksTimer.Stop();
    
  mSortLinksTimer1.Start();
  sort( mInterentityLinks.begin(), mInterentityLinks.end(),
        order_links_by_contigs_and_reads() );
  mSortLinksTimer1.Stop();

  mEraseLinksTimer.Start();
  mInterentityLinks.erase( unique( mInterentityLinks.begin(),
                                   mInterentityLinks.end() ),
                           mInterentityLinks.end() );
  mEraseLinksTimer.Stop();

  mSortLinksTimer2.Start();
  sort( mInterentityLinks.begin(), mInterentityLinks.end(),
        order_links_by_contigs_and_orient() );
  mSortLinksTimer2.Stop();
  
  vec<WeightedContigOffset> newOffsets;

  mBuildTimer.Start();
  pair<vec<Link>::iterator,vec<Link>::iterator> range;
  range.first = range.second = mInterentityLinks.begin();
  for ( ; range.first != mInterentityLinks.end(); range.first = range.second )
  {
    range = equal_range( range.first, mInterentityLinks.end(),
                         *(range.first),
                         order_links_by_contigs_and_orient() );

    this->BuildFromLinks( range.first, range.second, newOffsets );
  }
  mBuildTimer.Stop();

  if ( offsets.empty() )
    offsets.reserve( newOffsets.size() );

  mPushTimer.Start();
  for ( vec<WeightedContigOffset>::iterator newOffsetIter = newOffsets.begin();
        newOffsetIter != newOffsets.end(); ++newOffsetIter )
  {
    offsets.push_back( *newOffsetIter );

    push_heap( offsets.begin(), offsets.end() );
  }
  mPushTimer.Stop();

  if ( mpLog && mGetLinksTimer.GetTasksTimed() == 100 )
  {
    PRINT_TO( *mpLog, mGetLinksTimer );
    PRINT_TO( *mpLog, mSortLinksTimer1 );
    PRINT_TO( *mpLog, mEraseLinksTimer );
    PRINT_TO( *mpLog, mSortLinksTimer2 );
    PRINT_TO( *mpLog, mBuildTimer );
    PRINT_TO( *mpLog, mPushTimer );

    mGetLinksTimer.Reset();
    mSortLinksTimer1.Reset();
    mEraseLinksTimer.Reset();
    mSortLinksTimer2.Reset();
    mBuildTimer.Reset();
    mPushTimer.Reset();
  }
}
