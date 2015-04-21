// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#ifndef ASSEMBLY_OFFSET_H
#define ASSEMBLY_OFFSET_H

#include "assembly/Assembly.h"

#include "TaskTimer.h"

#include <set>

// This class encapsulates a potential positioning of one super or
// contig relative to another according to one link.

template <class Entity>
class Offset 
{
 public:
  Offset() {}
  Offset( const Link &theLink );

  Entity GetAnchored() const { return mAnchored; }
  Entity GetFloating() const { return mFloating; }
  
  void SetAnchored( const Entity &anchored ) { mAnchored = anchored; }
  void SetFloating( const Entity &floating ) { mFloating = floating; }

  int         GetOffsetAmount() const { return mOffset; }
  int         GetDev() const          { return mDev; }
  Orientation GetOrientation() const  { return mOrient; }

  void SetOffsetAmount( const int offsetAmount )   { mOffset = offsetAmount; }
  void SetDev         ( const int dev )            { mDev = dev; }
  void SetOrientation( const Orientation &orient ) { mOrient = orient; }

  ReadPair    GetPair() const         { return mPair; }
  int         GetInsertSize() const   { return mInsertSize; }

  bool IsForward() const  { return ( mOrient == orient_FW ); }
  bool IsReverse() const  { return ( mOrient == orient_RC ); }

  Interval GetIntervalOnAnchored() const { return mIntervalOnAnchored; }
  Interval GetIntervalOnFloating() const { return mIntervalOnFloating; }
  
  bool operator< ( const Offset &other ) const
  {
    return ( mOffset < other.mOffset ||
             mOffset == other.mOffset &&
             mPair.GetId() < other.mPair.GetId() );
  }
  
  friend 
  ostream & operator<< ( ostream & out, const Offset &so )
  {
    return out << so.mOffset << "\t" 
               << "+/- " << so.mDev 
               << "\t" << so.mOrient;
  }

  Offset<Entity> 
  GetSwapped() const
  {
    Offset<Entity> swapped( *this );

    swap( swapped.mAnchored, swapped.mFloating );

    if ( mOrient == orient_FW )
      swapped.mOffset = -mOffset;
    else
      swapped.mOffset = mOffset + mAnchored.GetLength() - mFloating.GetLength();

    return swapped;
  }
  
 private:
  Entity mAnchored;
  Entity mFloating;

  Interval mIntervalOnAnchored, mIntervalOnFloating;

  int mOffset;
  int mDev;
  Orientation mOrient;

  int mInsertSize;
  ReadPair mPair;
};

template <typename Entity> 
struct OrderOffsetByFloating {
  bool operator() ( const Offset<Entity>& lhs, const Offset<Entity>& rhs ) const {
    return lhs.GetFloating() < rhs.GetFloating();
  }
};

typedef Offset<Super>  SuperOffset;
typedef Offset<Contig> ContigOffset;


// This class encapsulates a potential positioning of one super or
// contig relative to another according to a weighted set of links.

template <class Entity>
class WeightedOffset 
{
 public:
  WeightedOffset() 
    : mWeight( -1.0 ) {}

  // The Offsets in the range [begin,end) must all be between the same
  // two Entities (which we will here call s1 and s2) and have the
  // same orientation.  These Offsets are used to determine the offset
  // of s2 relative to s1 by finding a weighted average of the
  // remaining Links, where the weight of each Link is inversely
  // proportional to its variance.
  WeightedOffset( typename set< Offset<Entity> >::iterator begin, 
                  typename set< Offset<Entity> >::iterator end );
  
  Entity GetAnchored() const { return mAnchored; }
  Entity GetFloating() const { return mFloating; }
  
  float GetWeight()       const { return mWeight; }
  int   GetNumLinks()     const { return mPairs.size(); }
  int   GetOffsetAmount() const { return mOffset; }
  int   GetStdev()        const { return mStdev; }
  Orientation GetOrient() const { return mOrient; }

  int GetSpread() const { return mSpread; }

  bool  IsForward() const  { return ( mOrient == orient_FW ); }
  bool  IsReverse() const  { return ( mOrient == orient_RC ); }
  
  void  GetPairs( vec<ReadPair> &pairs ) const { pairs = mPairs; }

  bool IsInvalid() const { return mWeight < 0.0; }
  void Invalidate() { mWeight = -1.0; }
  
  WeightedOffset<Entity> 
  GetSwapped() const
  {
    WeightedOffset<Entity> swapped( *this );

    swap( swapped.mAnchored, swapped.mFloating );

    if ( mOrient == orient_FW )
      swapped.mOffset = -mOffset;
    else
      swapped.mOffset = mOffset + mAnchored.GetLength() - mFloating.GetLength();

    return swapped;
  }
  
  bool operator< ( const WeightedOffset &other ) const
  {
    return mWeight < other.mWeight; 
  }
  
  friend
  ostream & operator<< ( ostream &out, const WeightedOffset<Entity> &obj )
  {
    return out << obj.mAnchored 
               << " (l=" << obj.mAnchored.GetLength() << ")"
               << " < " << obj.mOrient << " > " 
               << obj.mFloating
               << " (l=" << obj.mFloating.GetLength() << ")"
               << ": " << obj.mOffset
               << " +/- " << obj.mStdev 
               << " [" << obj.mPairs.size() << "]"
               << " (" << obj.mWeight << ")";
  }
 
  struct Filter : public unary_function<WeightedOffset<Entity>,bool>
  {
    Filter( const int minSize, const int maxSize,
            const int minLinks,
    	    const int medStdDev, const int minMedStdDevLinks,
	    const int highStdDev, const int minHighStdDevLinks,
	    const int maxStdDev,
	    const int minSpread=0 )
      : mMinSize( minSize ), mMaxSize( maxSize ),
        mMinLinks( minLinks ),
        mMedStdDev( medStdDev ), mMinMedStdDevLinks( minMedStdDevLinks ),
        mHighStdDev( highStdDev ), mMinHighStdDevLinks( minHighStdDevLinks ),
        mMaxStdDev( maxStdDev ), mMinSpread(minSpread)
    {
    }

    int GetMinLinks() const { return mMinLinks; }
    
    // Should we use this offset?
    bool operator() ( const WeightedOffset &theOffset ) const
    {
      // Is its deviation too high?
      if ( theOffset.GetStdev() > mMaxStdDev )
        return false;

      // Are there enough links in general?
      if ( theOffset.GetNumLinks() < mMinLinks )
        return false;

      // For high (but acceptable) deviation links, are there enough
      // links?  (Could be more restrictive than minLinks.)
      if ( theOffset.GetStdev() >= mHighStdDev &&
           theOffset.GetNumLinks() < mMinHighStdDevLinks )
        return false;

      // For medium deviation links, are there enough links?  (Could
      // be more restrictive than minLinks.)
      if ( theOffset.GetStdev() < mHighStdDev &&
           theOffset.GetStdev() >= mMedStdDev &&
           theOffset.GetNumLinks() < mMinMedStdDevLinks )
        return false;

      // If we've put a cap on the size, don't pass offsets between
      // entities larger than that limit.
      if ( mMaxSize > 0 && 
           theOffset.GetAnchored().GetLength() > mMaxSize &&
           theOffset.GetFloating().GetLength() > mMaxSize )
        return false;
      
      // If we've got a minimum size, don't pass offsets involving
      // entities smaller than that limit.
      if ( mMinSize > 0 && 
           ( theOffset.GetAnchored().GetLength() < mMinSize ||
             theOffset.GetFloating().GetLength() < mMinSize ) )
        return false;
      
      // If we've got a minimum spread, don't pass offsets which 
      // don't spread out enough
      if ( mMinSpread > 0 &&
	   theOffset.GetSpread() < mMinSpread )
	return false;
 
      return true;
    }

   private:
    int mMinSize, mMaxSize;
    int mMinLinks;
    int mMedStdDev, mMinMedStdDevLinks;
    int mHighStdDev, mMinHighStdDevLinks;
    int mMaxStdDev;
    int mMinSpread;
  };
  
 private:

  // Sets mWeight, mPairs, mOffset, and mStdDev by performing a weighted
  // average of the Offsets.
  void GetWeightedMeanLoc( typename set< Offset<Entity> >::iterator begin,
                           typename set< Offset<Entity> >::iterator end );
  

// Sets mSpread, which is the spread on the Entity of the reads making up 
// the offsets 
  void Spread( typename set< Offset<Entity> >::iterator begin,
	       typename set< Offset<Entity> >::iterator end );


  Entity mAnchored, mFloating;
  float mWeight;
  int mOffset;
  int mStdev;
  Orientation mOrient;
  vec<ReadPair> mPairs;  

  int mSpread;
};

typedef WeightedOffset<Super>  WeightedSuperOffset;
typedef WeightedOffset<Contig> WeightedContigOffset;




// This class builds WeightedOffsets from an assembly.

template <class Entity>
class WeightedOffsetBuilder
{
 public:
  WeightedOffsetBuilder( const float numDevs,
                         const int minGoodInternalLinks = 0,
                         ostream *pLog = 0 )
    : mpFilter( 0 ), 
      mNumDevs( numDevs ), 
      mUseInternalLinks( false ),
      mUseRepetitiveReads( true ),
      mMinGoodInternalLinks( minGoodInternalLinks ),
      mpLog( pLog ) {}

  ~WeightedOffsetBuilder()
  {
    delete mpFilter;
  }

  // This function builds offsets for the given contigsToUse in
  // theAssembly.  If expand is false, only links between contigs in
  // contigsToUse are used.  If expand is true, links between any
  // contig and a contig in contigsToUse are used.

  void Build( const set<Contig> &contigsToUse,
              vec< WeightedOffset<Entity> > &offsets,
              bool expand = false );

  void SetUseInternalLinks( const bool use ) { mUseInternalLinks = use; }

  void SetUseRepetitiveReads( const bool use ) { mUseRepetitiveReads = use; }

  void SetFilter( const typename WeightedOffset<Entity>::Filter &theFilter )
  {
    delete mpFilter;
    mpFilter = new typename WeightedOffset<Entity>::Filter( theFilter );
  }

 private:

  void AddWeightedOffsetFromOffsets( vec< WeightedOffset<Entity> > &weightedOffsets,
                                     set< Offset<Entity> > &theSet );

  void BuildFromLinks( vec<Link>::iterator begin, vec<Link>::iterator end,
                       vec< WeightedOffset<Entity> > &weightedOffsets );

  bool IsGoodContig( const Contig &theContig );

  vec<Link> mInterentityLinks, mContigLinks;
  typename WeightedOffset<Entity>::Filter *mpFilter;

  float    mNumDevs;
  bool     mUseInternalLinks;
  bool     mUseRepetitiveReads;
  int      mMinGoodInternalLinks;
  ostream *mpLog;

  set<Contig> mGoodContigs;
  set<Contig> mBadContigs;

  TaskTimer mGetLinksTimer;
  TaskTimer mSortLinksTimer1;
  TaskTimer mEraseLinksTimer;
  TaskTimer mSortLinksTimer2;
  TaskTimer mBuildTimer;
  TaskTimer mPushTimer;
};

typedef WeightedOffsetBuilder<Super>  WeightedSuperOffsetBuilder;
typedef WeightedOffsetBuilder<Contig> WeightedContigOffsetBuilder;

#endif
