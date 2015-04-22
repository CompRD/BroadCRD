// Copyright (c) 2003 Broad Institute/Massachusetts Institute of Technology

// A suspectAlignmentMgr manages a list of alignments such that 
// reads which appear to be problematic can be identified and interesting statistics 
// can be generated.
//

#ifndef SUSPECT_READ_ALIGNMENT_MGR
#define SUSPECT_READ_ALIGNMENT_MGR

#include "CoreTools.h"
#include "Qualvector.h"
#include "ReadPairing.h"

#include "lookup/LookAlign.h"

#include <map>

class suspectReadAlignmentMgr
{
 public:

  suspectReadAlignmentMgr() 
  {
      mpVecAligns = 0; 
  }

  suspectReadAlignmentMgr( vec<look_align_plus> *vecAligns );

 private:
  
  vec<look_align_plus> *mpVecAligns;
  map<int, pair<int,int> > mMapIdToAlignIndx;


 public:

  // which 'side' of the (fw-oriented) read aligns
  enum alignSide { left, right, center, full };

 
  // count the number of alignments a read has, or which a read has to a side ----------

  int numberAligns( int readId )
  { 
    pair<int,int> indxs = readIndicies( readId );
    return indxs.second - indxs.first;
  }

 
  int numberAligns( int readId, alignSide sideAligning )
  {
    int numaligns(0);
    pair<int,int> indxs = readIndicies( readId );
    for ( int i=indxs.first; i<indxs.second; ++i )
    {
      if ( IsLeft( (*mpVecAligns)[i] ) && sideAligning == left )
	++numaligns;
      else if ( IsRight( (*mpVecAligns)[i] ) && sideAligning == right )
	++numaligns;
    }
    
    return numaligns;
    
  }


  // if a read has any alignments, it counts.
  int numberReadsWithAligns()
  {
    return mMapIdToAlignIndx.size();
  }


  int readLength( int readId )
  {
    return (*mpVecAligns)[ this->readIndicies( readId ).first ].query_length;
  }

  // get the ids of each read appearing in aligns list
  vec<int> readIds()
  {
    vec<int> vecReadIds;
    int numAligns = mpVecAligns->size();
    for ( int i=0; i<numAligns; ++i )
    {
      vecReadIds.push_back( (*mpVecAligns)[i].query_id );
    }
    UniqueSort( vecReadIds );

    return vecReadIds;
  }

  
  // make sure this read is in the list of alignments
  bool isValidId( int readId )
  { 
    map<int, pair<int,int> >::iterator mapIter = mMapIdToAlignIndx.find( readId );
    if ( mapIter == mMapIdToAlignIndx.end() )
    {
      cout << "Read " << readId << " is not in the list of alignments." << endl;
      return false;
    }

    return true;
  }


  // return begin, end indicies into alignment vector for a given read id
  pair<int, int> readIndicies( int readId )
  {
    if ( isValidId( readId ) )
    {
      map<int, pair<int,int> >::iterator mapIter = mMapIdToAlignIndx.find( readId );
      return make_pair( mapIter->second.first, mapIter->second.second );
    }
    else
      return make_pair( -1, -1 );
  }


  // partner info ------

  bool partnerAlignedWell( int readId, vec<read_pairing> &pairs, vec<int> &pairs_index,
			   int max_link_stretch );

  bool partnerAlignedWell( look_align_plus &align, 
			   vec<read_pairing> &pairs, 
			   vec<int> &pairs_index,
			   int max_link_stretch );

  pair<int,int>  partnerIndicies( int readId, 
				  vec<read_pairing> &pairs, 
				  vec<int> &pairs_index );


  // a read has at least one alignment satisfying left,right, or center criteria -----

  bool hasLeftAlign( int readId );
  bool hasRightAlign( int readId );
  bool hasCenterAlign( int readId );
  bool hasFullLengthAlign( int readId );

  bool hasQualityFullLengthAlign( int readId, float maxerrRate );

  // functions for calculating stats based on quality scores or alignment errors -----

  int maxAlignedBases( int readId );

  int numAlignedBases( int readId, alignSide sideAligning );
  int sumAlignedQuality( int readId, alignSide sideAligning, vecqualvector &quals );
 
  int numUnalignedBases( int readId, alignSide sideAligning );
  int sumUnalignedQuality( int readId, alignSide sideAligning, vecqualvector &quals );

  int numErrors( int readId, alignSide sideAligning );

  float aveQualityAligned( vec<int> &readIds, 
			   alignSide sideAligning, 
			   vecqualvector &quals );
  float aveLengthAligned( vec<int> &readIds, alignSide sideAligning );
  float aveErrorRate(  vec<int> &readIds, alignSide sideAligning );

  float aveQualityUnaligned(  vec<int> &readIds,
			      alignSide sideAligning, 
			      vecqualvector &quals );
  float aveLengthUnaligned(  vec<int> &readIds,alignSide sideAligning );


  float aveQualityInRegion( int readId, int beginRegion, int endRegion, vecqualvector &quals );


  // functions to refine the list of alignments ------------------------

  // Do both the left side of the read and the right side align to the 
  // same chr, not just that one side aligns multiple times.
  bool placedOnOneChr( int readId );

  //  is the hit to a real human chr (1-24)?
  bool onTarget( int readId );

 
  // chimeric specific -------------------------------------------------

  // Is this read chimeric? 
  bool IsChimericCandidate( int readId );

  // return left, right estimates of chimeric point
  pair<int,int> ligationRegion( int readId );

  // does the ligation region involve an overlap, a gap, or do the lhs and rhs of the
  // chimeric read abut?

  bool LigationOverlap( int readId );
  bool LigationGap( int readId );
  bool LigationAbut( int readId );

  // counts of above

  int numLigationOverlaps();
  int numLigationGaps();
  int numLigationAbuts();

  int numBasesLigationOverlaps();
  int numBasesLigationGaps();
  int numBasesLigationAbuts();


  float aveQualityLigationOverlaps( vecqualvector &qual );
  float aveQualityLigationGaps( vecqualvector &qual );
  float aveQualityLigationAbuts( vecqualvector &qual );

  float aveQualityLigationOverlaps( int readId, vecqualvector &qual );
  float aveQualityLigationGaps(  int readId, vecqualvector &qual );
  float aveQualityLigationAbuts(  int readId, vecqualvector &qual );

 
};

#endif
