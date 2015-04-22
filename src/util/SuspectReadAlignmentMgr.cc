// Copyright (c) 2003 Broad Institute/Massachusetts Institute of Technology

// A suspectAlignmentMgr manages a list of alignments such that 
// reads which appear to be problematic can be identified and interesting statistics 
// can be generated.
//

#include "util/SuspectReadAlignmentMgr.h"

suspectReadAlignmentMgr::suspectReadAlignmentMgr( vec<look_align_plus> *vecAligns ) : 
  mpVecAligns( vecAligns ) 
{

  vec<look_align_plus>::iterator alignIter = mpVecAligns->begin();
  
  while ( alignIter != mpVecAligns->end() ) 
  {
    int queryId = alignIter->query_id;
    pair< vec<look_align_plus>::iterator, vec<look_align_plus>::iterator > range;
      
    range = equal_range( mpVecAligns->begin(), mpVecAligns->end(), *alignIter,
			 order_lookalign_Query() );

    int beginIndx = distance( mpVecAligns->begin(), range.first );
    int endIndx = distance( mpVecAligns->begin(), range.second );

    if ( !mMapIdToAlignIndx.insert( make_pair( queryId,
					      make_pair(beginIndx, endIndx)) ).second )
    {
      cout << "Already have " << queryId << endl;
    }
    
    alignIter = range.second;

  }
}


float suspectReadAlignmentMgr::aveQualityAligned( vec<int> &readIds, 
						  alignSide sideAligning, 
						  vecqualvector &quals )
{
  int numReads = readIds.size();
  int sumQual(0), sumBases(0);
  for ( int i=0; i<numReads; ++i )
  {
    sumBases += numAlignedBases( readIds[i], sideAligning );
    sumQual += sumAlignedQuality( readIds[i], sideAligning, quals );
  }

  if ( sumBases == 0 )
    return -1.0;

  return (float) sumQual/sumBases;

}

float suspectReadAlignmentMgr::aveLengthAligned( vec<int> &readIds, alignSide sideAligning )
{
  int numReads = readIds.size();
  if ( numReads == 0 )
    return -1.0;

  int sumBases(0), numAligns(0);
  for ( int i=0; i<numReads; ++i )
  {
    sumBases += numAlignedBases( readIds[i], sideAligning );
    numAligns += numberAligns( readIds[i], sideAligning );
  }

  return (float) sumBases/numAligns;
}

float suspectReadAlignmentMgr::aveErrorRate( vec<int> &readIds, alignSide sideAligning )
{
  int numReads = readIds.size();
  if ( numReads == 0 )
    return -1.0;

  int sumErrors(0), sumBases(0);
  for ( int i=0; i<numReads; ++i )
  {
    sumBases += numAlignedBases( readIds[i], sideAligning );
    sumErrors += numErrors( readIds[i], sideAligning );
  }

  return (float) sumErrors/sumBases;
}

float suspectReadAlignmentMgr::aveQualityUnaligned( vec<int> &readIds, 
						    alignSide sideAligning, 
						    vecqualvector &quals )
{
  int numReads = readIds.size();
  int sumQual(0), sumBases(0);
  for ( int i=0; i<numReads; ++i )
  {
    sumBases += numUnalignedBases( readIds[i], sideAligning );
    sumQual += sumUnalignedQuality( readIds[i], sideAligning, quals );
  }

  if ( sumBases == 0 )
    return -1.0;

  return (float) sumQual/sumBases;

}

float suspectReadAlignmentMgr::aveLengthUnaligned( vec<int> &readIds, 
						   alignSide sideAligning )
{
  int numReads = readIds.size();
  if ( numReads == 0 )
    return -1.0;


  int sumBases(0), numAligns(0);
  for ( int i=0; i<numReads; ++i )
  {
    sumBases += numUnalignedBases( readIds[i], sideAligning );
    numAligns += numberAligns( readIds[i], sideAligning );
  }

  return (float) sumBases/numAligns;
}

bool suspectReadAlignmentMgr::hasLeftAlign( int readId )
{
  pair<int,int> indxs = readIndicies( readId );
  
  for ( int i=indxs.first; i<indxs.second; ++i )
  {
    if ( IsLeft( (*mpVecAligns)[i] ) )
      return true;
  }

  return false;
}
  
bool suspectReadAlignmentMgr::hasRightAlign( int readId )
{
  pair<int,int> indxs = readIndicies( readId );
  
  for ( int i=indxs.first; i<indxs.second; ++i )
  {
    if ( IsRight( (*mpVecAligns)[i] ) )
      return true;
  }

  return false;
}
    
bool suspectReadAlignmentMgr::hasCenterAlign( int readId )
{
  pair<int,int> indxs = readIndicies( readId );
  
  for ( int i=indxs.first; i<indxs.second; ++i )
  {
    if ( IsCenter( (*mpVecAligns)[i] ) )
      return true;
  }

  return false;
}

int suspectReadAlignmentMgr::numErrors( int readId, alignSide sideAligning )
{
  int numErrors(0);
  
  pair<int,int> indxs = readIndicies( readId );
  for ( int i=indxs.first; i<indxs.second; ++i )
  {
    if ( IsLeft( (*mpVecAligns)[i] ) && sideAligning == left )
      numErrors +=  (*mpVecAligns)[i].Errors();
    else if ( IsRight( (*mpVecAligns)[i] ) && sideAligning == right )
      numErrors +=  (*mpVecAligns)[i].Errors();
    else
      numErrors += (*mpVecAligns)[i].Errors();
  }

  return numErrors;
}
  

int suspectReadAlignmentMgr::numAlignedBases( int readId, alignSide sideAligning )
{
  int numBases(0);

  pair<int,int> indxs = readIndicies( readId );
  for ( int i=indxs.first; i<indxs.second; ++i )
  {
    int beginAlign(0), endAlign(0);

    look_align_plus &thisAlign = (*mpVecAligns)[i];
    if ( !thisAlign.rc1 )
    { 
      beginAlign = thisAlign.a.pos1();
      endAlign = thisAlign.a.Pos1();
    }
    else
    {
      beginAlign = thisAlign.query_length-thisAlign.a.Pos1();
      endAlign = thisAlign.query_length-thisAlign.a.pos1();
    }

    if ( sideAligning == left && IsLeft( thisAlign ) )
    {
      numBases += endAlign-beginAlign;
    }
    else if ( sideAligning == right && IsRight( thisAlign ) )
    {
      numBases += endAlign-beginAlign;
    }
  }
  
  return numBases;
}

int suspectReadAlignmentMgr::sumAlignedQuality( int readId, 
						alignSide sideAligning, 
						vecqualvector &quals )
{
  int sumQuals(0);

  pair<int,int> indxs = readIndicies( readId );
  for ( int i=indxs.first; i<indxs.second; ++i )
  {
    look_align_plus &thisAlign = (*mpVecAligns)[i];

    int beginAlign(0), endAlign(0);
    if ( !thisAlign.rc1 )
    { 
      beginAlign = thisAlign.a.pos1();
      endAlign = thisAlign.a.Pos1();
    }
    else
    {
      beginAlign = thisAlign.query_length-thisAlign.a.Pos1();
      endAlign = thisAlign.query_length-thisAlign.a.pos1();
    }

    if ( sideAligning == left && IsLeft( thisAlign ) )
    {
      for ( int j=beginAlign; j<endAlign; ++j )
      {
	sumQuals += (int) quals[ readId ][j];
      }
    }
    else if ( sideAligning == right && IsRight( thisAlign ) )
    {
      for ( int j=beginAlign; j<endAlign; ++j )
      {
	sumQuals += (int) quals[ readId ][j];
      }
    }
  }
  
  return sumQuals;  
}


int suspectReadAlignmentMgr::numUnalignedBases( int readId, alignSide sideAligning )
{
  int numBases(0);

  pair<int,int> indxs = readIndicies( readId );
  //  cout << readId << " has indicies " << indxs.first <<" "<<indxs.second << endl;
  for ( int i=indxs.first; i<indxs.second; ++i )
  {
    look_align_plus &thisAlign = (*mpVecAligns)[i];
    
/*     cout << "\t" << thisAlign.query_id <<" "<<thisAlign.target_id <<" "<< thisAlign.a.pos1() */
/* 	 << " "<< thisAlign.a.Pos1() << " "<< thisAlign.query_length <<" "; */
	 
/*     if ( thisAlign.rc1 ) */
/*       cout << "rc " << endl; */
/*     else */
/*       cout << "fw " << endl; */

    pair<int,int> hanging = hangingEnd( thisAlign );
    //cout << "\thanging: " << hanging.first <<" "<<hanging.second << endl;

    if ( sideAligning == left && IsLeft( thisAlign ) )
    {
      numBases += (hanging.second-hanging.first);
      //cout << "\tleft: " << numBases <<endl;
    }
    else if ( sideAligning == right && IsRight( thisAlign ) )
    {
      numBases += (hanging.second-hanging.first);
      //cout << "\tright: " << numBases << endl;
    }
  }
  
  return numBases;
}


int suspectReadAlignmentMgr::sumUnalignedQuality( int readId, 
						  alignSide sideAligning, 
						  vecqualvector &quals )
{
  int sumQuals(0);

  pair<int,int> indxs = readIndicies( readId );
  for ( int i=indxs.first; i<indxs.second; ++i )
  {
    look_align_plus &thisAlign = (*mpVecAligns)[i];
    
    pair<int,int> hanging = hangingEnd( thisAlign );

    if ( sideAligning == left && IsLeft( thisAlign ) )
    {
      for ( int j=hanging.first; j<hanging.second; ++j )
      {
	sumQuals += (int) quals[ readId ][j];
      }
    }
    else if ( sideAligning == right && IsRight( thisAlign ) )
    {
      for ( int j=hanging.first; j<hanging.second; ++j )
      {
	sumQuals += (int) quals[ readId ][j];
      }
    }
  }
  
  return sumQuals;
}

// TODO: potentially dangerous truncation of indices by int args
float suspectReadAlignmentMgr::aveQualityInRegion( int readId, 
						   int beginRegion, 
						   int endRegion, 
						   vecqualvector &quals )
{
  if ( endRegion-beginRegion <= 0 )
    return -1.0;
  if ( static_cast<unsigned>(endRegion) > quals[readId].size() )
    endRegion = quals[readId].size();
  if ( beginRegion < 0 )
    beginRegion = 0;

  int sumQual(0);
  for ( int i=beginRegion; i<endRegion; ++i )
    sumQual += (int) quals[readId][i];

  return sumQual/(endRegion-beginRegion);
}



// Do both the left side of the read and the right side align to the 
// same chr, not just that one side aligns multiple times.
bool suspectReadAlignmentMgr::placedOnOneChr( int readId )
{
  pair<int,int> indxs = readIndicies( readId );
 
  int begin_index = indxs.first;
  int end_index = indxs.second;
    
  for ( int i=begin_index; i<end_index-1; ++i )
  {
    look_align_plus &thisAlign = (*mpVecAligns)[i];
    
    for ( int j=begin_index+1; j<end_index; ++j )
    {
      look_align_plus &nextAlign = (*mpVecAligns)[j];
      
      if ( thisAlign.target_id != nextAlign.target_id )
	continue;
      
      if ( IsLeft( thisAlign ) && IsRight( nextAlign ) || 
	   IsRight( thisAlign ) && IsLeft( nextAlign ) )
      {
	return true;
      }
    }
  }
	  
  return false;

}

// Does the read has a full length alignment?
bool suspectReadAlignmentMgr::hasFullLengthAlign( int readId )
{
  pair<int,int> indxs = readIndicies( readId );
 
  int begin_index = indxs.first;
  int end_index = indxs.second;

  for ( int i=begin_index; i<end_index; ++i )
  {
    look_align_plus &thisAlign = (*mpVecAligns)[i];
    
    if ( thisAlign.FullLength() )
      return true;
  }

  return false;
}

//  is the hit to a real human chr (1-24)?
bool suspectReadAlignmentMgr::onTarget( int readId )
{
  pair<int,int> indxs = readIndicies( readId );
 
  int begin_index = indxs.first;
  int end_index = indxs.second;

  bool isCandidate(false);
  for ( int i=begin_index; i<end_index; ++i )
  {
    look_align_plus &thisAlign = (*mpVecAligns)[i];

    if ( thisAlign.target_id < 1 ||  thisAlign.target_id > 24 )
      continue;
      
    isCandidate = true;
    break;
  }
    
  return isCandidate;
}


//  chimeric functions


// Is this read chimeric? Assumes we know that it aligns to more than one 
// chromosome, but not how it aligns.
bool suspectReadAlignmentMgr::IsChimericCandidate( int readId )
{
  if ( !isValidId( readId ) )
    return false;

  if ( !onTarget( readId ) )
    return false;
  
  if ( hasFullLengthAlign( readId ) )
    return false;
    
  if ( placedOnOneChr( readId ) )
    return false;

  return true;

}

// return left, right estimates of chimeric point
pair<int,int> suspectReadAlignmentMgr::ligationRegion( int readId )
{
  pair<int,int> indxs = readIndicies( readId );

  int fromLeft(0), fromRight(1000000);
  for ( int i=indxs.first; i<indxs.second; ++i )
  {
    look_align_plus &aligner = (*mpVecAligns)[i];
//     cout << "\t" << aligner.query_id <<" "<<aligner.target_id <<" "<< aligner.a.pos1()
// 	 << " "<< aligner.a.Pos1() << " "<< aligner.query_length <<" "; 
	 
//     if ( aligner.rc1 )
//       cout << "rc " << endl; 
//     else 
//        cout << "fw " << endl;

    pair<int,int> hanging = hangingEnd( aligner );
//     cout << "\thanging: " << hanging.first <<" "<<hanging.second << endl;

    if ( IsLeft( aligner ) )
      fromLeft = Max( hanging.first, fromLeft );
    else if ( IsRight( aligner ) )
      fromRight = Min( hanging.second, fromRight );
    else
      return make_pair( -1, -1 );
  }
    
  return make_pair( fromLeft, fromRight );
}

bool suspectReadAlignmentMgr::LigationOverlap( int readId )
{
  pair<int,int> lig_reg = ligationRegion( readId );
  return ( lig_reg.first > lig_reg.second );
}
  
bool suspectReadAlignmentMgr::LigationGap( int readId )
{
  pair<int,int> lig_reg = ligationRegion( readId );
  return ( lig_reg.first < lig_reg.second );
}

bool suspectReadAlignmentMgr::LigationAbut( int readId )
{
  pair<int,int> lig_reg = ligationRegion( readId );
  if ( lig_reg.first == -1 && lig_reg.second == -1 )
    return false;

  return ( lig_reg.first == lig_reg.second );
}

int suspectReadAlignmentMgr::numLigationOverlaps()
{
  int numOvlps(0);

  vec<int> allReads = readIds();
  for ( int i=0; i<(int) allReads.size(); ++i )
  {
    if ( IsChimericCandidate( allReads[i] ) )
    {
      if ( LigationOverlap( allReads[i] ) )
	++numOvlps;
    }
  }

  return numOvlps;
}

int suspectReadAlignmentMgr::numLigationGaps()
{
  int numGaps(0);

  vec<int> allReads = readIds();
  for ( int i=0; i<(int) allReads.size(); ++i )
  {
    if ( IsChimericCandidate( allReads[i] ) )
    {
      if ( LigationGap( allReads[i] ) )
	++numGaps;
    }
  }

  return numGaps;
}

int suspectReadAlignmentMgr::numLigationAbuts()
{
  int numbs(0);

  vec<int> allReads = readIds();
  for ( int i=0; i<(int) allReads.size(); ++i )
  {
    if ( IsChimericCandidate( allReads[i] ) )
    {
      if ( LigationAbut( allReads[i] ) )
	++numbs;
    }
  }

  return numbs;
}

int suspectReadAlignmentMgr::numBasesLigationOverlaps()
{
  int numOvlps(0);

  vec<int> allReads = readIds();
  for ( int i=0; i<(int) allReads.size(); ++i )
  {
    if ( IsChimericCandidate( allReads[i] ) )
    {
      if ( LigationOverlap( allReads[i] ) )
      {
	pair<int,int> ligreg = ligationRegion( allReads[i] );
	numOvlps += (ligreg.first-ligreg.second);
      }
    }
  }

  return numOvlps;
}

int suspectReadAlignmentMgr::numBasesLigationGaps()
{
  int numGaps(0);

  vec<int> allReads = readIds();
  for ( int i=0; i<(int) allReads.size(); ++i )
  {
    if ( IsChimericCandidate( allReads[i] ) )
    {
      if ( LigationGap( allReads[i] ) )
      {
	pair<int,int> ligreg = ligationRegion( allReads[i] );
	numGaps += (ligreg.second-ligreg.first);
      }
    }
  }

  return numGaps;
}

int suspectReadAlignmentMgr::numBasesLigationAbuts()
{
  int numbs(0);

  vec<int> allReads = readIds();
  for ( int i=0; i<(int) allReads.size(); ++i )
  {
    if ( IsChimericCandidate( allReads[i] ) )
    {
      if ( LigationAbut( allReads[i] ) )
	++numbs;
    }
  }

  return numbs;
}


float suspectReadAlignmentMgr::aveQualityLigationOverlaps( vecqualvector &qual )
{

  int numBases(0);
  int sumQuals(0);

  vec<int> allReads = readIds();
  for ( int i=0; i<(int) allReads.size(); ++i )
  {
    if ( IsChimericCandidate( allReads[i] ) )
    {
      if ( LigationOverlap( allReads[i] ) )
      {
	pair<int,int> ligreg = ligationRegion( allReads[i] );
	numBases += (ligreg.first-ligreg.second);
	
	for ( int j=ligreg.second; j<ligreg.first; ++j )
	{
	  sumQuals += qual[allReads[i]][j];
	}
      }
    }
  }

  if ( numBases == 0 )
    return -1.0;

  return (float) sumQuals/numBases;
}


float suspectReadAlignmentMgr::aveQualityLigationGaps( vecqualvector &qual )
{
  int numBases(0);
  int sumQuals(0);

  vec<int> allReads = readIds();
  for ( int i=0; i<(int) allReads.size(); ++i )
  {
    if ( IsChimericCandidate( allReads[i] ) )
    {
      if ( LigationGap( allReads[i] ) )
      {
	pair<int,int> ligreg = ligationRegion( allReads[i] );
	numBases += (ligreg.second-ligreg.first);

	for ( int j=ligreg.first; j<ligreg.second; ++j )
	{
	  sumQuals += qual[allReads[i]][j];
	}
      }
    }
  }

  if ( numBases == 0 )
    return -1.0;

  return (float) sumQuals/numBases;
}



float suspectReadAlignmentMgr::aveQualityLigationAbuts( vecqualvector &qual )
{
  int numBases(0);
  int sumQuals(0);

  vec<int> allReads = readIds();
  for ( int i=0; i<(int) allReads.size(); ++i )
  {
    if ( IsChimericCandidate( allReads[i] ) )
    {
      if ( LigationAbut( allReads[i] ) )
      {
	pair<int,int> ligreg = ligationRegion( allReads[i] );
	++numBases;

	sumQuals += qual[allReads[i]][ligreg.first];

      }
    }
  }

  if ( numBases == 0 )
    return -1.0;

  return (float) sumQuals/numBases;
}
  
float suspectReadAlignmentMgr::aveQualityLigationOverlaps( int readId, vecqualvector &qual )
{

  int numBases(0);
  int sumQuals(0);

  
  if ( IsChimericCandidate( readId ) )
  {
    if ( LigationOverlap( readId ) )
    {
      pair<int,int> ligreg = ligationRegion( readId );
      numBases += (ligreg.first-ligreg.second);
      
      for ( int j=ligreg.second; j<ligreg.first; ++j )
	{
	  sumQuals += qual[readId][j];
	}
    }
  }
  else
    return -1.0;

  return (float) sumQuals/numBases;
}


float suspectReadAlignmentMgr::aveQualityLigationGaps( int readId, vecqualvector &qual )
{
  int numBases(0);
  int sumQuals(0);

  if ( IsChimericCandidate( readId ) )
  {
    if ( LigationGap( readId ) )
    {
      pair<int,int> ligreg = ligationRegion( readId );
      numBases += (ligreg.second-ligreg.first);

      for ( int j=ligreg.first; j<ligreg.second; ++j )
      {
	sumQuals += qual[readId][j];
      }
    }
  }
  else
    return -1.0;

  return (float) sumQuals/numBases;
}



float suspectReadAlignmentMgr::aveQualityLigationAbuts( int readId, vecqualvector &qual )
{
  int numBases(0);
  int sumQuals(0);

  if ( IsChimericCandidate( readId ) )
  {
    if ( LigationAbut( readId ) )
    {
      pair<int,int> ligreg = ligationRegion( readId );
      ++numBases;
      
      sumQuals += qual[readId][ligreg.first];

    }
  }
  else
    return -1.0;

  return (float) sumQuals/numBases;
}
  
