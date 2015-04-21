#include "tiled/BaggedReads.h"

bagged_reads::bagged_reads()
{
}

bagged_reads::bagged_reads( const String &strReadLocsFilename )
    : mbIndexed( false )
{
    READX( strReadLocsFilename, mVecLocs );
}

bagged_reads::bagged_reads( const vec<read_location> &vecReadLocs )
    : mVecLocs( vecReadLocs ),
      mbIndexed( false )
{
}

void 
bagged_reads::Load( const String &strReadLocsFilename )
{
    mVecLocs.clear();
    READX( strReadLocsFilename, mVecLocs );
    mbIndexed = false;
}

void 
bagged_reads::CleanUp( )
{
    mVecLocs.erase( remove_if( mVecLocs.begin(), mVecLocs.end(),
                               mem_fun_ref( &read_location::IsDead ) ),
                    mVecLocs.end() );

    sort( mVecLocs.begin(), mVecLocs.end() );
    mbIndexed = false;
}

void 
bagged_reads::SortAndSave( const String &strReadLocsFilename )
{
    this->CleanUp( );
    WriteLocs( strReadLocsFilename, mVecLocs );
}

int
bagged_reads::GetNumberOfContigs( ) const
{
    if ( ! mbIndexed )
        this->Index();

    return mContigToLocIdx.size();
}

void 
bagged_reads::GetIndicesForContig( const int contigId,
                                   vec<int> &vecLocIndices ) const
{
    if ( ! mbIndexed )
        this->Index();

    ForceAssertGe( contigId, 0 );
    ForceAssertLt( contigId, (int)mContigToLocIdx.size() );

    vecLocIndices = mContigToLocIdx[ contigId ];
}
        
void 
bagged_reads::GetIndicesForIntervalOnContig( const int contigId,
                                             const seq_interval &theInterval,
                                             vec<int> &vecLocIndices ) const
{
    if ( ! mbIndexed )
        this->Index();

    ForceAssertGe( contigId, 0 );
    ForceAssertLt( contigId, (int)mContigToLocIdx.size() );

    vecLocIndices.clear();

    const vec<int> &vecAllIndices = mContigToLocIdx[ contigId ];
    
    vec<int>::const_iterator idxIter = vecAllIndices.begin();
    while ( idxIter != vecAllIndices.end() &&
            mVecLocs[ *idxIter ].StartOnContig() < theInterval.End() )
    {
        if ( mVecLocs[ *idxIter ].StopOnContig() >= theInterval.Begin() )
            vecLocIndices.push_back( *idxIter );
        ++idxIter;
    }
}
        
void
bagged_reads::GetIndicesForRead( const int readId,
                                 vec<int> &vecLocIndices ) const
{
    if ( ! mbIndexed )
        this->Index();

    ForceAssertGe( readId, 0 );
    ForceAssertLt( readId, (int) mContigToLocIdx.size() );

    vecLocIndices = mReadToLocIdx[ readId ];
}


struct order_read_location_indices : public binary_function<int,int,bool>
{
    order_read_location_indices( const vec<read_location> &vecReadLoc )
        : mrVecReadLocs( vecReadLoc ) {}
    
    bool operator() ( const int lhs, const int rhs )
    {
        return ( mrVecReadLocs[ lhs ] < mrVecReadLocs[ rhs ] );
    }
    
  private:
    const vec<read_location> &mrVecReadLocs;
};


int
bagged_reads::AddLocation( const read_location &newLoc )
{
    ForceAssertGe( newLoc.Contig(), 0 );
    
    int newLocIdx = mVecLocs.size();
    mVecLocs.push_back( newLoc );

    if ( mbIndexed )
    {
        if ( (int)mContigToLocIdx.size() <= newLoc.Contig() )
            mContigToLocIdx.resize( newLoc.Contig() + 1 );
        mContigToLocIdx[ newLoc.Contig() ].push_back( newLocIdx );

        vec<int> &thisContigIdxs = mContigToLocIdx[ newLoc.Contig() ];
        sort( thisContigIdxs.begin(), thisContigIdxs.end(),
              order_read_location_indices( mVecLocs ) );
        
        if ( (int)mReadToLocIdx.size() <= newLoc.ReadId() )
            mReadToLocIdx.resize( newLoc.ReadId() + 1 );
        mReadToLocIdx[ newLoc.ReadId() ].push_back( newLocIdx );
    }
    
    return newLocIdx;
}

void
bagged_reads::RemoveLocationAt( const int locIndex )
{
    read_location &theLoc = mVecLocs[ locIndex ];
    
    if ( mbIndexed )
    {
        const int contigId = theLoc.Contig();
        if ( contigId >= 0 )
        {
            vec<int> &thisMap = mContigToLocIdx[ contigId ];
            
            vec<int>::iterator locIdxIter = find( thisMap.begin(), thisMap.end(), 
                                                  locIndex );
            if ( locIdxIter != thisMap.end() )
                thisMap.erase( locIdxIter );
        }

        const int readId = theLoc.ReadId();
        if ( readId >= 0 ) 
        {
            vec<int> &thisMap = mReadToLocIdx[ readId ];
            
            vec<int>::iterator locIdxIter = find( thisMap.begin(), thisMap.end(), 
                                                  locIndex );
            if ( locIdxIter != thisMap.end() )
                thisMap.erase( locIdxIter );
        }
    }

    theLoc.Kill();
}

void
bagged_reads::RemoveContig( const int contigId )
{
    if ( ! mbIndexed )
        this->Index();
    
    vec<int> &locsToKill = mContigToLocIdx[ contigId ];
    
    for ( vec<int>::iterator locToKillIter = locsToKill.begin();
          locToKillIter != locsToKill.end(); ++locToKillIter )
    {
        read_location &theLoc = mVecLocs[ *locToKillIter ];
        theLoc.Kill();
    }

    locsToKill.clear();
}

void
bagged_reads::JoinContigs( const int firstContigId,
                           const int secondContigId,
                           const int offsetOfSecondFromFirst )
{
    if ( ! mbIndexed )
        this->Index();

    ForceAssertGe( firstContigId, 0 );
    ForceAssertLt( firstContigId, (int)mContigToLocIdx.size() );

    ForceAssertGe( secondContigId, 0 );
    ForceAssertLt( secondContigId, (int)mContigToLocIdx.size() );

    vec<int> &firstContigMap =  mContigToLocIdx[ firstContigId ];
    vec<int> &secondContigMap = mContigToLocIdx[ secondContigId ];

    copy( secondContigMap.begin(), secondContigMap.end(),
          back_inserter( firstContigMap ) );

    vec<int>::iterator secondMapIter = secondContigMap.begin();
    for ( ; secondMapIter != secondContigMap.end(); ++secondMapIter )
    {
        read_location &aLoc = mVecLocs[ *secondMapIter ];
        aLoc.SetContig( firstContigId );
        aLoc.SetStartOnContig( aLoc.StartOnContig() + offsetOfSecondFromFirst );
    }

    secondContigMap.clear();

    sort( firstContigMap.begin(), firstContigMap.end(),
          order_read_location_indices( mVecLocs ) );

    this->NormalizeLocs( firstContigId );
}

int
bagged_reads::SplitContig( const int contigId,
                           const vec<int> &locsToKeep )
{
    if ( ! mbIndexed )
        this->Index();

    set<int> sortedLocsToKeep( locsToKeep.begin(), locsToKeep.end() );
    
    ForceAssertGe( contigId, 0 );
    ForceAssertLt( contigId, (int)mContigToLocIdx.size() );

    int newContigId = mContigToLocIdx.size();
    mContigToLocIdx.push_back( vec<int>() );

    // We have to get this reference after the above push_back() to
    // avoid possibly invalidating the reference.
    vec<int> &origContigMap =  mContigToLocIdx[ contigId ];

    sort( origContigMap.begin(), origContigMap.end() );

    vec<int> &newContigMap = mContigToLocIdx.back();

    set_difference( origContigMap.begin(), origContigMap.end(),
                    sortedLocsToKeep.begin(), sortedLocsToKeep.end(),
                    back_inserter( newContigMap ) );

    vec<int> oldContigMap;
    
    set_intersection( origContigMap.begin(), origContigMap.end(),
                      sortedLocsToKeep.begin(), sortedLocsToKeep.end(),
                      back_inserter( oldContigMap ) );

    vec<int>::iterator newMapIter;
    for ( newMapIter = newContigMap.begin(); newMapIter != newContigMap.end(); ++newMapIter )
    {
        read_location &aLoc = mVecLocs[ *newMapIter ];
        aLoc.SetContig( newContigId );
    }

    origContigMap.swap( oldContigMap );
    
    sort( origContigMap.begin(), origContigMap.end(),
          order_read_location_indices( mVecLocs ) );

    sort( newContigMap.begin(), newContigMap.end(),
          order_read_location_indices( mVecLocs ) );

    this->NormalizeLocs( contigId );
    this->NormalizeLocs( newContigId );

    return newContigId;
}


int
bagged_reads::SplitNonContigs( )
{
    int numSplits = 0;

    for ( unsigned int contigId = 0; contigId < mContigToLocIdx.size(); ++contigId )
    {
        vec<int> &locIndices = mContigToLocIdx[ contigId ];

        if ( locIndices.empty() ) 
            continue;

        vec<int>::iterator locIndicesIter = locIndices.begin();
        
        int maxStop = mVecLocs[ *locIndicesIter ].StopOnContig();

        bool wasSplit = false;

        while ( locIndicesIter != locIndices.end() )
            if ( mVecLocs[ *locIndicesIter ].StartOnContig() > maxStop )
            {
                vec<int> vecLocsToKeep( locIndices.begin(), locIndicesIter );
                this->SplitContig( contigId, vecLocsToKeep ); // invalidates reference locIndices
                wasSplit = true;
                break;
            }
            else
                ++locIndicesIter;

        if ( wasSplit )
        {
            ++numSplits;
            break;
        }
    }

    return numSplits;
}


void
bagged_reads::RemoveEmptyContigs( )
{
    unsigned int newContigId = 0;
    unsigned int oldContigId = 0;

    while ( oldContigId < mContigToLocIdx.size() )
    {
        if ( mContigToLocIdx[ oldContigId ].empty() )
            ++oldContigId;
        else
        {
            if ( oldContigId != newContigId )
            {
                vec<int> &vecOldContigLocIdxs = mContigToLocIdx[ oldContigId ];
                for ( vec<int>::iterator locIdxIter = vecOldContigLocIdxs.begin();
                      locIdxIter != vecOldContigLocIdxs.end(); ++locIdxIter )
                    mVecLocs[ *locIdxIter ].SetContig( newContigId );

                mContigToLocIdx[ newContigId ] = vecOldContigLocIdxs;
            }
            ++newContigId, ++oldContigId;
        }
    }
    mContigToLocIdx.resize( newContigId );
}


void
bagged_reads::NormalizeLocs( const int contigId )
{
    if ( ! mbIndexed )
        this->Index();

    ForceAssertGe( contigId, 0 );
    ForceAssertLt( contigId, (int)mContigToLocIdx.size() );
    
    vec<int> &contigMap =  mContigToLocIdx[ contigId ];

    if ( contigMap.empty() )
        return;

    int firstLocIdx = *(contigMap.begin());
    read_location &firstLoc = mVecLocs[ firstLocIdx ];

    int leftmostStart = firstLoc.StartOnContig();
    int rightmostStop = firstLoc.StopOnContig();

    vec<int>::iterator mapIter;
    for ( mapIter = contigMap.begin(); mapIter != contigMap.end(); ++mapIter )
    {
        read_location &aLoc = mVecLocs[ *mapIter ];
        leftmostStart = min( leftmostStart, aLoc.StartOnContig() );
        rightmostStop = max( rightmostStop, aLoc.StopOnContig() );
    }

    int newLengthOfContig = rightmostStop - leftmostStart + 1;

    for ( mapIter = contigMap.begin(); mapIter != contigMap.end(); ++mapIter )
    {
        read_location &aLoc = mVecLocs[ *mapIter ];
        aLoc.SetStartOnContig( aLoc.StartOnContig() - leftmostStart );
        aLoc.SetLengthOfContig( newLengthOfContig );
    }
}
    

void 
bagged_reads::Index() const
{
    int maxContigId = -1;
    int maxReadId = -1;
    for ( unsigned int locIdx = 0; locIdx < mVecLocs.size(); ++locIdx )
    {
        const read_location &aLoc = mVecLocs[ locIdx ];

        maxContigId = max( maxContigId, aLoc.Contig() );
        maxReadId   = max( maxReadId,   aLoc.ReadId() );
    }

    mContigToLocIdx.resize( maxContigId + 1 );
    for_each( mContigToLocIdx.begin(), 
              mContigToLocIdx.end(),
              mem_fun_ref( &vec<int>::clear ) );

    mReadToLocIdx.resize( maxReadId + 1 );
    for_each( mReadToLocIdx.begin(), 
              mReadToLocIdx.end(),
              mem_fun_ref( &vec<int>::clear ) );

    for ( unsigned int locIdx = 0; locIdx < mVecLocs.size(); ++locIdx )
    {
        const read_location &aLoc = mVecLocs[ locIdx ];
        
        const int contigId = aLoc.Contig();
        if ( contigId < 0 ) continue;

        const int readId   = aLoc.ReadId();
        if ( readId < 0 ) continue;

        mContigToLocIdx[ contigId ].push_back( locIdx );
        mReadToLocIdx[ readId ].push_back( locIdx );
    }

    for ( unsigned int contigId = 0; contigId < mContigToLocIdx.size(); ++contigId )
        sort( mContigToLocIdx[ contigId ].begin(),
              mContigToLocIdx[ contigId ].end(),
              order_read_location_indices( mVecLocs ) );

    mbIndexed = true;
}
