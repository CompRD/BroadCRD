///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#include "assembly/PrefetchStrategy.h"

#include "assembly/AssemblyContig.h"
#include "assembly/ContigLocationManager.h"
#include "assembly/AssemblyRead.h"
#include "assembly/ReadLocationManager.h"
#include "assembly/Super.h"


void PrefetchReadsByContigStrategy::GetIdsToLoad( const int seedId, vec<int> &vecIdsToLoad )
{
    vecIdsToLoad.clear();

    vec<ReadLocation> vecLocationsOfRead;
    mpReadLocationMgr->GetByRead( seedId, vecLocationsOfRead );

    for ( unsigned int locIdx = 0; locIdx < vecLocationsOfRead.size(); ++locIdx )
    {
        int contigId = vecLocationsOfRead[ locIdx ].GetContig().GetId();
        vec<ReadLocation> vecLocationsInContig;
        mpReadLocationMgr->GetByContig( contigId, vecLocationsInContig );
        
        vec<Read> vecReads( vecLocationsInContig.size() );
        transform( vecLocationsInContig.begin(), vecLocationsInContig.end(),
                   vecReads.begin(),
                   mem_fun_ref( &ReadLocation::GetRead ) );

        transform( vecReads.begin(), vecReads.end(),
                   back_inserter( vecIdsToLoad ),
                   mem_fun_ref( &Read::GetId ) );

        vec<Read> vecPartners( vecReads.size() );
        transform( vecReads.begin(), vecReads.end(),
                   vecPartners.begin(),
                   mem_fun_ref( &Read::GetPartner ) );

        vecPartners.erase( remove_if( vecPartners.begin(), vecPartners.end(),
                                      not1( mem_fun_ref( &Read::IsValid ) ) ),
                           vecPartners.end() );

        transform( vecPartners.begin(), vecPartners.end(),
                   back_inserter( vecIdsToLoad ),
                   mem_fun_ref( &Read::GetId ) );
    }
    
    if ( vecIdsToLoad.empty() )
    {
        vecIdsToLoad.push_back( seedId );
        for ( int distance = 1; distance <= (int) mNeighborhood; ++distance )
            if ( seedId - distance >= 0 )
                vecIdsToLoad.push_back( seedId - distance );
    }

    sort( vecIdsToLoad.begin(), vecIdsToLoad.end() );
    vecIdsToLoad.erase( unique( vecIdsToLoad.begin(), vecIdsToLoad.end() ),
                        vecIdsToLoad.end() );

    // cout << "Found " << vecIdsToLoad.size() << " ids from seed id " << seedId << endl;
}
        

void PrefetchReadsBySuperStrategy::GetIdsToLoad( const int seedId, vec<int> &vecIdsToLoad )
{
    vecIdsToLoad.clear();
    vecIdsToLoad.push_back( seedId );

    vec<ReadLocation> vecReadLocations;
    mpReadLocationMgr->GetByRead( seedId, vecReadLocations );

    vec<ContigLocation> vecReadContigLocations;
    for ( unsigned int locIdx = 0; locIdx < vecReadLocations.size(); ++locIdx )
    {
        Contig theContig = vecReadLocations[ locIdx ].GetContig();

        vec<ContigLocation> vecThisContigLocations;
        theContig.GetSelfLocations( vecThisContigLocations );

        copy( vecThisContigLocations.begin(), vecThisContigLocations.end(),
              back_inserter( vecReadContigLocations ) );
    }

    vec<ContigLocation> vecAllContigLocations;
    for ( unsigned int locIdx = 0; locIdx < vecReadContigLocations.size(); ++locIdx )
    {
        Super theSuper = vecReadContigLocations[ locIdx ].GetSuper();

        vec<ContigLocation> vecThisSuperContigLocations;
        theSuper.GetContigLocations( vecThisSuperContigLocations );

        copy( vecThisSuperContigLocations.begin(), vecThisSuperContigLocations.end(),
              back_inserter( vecAllContigLocations ) );
    }

    for ( unsigned int locIdx = 0; locIdx < vecAllContigLocations.size(); ++locIdx )
    {
        int contigId = vecAllContigLocations[ locIdx ].GetContig().GetId();

        vec<ReadLocation> vecLocationsInContig;
        mpReadLocationMgr->GetByContig( contigId, vecLocationsInContig );
        
        vec<ReadToken> vecReads( vecLocationsInContig.size() );
        transform( vecLocationsInContig.begin(), vecLocationsInContig.end(),
                   vecReads.begin(),
                   mem_fun_ref( &ReadLocation::GetRead ) );

        transform( vecReads.begin(), vecReads.end(),
                   back_inserter( vecIdsToLoad ),
                   mem_fun_ref( &ReadToken::GetId ) );
    }

    if ( vecIdsToLoad.empty() )
        for ( int distance = 0; distance <= (int) mNeighborhood; ++distance )
            if ( seedId - distance >= 0 )
                vecIdsToLoad.push_back( seedId - distance );

    sort( vecIdsToLoad.begin(), vecIdsToLoad.end() );
    vecIdsToLoad.erase( unique( vecIdsToLoad.begin(), vecIdsToLoad.end() ),
                        vecIdsToLoad.end() );

    cout << "Found " << vecIdsToLoad.size() << " ids from seed id " << seedId << endl;
}



void PrefetchContigsBySuperStrategy::GetIdsToLoad( const int seedId, vec<int> &vecIdsToLoad )
{
    vecIdsToLoad.clear();
    vecIdsToLoad.push_back( seedId );

    vec<ContigLocation> vecContigLocations;
    mpContigLocationMgr->GetByContig( seedId, vecContigLocations );

    vec<ContigLocation> vecAllContigLocations;
    for ( unsigned int locIdx = 0; locIdx < vecContigLocations.size(); ++locIdx )
    {
        Super theSuper = vecContigLocations[ locIdx ].GetSuper();

        vec<ContigLocation> vecThisSuperContigLocations;
        theSuper.GetContigLocations( vecThisSuperContigLocations );

        copy( vecThisSuperContigLocations.begin(), vecThisSuperContigLocations.end(),
              back_inserter( vecAllContigLocations ) );
    }

    vec<ContigToken> vecContigs( vecAllContigLocations.size() );
    transform( vecAllContigLocations.begin(), vecAllContigLocations.end(),
               vecContigs.begin(),
               mem_fun_ref( &ContigLocation::GetContig ) );

    transform( vecContigs.begin(), vecContigs.end(),
               back_inserter( vecIdsToLoad ),
               mem_fun_ref( &ContigToken::GetId ) );

    sort( vecIdsToLoad.begin(), vecIdsToLoad.end() );
    vecIdsToLoad.erase( unique( vecIdsToLoad.begin(), vecIdsToLoad.end() ),
                        vecIdsToLoad.end() );

    cout << "Found " << vecIdsToLoad.size() << " ids from seed id " << seedId << endl;
}
