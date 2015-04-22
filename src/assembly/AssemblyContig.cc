///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#include "assembly/AssemblyContig.h"
#include "assembly/ContigDataManager.h"
#include "assembly/AssemblyRead.h"

int
Contig::GetLength( ) const
{
    return GetMgrPtr()->GetLength( GetId() );
}

const basevector &
Contig::GetBases( ) const
{
    return GetMgrPtr()->GetBases( GetId() );
}

const qualvector &
Contig::GetQuals( ) const
{
    return GetMgrPtr()->GetQuals( GetId() );
}

int
Contig::GetNumReadLocations( ) const
{
    return GetMgrPtr()->GetNumReadLocations( GetId() );
}

void
Contig::GetReadLocations( vec<ReadLocation> &vecRLs ) const
{
    GetMgrPtr()->GetReadLocations( GetId(), vecRLs ); 
}

void
Contig::GetSelfLocations( vec<ContigLocation> &vecCLs ) const
{
    GetMgrPtr()->GetContigLocations( GetId(), vecCLs ); 
}

#if __GNUC__ > 2
#include <ext/hash_set>
using __gnu_cxx::hash_set;
#else
#include <hash_set>
#endif

void
Contig::GetLinks( vec<Link> &vecLinks, bool externalOnly ) const
{
    vecLinks.clear();

    vec<ContigLocation> vecContigLocs;
    vec<ReadLocation>   vecReadLocs;
    vec<ReadLocation>   vecPartnerLocs;
    vec<ContigLocation> vecPartnerContigLocs;

    hash_set<ReadToken> setUsedPartners;
    
    this->GetSelfLocations( vecContigLocs );

    // If this contig is not placed, we still want to get its links.
    if ( vecContigLocs.empty() )
        vecContigLocs.push_back( ContigLocation( SuperToken(), *this, 
                                                 Interval( 0, this->GetLength() ),
                                                 orient_FW ) );

    vec<ContigLocation>::iterator contigLocIter = vecContigLocs.begin();
    for ( ; contigLocIter != vecContigLocs.end(); ++contigLocIter )
    {
        Contig theContig = contigLocIter->GetContig();

        theContig.GetReadLocations( vecReadLocs );

        vec<ReadLocation>::iterator readLocIter = vecReadLocs.begin();
        for ( ; readLocIter != vecReadLocs.end(); ++readLocIter )
        {
            Read theRead = readLocIter->GetRead();

            if ( setUsedPartners.count( theRead ) )
                continue;

            ReadPair thePair = theRead.GetPair();
            if ( ! thePair.IsValid() )
                continue;

            Read thePartner = thePair.GetOtherRead( theRead );

            thePartner.GetLocations( vecPartnerLocs );

            if ( ! vecPartnerLocs.size() )
                continue;

            vec<ReadLocation>::iterator partnerLocIter = vecPartnerLocs.begin();
            for( ; partnerLocIter != vecPartnerLocs.end(); ++partnerLocIter )
            {
                Contig partnerContig = partnerLocIter->GetContig();

                if ( externalOnly && partnerContig == theContig )
                  continue;

                partnerContig.GetSelfLocations( vecPartnerContigLocs );

                // If the partner contig is not placed, we still want to get its links.
                if ( vecPartnerContigLocs.empty() )
                    vecPartnerContigLocs.push_back( ContigLocation( SuperToken(), partnerContig, 
                                                                    Interval( 0, partnerContig.GetLength() ),
                                                                    orient_FW ) );

                vec<ContigLocation>::iterator partnerContigLocIter = vecPartnerContigLocs.begin();
                for ( ; partnerContigLocIter != vecPartnerContigLocs.end(); ++partnerContigLocIter )
                    vecLinks.push_back( Link( *contigLocIter, *readLocIter,
                                              *partnerContigLocIter, *partnerLocIter,
                                              thePair ) );
            }

            setUsedPartners.insert( thePartner );
        }
    }
}
                

void
Contig::GetLinks( vec<Link> &vecLinks ) const
{
  this->GetLinks( vecLinks, false );
}


void
Contig::GetExternalLinks( vec<Link> &vecLinks ) const
{
  this->GetLinks( vecLinks, true );
}
