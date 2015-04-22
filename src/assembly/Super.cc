///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "assembly/Super.h"

#include "assembly/AssemblyContig.h"
#include "assembly/ContigLocation.h"
#include "assembly/Interval.h"
#include "assembly/AssemblyRead.h"
#include "assembly/SuperDataManager.h"

void
Super::GetContigLocations( vec<ContigLocation> &vecCLs ) const
{
    GetMgrPtr()->GetContigLocations( GetId(), vecCLs );
}

int
Super::GetNumContigLocations( ) const
{
    return GetMgrPtr()->GetNumContigLocations( GetId() );
}

int 
Super::GetLength( ) const
{
    return GetMgrPtr()->GetLength( GetId() );
}

longlong
Super::GetSumOfContigLengths( ) const
{
    return GetMgrPtr()->GetSumOfContigLengths( GetId() );
}

int
Super::GetNumReadLocations( ) const
{
    vec<ContigLocation> contigLocs;
    this->GetContigLocations( contigLocs );
    
    longlong numLocs = 0;
    for ( vec<ContigLocation>::iterator iLoc = contigLocs.begin();
          iLoc != contigLocs.end(); ++iLoc )
      numLocs += Contig(iLoc->GetContig()).GetNumReadLocations();

    return numLocs;
}

void
Super::Print( ostream &out ) const
{
    vec<ContigLocation> vecContigLocs;
    this->GetContigLocations( vecContigLocs );
    sort( vecContigLocs.begin(), vecContigLocs.end() );
    
    out << *this << " (l=" << this->GetLength() << ")" << endl;
    
    copy( vecContigLocs.begin(), vecContigLocs.end(),
          ostream_iterator<ContigLocation>( out, "\n" ) );
}

#if __GNUC__ > 2
#include <ext/hash_set>
using __gnu_cxx::hash_set;
#else
#include <hash_set>
#endif

void
Super::GetLinks( vec<Link> &vecLinks ) const
{
    vecLinks.clear();

    vec<ContigLocation> vecContigLocs;
    vec<ReadLocation>   vecReadLocs;
    vec<ReadLocation>   vecPartnerLocs;
    vec<ContigLocation> vecPartnerContigLocs;

    hash_set<ReadToken> setUsedPartners;

    this->GetContigLocations( vecContigLocs );

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

            // if ( vecPartnerLocs.size() > 1 )
            //     cout << thePartner << " has " << vecPartnerLocs.size() << " locations!" << endl;

            vec<ReadLocation>::iterator partnerLocIter = vecPartnerLocs.begin();
            for( ; partnerLocIter != vecPartnerLocs.end(); ++partnerLocIter )
            {
                Contig partnerContig = partnerLocIter->GetContig();

                partnerContig.GetSelfLocations( vecPartnerContigLocs );

                // if ( vecPartnerContigLocs.size() > 1 )
                //     cout << partnerContig << " has multiple locations!" << endl;

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
                
