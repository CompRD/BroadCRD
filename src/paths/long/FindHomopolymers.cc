///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * FindHomopolymers.cc
 *
 *  Created on: Mar 27, 2013
 *      Author: tsharpe
 */

#include "paths/long/FindHomopolymers.h"

void HomopolymerFinder::explore( int vvv )
{
    vec<int> path;
    vec<int> const& froms = mpHBV->From(vvv);
    int nVs = froms.size();
    for ( int iii = 0; iii != nVs; ++iii )
    {
        int edgeId = mpHBV->EdgeObjectIndexByIndexFrom(vvv, iii);
        bvec const& edge = mpHBV->EdgeObject(edgeId);
        auto beg = edge.begin(), end = edge.end();
        auto itr = findStartPos(vvv, edge);
        while ( itr != end )
        {
            auto itr2 = std::find_if(itr, end, getNEPred(*itr));
            unsigned runLen = std::distance(itr, itr2);
            if ( itr2 == end )
            {
                path.push_back(edgeId);
                unsigned off = std::distance(beg, itr);
                extend(froms[iii], *itr, off, runLen, false, path);
                path.clear();
            }
            else if ( runLen >= mMinRunLen )
                mpOut->push_back(
                        HomopolymerRun(edgeId, std::distance(beg,itr),
                                        std::distance(itr2,end), runLen,
                                        false, *itr));
            itr = itr2;
        }
    }
}

bool HomopolymerFinder::extend( int vvv, unsigned char base, unsigned off,
                                  unsigned curLen, bool cyclic, vec<int>& path )
{
    bool productive = false;
    Pred_t pred = getNEPred(base);
    vec<int> const& froms = mpHBV->From(vvv);
    int nVs = froms.size();
    for ( int iii = 0; iii != nVs; ++iii )
    {
        int edgeId = mpHBV->EdgeObjectIndexByIndexFrom(vvv, iii);

        // cycle detection
        auto pathBeg = path.begin(), pathEnd = path.end();
        auto pathItr = std::find(pathBeg, pathEnd, edgeId);
        bool isCyclic = false;
        if ( pathItr == pathBeg )
            isCyclic = !off;
        else
            isCyclic = pathItr != pathEnd;
        if ( isCyclic )
        {
            cyclic = true;
            continue;
        }

        bvec const& edge = mpHBV->EdgeObject(edgeId);
        if ( edge.front() != base )
            continue;
        auto beg = edge.begin(), end = edge.end();
        auto itr = std::find_if(beg, end, pred);
        unsigned runLen = curLen + std::distance(beg, itr);
        if ( itr == end )
        {
            path.push_back(edgeId);
            productive = extend(froms[iii], base, off, runLen, cyclic, path);
            path.pop_back();
        }
        else if ( runLen >= mMinRunLen || cyclic )
        {
            path.push_back(edgeId);
            unsigned lastOff = std::distance(itr,end);
            mpOut->push_back(
                    HomopolymerRun(path, off, lastOff, runLen, cyclic, base));
            productive = true;
            path.pop_back();
        }
    }
    if ( !productive && (curLen >= mMinRunLen || cyclic) )
    {
        mpOut->push_back(
                HomopolymerRun(path, off, 0, curLen, cyclic, base));
        productive = true;
    }
    return productive;
}

void HomopolymerFinder::cluster( VecHomopolymerRunVec* pOut,
                                    VecHomopolymerRun const& runs,
                                    unsigned minClusterSize )
{
    if ( minClusterSize < 2 )
        pOut->reserve(runs.size());
    VecHomopolymerRun singleton(1);
    VecHomopolymerRun scratch;
    std::vector<bool> used(runs.size(),false);
    auto beg=runs.begin(), end=runs.end();
    auto pUsed = used.begin();
    for ( auto itr=beg; itr != end; ++itr,++pUsed )
    {
        if ( *pUsed ) continue;

        if ( itr->isSolo() )
        {
            if ( minClusterSize <= 1 )
            {
                singleton.back() = *itr;
                pOut->push_back(singleton);
            }
        }
        else
        {
            scratch.push_back(*itr);
            auto pUsed2 = pUsed;
            auto itr2 = itr;
            while ( ++itr2 != end )
            {
                if ( *++pUsed2 ) continue;
                for ( HomopolymerRun const& tst : scratch )
                {
                    if ( overlaps(tst,*itr2) )
                    {
                        scratch.push_back(*itr2);
                        *pUsed2 = true;
                        break;
                    }
                }
            }
            if ( scratch.size() >= minClusterSize )
                pOut->push_back(scratch);
            scratch.clear();
        }
    }
}

namespace
{

bool checkForRun( HyperBasevector const& hbv, int vvv, int idx, int v2,
                    vec<int>& path, unsigned char baseCode )
{
    bvec const& edge = hbv.EdgeObjectByIndexFrom(vvv,idx);
    auto end = edge.end();
    if ( std::find_if(edge.begin(hbv.K()-1),end,
            [baseCode]( unsigned char code )
            { return code != baseCode; }) != end )
        return false;
    vvv = hbv.From(vvv)[idx];
    if ( vvv == v2 || std::find(path.begin(),path.end(),vvv) != path.end() )
        return true;
    path.push_back(vvv);
    bool result = true;
    vec<int> downstreamNodes = hbv.From(vvv);
    int nNodes = downstreamNodes.size();
    for ( idx = 0; result && idx != nNodes; ++idx )
        result = checkForRun(hbv,vvv,idx,v2,path,baseCode);
    path.pop_back();
    return result;
}

bool checkForRun( HyperBasevector const& hbv, int vvv, int idx, int v2,
                    vec<int>& path, unsigned char* pBaseCode )
{
    bvec const& edge = hbv.EdgeObjectByIndexFrom(vvv,idx);
    unsigned char baseCode = edge[hbv.K()-1];
    *pBaseCode = baseCode;
    bvec::const_iterator end = edge.end();
    if ( std::find_if(edge.begin(hbv.K()),end,
            [baseCode]( unsigned char code )
            { return code != baseCode; }) != end )
        return false;
    vvv = hbv.From(vvv)[idx];
    if ( vvv == v2 || std::find(path.begin(),path.end(),vvv) != path.end() )
        return true;
    path.push_back(vvv);
    bool result = true;
    vec<int> downstreamNodes = hbv.From(vvv);
    int nNodes = downstreamNodes.size();
    for ( idx = 0; result && idx != nNodes; ++idx )
        result = checkForRun(hbv,vvv,idx,v2,path,baseCode);
    path.pop_back();
    return result;
}

}

bool HomopolymerFinder::isHomopolymerCell( HyperBasevector const& hbv,
                                            int v1, int v2 )
{
    unsigned K = hbv.K();
    vec<int> path(1,v1);
    vec<int> const& downstreamNodes = hbv.From(v1);
    int nNodes = downstreamNodes.size();
    unsigned char baseCode = 0;
    ForceAssertGt(nNodes,0);
    if ( nNodes == 1 )
        return checkForRun(hbv,v1,0,v2,path,&baseCode);
    else
    {
        bool result = checkForRun(hbv,v1,0,v2,path,&baseCode);
        for ( int idx = 1; result && idx != nNodes; ++idx )
            result = checkForRun(hbv,v1,idx,v2,path,baseCode);
        return result;
    }
}

#include "feudal/SmallVecDefs.h"
template class SmallVec<HomopolymerRun,MempoolAllocator<HomopolymerRun>>;

#include "feudal/OuterVecDefs.h"
template class OuterVec<VecHomopolymerRun>;
