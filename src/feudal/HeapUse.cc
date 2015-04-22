///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * HeapUse.cc
 *
 *  Created on: Apr 22, 2014
 *      Author: tsharpe
 */
#if 0
#include "feudal/HeapUse.h"
#include "feudal/HashSet.h"
#include "math/Hash.h"
#include <cstddef>
#include <new>

namespace
{

class HeapUse
{
public:
    HeapUse()=default;
    HeapUse( void* mem, size_t siz=0 )
    : mMem(reinterpret_cast<size_t>(mem)), mSiz(siz) {}

    size_t getAddr() const { return mMem; }
    size_t getSize() const { return mSiz; }

    friend bool operator==( HeapUse const& hu1, HeapUse const& hu2 )
    { return hu1.getAddr()==hu2.getAddr(); }

    struct Hasher
    { typedef HeapUse argument_type;
      size_t operator()( HeapUse const& hu ) const
      { char const* ppp = reinterpret_cast<char const*>(&hu.mMem);
        return FNV1a(ppp,ppp+sizeof(size_t)); } };
    friend struct Hasher;

private:
    size_t mMem;
    size_t mSiz;
};

typedef HashSet<HeapUse,HeapUse::Hasher> HeapUseSet;
HeapUseSet gHeapUseSet(5000000);
class HUSPtr
{
public:
    HUSPtr() { mpHUS = &gHeapUseSet; }
    ~HUSPtr() { mpHUS = nullptr; }
    HeapUseSet* mpHUS;
};
HUSPtr gPtr;

}

void* operator new( size_t siz )
{
    if ( !siz ) siz = 1ul;
    void* mem;
    while ( !(mem  = malloc(siz)) )
    {
        // std::get_new_handler is not yet defined by gcc 4.7.2
        std::new_handler handler = std::set_new_handler(nullptr);
        std::set_new_handler(handler);
        if ( !handler )
            throw std::bad_alloc();
        (*handler)();
    }
    if ( gPtr.mpHUS )
        gPtr.mpHUS->insertUniqueValue(HeapUse(mem,siz));
    return mem;
}

void operator delete( void* mem )
{
    if ( !mem ) return;
    if ( gPtr.mpHUS )
        gPtr.mpHUS->remove(HeapUse(mem));
    free(mem);
}

void reportMemoryInUse()
{
    std::cout << "Memory In Use:" << '\n';
    for ( auto const& hhs : *gPtr.mpHUS )
        for ( HeapUse const& hu : hhs )
            std::cout << hu.getAddr() << '\t' << hu.getSize() << '\n';
    std::cout << std::endl;
}
#endif
