///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file PriorityQueue.h
 * \author tsharpe
 * \date Jul 21, 2011
 *
 * \brief
 */
#ifndef PRIORITYQUEUE_H_
#define PRIORITYQUEUE_H_

#include "system/Assert.h"
#include <algorithm>
#include <cstddef>
#include <functional>

// T must be copyable and default-constructible, and must be capable of being
// compared by the comparator.
template <class T, class Comp = std::less<T> >
class PriorityQueue
{
public:
    typedef size_t size_type;
    typedef T value_type;
    typedef T const& const_reference;

    explicit PriorityQueue( size_type maxCap, Comp const& comp = Comp() )
    : mComp(comp), mSize(0), mCap(maxCap), mArr(new T[maxCap]) {}

    PriorityQueue( PriorityQueue const& that ) : mArr(0) { *this = that; }
    ~PriorityQueue() { delete [] mArr; }

    PriorityQueue& operator=( PriorityQueue const& that )
    { delete [] mArr;
      mComp = that.mComp; mSize = that.mSize; mCap = that.mCap;
      mArr = new T[mCap]; }

    bool empty() const { return !mSize; }
    bool full() const { return mSize == mCap; }
    size_type size() const { return mSize; }

    const_reference top() const { Assert(!empty()); return mArr[0]; }

    PriorityQueue& push( const_reference val )
    { Assert(!full());
      mArr[mSize++] = val;
      keyUp();
      return *this; }

    PriorityQueue& pop()
    { Assert(!empty());
      if ( --mSize )
      { mArr[0]=mArr[mSize];
        keyDown(); }
      return *this; }

    PriorityQueue& replaceTop( const_reference val )
    { Assert(!empty());
      mArr[0] = val;
      keyDown();
      return *this; }

private:
    void keyUp()
    { T* arr = mArr - 1; // using 1-based indices
      size_type idx = mSize;
      size_type parent;
      while ( (parent = idx>>1) && mComp(arr[idx],arr[parent]) )
      { using std::swap; swap(arr[idx],arr[parent]); idx = parent; } }

    void keyDown()
    { T* arr = mArr - 1; // using 1-based indices
      size_type idx = 1;
      size_type child;
      while ( (child = idx<<1) <= mSize )
      { size_type idx2 = idx;
        if ( mComp(arr[child],arr[idx]) ) idx2 = child;
        if ( ++child <= mSize && mComp(arr[child],arr[idx2]) ) idx2 = child;
        if ( idx == idx2 ) break;
        using std::swap; swap(arr[idx],arr[idx2]);
        idx = idx2; } }

    Comp mComp;
    size_type mSize;
    size_type mCap;
    T* mArr;
};
#endif /* PRIORITYQUEUE_H_ */
