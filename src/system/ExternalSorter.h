///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file ExternalSorter.h
 * \author tsharpe
 * \date Feb 9, 2010
 *
 * \brief Sorts more stuff than can be held in memory.
 * Implementation of an external sort.
 * Sorts more stuff than can be held all at once in memory.
 */
#ifndef EXTERNALSORTER_H_
#define EXTERNALSORTER_H_

#include "String.h"
#include "feudal/BinaryStream.h"
#include "feudal/PriorityQueue.h"
#include "system/ErrNo.h"
#include "system/SpinLockedData.h"
#include "system/System.h"
#include "system/file/TempFile.h"
#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <vector>

/// External sorter.
/// Insert blocks of sorted elements.
/// When done inserting, instantiate a Merger object, passing it the sorter.
/// The Merger object will present you your sorted results, one by one.
/// The type T must be binary read/writable.  The type C compares two T's with
/// the usual strict, weak ordering.
template <class T>
class ExternalSorter : private SpinLockedData
{
public:
    template < class Container = std::vector<T>, class Comp = std::less<T> >
    class Merger
    {
        class Reader
        {
        public:
            Reader( String const& filename )
            : mFilename(filename), mReader(filename.c_str())
            { mReader.read(&mRemaining); }

            ~Reader()
            { if ( unlink(mFilename.c_str()) )
              { ErrNo err;
                FatalErr("Unable to unlink " << mFilename << err); } }

            bool getNext( T& val )
            { bool result = false;
              if ( mRemaining )
              { mReader.read(&val); --mRemaining; result = true; }
              return result; }

        private:
            Reader( Reader const& ); // unimplemented -- no copying
            Reader& operator=( Reader const& ); // unimplemented -- no copying

            String mFilename;
            BinaryReader mReader;
            size_t mRemaining;
        };

        class IndirectComparator
        {
        public:
            IndirectComparator( Comp const& comp, Container const& vals )
            : mComp(comp), mVals(vals) {}

            // compiler-supplied copy constructor and destructor are OK

            bool operator()( size_t idx1, size_t idx2 )
            { return mComp(mVals[idx1],mVals[idx2]); }

        private:
            Comp mComp;
            Container const& mVals;
        };

    public:
        Merger( ExternalSorter& sorter,
                Container const& container = Container(),
                Comp const& comp=Comp() )
        : mNRecs(sorter.size()), mVals(container),
          mPQ(sorter.mFiles.size(),IndirectComparator(comp,mVals))
        { SpinLocker locker(sorter);
          std::vector<String> const& files = sorter.mFiles;
          size_t nFiles = files.size();
          mRdrs.reserve(nFiles);
          mVals.resize(nFiles);
          T* pVals = &mVals.front();
          typedef std::vector<String>::const_iterator Itr;
          for ( Itr itr(files.begin()), end(files.end()); itr != end; ++itr )
          { Reader* pReader = new Reader(*itr);
            if ( !pReader->getNext(*pVals) )
                delete pReader;
            else
            { size_t idx = mRdrs.size();
              mRdrs.push_back(pReader);
              pVals += 1;
              mPQ.push(idx); } }
          sorter.clear(); }

        ~Merger()
        { typedef typename std::vector<Reader*>::iterator Itr;
          for ( Itr itr(mRdrs.begin()), end(mRdrs.end()); itr != end; ++itr )
              delete *itr; }

        size_t size() const { return mNRecs; }

        bool getNext( T& val )
        { bool result = !mPQ.empty();
          if ( result )
          { size_t idx = mPQ.top();
            T& vvv = mVals[idx];
            val = vvv;
            if ( mRdrs[idx]->getNext(vvv) ) mPQ.replaceTop(idx);
            else mPQ.pop(); }
          return result; }

    private:
        Merger( Merger const& ); // unimplemented -- no copying
        Merger& operator=( Merger const& ); // unimplemented -- no copying

        size_t mNRecs;
        Container mVals;
        std::vector<Reader*> mRdrs;
        PriorityQueue<size_t,IndirectComparator> mPQ;
    };

    /// Make one.
    /// Temp files are named with your specified path prefix, with six chars
    /// tacked onto the end for uniqueness.
    ExternalSorter( String const& tmpFilePrefix )
    : mTmpFileTemplate(tmpFilePrefix+"XXXXXX"), mNRecs(0)
    {}

    ~ExternalSorter()
    { typedef std::vector<String>::const_iterator Itr;
      for ( Itr itr(mFiles.begin()), end(mFiles.end()); itr != end; ++itr )
          unlink(itr->c_str()); }

    size_t size() const { return mNRecs; }

    /// Insert a sorted sequence of elements.
    template <class Itr>
    void insert( Itr start, Itr end )
    { insert(start,end,static_cast<typename std::iterator_traits<Itr>::value_type*>(0)); }

    /// Insert a sorted sequence of elements.
    void insert( T* start, T* end )
    { size_t sz = end - start;
      if ( sz )
      { BinaryWriter writer(newFile(sz));
        writer.write(sz); writer.write(start,end); writer.close(); } }

private:
    ExternalSorter( ExternalSorter const& ); // unimplemented -- no copying
    ExternalSorter& operator=( ExternalSorter const& ); // unimplemented -- no copying

    template <class Itr>
    void insert( Itr start, Itr end, T* )
    { using std::distance; size_t sz = distance(start,end);
      if ( sz )
      { BinaryWriter writer(newFile(sz));
        writer.write(sz); writer.writeItr(start,end); writer.close(); } }

    char const* newFile( size_t sz )
    { SpinLocker locker(*this);
      mNRecs += sz;
      mFiles.push_back(temp_file::generateName(mTmpFileTemplate.c_str()));
      return mFiles.back().c_str(); }

    void clear()
    { mFiles.clear();
      mNRecs = 0; }

    template <class Container, class Comp> friend class Merger;

    String mTmpFileTemplate;
    std::vector<String> mFiles;
    size_t mNRecs;
};

#endif /* EXTERNALSORTER_H_ */
