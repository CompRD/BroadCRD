//Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#ifndef COMP_MASTER_VEC_TEMPLATE_H
#define COMP_MASTER_VEC_TEMPLATE_H

#include <fcntl.h>
#include <sys/types.h>

#include "CoreTools.h"
#include "Intvector.h"
#include "CompMasterVec.h"
#include "HashLinearRemoteStorage.h"
#include "system/file/FileReader.h"
#include "system/file/FileWriter.h"

using std::cout;
using std::endl;

template<class X, class H>
const size_t compmastervec<X,H>::DYNAMIC_DATA_SIZE = 10000;
template<class X, class H>
const size_t compmastervec<X,H>::OBJECTS = 10000;
template<class X, class H>
const size_t compmastervec<X,H>::UNIQUE_OBJECTS = 500;


template<class X, class H>
compmastervec<X,H>::compmastervec( ):
  m_mvec(),
  m_hashtable(&m_mvec, UNIQUE_OBJECTS),
  m_indices()
{
  Reserve(DYNAMIC_DATA_SIZE, UNIQUE_OBJECTS,OBJECTS);
}

template<class X, class H>
compmastervec<X,H>::compmastervec( size_t n ):
  m_mvec(),
  m_hashtable(&m_mvec, UNIQUE_OBJECTS),
  m_indices()
{
  Reserve(DYNAMIC_DATA_SIZE, UNIQUE_OBJECTS, n);
}


template<class X, class H>
compmastervec<X,H>::compmastervec( const String &filename ):
  m_mvec(),
  m_hashtable(&m_mvec, UNIQUE_OBJECTS),
  m_indices()
{
  ReadAll( filename );
}


template<class X, class H>
compmastervec<X,H>::compmastervec( const MasterVec<X> & mvec ):
  m_mvec(),
  m_hashtable(&m_mvec, UNIQUE_OBJECTS),
  m_indices()
{
  Append(mvec,0,mvec.size());
}


template<class X, class H>
compmastervec<X,H>::~compmastervec( )
{
}

template<class X, class H>
void
compmastervec<X,H>::clear( )
{
  m_mvec.clear();
  m_indices.clear();
  m_hashtable.clear();
}

template<class X, class H>
void
compmastervec<X,H>::destroy( )
{
  m_mvec.destroy();
  m_indices.clear();
  m_hashtable.clear();
}

template<class X, class H>
void
compmastervec<X,H>::resize( size_t n )
{
  m_indices.resize(n, HashLinearRemoteStorage<X,H>::EMPTY);
}


template<class X, class H>
void
compmastervec<X,H>::reserve( size_t n )
{
  m_mvec.reserve(n);
}


template<class X, class H>
void
compmastervec<X,H>::Reserve( size_t, size_t unique_objects,
                               size_t total_objects )
{
    m_mvec.reserve(unique_objects);
    m_hashtable.Reserve(unique_objects);
    m_indices.reserve(total_objects);
}

template<class X, typename H>
void
compmastervec<X,H>::CreateReverseIndex(VecIntVec& rindex) {
  vec<IntVec> v(UniqueSize());
  const size_t S = size();
  const size_t U = UniqueSize();
  for (size_t i = 0; i != S; ++i ) {
    v[m_indices[i]].push_back(i);
  }

  rindex.clear();
  rindex.Reserve(S, U);
  for (size_t i = 0; i != U; ++i) {
    rindex.push_back(v[i]);
  }
}

template<class X, typename H>
void
compmastervec<X,H>::Append( const compmastervec<X,H> & orig, size_t from,
                              size_t to) {
  //AssertGe(from,0);
  AssertGe(to,from);
  AssertGe(orig.size(),to);
  vec<int> entries(to-from);
  for (size_t i = from; i != to; ++i) {
    entries[i-from] = static_cast<int>(i);
  }
  Append(orig, entries);
}

template<class X, typename H>
void
compmastervec<X,H>::Append( const compmastervec<X,H> & orig,
                              const vec<int> & entries ) {
  for (size_t i = 0; i != entries.size(); ++i) {
    push_back(orig[entries[i]]);
  }
  //We don't know how many of the objects being appended are unique to us,
  // so it does
  //not make much sense to resize, given that compmastervec resizing is, in
  //general, efficient. So we don't add any space for orig objects.
}

template<class X, typename H>
void
compmastervec<X,H>::Append(const MasterVec<X> & orig, size_t from, size_t to){
  //AssertGe(from,0);
  AssertGe(to,from);
  AssertGe(orig.size(),to);
  vec<int> entries(to-from);
  for (size_t i = from; i != to; ++i) {
    entries[i-from] = static_cast<int>(i);
  }
  Append(orig, entries);
}

template<class X, typename H>
void
compmastervec<X,H>::Append( const MasterVec<X> & orig,
                              const vec<int> & entries) {
  for (size_t i = 0; i != entries.size(); ++i) {
    push_back(orig[entries[i]]);
  }
  //We don't know how many of the objects being appended are unique to us,
  // so it does
  //not make much sense to resize, given that compmastervec resizing is, in
  //general, efficient. So we don't add any space for orig objects.
}

///Need to deal with appending to an existing compmastervec!
template<class X, class H>
void compmastervec<X,H>::Write( const String& filename,
                                  size_t from,
                                  size_t to ) const
{
  AssertGe(to,from);
  AssertGt(size(),to);
  vec<int> entries(to - from);
  for (size_t i = from; i != to; ++i) {
    entries[i] = static_cast<int>(i);
  }
  Write(filename, entries);
}

/**
 * \todo find a more economical way to append a compmastervec to a file.
 */
template<class X, class H>
void compmastervec<X,H>::Write( const String& filename,
                                  const vec<int> & entries) const
{
  compmastervec<X, H> * temp = NULL;
  if (IsSomeSortOfFile(filename)) {
    temp = new compmastervec<X, H>(filename);
  } else {
    temp = new compmastervec<X, H>;
  }
  temp->Append(*this, entries);
  temp->WriteAll(filename);
  delete temp;
}

template<class X, class H>
void
compmastervec<X,H>::WriteAll( const String& filename ) const
{
  Remove( filename );
  Remove( filename + ".gz" );
  m_mvec.WriteAll(filename);
  FileWriter fw(filename, true);
  size_t isize = m_indices.size();
  fw.write(&isize, sizeof(isize));
  fw.write(&m_indices[0], isize * sizeof(isize));
  fw.close();
  SetFileCompBit(filename);
}

template<class X, class H>
void compmastervec<X,H>::SetFileCompBit(const String & filename,
                                          bool on) const {
  //Label the file as a CompMasterVec file
  FeudalControlBlock fcb(filename.c_str(),false,0);
  fcb.setCompressed(on);

  FileWriter fw(filename, true);
  fw.seek(0);
  fw.write(&fcb, sizeof(fcb));
  fw.close();
}

template<class X, class H>
void
compmastervec<X,H>::Read( const String& filename,
                            const vec<int>& entries,
                            size_t extra,
                            bool pre_reserved,
                            bool pre_resized )
{
  cerr << "not yet implemented " << endl;
  Assert(0);
}

template<class X, class H>
void
compmastervec<X,H>::ReadAll( const String& filename,
                               bool append )
{
  //This compmastervec must be empty! No appending allowed for now!
  AssertEq(append, false);
  AssertEq(m_indices.size(), 0u);
  //this call to ReadAll will ensure the file is mostly correct.
  m_mvec.ReadAll(filename);
  //Now get the m_indices.
  FileReader fr(filename.c_str());
  size_t fileSize;
  FeudalControlBlock fcb(fr,false,&fileSize);

  //Find how far into the file we need to read to get the m_indices.
  size_t vecOffset = fcb.getFixedOffset() + fcb.getNElements()*X::fixedDataLen();

  // Go to the beginning of the m_indices vec.
  size_t vecSize = 0;
  unsigned int wordLen = 4;
  // wordLength detection
  ConfigByteLength:
  fr.seek(vecOffset);
  vecSize = 0;
  fr.read(&vecSize, wordLen);
      if ((vecSize + 1) * wordLen + vecOffset != fileSize)
      {
        if (wordLen == 4) { wordLen = 8; goto ConfigByteLength; }
        if (wordLen == 8)
        {
            cout << "Can't figure out word length in CompMasterVec file: " << filename << endl;
            CRD::exit(1);
        }
      }
  // Reality Checks
  AssertGe(vecSize, m_mvec.size());
  AssertEq((vecSize + 1) * wordLen + vecOffset, fileSize);
  m_indices.resize(vecSize);
  // If our Word Length is four, the alignment of the data we are copying
  // is not of the same alignment as what we are copying it into.  So,
  // we need to copy value by value.
  if (wordLen == 4)
  {
    for(size_t idx = 0; idx < vecSize; ++idx)
      fr.read(&m_indices[idx], wordLen); // Little Endian assumed
  } else 
  {
      fr.read(&m_indices[0], vecSize * wordLen);
  }
  //OK, now we have all the data, we need to recreate the hash_table
  //for all the objects in the mastervec. The problem is that our hashtable
  //only knows how to insert an object by putting it into the mastervec...
  //what we need is a Load(const X & object, int index) method that will
  //act like Insert but not actually put the item in the mastervec, just
  //save the index we give it into the table. Ah ha!

  //Now create the hashtable
  m_hashtable.storage = &m_mvec;
  m_hashtable.Rehash();
}

#endif //COMPMASTERVECTEMPLATE_H
