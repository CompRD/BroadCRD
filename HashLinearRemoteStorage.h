// Copyright (c) 2000-2004 Whitehead Institute for Biomedical Research

#ifndef HASH_LINEAR_REMOTE_STORAGE_H
#define HASH_LINEAR_REMOTE_STORAGE_H

// =======================================================================
/** HashLinearRemoteStorage.h Hash table designed to save mallocs

@class HashLinearRemoteStorage
\author Pablo Alvarez

This hash table is designed to minimize the number of memory allocations,
and thereby work reasonably well even with a very large number of items.
The way it does this is by storing the actual items in a mastervec, and
keeping the indices of the items in its hash table. In addition, the hash
table uses linear probing to avoid any additional mallocs for individual 
items (as would be the case in separate chaining hashing). It will rehash
itself when there are too many items for its size, and it also checks to
see that the mastervec has enough space, and calls mastervec::Reserve() to 
increase the mastervec's size exponentially as needed.

A problem with rehashing is that we are providing indices into the 
hash table itself to the user when we Insert(), and those indices are used
in GetItem(). So those indices should not change even when we rehash. The 
solution to this is to simply return to the user the index into the mastervec,
not the index into the table! Indices into the mastervec never change,
because a mastervec is only capable of push_back(). The key insight is 
that the only role of the table is for searching during Insert(), but 
once we know whether we have the X or not, we can get it directly from 
the mastervec.

Note that this class knows nothing about where the mastervec comes from. 
That is supplied by whoever is using the hashtable. Ideally that 
code has some idea of how big the mastervec should be, but the hashtable
is capable of resizing it economically as needed.

*/

#include "Vec.h"
#include "feudal/MasterVec.h"

template<class X, typename H> class compmastervec;
/*The template types X and T are the same as for mastervec.
  Type X must provide operator==.
  Template type H needs to be a hash function that accepts Xs and a 
  table size M as input and returns an int that is less than M.
*/

template<class X, typename H, typename Eq = std::equal_to<X> >
class HashLinearRemoteStorage {
  friend class compmastervec<X,H>;

 public:
  static const size_t EMPTY; ///<marker value for empty cells in the table
  static const double FILL_FACTOR; ///< between 0 and 1, trigger for rehashing.
  static const double RESIZE_FACTOR; ///how much bigger to make the table.

  /**Constructor
   * @param mastervec: this is the location that will hold all the dynamic
   * data passed in. It is not owned by the hash table.
   * @param tablesize: initial size of the hash table. If you know
   * approximately how many items you are likely to need to hash (M), you can
   * set this to 2*M for maximum efficiency, but it will work efficiently
   * even if it has to expand and rehash.
   */
  HashLinearRemoteStorage(MasterVec<X> * mvec, size_t tableSize = 10000):
    storage(mvec),
    table(tableSize, EMPTY),
    m_hash(),
    m_eq()
  {
    Assert(mvec);
    if ((longlong) storage->capacity() < (longlong) table.size()) {
      storage->reserve(table.size());
    }
  }

  void Reserve(size_t size) {
    if (size > table.size()) {
      table.reserve(size);
    }
  }

  size_t size() {
    return table.size();
  }

  void clear() {
    table.clear();
  }

  /**Insert an item in the table, return the index within the mastervec.
   * Note that the table will actually contain, not the item, but the 
   * index to the item in the "storage" mastervec. It is this index that we
   * return, and it will always stay the same even if the table is rehashed.
   * @param X: object to be inserted
   * @param extra_space: extra space to be reserved for that object to grow
   * in the mastervec (as in mastervec::push_back()).
   * @param external: if true and if the object has never been seen before, 
   * the object's memory is self-owned.
   */ 
  size_t Insert(const X & x, int extra_space = 0, Bool external = False);
    
  pair<size_t, bool> InsertCheck(const X & x, int extra_space = 0, 
			      Bool external = False);

  /**Find the index for an item in the mastervec. Returns EMPTY if not found.
   */
  size_t Find(const X & x);

  /**get the item stored at index i in the mastervec.
   * This is the index that we returned to the user in Insert() or Find().
   */
  const X & GetItem(size_t i) {
    AssertLt(i, storage->size());
    return (*storage)[i]; 
  }

  /** We have run out of space, so we rehash every element into a bigger table.
   * Because we return mastervec indices to the user, and those do not change,
   * we don't need to worry about the internal structure of the table. This is
   * also used to create a new hashtable from an already existing mastervec,
   * for example when we are copying a complete mastervec or reading one
   * in from a file. 
   */
  void Rehash();

 private:
  ///Forbid copy 
  HashLinearRemoteStorage(const HashLinearRemoteStorage &);
  ///Forbid assignment
  HashLinearRemoteStorage& operator=(const HashLinearRemoteStorage &);
    
 private:
  MasterVec<X> * storage;///< this is where we store the actual data to
  ///avoid having to do a lot of mallocs

  vec<size_t> table; ///< this is the actual hash table where we store the indices into storage
  H m_hash; ///< hash object to use, declare here to construct only once.
  Eq m_eq; ///< equality object to use, declare here to construct only once.
};

template<class X, typename H, typename Eq>
const size_t HashLinearRemoteStorage<X,H,Eq>::EMPTY = -1; 
template<class X, typename H, typename Eq>
 const double HashLinearRemoteStorage<X,H,Eq>::FILL_FACTOR = 0.666; 
template<class X, typename H, typename Eq>
 const double HashLinearRemoteStorage<X,H,Eq>::RESIZE_FACTOR = 4.0; 

template<class X, typename H, typename Eq>
size_t HashLinearRemoteStorage<X,H,Eq>::Insert(const X & x, int extra_space,
                                           Bool external) {
  return (InsertCheck(x, extra_space, external)).first;
}

template<class X, typename H, typename Eq>
pair<size_t,bool> HashLinearRemoteStorage<X,H,Eq>::
InsertCheck(const X & x, int extra_space, Bool external) {
  const size_t M = table.size();
  size_t i = m_hash(x, M);
  while (table[i] != EMPTY) {
    if (m_eq(x, (*storage)[table[i]])) { 
      return make_pair(table[i], true); 
    } //found!
    if (M == ++i) i = 0; //avoid using %
  }

  storage->push_back(x);

  //save the position before we rehash, because it will change!
  size_t ret = table[i] = storage->size() - 1;
  if (ret > (M * FILL_FACTOR)) {
    Rehash(); 
  }
  return make_pair(ret, false);
}

template<class X, typename H, typename Eq>
void HashLinearRemoteStorage<X,H,Eq>::Rehash() {
  const size_t K = storage->size();
  table.assign(std::max(K, table.size()) * RESIZE_FACTOR, EMPTY);
  size_t M = table.size();
  if (storage->capacity() < M) {
    storage->reserve(M);
  }
  for (size_t i = 0; i != K; ++i) {
    size_t j = m_hash((*storage)[i], M);
    while (table[j] != EMPTY) { //linear probing
      if (M == ++j) j = 0; //avoid using %
    }
    table[j] = i;
  }
  //gettimeofday(&end, NULL);
  //double t = end.tv_sec - start.tv_sec + 
  //  (end.tv_usec - start.tv_usec) /1000000.0;
  //cout << "Rehash to size " << M << " took " << t << " seconds" << endl;
}

template<class X, typename H, typename Eq>
size_t HashLinearRemoteStorage<X,H,Eq>::Find(const X & x) {
  const size_t M = table.size();
  size_t i = m_hash(x, M);
  while (table[i] != EMPTY) {
    if (m_eq(x, (*storage)[table[i]]) ) { break; }
    if (M == ++i) i = 0; //avoid using %
  }
  return table[i]; //EMPTY if we never found a match.
}

#endif //HASH_LINEAR_REMOTE_STORAGE_H
