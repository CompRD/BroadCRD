//Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#ifndef COMPMASTERVEC_H
#define COMPMASTERVEC_H

// =======================================================================
/** CompMasterVec.h mastervec-like class for data that is very redundant.

@class compmastervec
\author Pablo Alvarez

The compmastervec class works like mastervec,
but is intended to be used when the data to be stored are very redundant. 
For example, we want to store 40 M Strings, but we know that in many cases
30000 instances of the String will actually be the same String.

Note: mastervec is a vector whose elements are dynamically allocated. The
trick is that mastervec initially allocates and manages all the memory,
so we can have for instance 40 M different Strings without doing 40M 
mallocs which would really mess up the system. 

This class only uses single file storage, unlike mastervec
which has the option of using one or three files for storage.

Performance:
  Insertion of N objects, of which K are unique, is of O(N + KlogK).
  If you have a good guess as to K, you can reduce that to O(N+K).
  Reading an object is O(1);

Template types: X and A are the same as for mastervec and must meet the 

  In addition, type X must provide operator== for hashing.

  Template type H needs to be a hash function that accepts Xs and a 
  table size M as input and returns an int that is less than M, such
  as the default CMVHash struct below. The default hash should work fine for
  most things.

Additional use of compmastervec as a memory efficient hash table:
You can use compmastervec even if you don't want to compress information, 
 (that is, if all objects are unique) but just want to be able to find it 
quickly and avoid multiple mallocs. That is what the FindUnique method is 
for, but it requires that there be no duplicates.

There is a test program associated with this class, 
testing/CompMasterVecTest.cc. Please run it if you change anything, and add
tests to it as you add functionality.


More detailed notes are in the .h file.

  \todo: design a readonly compmastervec, for use when 
  reading from files mainly, which does not contain a hashtable and does not
  Insert().


*/


/* Detailed notes on class design:
-Here's where the class comes from at the algorithmic level. Assume we have 
N objects coming in (e.g. 10^9) of which K (e.g. 10^8) are unique. What's the 
best way of storing them?

There are two problems here:

First, we need some form of unique storage for the K unique objects 
 (collection U), and we also need a structure of size N (collection I) that 
will hold indices or pointers to the appropriate storage location inside U.

Second, we want to avoid doing holding K separate mallocs in memory, if 
possible. This could be done if we have some idea of the value of K and the
size of each unique item, by mallocing a large area M of (K * individual 
size) and storing them in there sequentially, like mastervec does.

I was trying to solve both problems at the same time, but in fact they can be 
separated by having structure U hold not the objects themselves, but rather 
indices into the large area M. This assumes that structure U does not require 
a separate malloc for each item in the structure. Unfortunately, this is not 
true of trees, and is also not true of most hash tables. It is true of a hash 
table with linear probing or double hashing, however!

Ok, now some numbers. 

Tree: we have N insertions into a tree of final size K, so the total number 
of insertions is O(NlogK), and retrieval takes O(logK) time. We need O(K)
separate mallocs.

Hash table: we have N insertions in constant time, so insertion is O(N). If 
we do not know K a priori, we have O(logK) rehashings needed, each of which 
involves:
-the actual rehashing, which is O(K).
-reindexing the vector of indices, which is O(N).
Thus, the rehashings will take O((N+K)logK) = O(NlogK). However, the constant
factor for this will be much lower than the constant factor for the insertions
into the tree.

it is also possible that we would not need to touch the O(N) vector of indices
if we have an intermediate vector of size K which is the one that gets
reordered, in which case the rehashings will only take O(KlogK) instead of
O(NlogK). That leaves us with a hashtable method of O(N + KlogK) which is almost
certainly better than the tree method if K is very big. Of course, if K is
small, then the tree method is much simpler to implement and faster.

Another idea from Michael: there are really two states of this object, writing
and reading. and it turns out that if you know you are in a read-only state, you
can get rid of the intermediate structure (U) and make the indices in I index
directly into M. But having two-state objects is a problem, so the right way to
do it is to have two different objects: compmastervec and readonlycompmastervec.
The file format should be such that readonlycompmastervec can read it and create
itself as a readonly, and that compmastervec can also recreate its complete
internals easily (obviously, we want K in there, so compmastervec will know how
much space to allocate for the U structure!).
The key difference between readonly and compmastervec is that readonly does not
support any form of push_back(), resize(), reserve(), any constructor other than
from another compmastervec or a file, clear(), destroy(), SwapElements(), or
X& operator[], but rather only const X & operator[];

Note that the transformation from a compmastervec to a readonlycompmastervec
should take only O(N) whether we are using a hashtable or a tree, because we do
not have to search, we only have to index through to get the index into the
structure M.

having said all this, here is my current plan:
I is a vec of all the indices. We can calculate N by counting and preallocate if we want.
U is a hashtable with linear probing, that resizes itself if it gets more than
    2/3 full, and resizes x 4 for starters. U is actually a little more complex
    than that: it is first a vec of size K and then a hashtable, with the vec
    containing the index into the table and I containing the index into the vec.
    This way, when we rehash, we can extend the vec, adjust its existing
    indices, and not need to touch I.
M is a mastervec. If M gets full, we resize it with the same factor as U.

for readonlycompmastervec, I can contain either an index or directly the offset
of the data. The second would be more efficient but would require more rewriting.

  So how do we do this?
  -Keep two vectors: one of them contains the actual dynamic data, while the
  other is a vector of indices to the data. For example, if the input data has
  data[13333] = "string1", data[26666] = "string1", data[2222] = "string2", then
  we have:
  m_mvec[0] = "string1"
  m_mvec[1] = "string2"
  indices[13333] = 0;
  indices[26666] = 0;
  indices[2222] = 1;

  in this way the total storage is Nx4 + (unique_data), where N is the number of 
  data items. This is likely to be much shorter than 
  N*(each dynamic data length).

  -The dynamic data are kept in a mastervec. The mastervec is indexed into
  by a hashtable, which is how we handle the searches needed to check for 
  repetitions when we put in new data with push_back() or non-const operator[]

  -We need to have smart iterator and reference classes (see below).

  -problem with erasing or replacing data. if we want to replace the data
  at index i, in theory we should check to see if any other element in the 
  indices vector is interested in the same data, so we can remove it from the 
  mastervec. But you cannot remove data from mastervecs, you can only 
  push_back. This means that if you overwrite something that was unique, 
  you now have "dead" data in your mastervec, and there is no way to clean 
  it up. 

*/

#include "feudal/MasterVec.h"
#include "Intvector.h"
#include "Vec.h"
#include "String.h"
#include "system/Assert.h"
#include "HashLinearRemoteStorage.h"
#include <cstddef>

#include <ext/hash_set>


/** CMVHash: Default hash class for compmastervec.
 *
 * @class CMVHash
 * \author Pablo Alvarez
 *
 * It will treat all Xs as vecs of char, glom each group of 4 chars
 * into an int, xor all the ints together, call gnu's hash<int>
 * method on that, and return that result % M.
 * There is a specialization for String that uses the djb2 hash.
 */
template<class X>
struct CMVHash: public binary_function<const X &, int, int> { 
  ///Returns an int in the interval [0,M)
  unsigned int operator()(const X & x, int M) const {
    return operator()(x) % M;
  }
 
  ///Returns an int from 0 to UINT_MAX.
  /// xors together all the ints and then calls hash<int>() on the result.
  /// \todo use extra bytes for things longer than an int (e.g. with 7 bytes).
  unsigned int operator()(const X & x) const {
      int hash = 5381;
      char const* end = reinterpret_cast<char const*>(&*x.end());
      char const* begin = reinterpret_cast<char const*>(&*x.begin());
      while ( begin != end ) {
        hash = ((hash << 5) + hash) + *begin++;
      }
      return hash;
  }
};

template<class X, class H = CMVHash<X> >
class compmastervec {
  
public:

  /** Reference: reference class for compmastervec operator[] non-const.
 
  @class Reference

  \author Pablo Alvarez.


  This class is needed so we can use operator[] as an l-value with 
  compmastervec. It is an inner class of compmastervec.

  Operator= can be overloaded to do the following:
  -figure out if we are changing the value.
  -if so, check whether X is already in the mastervec, and add it if 
  necessary, and set the appropriate value in the indices vector.

  The class does not need to have an X & as a member, just the index, and 
  can  use the index to access the value as needed. That way the Xs are 
  always contained within the compmastervec and we cannot mess things up 
  from here. 
  */ 
  class Reference {
   public:
    Reference(compmastervec<X,H> * o, size_t i):
      owner(o), index(i) {
      AssertLt(index, owner->size());
    }

    ///Get the item we refer to.
    X& theX() const {
      return owner->m_mvec[owner->m_indices[index]];
    }

    ///Automatic conversion to const X&
    operator X& () const { 
      return theX(); 
    }

    bool operator==(const X& rhs) const { return theX() == rhs; }
    bool operator!=(const X& rhs) const { return theX() != rhs; }
    bool operator<=(const X& rhs) const { return theX() <= rhs; }
    bool operator>=(const X& rhs) const { return theX() >= rhs; }
    bool operator<(const X& rhs) const { return theX() < rhs; }
    bool operator>(const X& rhs) const { return theX() > rhs; }

    ///Change a value inside the compmastervec.
    const X& operator=(const X & x) 
    {
        if ((owner->m_indices[index] == HashLinearRemoteStorage<X,H>::EMPTY) 
            || (x != theX()))
        {
            size_t i = owner->m_hashtable.Insert(x);
            owner->m_indices[index] = i;
        }
        return x;
    }

    ///Change a value inside the compmastervec.
    const X& operator=(const Reference & rhs) {
      return &rhs == this ? rhs : operator=((const X &) rhs);
    }

    ///Required to fix problem in gcc 4.1.0
    ///because it cannot find the operator const X & () when sending to
    ///an ostream
    friend ostream & operator<<(ostream & os, const Reference & r) {
      os << r.theX();
      return os;
    }

    /*
      And I think those are really the only two operations that we need. 
      Autogenerated destructor and copy constructor are fine. 
    */
   private:
    compmastervec<X,H> * owner;
    size_t index; ///< the index of our X object within the compmastervec.
  };


  /** Iterator: a random-access iterator for compmastervec
      @class Iterator
      \author Pablo Alvarez
    
      Ok, now let's think about iterators. They need to contain the current 
      index, increment it and decrement it as needed, and return either a 
      const X & or a Reference when dereferenced. That's not too hard...

      This class keeps track of the owning compmastervec and of its
      current index, and uses them to provide all the needed iterator
      facilities. Note that to get a const_iterator, all you need to do
      is use a const Iterator. The only function affected by its being const
      is operator *().

      To do: figure out if we need to assert that the owner is not empty.
      It sort of seems superfluous because otherwise we would not have been
      able to create the Iterator in the first place, but on the other hand
      the size can be changed out from under us...
      One possible principle is this: we make sure everything is OK whenever 
      we dereference or when we create a new iterator or modify this one.
  */
  class Iterator {
   public:
    Iterator(compmastervec<X,H> * o, size_t i): owner(o), index(i) {
      //Need to allow for index = size, because that is end().
      AssertLe(index, owner->size());
    }

    bool operator==(const Iterator & rhs) { return index == rhs.index; }
    bool operator!=(const Iterator & rhs) { return index != rhs.index; }
    bool operator<(const Iterator & rhs) { return index < rhs.index; }


    ///Addition operator, based on operator+=
    ///See More Effective C++, p. 108 for reasons.
    Iterator operator+(size_t i) {
      return (Iterator(*this)+= i);
    }

    ///Subtraction operator, based on operator-=
    ///See More Effective C++, p. 108 for reasons.
    Iterator operator-(size_t i) {
      return (Iterator(*this)-= i);
    }

    ///Assignment addition operator, asserts that the return value is valid. 
   Iterator & operator+=(size_t i) {
      index += i;
      AssertLe(index, owner->size());
      AssertGt(owner->size(), 0u);
      return *this;
    }
    ///Assignment subtraction operator, asserts that the return value is valid.
    Iterator & operator-=(size_t i) {
      index -= i;
      AssertLe(index, owner->size());
      AssertGt(owner->size(), 0u);
      return *this;
    }

    ///Pre-increment, asserts that the result iterator is valid.
    Iterator& operator++() { 
      AssertLe(index, owner->size()); 
      ++index;
      return *this;
    } 

    ///Post-increment, asserts that the result iterator is valid.
    Iterator operator++(int) { 
      AssertLe(index, owner->size()); 
      Iterator temp = *this;
      ++index; 
      return temp;
    } 
    ///Pre-decrement, asserts that the result iterator is valid.
    Iterator& operator--() { 
      AssertGt(index, 0u);
      --index;
      return *this;
    } 
    ///Post-decrement, asserts that the result iterator is valid.
    Iterator operator--(int) { 
      AssertGt(index, 0u);
      Iterator temp = *this;
      --index; 
      return temp;
    } 

    ///Return a Reference that can be used as an l-value
    Reference operator*() { return Reference(owner, index); }

    ///Return a const X & that can only be an r-value.
    const X & operator*() const {
      return owner->m_mvec[owner->m_indices[index]];
    }

    ///Subtract two Iterators to obtain their distance.
    friend 
    std::ptrdiff_t operator-(const Iterator & lhs, const Iterator & rhs) {
      return lhs.index - rhs.index; 
    }

    ///Specialization of std::distance()
    friend 
    std::ptrdiff_t distance(const Iterator & first,
                       const Iterator & last) {
      return first.index - last.index; 
    }

   private:
    compmastervec<X,H> * owner;
    size_t index; ///< the index of our X object within the compmastervec.
  };
  
 public:

  typedef Iterator                      iterator;
  typedef Reference                     reference;
  typedef size_t                        size_type;
  typedef std::ptrdiff_t                difference_type;

  static const size_t DYNAMIC_DATA_SIZE;///<Default value
  static const size_t OBJECTS; ///<Default value
  static const size_t UNIQUE_OBJECTS; ///<Default value

  compmastervec( );

  ///Constructor reserves space for total number of objects.
  compmastervec( size_t total_objects );

  ///Constructor from a datafile.
  compmastervec( const String &filename ); 

  ///Constructor from a mastervec.
  compmastervec( const MasterVec<X> & mvec);

  ~compmastervec();
     
  iterator       begin()       { return iterator(this, 0); }
  iterator       end()         { return iterator(this, size()); }

  reference       front()       { return reference(*begin()); }
  reference       back()        { return reference(*(end()-1)); }


  /** Reserve memory in the different components
   * @param raw_mem_size: total size of dynamic data of unique objects, 
   measured in units of sizeof(A), not char
   * If any of the parameters are less than their respective current values,
   * no change takes place.
   */
  void Reserve( size_t raw_mem_size, size_t unique_objects, 
                size_t total_objects );

  /// Change the total number of objects
  void resize( size_t total_objects );
     
  /// Add capacity for the total number of objects
  void reserve( size_t total_objects );
     
  ///Clear out all objects, but do not reduce memory usage
  void clear( );
     
  ///Clear out all objects and deallocate all memory.
  void destroy( );

  ///True if object is in the compmastervec.
  bool Contains(const X & obj) { 
    return m_hashtable.Find(obj) != HashLinearRemoteStorage<X,H>::EMPTY; 
  }

  ///Return the unique index for this object in the mastervec.
  size_t UniqueIndex(const X & obj) {
    return m_hashtable.Find(obj);
  }

  ///Return the unique index for the object pointed to by this index.
  size_t UniqueIndex(size_t index) { return m_indices[index]; }

  ///Find first index for an object. Inefficient, method is O(total_objects)!
  /// returns -1 if not found.
  size_t Find (const X & obj) { 
    size_t i = m_hashtable.Find(obj);
    if (i == HashLinearRemoteStorage<X,H>::EMPTY) return -1;
    return (find(m_indices.begin(), m_indices.end(), i) - m_indices.begin());
  }

  /**Only works if we use this class for hashing, not for compression.
     To be precise, every object inserted must be unique, and there can have
been no swaps.

     If the object is not found, if there are any duplicates, or if there
have been any swaps that affect this object, returns -1.

     So in most cases, you need to create a reverse index, and then
     use UniqueIndex and look the result up in the reverse index.
  */
  size_t FindUnique (const X & obj) { 
    if (m_indices.size() != m_mvec.size()) {
      return -1;//there are duplicates
    }
    size_t i = m_hashtable.Find(obj);
    if (i == HashLinearRemoteStorage<X,H>::EMPTY) return -1;
    if (m_indices[i] == i) return i;
    else return -1;
  }

  /**Create a mastervec that is a reverse index for this compmastervec.
     \param rindex: MasterVec of size UniqueSize(). Each of its component
     SerfVecs of index i contains all the indices in m_indices that correspond
     to the same underlying object whose UniqueIndex() is i.
     This reverse index can be useful in a variety of situations, for 
     example when finding reads that have the same template_id.
   */
  void CreateReverseIndex(VecIntVec& rindex);

  /**Find a permutation that will make *this match to.
     Preconditions:
     - all elements of *this must be unique.
     - class V must provide size() and operator[].
     - to.size() < this->size()
     Any items in *this that are not found in to will be given an index of
     -1 in the permutation vector.
     The indices for any elements in to that are not found in *this will be 
     placed in the notFound vector
  */
  template<class V>
  void FindPermutationTo(const V & to, vec<int> & perm, vec<int> & notFound);

  /** Add a new item to the compmastervec.
   * @param extra_space space for the object X to expand measuredin sizeof(A).
   * We store the item in the hashtable, and store the index returned by 
   * the hashtable into m_indices.
   */
  void push_back( const X& obj, size_t extra_space = 0 ) {
    size_t i = m_hashtable.Insert( obj, extra_space);
    m_indices.push_back(i);
  }

  ///Add a new item that owns its own dynamic memory.
  void push_back_external( const X& obj ) 
  {
    size_t i = m_hashtable.Insert( obj, 0, true);
    m_indices.push_back(i);
  }

  /**Append elements in range [from,to) from orig to this compmastervec.
   * Precondition: 0 <= from <= to <= orig.size().
   * @param extraRaw: additional dynamic memory to be Reserved for future use.
   * @param extraSize: additional space for unique objects to be Reserved for future use.
   * Since append is likely to cause a memory reallocation, the two default
   * parameters allow the programmer to make that reallocation larger if they
   * know they are going to need it later.
   */
  void Append( const compmastervec & orig, size_t from, size_t to);

  /**Append all elements indexed by entries to this mastervec.
   * Precondition: all indices in entries are valid indices for orig.
   * See notes for append(orig, from, to).
   */
  void Append( const compmastervec & orig, const vec<int> & entries);

  /**Append elements in range [from,to) from orig to this compmastervec.
   * Precondition: 0 <= from <= to <= orig.size().
   * @param extraRaw: additional dynamic memory to be Reserved for future use.
   * @param extraSize: additional space for unique objects to be Reserved for future use.
   * Since append is likely to cause a memory reallocation, the two default
   * parameters allow the programmer to make that reallocation larger if they
   * know they are going to need it later.
   */
  void Append( const MasterVec<X> & orig, size_t from, size_t to);

  /**Append all elements indexed by entries to this mastervec.
   * Precondition: all indices in entries are valid indices for orig.
   * See notes for append(orig, from, to).
   */
  void Append( const MasterVec<X> & orig, const vec<int> & entries);

  ///Delegate to internal mastervec
  size_t sumSizes( ) const { return m_mvec.sumSizes(); }

  ///Number of unique items.
  size_t UniqueSize() const { return m_mvec.size(); }

  ///Capacity for unique items.
  size_t UniqueCapacity() const { return m_mvec.capacity(); }

  ///Allow direct reading from our mastervec.
  const MasterVec<X> & UniqueVec() const { return m_mvec; }
  
  ///Size of m_indices array
  size_t      size( )        const { return m_indices.size(); }

  bool empty() const { return m_indices.empty(); }
  

  ///Capacity of m_indices array
  size_t      capacity( )    const { return m_indices.capacity(); }

  ///Delegate to internal mastervec
  //const A* rawdata( )     const { return m_mvec.rawdata(); }

  ///Swap the elements pointed to by i and j, without actually moving the data.
  void SwapElements( size_t i, size_t j ) 
  {
    if ( i == j || m_indices[i] == m_indices[j] ) return;
    size_t temp = m_indices[i];
    m_indices[i] = m_indices[j];
    m_indices[j] = temp;
  }

  /**This indexing operator can be used to get an l-value.
   * Thus, you can change the mastervec's contents. The Reference class is 
   * smart about that, and has an automatic cast to const X &. If the item
   * put into the l-value is new, it is added to the list of unique items,
   * and the old item is not removed (it may be pointed to by other m_indices,
   * and besides mastervec does not support erase() ).
   */
  Reference operator[] (size_t i) 
  {
    AssertLt( i, size() );
    return Reference(this, i);   
  }

  const X& operator[] (size_t i) const
  {
    AssertLt( i, size() );
    return m_mvec[m_indices[i]];   
  }

  // The Sort() and Unique() methods are not implemented for compmastervecs
  // because they do not match the way the class is meant to be used.

  /** Write complete compmastervec to a file. Efficient.
   */
  void WriteAll( const String& filename ) const;

  /**Write all elements between from and to to a file.
   * \todo: write it so it is economical and does not use 
   * Write(String, vec).
   */
  void Write( const String& filename,
	      size_t from,
	      size_t to ) const;
  
  /** Write all elements with m_indices in entries to a file
   * This function first creates a new compmastervec with all the entries
   * using the partial copy constructor and then writes that to file.
   */   
  void Write( const String& filename,
              const vec<int> & entries) const;
  
  ///Not implemented yet   
  void Read( const String& filename,
	     const vec<int>& entries,
	     size_t extra = 0,
	     bool pre_reserved = false,
	     bool pre_resized = false );

  ///Not implemented yet
  void ReadRange( const String& filename,
		  size_t from,
		  size_t to,
		  size_t extra = 0,
		  bool pre_reserved = false,
		  bool pre_resized = false );

  ///Not implemented yet
  void ReadOne( const String& filename,
		size_t it,
		size_t extra = 0,
		bool pre_reserved = false,
		bool pre_resized = false )
  {
    // ReadRange( filename, it, it + 1, extra, pre_reserved, pre_resized );   
  }
     
  ///Not implemented yet
  void SparseRead( const String& filename,
		   const vec<int>& entries,
		   size_t extra,
		   bool pre_reserved = false );

  ///Not implemented yet
  void SparseReadRange( const String& filename,
			size_t from,
			size_t to,
			size_t extra,
			bool pre_reserved = false );

  /**Read a complete compmastervec from a file.
   * This is efficient.
   */    
  void ReadAll( const String& filename, bool append = false );

  /**Assignment operator makes a deep copy.
   * It copies both m_mvec and m_indices vectors.
   * Also recreates the hash table with Rehash().
   */
  compmastervec& operator= ( const compmastervec& original ) {
    if (this == &original) return *this;
    m_mvec = original.m_mvec;
    m_indices = original.m_indices;
    m_hashtable.Rehash();
    return *this;
  }

  /** Returns true if all m_indices and elements are equal.
   * \todo: change it so that if the insertion/hashing order was different
   * but all m_indices still point to the same elements it also returns true.
   */
  bool operator==(const compmastervec & lhs) {
    if (m_indices != lhs.m_indices) return false;
    if (m_mvec.size() != lhs.m_mvec.size()) return false;
    for (size_t i = 0; i != m_mvec.size(); ++i) {
      //Use operator== because that is the only one guaranteed for type X.
      if (!(m_mvec[i] == lhs.m_mvec[i])) return false;
    }
    return true;
  }

  bool operator!= (const compmastervec & lhs) {
    return !(operator==(lhs));
  }
    
  friend
  std::ostream & operator<<(std::ostream & os, const compmastervec & cmv) {
    os << "Unique Objects: \n" << cmv.m_mvec.size() << endl;
    os << "Indices, size: " << cmv.m_indices << endl;
    return os;
  }

 private:

  ///Set the isCompFile bit in the file header.
  void SetFileCompBit(const String & filename, bool on = true) const ;

  /// The copy constructor is declared private to prevent accidental passing of 
  /// compmastervec arguments or return values by value.  

  compmastervec( const compmastervec& original );

 protected:
  MasterVec<X> m_mvec;
  HashLinearRemoteStorage<X,H> m_hashtable;
  vec<size_t> m_indices;
};
extern template class compmastervec<String>;

///Total size of unique dynamic data in the file
size_t CompmastervecFileRawCount( const String& filename );

///Total number of unique objects in the file
size_t CompmastervecFileUniqueObjectCount( const String& filename );

///Calls m.destroy()
template<class X, class H> 
inline void Destroy( compmastervec<X, H>& m )
{
  m.destroy( );    
}

template<class X, class H>
template <class V>
void compmastervec<X,H>::FindPermutationTo(const V & to, vec<int> & perm,
                                             vec<int> & notFound) {
  const size_t N = size();
  AssertLe(to.size(), N);
  perm.assign(N, -1);
  notFound.clear();
  size_t K = to.size();
  for (size_t i = 0; i != K; ++i) {
    long pos = FindUnique(to[i]);
    if (-1 == pos) {
      notFound.push_back(i);
    } else {
      perm[pos] = i;
    }
  }
}

///Do a WEAK test to see if a file is in the compmastervec format.
bool IsGoodCompmastervecFile( const String& filename );

#endif //COMPMASTERVEC_H
