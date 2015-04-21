#ifndef __SPARSE_ARRAY_H
#define __SPARSE_ARRAY_H

#include "system/Assert.h"
#include <vector>
#include <map>


#ifdef __GNUC__
#include <ext/hash_map>
using __gnu_cxx::hash_map;
#endif 

/// Wrapper. Our SparseArray class template takes a template template
/// argument container<I,V> that specifies the container type to 
/// store the data internally; this associative container must have
/// two template arguments - key and value (cf. map). Thus we use this
/// wrapper to explicitly define vector of (index, value) pairs as
/// a special template with two template parameters, so now it can be passed
/// to SparseArray as backend storage.
template <class I, class V> class PairVector : public vector< pair<I,V> > {};

// forward declaration
template <class T, class R> class SparseArray;

/// Auxiliary class used to return indexed values from sparse array
/// (as in my_sparse_array[i]). This "reference class" trick is needed
/// because the storage for most of the sparse array elelments is not
/// allocated (this is what makes the array sparse and helps save
/// huge amounts of memory). At the moment when the array element is 
/// accessed, there is no way of knowing how it is going to be used:
/// one may simply read the value from the array (e.g.
/// v=my_sparse_array[i]), or one may want to modify that value
/// (e.g. my_sparse_array[i]=v). In the latter case, the storage might
/// need to be allocated if the new value is not the default value
/// assigned to the majority of (unallocated) elements of the array.
/// Using this auxliary lvalue class we can defer resolution of 
/// this allocation issue to the point when the context is known to
/// the compiler. For instance, in both cases of 1) assignment 
/// (my_sparse_array[i]=v) or 2) reading the value (v=my_sparse_array[i]),
/// sparse array's indexing operator
/// returns lvalue wrapper class (no allocation performed yet), then
/// lvalue's assignment operator (case 1, allocation will be performed)
/// or type conversion operator (case 2, read the value off, no allocation)
/// will be invoked. NOTE: this trick allows using indexing operators
/// on sparse arrays trasparently, which helps with implementing some
/// algorithms generically, but this also incurs some performance penalty.
/// See SparseArray class and its Set and Get methods for faster access.

template <class T, class R> class _SparseArray_lvalue {

	/// private constructor: only friend class (SparseArray) 
    /// can instantiate this wrapper
	_SparseArray_lvalue(const SparseArray<T,R> * a) :
	  m_array( *( const_cast<SparseArray<T,R> *>(a) ) ) {}
	/// private assignment: lvalue wrapper class instances 
	/// can not be assigned by clients

        const _SparseArray_lvalue & operator=(const _SparseArray_lvalue &) { return *this; }

	unsigned int m_ind;
	SparseArray< T, R > & m_array;
public:  
	void SetIndex(unsigned int i) { m_ind = i; }

	operator T () const { return m_array.Get(m_ind); }

	const T & operator=(const T & val) {
		return m_array.Set(m_ind, val);
	}

	const T operator++() {
	    typename SparseArray<T,R>::t_iter it = m_array.FindIndex(m_ind); 
	    T val;
	    if ( it == m_array.StorageEnd() ) {
	      val = m_array.GetDefault();
	      ++val;
	      m_array.PushBack(m_ind, val);
	    } else {
	      val = (++(it->second));
	      if ( val == m_array.GetDefault() ) m_array.Erase(it);
	    }
	    return val;
	}


	const T operator++(int unused) {
	    typename SparseArray<T,R>::t_iter it = m_array.FindIndex(m_ind); 
	    T ret;
	    if ( it == m_array.StorageEnd() ) {
	      T val = ret = m_array.GetDefault();
	      ++val;
	      m_array.PushBack(m_ind, val);
	    } else {
	      ret = it->second;
	      ++(it->second);
	      if ( it->second == m_array.GetDefault() ) m_array.Erase(it);
	    }
	    return ret;
	}

	T operator+=(T add_val ) {
	    typename SparseArray<T,R>::t_iter it = m_array.FindIndex(m_ind); 
	    T val;
	    if ( it == m_array.StorageEnd() ) {
	      val = m_array.GetDefault();
	      val += add_val;
	      if ( val != m_array.GetDefault() ) m_array.PushBack(m_ind, val);
	    } else {
	      it->second += add_val;
	      val = it->second;
	      if ( val == m_array.GetDefault() ) m_array.Erase(it);
	    }
	    return val;
	}

	friend class SparseArray<T,R>;
};


/// This is an auxiliary class not intended for direct use by the client 
/// code. This class provides basic storage functionality for (index, value)
/// pairs backed up by a vector-like container (see also specializations for 
/// map-like containers).
template <class T, class CONTAINER> 
class IndexedStoragePolicy {
 private:
  CONTAINER m_elements;
  bool m_sorted;
 public:
  typedef typename CONTAINER::iterator t_iter;
  typedef typename CONTAINER::const_iterator t_const_iter;


  IndexedStoragePolicy() : m_sorted(true) {}

  bool IsSorted() const { return m_sorted; }

/// Looks up the specified index in the sparse array.
/// If this index is already used (i.e. non-default value is set
/// at that position in the sparse array), the returned iterator into the
/// container of \em indices points to the element where the index \c i
/// is stored. Otherwise the returned iterator points immediately
/// after the end of the internal container of indices.
/// Note: this method checks whether the sparse array's internal
/// representation is sorted or not and acts accordingly - finding
/// indices in sorted instances is faster (log n), finding indices
/// in unsorted instances is linear (n), where n is the number of
/// currently held non-default values in the sparse array (\em not
/// the total size of the original, unpacked array it represents!)

  t_const_iter FindIndex(unsigned int i)  const {
 
    t_const_iter it_end = m_elements.end();

    if ( m_sorted ) {
        // do binary search
      t_const_iter it_begin = m_elements.begin();
      while ( it_begin != it_end ) {
	t_const_iter it_mid = it_begin;
	advance( it_mid, distance(it_begin, it_end)/2 );
	unsigned int index = it_mid -> first;
	if ( index == i ) return it_mid; // index found!
	if ( index > i ) {
	  it_end = it_mid;
	} else {
	  it_begin = it_mid;
	  it_begin++;
	}
      }
      return m_elements.end();
    } else {	
      for ( t_const_iter it = m_elements.begin() ; it != it_end ; it++ ) {
	if ( (it->first) == i ) return it;
      }
      return it_end;
    }

  }



  t_iter FindIndex(unsigned int i)  {
 
    t_iter it_end = m_elements.end();

    if ( m_sorted ) {
        // do binary search
      t_iter it_begin = m_elements.begin();
      while ( it_begin != it_end ) {
	t_iter it_mid = it_begin;
	advance( it_mid, distance(it_begin, it_end)/2 );
	unsigned int index = it_mid -> first;
	if ( index == i ) return it_mid; // index found!
	if ( index > i ) {
	  it_end = it_mid;
	} else {
	  it_begin = it_mid;
	  it_begin++;
	}
      }
      return m_elements.end();
    } else {	
      for ( t_iter it = m_elements.begin() ; it != it_end ; it++ ) {
	if ( (it->first) == i ) return it;
      }
      return it_end;
    }
  }

  void PushBack(unsigned int i, const T & val ) {
    // if the new index being added is not in order,
    // set "sorted" flag to false
    if ( m_sorted && 
	 ! m_elements.empty() && 
	   i < m_elements.back().first ) m_sorted = false;

    m_elements.push_back(pair<unsigned int, T>(i, val));
  }

  t_const_iter StorageBegin() const { return m_elements.begin(); }
  t_const_iter StorageEnd() const { return m_elements.end() ; }
  t_iter StorageBegin() { return m_elements.begin(); }
  t_iter StorageEnd() { return m_elements.end() ; }
  void Erase(t_iter pos) { 
    m_elements.erase(pos); 
  }
  unsigned int StorageSize() { return m_elements.size(); }
  void ClearStorage() { m_elements.clear(); m_sorted = true; }

  /// Reserves the storage for \c n \em non-trivial elements
  /// (i.e. the elements with non-default values that will be
  /// physically stored in the sparse array representation).
  /// This is \em not the total size of the sparse array!
  void ReserveStorage(unsigned int n) { m_elements.reserve(n); }

  void RollbackIterator ( t_iter & i ) const  {--i; }
  void RollbackIterator ( t_const_iter & i ) const { --i; }
};

/// Specialization for (index, value) pair storage backed by a map
template <class T> 
class IndexedStoragePolicy<T, map<unsigned int, T> > {
 private:
  map<unsigned int, T> m_elements;
 public:
  bool IsSorted() const { return true; }
  typedef typename map<unsigned int, T >::iterator t_iter;
  typedef typename map<unsigned int, T >::const_iterator t_const_iter;

  t_iter FindIndex(unsigned int i) { return m_elements.find(i); }
  t_const_iter FindIndex(unsigned int i) const { return m_elements.find(i); }
  void PushBack(unsigned int i, const T & val ) {
    m_elements.insert(pair<unsigned int , T>(i, val));
  }
  t_const_iter StorageBegin() const { return m_elements.begin(); }
  t_const_iter StorageEnd() const { return m_elements.end() ; }
  t_iter StorageBegin() { return m_elements.begin(); }
  t_iter StorageEnd()  { return m_elements.end() ; }
  void Erase(t_iter pos) { 
    m_elements.erase(pos); 
  }
  unsigned int StorageSize() { return m_elements.size(); }
  void ClearStorage() { m_elements.clear(); }

  /// Required by the interface, does nothing. There is no way to reserve
  /// storage for the map
  void ReserveStorage(unsigned int n) {  }
  void RollbackIterator ( t_iter & i ) const  {--i; }
  void RollbackIterator ( t_const_iter & i ) const { --i; }
};


#ifdef __GNUC__

/// Specialization for (index, value) pair storage backed by a hash_map
template <class T> 
class IndexedStoragePolicy<T,hash_map<unsigned int,T> > {
 private:
  hash_map<unsigned int, T> m_elements;
 public:
  bool IsSorted() const { return false; }
  typedef typename __gnu_cxx::hash_map<unsigned int, T >::iterator t_iter;
  typedef typename __gnu_cxx::hash_map<unsigned int, T >::const_iterator t_const_iter;

  t_iter FindIndex(unsigned int i) { return m_elements.find(i); }
  t_const_iter FindIndex(unsigned int i) const { return m_elements.find(i); }
  void PushBack(unsigned int i, const T & val ) {
    m_elements.insert(pair<unsigned int , T>(i, val));
  }
  t_const_iter StorageBegin() const { return m_elements.begin(); }
  t_const_iter StorageEnd() const { return m_elements.end() ; }
  t_iter StorageBegin() { return m_elements.begin(); }
  t_iter StorageEnd()  { return m_elements.end() ; }
  void Erase(t_iter pos) { 
    m_elements.erase(pos); 
  }
  unsigned int StorageSize() { return m_elements.size(); }
  void ClearStorage() { m_elements.clear(); }


  void ReserveStorage(unsigned int n) { m_elements.resize(n); }
  void RollbackIterator ( t_iter & i ) const  { }
  void RollbackIterator ( t_const_iter & i ) const { }
};

#endif

template <class T, class CONTAINER=PairVector<unsigned int,T> >
  class SparseArray : protected IndexedStoragePolicy<T, CONTAINER> {
private:
  using IndexedStoragePolicy<T, CONTAINER>::StorageBegin;  
  using IndexedStoragePolicy<T,CONTAINER>::StorageEnd;  
  using IndexedStoragePolicy<T,CONTAINER>::ClearStorage;  
  using IndexedStoragePolicy<T,CONTAINER>::FindIndex;  
  using IndexedStoragePolicy<T,CONTAINER>::Erase;  

public:
  using IndexedStoragePolicy<T,CONTAINER>::IsSorted;  
  using IndexedStoragePolicy<T,CONTAINER>::StorageSize;  

  
// have to declare as a friend to be able to instantiate!
// lvalue is a utility class completely hidden from applications!
friend class _SparseArray_lvalue<T, CONTAINER>;

SparseArray() : m_lval(this),
	        m_default(0),
		m_size(0) {};

/// copy constructor 
SparseArray(const SparseArray & s) : 
    IndexedStoragePolicy<T, CONTAINER>(s), /// default copy ctor for the storage
    m_lval(this), ///< CRITICAL! if copy constructor were not provided, this 
                  /// would be copied from s!!!
    m_default(s.m_default),
    m_size(s.m_size) { }


explicit SparseArray(unsigned int size, const T & def) : 
  m_lval(this),
  m_default(def),
  m_size(size) {}

/// Returns the size of the vector; \em not the size of the
/// non-trivial representation (i.e. the number of actually
/// stored non-default elements
unsigned int size() const { return m_size; }

void resize(unsigned int new_size, T val) {
  if ( new_size >= m_size ) {
    if ( val != GetDefault() ) {
      // inserting non default values:
      for ( unsigned int i = m_size; i < new_size ; i++ ) {
	this->PushBack(i, val);
      }
    } // if val is actually equal to default, there was nothing
      // to physically insert!
    m_size = new_size; // reset the size
  } else { 
    // if we are shrinking the array, we may need to erase some elements
    // (those with indices i >= new_size)
    bool done = false;
    while ( ! done ) {
      t_iter it = StorageBegin();
      for ( ; it != StorageEnd() ; it++ ) {
	if ( it->first >= new_size ) {
	  Erase(it);
	  break;
	}
      }
      if ( it == StorageEnd() ) done = true;
    }
  }
}
  

void resize(unsigned int new_size) { 
  if ( new_size >= m_size ) {
    m_size = new_size;
    return;
  } else {
    bool done = false;
    while ( ! done ) {
      t_iter it = StorageBegin();
      for ( ; it != StorageEnd() ; it++ ) {
	if ( it->first >= new_size ) {
	  Erase(it);
	  break;
	}
      }
      if ( it == StorageEnd() ) done = true;
    }
  }
}

using IndexedStoragePolicy<T, CONTAINER>::ReserveStorage;

/// Resets all elements of the sparse array to the default value.
/// Does \em not change the size of the vector!
void ResetToDefault(void) { ClearStorage(); }

/// Returns default value (the value, for which this sparse
/// array instance does not allocated physical storage).
T GetDefault(void) const { return m_default; }


/// Returns the amount of memory (in bytes) actually used 
/// by the instance of SparseArray. The memory used does include
/// dynamically allocated \em and used memory (such as vector::size()), 
/// but does not include preallocated and yet unused
/// memory (such as vector::capacity())
  unsigned int MemoryUsage() {
    return (sizeof(SparseArray<T,CONTAINER>) + 
	    sizeof(pair<unsigned int, T>)*StorageSize());
  }

  /// Sets the \em i-th element of the array to \c val.
  /// This method checks the index at which the value is being added
  /// and monitors the \c sorted status of the sparse array's internal
  /// representation. If a number of elements are set in the order of
  /// increasing indices, the most efficient internal representation 
  /// will be kept and subsequent setting/retrieval operations will be
  /// much faster. This method is checked (breaks execution if the index
  /// is out of range).
  const T & Set(unsigned int i, const T & val ) {
    //    ForceAssertLt(i, m_size);
    t_iter it = FindIndex(i);
    if ( it != StorageEnd() ) {
      // element with this index is already explicitly
      // stored by the sparse array 
		if ( val == m_default ) {
			// we are setting a default value! the explicit
			// storage for element i must be destroyed!
		  	Erase(it);
		} else {
			// we are setting non-default value: need to update;
			(it->second) = val;
		}
    } else {
        // the storage for element i did not exist. we need to create it
        // only if we are setting a non-default value! otherwise there is
        // nothing to do
        if ( val == m_default ) return val;


	// Use the storage policy's PushBack - it might need
	// to take care of fine details, such as tracking
	// the "sorted" status (if we are backed by an array rather than a map)

	this->PushBack(i, val);
    }
    return val;
  }

  /// Returns the value stored in the array at position \c i.
  /// This method is checked and breaks execution if index is out
  /// of range.
  T Get(unsigned int i) const {
    ForceAssertLt(i, m_size);
    t_const_iter it = FindIndex(i);
    // if the specified index is not stored explicitly, then
    // this means that the value at that position is default one:
    return ( ( it == StorageEnd() ) ? m_default : it->second );
  }

  /// indexing operator for const sparse arrays. values can not
  /// be assigned to const sparse arrays so it's safe to return 
  /// the stored value directly as this method does
  T operator [] (unsigned int i) const { return Get(i); }

  /// indexing operator for non-const sparse arrays. Indexed elements
  /// of non-const sparse arrays must be lvalues - so that e.g. a new value
  /// could be assigned to them with array[i]=v syntax. In a sparse array,
  /// the actual storage at index i might not even exist at the moment of 
  /// invocation of array[i]. It is also not known at that moment whether
  /// we are simply accessing the stored value
  /// so that there is no need to create explicit storage for the position i in 
  /// the array, or trying to assign a new value (via array[i]=v or some other 
  /// call, e.g. array[i]++), in which case explicit storage must be created 
  /// if it did not exist before. To achieve this level of flexibility, the 
  /// evaluation should be deferred till the moment when the actual operation
  /// on array[i] is performed. This is why a trick with an auxilliary 'lvalue'
  /// class is used here. Depending on the context, different operators of that 
  /// lvalue class will either simply return the stored value (as in v=array[i]) or 
  /// allocate explicit storage (as in array[i]=v).
  _SparseArray_lvalue<T,CONTAINER> & operator[] (unsigned int i) {
    ForceAssertLt(i, m_size);
    m_lval.SetIndex(i);
    return m_lval;
  }


  class const_iterator {

    friend class SparseArray<T,CONTAINER>;
  private:
    const SparseArray<T,CONTAINER> & m_array; // who am i?
    typename SparseArray<T,CONTAINER>::t_const_iter m_iter; // where am i?
    unsigned int m_i; // index into the (unfolded) array
   
    /// private constructor; only friends can use. This constructor
    /// is used internally to initialize new iterator to point to the
    /// beginning of the sparse array
    const_iterator(const SparseArray<T,CONTAINER> *a) : 
       m_array(*a),
       m_i(0)
    { 
      if ( m_array.IsSorted() ) {
	// if the physical storage array is sorted,
	// m_iter always points either to the storage location
	// with the smallest index I such that I >= m_i or to the location
	// with the largest I such that I <= m_i;
	// this is more efficient since we do not have to 
	// search for each index value in the physical storage
	// (see also iterator increment/decrement operators)
	m_iter = m_array.StorageBegin();
      } else {
	// otherwise we do not know where (and if) individual index
	// is stored, so we attempt to find every index directly:
	m_iter = m_array.FindIndex(m_i);
      }
    }

    /// private constructor; only friends can use. This constructor
    /// is used internally to initialize new iterator to point to the
    /// end of the sparse array. Dereferencing an iterator constructed
    /// by this method is undefined (and will probably result in a segmentation
    /// fault; but may also return array's default value), 
    /// just like dereferencing of any iterator returned by container's
    /// ::end() method (this is what this constructor is used for). However,
    /// this iterator can be safely decremented (again, just like any other
    /// ::end() iterator). All this
    /// extra machinery is provided only to allow STL-style manipulation of 
    /// sparse arrays using iterators. 
    const_iterator(const SparseArray<T,CONTAINER> *a, unsigned int size) : 
       m_array(*a),
       m_iter( a->StorageEnd() ),
       m_i(size)
    {
      if ( m_array.IsSorted() ) {
	// if the physical storage array is sorted,
	// m_iter always points either to the storage location
	// with the smallest index I such that I >= m_i or to the location
	// with the largest I such that I <= m_i;
	// this is more efficient since we do not have to 
	// search for each index value in the physical storage
	// (see also iterator increment/decrement operators)
	m_array.RollbackIterator(m_iter);
      } // otherwise, leave m_iter pointing to StorageEnd() since
        // index value m_i=size is not in the array [0...size-1]
    }
  public:
    /// copy constructor:
    const_iterator(const const_iterator & it) :
      m_array(it.m_array),
      m_iter(it.m_iter),
      m_i(it.m_i)
    {}

    /// prefix increment
    const_iterator & operator++() {

      if ( m_array.IsSorted() && m_iter != m_array.StorageEnd()
	   && m_i >= m_iter->first ) {
	// we do need to perform this separate check here:
	// if ++ and -- operations are intertwined, we are
	// only guaranteed that m_iter points either to smallest
	// index I >= m_i or to largest I <= m_i; here we reduce
	// the representation: forward movement requires I >= m_i:
	++m_iter;
      } // now we are guaranteed that m_iter->first >= m_i

      ++m_i;

      if ( m_array.IsSorted() ) {
	
	// if we just crossed a physically stored element,
	// advance the iterator into the stored elements:
	if ( m_iter !=  m_array.StorageEnd() &&
	     m_i > m_iter->first ) {
	  ++m_iter;
	}
      } else {
	m_iter = m_array.FindIndex(m_i);
      }
      return *this;
    }

    /// postfix increment: SLOWER because of extra copy constructor!
    const_iterator operator++(int unused) {
      const_iterator tmp_it(*this);
      if ( m_array.IsSorted() && m_iter != m_array.StorageEnd()
	   && m_i >= m_iter->first ) {
	++m_iter;
      } 

      ++m_i;
      // if we just crossed a physically stored element,
      // advance the iterator into the stored elements:
      if ( m_array.IsSorted() ) {
	if ( m_iter != m_array.StorageEnd() &&
	     m_i > m_iter->first ) {
	  ++m_iter;
	}
      } else {
	m_iter = m_array.FindIndex(m_i);
      }
      return tmp_it;
    }

    /// prefix decrement
    const_iterator & operator--() {
      if ( m_array.IsSorted() ) {
	if ( m_iter == m_array.StorageEnd() ) m_array.RollbackIterator(m_iter);
	else {
	  if  ( m_i <= m_iter->first ) {
	    // we do need to perform this separate check here:
	    // if ++ and -- operations are intertwined, we are
	    // only guaranteed that m_iter points either to smallest
	    // index I >= m_i or to largest I <= m_i; here we reduce
	    // the representation: reverse movement requires I <= m_i:
	    m_array.RollbackIterator(m_iter);
	  }
	}
      } // now we are guaranteed that m_iter->first <= m_i
      --m_i;
      // if we just crossed a physically stored element,
      // decrease the iterator into the stored elements:
      if ( m_array.IsSorted() ) {
	if ( m_iter != m_array.StorageBegin() &&
	     m_i < m_iter->first ) m_array.RollbackIterator(m_iter);
      } else {
	m_iter = m_array.FindIndex(m_i);
      }
      return *this;
    }

    /// postfix increment: SLOWER because of extra copy constructor!
    const_iterator operator--(int unused) {
      const_iterator tmp_it(*this);
      if ( m_array.IsSorted() ) {
	if ( m_iter == m_array.StorageEnd() ) m_array.RollbackIterator(m_iter);
	else {
	  if  ( m_i <= m_iter->first ) {
	    m_array.RollbackIterator(m_iter);
	  }
	}
      } 
      --m_i;
      // if we just crossed a physically stored element,
      // decresase the iterator into the stored elements:
      if ( m_array.IsSorted() ) {
	if ( m_iter !=  m_array.StorageBegin() &&
	   m_i < m_iter->first ) --m_iter;
      } else {
	m_iter = m_array.FindIndex(m_i);
      }
      return tmp_it;
    }

    /// Dereferencing operator: returns a sparse array's element pointed
    /// to by this iterator; TODO: this is not an lvalue! we implement
    /// only constant (immutable) iterator so far!
    T operator*() {
      if ( m_array.IsSorted() ) {
	// we already keep the valid iterator; it points either to 
	// storage bin with the correct index (if the current index is stored
	// at all) or to a bin with the next or prev stored index (if the value at the 
	// current index is a default one, i.e. the index is not stored)
	return ( ( m_i == m_iter->first ) ? m_iter->second : m_array.GetDefault() );
      } else {
	// we attempted to retrieve m_iter on each iterator increment, but it can 
	// point to StorageEnd if the current index is not stored
	return ( m_iter == m_array.StorageEnd() ? m_array.GetDefault() : m_iter->second );
      }
      return m_array.GetDefault(); // should never get here, but compiler complains...
    }

    friend bool operator == ( const const_iterator & I1, const const_iterator & I2) {
      return I1.m_i == I2.m_i;
    }
    friend bool operator != ( const const_iterator & I1, const const_iterator & I2) {
      return I1.m_i != I2.m_i;
    }

  };

  /// Returns an iterator pointing to the beginning of the sparse array.
  /// This iterator will traverse over \em all elements of the array,
  /// not just the stored (non-default values) ones. Dereferencing this
  /// iterator at positions in the sparse array where no data is physically
  /// stored (i.e. array element has, by convention, a default value) will
  /// correctly return the default value associated with this array.
  const_iterator begin() const {
    return const_iterator(this);
  }

  /// Returns an iterator pointing right after the last element of the 
  /// sparse array. This iterator will traverse over \em all elements of 
  /// the array, not just the stored (non-default values) one. Dereferencing 
  /// this iterator at positions in the sparse array where no data is physically
  /// stored (i.e. array element has, by convention, a default value) will
  /// correctly return the default value associated with this array.
  const_iterator end() const {
    return const_iterator(this,size());
  }


    //  friend class SparseArray<T,CONTAINER>::const_iterator;
  typedef typename SparseArray<T,CONTAINER>::const_iterator iterator;


private:
    typedef typename IndexedStoragePolicy<T,CONTAINER>::t_iter t_iter;
    typedef typename IndexedStoragePolicy<T,CONTAINER>::t_const_iter t_const_iter;



  _SparseArray_lvalue<T,CONTAINER> m_lval;
  T m_default;
  unsigned int m_size;
};


#endif
