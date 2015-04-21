///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file MultiVec.h
 * \author neilw
 * \date Oct 10, 2013
 *
 * \brief Wrap vec-like things into a single, virtual vec-like thing accessed with a linear index.
 *
 */
#ifndef MULTIVEC_H_
#define MULTIVEC_H_

template< class TVecLike, class TValueType = typename TVecLike::value_type >
class MultiVec
{
public:
    typedef TVecLike container_type;
    typedef TValueType value_type;

    MultiVec( std::vector<container_type> const& vlist ) : _vlist(vlist), _csize(), _vsize() {
    	build_csize_vsize(vlist);
    }

    MultiVec( MultiVec<TVecLike> const& that ) : _vlist( that._vlist ), _csize( that._csize ), _vsize( that._vsize ) {};

    // compiler-supplied destructor is OK

    MultiVec& operator=( MultiVec<TVecLike> const& that ) {
        if ( this != &that ) {
            _vlist = that._vlist;
            _csize = that._csize;
            _vsize = that._vsize;
        }
        return *this;
    }


    size_t size() const { return _csize.back(); }
    bool empty() const { return size() == 0; }

    value_type& front() { return obj(0); }
    value_type const& front() const { return obj(0); }
    value_type& back() { return obj(size()-1); }
    value_type const& back() const { return obj(size()-1); }
    value_type& operator[]( size_type idx ) { return obj(idx); }
    value_type const& operator[]( size_type idx ) const { return obj(idx); }
    value_type& at( size_type idx ) { return obj(idx); }
    value_type const& at( size_type idx ) const { return obj(idx); }


private:
    // care should be taken to keep these consistent
    std::vector<container_type>  _vlist;        // internal subcontainers
    std::vector<size_t>          _csize;        // cumulative sum of _vsize
    std::vector<size_t>          _vsize;        // size of each subcontainer

    value_type& obj( size_t idx ) {
        // use cumulative size array to figure out which subcontainer
        // holds what we're interested in and balk if we're out of
        // bounds.
    	auto itr = std::upper_bound( _csize.begin(), _csize.end(), idx );
    	if ( itr == _csize.end() ) OutOfBoundsReporter::oob("MultiVec", idx, size() );
        // iterator points to a _csize element, we translate it to an
        // offset so we can look at the corresponding element in the
        // other arrays
    	size_t vidx = itr - _csize.begin();
        // easiest to tell how far we are from the *end* of the subarray
    	size_t iidx = *itr - idx;
    	return _vlist[ vidx ][ _vsize[vidx] - iidx ];
    }

    void build_csize_vsize( std::vector<container_type> const& vlist ) {
    	_vsize.resize( vlist.size() );
    	_csize.resize( vlist.size() );

    	if ( vlist.size() > 0 ) {

    		_vsize[0] = _csize[0] = vlist[0].size();

                for ( size_t i = 1; i < vlist.size(); ++i ) {
                        _vsize[i] = vlist[i].size();
                        _csize[i] = _csize[i-1] + _vsize[i];
                }
    	}
    }
};

#endif // MULTIVEC_H_
