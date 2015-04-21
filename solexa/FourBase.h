/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef FOUR_BASE_H
#define FOUR_BASE_H

#include "Basevector.h"
#include "math/Matrix.h"
#include "solexa/FourBase2.h"
#include "feudal/MasterVec.h"
#include "feudal/SerfVec.h"

/*
   Class: four_base

   Represents our confidence, for each of the four possible bases, that a given
   position of a read is that base.

   Related to <quality scores>.
*/
class four_base {

  typedef short value_type;

public:
  typedef float intens_t;

  four_base( ) { }
  four_base( intens_t A, intens_t C, intens_t G, intens_t T )
  {
    Init(A,C,G,T);
  }

  four_base(const four_base2 & f) {
    if (!f.ShortOk()) cout << f << endl;
    Init(f.A(), f.C(), f.G(), f.T());
  }

  value_type Cast(intens_t f) {
    if (f >= numeric_limits<value_type>::max())
      return numeric_limits<value_type>::max();
    if (f <= numeric_limits<value_type>::min())
      return numeric_limits<value_type>::min();
    else return static_cast<value_type>(f);
  }

  void Init(intens_t A, intens_t C, intens_t G, intens_t T ) {
    base_[0] = Cast(A);
    base_[1] = Cast(C);
    base_[2] = Cast(G);
    base_[3] = Cast(T);
  }

  intens_t A( ) const { return intens_t(base_[0]); }
  intens_t C( ) const { return intens_t(base_[1]); }
  intens_t G( ) const { return intens_t(base_[2]); }
  intens_t T( ) const { return intens_t(base_[3]); }

  intens_t base(int i) const { return intens_t(base_[i]); }

  friend ostream& operator<<( ostream& out, const four_base& b );

  // Method: Call
  // Call a base.  This always returns an answer, in the case of equality
  // favoring A > C > G > T.  Return 0, 1, 2, or 3.
  // See also <CallQuality()>.
  int Call( ) const;

  // Method: CallQuality
  // Return a primitive metric with value between 1 and infinity
  // (1000000000).  It's the ratio of the value for the best base to the next
  // best base.  The worse value (1) should be interpreted for now as
  // "completely random".
  float CallQuality( ) const;

  intens_t MaxInt() const { return intens_t(*max_element(base_,base_+4)); }

  intens_t Sum() const { return accumulate(base_,base_+4, 0.f); }

  intens_t Max() const
  {
        intens_t max = -32768;
        for (int i = 0; i < 4; i++)
        {
            if (base_[i] > max) { max = base_[i]; }
        }
        return max;
  }

  intens_t Min() const
  {
        intens_t min = 32768;
        for (int i = 0; i < 4; i++)
        {
            if (base_[i] < min) { min = base_[i]; }
        }
        return min;
  }

  friend
  bool operator==(const four_base & lhs, const four_base & rhs) {
    return 0 == memcmp(lhs.base_, rhs.base_, 4 * sizeof(value_type));
  }

private:

  value_type base_[4];

};  // class four_base

inline bool operator !=(const four_base & lhs, const four_base & rhs) {
  return !(lhs == rhs);
}

TRIVIALLY_SERIALIZABLE(four_base);

typedef SerfVec<four_base> FourBaseVec;
typedef MasterVec<FourBaseVec> VecFourBaseVec;
extern template class SmallVec< four_base, MempoolAllocator<four_base> >;
extern template class OuterVec<FourBaseVec>;

/// Call bases for an entire intensity vector:
void Call( const VecFourBaseVec& I, vecbasevector& bases );

/// Decay in signal of top base.
/// Return sum(signal[first, first+window))/sum(signal[second,second+window)),
/// or 0 if any of the signals are 0.
/// The vectors must be of size nbases.
template<class V = FourBaseVec >
class Decay: public unary_function<V, double> {
private:
  // TODO: potentially dangerous truncation of indices
  int nbases_, first_, second_, window_;

public:

  Decay(int nbases, int first=0, int second=20, int window=5):
   nbases_(nbases), first_(first), second_(second), window_(window) {
     first_ = min(nbases_ - window_, first_);
     second_ = min(nbases_ - window_, second_);
     ForceAssertGt(second_, first_);
     ForceAssertGe(first_, 0);
     ForceAssertGt(second_, 0);
   }

  double operator()(const V & v) {
    float sum1=0, sum2=0;
    AssertEq(v.size(),static_cast<typename V::size_type>(nbases_));
    for (int i=0; i !=window_; ++i) {
      float f1 = v[i + first_].MaxInt();
      float f2 = v[i + second_].MaxInt();
      if (f1 <=0 || f2 <=0) return 0;
      sum1 += f1;
      sum2 += f2;
    }
    return sum2/sum1;
  }
};

/// Decay in signal of all bases together
/// Return sum(signal[first, first+window))/sum(signal[second,second+window)),
/// or 0 if any of the signals are 0.
/// The vectors must be of size nbases.
template<class V = FourBaseVec >
class DecayNoPhasing: public unary_function<V, double> {
private:
    // TODO: potentially dangerous truncation of indices
  int nbases_, first_, second_, window_;

public:

  DecayNoPhasing(int nbases, int first=0, int second=20, int window=5):
   nbases_(nbases), first_(first), second_(second), window_(window) {
     first_ = min(nbases_ - window_, first_);
     second_ = min(nbases_ - window_, second_);
     ForceAssertGt(second_, first_);
     ForceAssertGe(first_, 0);
     ForceAssertGt(second_, 0);
  }

  double operator()(const V & v) {
    float sum1=0, sum2=0;
    AssertEq(v.size(),static_cast<typename V::size_type>(nbases_));
    for (int i=0; i !=window_; ++i) {
      float f1 = v[i + first_].Sum();
      float f2 = v[i + second_].Sum();
      if (f1 <=0 || f2 <=0) return 0;
      sum1 += f1;
      sum2 += f2;
    }
    return sum2/sum1;
  }
};

/// Output the average intensity over the four bases for the first cycle.
/// Use this measure rather than height of highest base because this is
/// what Solexa uses and outputs in their metrics.
template<class V>
float FirstHeight(const V & v) {
  return v[0].Sum()/4.0;
}

/// Assumes the joint distribution has no covariance.
template<class V>
float ZscoreDecayHeight(const V & v, const NormalDistribution & decay,
			const NormalDistribution & height) {
  if (v.empty()) return 0;			
  Decay<V> dec(v.size());
  float zd = decay.Zscore(dec(v));
  float zh = height.Zscore(FirstHeight(v));
  return sqrt(zd * zd + zh * zh);
}

template<class V>
float MinQuality(const V & v, int end, int start=0) {
  float ret = 1e6;
  start = max(start, (int)0);
  end = min(end, (int) (v.size()) );
  for (int i=start; i != end; ++i) {
    ret = min(ret, v[i].CallQuality());
  }
  return ret;
}

/// Return the lowest MaxInt() in [start, end)
template<class V>
float LowestHeight(const V & v, int end, int start=0) {
  float ret = 1e6;
  start = max(start, (int)0);
  end = min(end, (int) (v.size()) );
  for (int i=start; i != end; ++i) {
    ret = min(ret, v[i].MaxInt());
  }
  return ret;
}

///Return the inverse of the greatest change between two successive cycles.
///Returns 0 if any cycle has a Sum() of 0.
template<class V>
float MinIncreaseOrDecrease(const V & v, int bases = -1) {
  float ret = 1;
  if (-1 == bases) bases = v.size();
  for (int i=0; i < bases-1; ++i) {
    if (0 == v[i].Sum() || 0 == v[i+1].Sum()) return 0;
    float r = v[i].Sum() / v[i+1].Sum();
    ret = std::min( {ret, r, 1/r} );
  }
  return ret;
}

///Return true if any cycle has Sum() <= 0
template<class V>
bool HasDot(const V & v, int start = 0, int end= -1) {
  if (-1 == end) end = v.size();
  AssertLe(start, end);
  for (int i=start; i != end; ++i) {
    if (0 >= v[i].Sum()) return true;
  }
  return false;
}

class DeconvolveIntensities {
  matrix<double> m;
  Permutation permut;
  double work[4];
public:
  DeconvolveIntensities(String filename); // Read crosstalk file from
					  // disk; if filename is
					  // empty, do nothing
  // Deconvolve intensities
  void operator()(double * v)
  {
    m.solve_from_lu(permut, v, work);
    for (int i=0; i<4; ++i)
      v[i] = work[i];
  }
  Bool empty() { return m.empty(); } // Whether we have a matrix.
};

#endif
