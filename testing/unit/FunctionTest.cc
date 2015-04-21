/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "math/Functions.h"
#include "TestAssert.h"

template<class T>
pair<typename vec<T>::const_iterator, typename vec<T>::const_iterator>
Neighborhood(const vec<T> & q, int pos, int nsize) {
  AssertLt(pos, q.isize());
  pair<typename vec<T>::const_iterator, typename vec<T>::const_iterator> ret;
  ret.first = max(q.begin(), q.begin() + pos - nsize);
  ret.second = min(q.end(), q.begin() + pos + nsize + 1);
  return ret;
}

template<class T>
T MinNeighborhood(const vec<T> & q, int pos, int nsize) {
  pair<typename vec<T>::const_iterator, typename vec<T>::const_iterator> n
    = Neighborhood(q,pos,nsize);
  return *min_element(n.first, n.second);
}

template<class T>
T MaxNeighborhood(const vec<T> & q, int pos, int nsize) {
  pair<typename vec<T>::const_iterator, typename vec<T>::const_iterator> n
    = Neighborhood(q,pos,nsize);
  return *max_element(n.first, n.second);
}

template<class T>
T SumNeighborhood(const vec<T> &q, int pos, int nsize) {
  pair<typename vec<T>::const_iterator, typename vec<T>::const_iterator> n
    = Neighborhood(q,pos,nsize);
  return accumulate(n.first, n.second, T(0));
}

template<class T>
double MeanNeighborhood(const vec<T> &q, int pos, int nsize) {
  pair<typename vec<T>::const_iterator, typename vec<T>::const_iterator> n
    = Neighborhood(q,pos,nsize);
  return accumulate(n.first, n.second, double(0))/distance(n.first, n.second);
}

class FunctionTester {
  int verbose;

public:

  FunctionTester(int verbose): verbose(verbose) { }
  
  int TestCosine() {
    vec<short> v1, v2, v3, v4;
    v1.push_back(1,1);
    v2.push_back(2,2,44);
    v3.push_back(-1,1);
    v4.push_back(32000,0);
    double c = Cosine(v1.begin(), v1.end(), v2.begin(), v2.end());
    cout << c << endl;
    TestAssertLt(c-1,0.001);
    c = Cosine(v1.begin(), v1.end(), v3.begin(), v3.end());
    cout << c << endl;
    TestAssertLt(c,0.001);
    c = Cosine(v2.begin(), v2.end(), v3.begin(), v3.end());
    cout << c << endl;
    TestAssertLt(c,0.001);
    c = Cosine(v1.begin(), v1.end(), v4.begin(), v4.end());
    cout << c << endl;
    TestAssertLt(c-1.0/sqrt(2),0.001);
    cout << "Cosine tests passed." << endl;
    return 0;
  }

  int TestQuadratic() {
    QuadraticFunction q(1,-2,1);
    double epsilon = 0.0001;
    TestAssert(abs(q(1) -0) < epsilon);
    pair<double,double> sol = q.solutions();

    TestAssert(abs(sol.first -sol.second) < epsilon);
    TestAssert(abs(sol.first - 1) < epsilon);

    q = QuadraticFunction(1,0,1);
    TestAssert(abs(q(1) -2) < epsilon);
    sol = q.solutions();
    TestAssert(!isfinite(sol.first));

    q = QuadraticFunction(2,0,-2);
    TestAssert(abs(q(2) -6) < epsilon);
    sol = q.solutions();
    TestAssert(abs(sol.first + 1) < epsilon);
    TestAssert(abs(sol.second - 1) < epsilon);

    q = QuadraticFunction(1,2,-1);
    TestAssert(abs(q(1) -2) < epsilon);
    sol = q.solutions();
    TestAssert(abs(sol.first + 2.4142) < epsilon);
    TestAssert(abs(sol.second - 0.4142) < epsilon);

    q = QuadraticFunction(.5,1,-0.5);
    TestAssert(abs(q(1) -1) < epsilon);
    sol = q.solutions();
    TestAssert(abs(sol.first + 2.4142) < epsilon);
    TestAssert(abs(sol.second - 0.4142) < epsilon);

    cout << "QuadraticFunction tests passed " << endl;
    return 0;
  }

  int TestOptimalCutoff() {
    NormalDistribution d1(0,1);
    NormalDistribution d2(1,1);

    TestClose(OptimalCutoff(d1,d2),0.5);
    TestClose(OptimalCutoff(d1,d2,1),0.5);
    TestClose(OptimalCutoff(d1,d2,1.001),0.5005);
    TestClose(OptimalCutoff(d1,d2,1.1), 0.547655);
    TestClose(OptimalCutoff(d1,d2,2), 0.846574);
 
    d2.sigma_=sqrt(2);
    TestClose(OptimalCutoff(d1,d2),0.414214);
    TestClose(OptimalCutoff(d1,d2,1),0.414214);
    TestClose(OptimalCutoff(d1,d2,1.001),0.41492);
    TestClose(OptimalCutoff(d1,d2,1.1),0.480074);
    TestClose(OptimalCutoff(d1,d2,2),0.840189);
    OptimalCutoff(d1,d2,100);

    d2.sigma_=1/sqrt(2);
    OptimalCutoff(d1,d2);
    OptimalCutoff(d1,d2,1);
    OptimalCutoff(d1,d2,1.001);
    OptimalCutoff(d1,d2,1.1);
    OptimalCutoff(d1,d2,2);
    OptimalCutoff(d1,d2,100);

    cout << "OptimalCutoff tests passed " << endl;

    return 0;
  }

  template <typename T>
  void check(const vec<int> &v, const vec<T> &m, const vec<T> &e, int n)
  {
    if (verbose) {
      cout << "\nInput:  ";
      for (int i=0; i<v.isize(); ++i)
	cout << v[i];
      cout << "\tn=" << n;
      cout << "\nOutput: ";
      for (int i=0; i<v.isize(); ++i)
	cout << m[i];
      cout << "\nExpect: ";
      for (int i=0; i<v.isize(); ++i)
	cout << e[i];
      cout << "\n";
    }
    TestAssert(equal(m.begin(), m.end(), e.begin()));
  }
  
  void checkMinWindow(const vec<int> &v)
  {
    vec<int> e(v.size()), m;
    for (int nbhd=0; nbhd<10; ++nbhd) {
      for (int i=0; i<v.isize(); ++i) {
	e[i] = MinNeighborhood(v,i,nbhd);
      }
      m.clear();
      MinWindow(v.begin(), v.end(), back_inserter(m), nbhd);
      check(v,m,e,nbhd);
    }
  }
  
  void checkMaxWindow(const vec<int> &v)
  {
    vec<int> e(v.size()), m;
    for (int nbhd=0; nbhd<10; ++nbhd) {
      for (int i=0; i<v.isize(); ++i) {
	e[i] = MaxNeighborhood(v,i,nbhd);
      }
      m.clear();
      MaxWindow(v.begin(), v.end(), back_inserter(m), nbhd);
      check(v,m,e,nbhd);
    }
  }

  void checkSumWindow(const vec<int> &v)
  {
    vec<int> e(v.size()), m;
    for (int nbhd=0; nbhd<10; ++nbhd) {
      for (int i=0; i<v.isize(); ++i) {
	e[i] = SumNeighborhood(v,i,nbhd);
      }
      m.clear();
      SumWindow(v.begin(), v.end(), back_inserter(m), nbhd);
      check(v,m,e,nbhd);
    }
  }

  void checkMeanWindow(const vec<int> &v)
  {
    vec<double> e(v.size()), m;
    for (int nbhd=0; nbhd<10; ++nbhd) {
      for (int i=0; i<v.isize(); ++i) {
	e[i] = MeanNeighborhood(v,i,nbhd);
      }
      m.clear();
      MeanWindow(v.begin(), v.end(), back_inserter(m), nbhd);
      check(v,m,e,nbhd);
    }
  }

  int TestMinWindow() {
    vec<int> v, m, e;

    v.push_back(7,6,5,4,3,2,1,0);
    checkMinWindow(v);

    v.clear();
    v.push_back(0,1,2,3,4,5,6,7);
    checkMinWindow(v);

    v.clear();
    v.push_back(3,1,4,1,5,9,2,6);
    checkMinWindow(v);

    v.clear();
    v.push_back(3);
    checkMinWindow(v);

    v.clear();
    checkMinWindow(v);

    cout << "MinWindow tests passed " << endl;
    return 0;
  }

  int TestMaxWindow() {
    vec<int> v, m, e;

    v.push_back(7,6,5,4,3,2,1,0);
    checkMaxWindow(v);

    v.clear();
    v.push_back(0,1,2,3,4,5,6,7);
    checkMaxWindow(v);

    v.clear();
    v.push_back(3,1,4,1,5,9,2,6);
    checkMaxWindow(v);

    v.clear();
    v.push_back(3);
    checkMaxWindow(v);

    v.clear();
    checkMaxWindow(v);

    cout << "MaxWindow tests passed " << endl;
    return 0;
  }

  int TestSumWindow() {
    vec<int> v, m, e;

    v.push_back(3,2,1,0,1,0,1,0);
    checkSumWindow(v);

    v.clear();
    v.push_back(0,1,2,3,1,0,1,1);
    checkSumWindow(v);

    v.clear();
    v.push_back(0,0,0,0,3,1,4,1);
    checkSumWindow(v);

    v.clear();
    v.push_back(3);
    checkSumWindow(v);

    v.clear();
    checkSumWindow(v);

    cout << "SumWindow tests passed " << endl;
    return 0;
  }

  int TestMeanWindow() {
    vec<int> v, m, e;

    v.push_back(3,2,1,0,1,0,1,0);
    checkMeanWindow(v);

    v.clear();
    v.push_back(0,1,2,3,1,0,1,1);
    checkMeanWindow(v);

    v.clear();
    v.push_back(0,0,0,0,3,1,4,1);
    checkMeanWindow(v);

    v.clear();
    v.push_back(3);
    checkMeanWindow(v);

    v.clear();
    checkMeanWindow(v);

    cout << "MeanWindow tests passed " << endl;
    return 0;
  }


  void TestInverseNormalCDF()
  {
    // Walk through the standard points on the standard normal distribution
    TestCloserEpsilon(-3.0, InverseNormalCDF(1.0 - 0.998650101968370), 1e-8);
    TestCloserEpsilon(-2.0, InverseNormalCDF(1.0 - 0.977249868051821), 1e-8);
    TestCloserEpsilon(-1.0, InverseNormalCDF(1.0 - 0.841344746068543), 1e-8);
    TestCloserEpsilon(-0.5, InverseNormalCDF(1.0 - 0.691462461274013), 1e-8);
    TestAssert(0.0==InverseNormalCDF(0.5));
    TestCloserEpsilon(0.5, InverseNormalCDF(0.691462461274013), 1e-8);
    TestCloserEpsilon(1.0, InverseNormalCDF(0.841344746068543), 1e-8);
    TestCloserEpsilon(2.0, InverseNormalCDF(0.977249868051821), 1e-8);
    TestCloserEpsilon(3.0, InverseNormalCDF(0.998650101968370), 1e-8);
    // A couple of checks that variant mu and sigma work correctly, at the z=1 point.
    TestCloserEpsilon(10.0, InverseNormalCDF(0.841344746068543, 0.0, 10.0), 1e-8);
    TestCloserEpsilon(10.0, InverseNormalCDF(0.841344746068543, 9.0, 1.0), 1e-8);
    cout << "InverseNormalCDF tests passed " << endl;
  }

  void TestBinomialCI()
  {
    pair<double, double> res;
    for (int n=1; n<10000; n *= 2) {
      res = BinomialConfidenceInterval(0, n);
      TestAssert(res.first==0.0);
      TestCloserEpsilon(res.second, 3.84/(n+3.84), 0.01);
    }
    cout << "BinomialConfidenceInterval tests passed " << endl;
  }

};

int main(int argc, char ** argv) {
  FunctionTester t(argc-1);
  t.TestCosine();
  t.TestQuadratic();
  t.TestOptimalCutoff();
  t.TestMinWindow();
  t.TestMaxWindow();
  t.TestSumWindow();
  t.TestMeanWindow();
  t.TestInverseNormalCDF();
  t.TestBinomialCI();
  return 0;
}
