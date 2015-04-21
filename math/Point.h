/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef POINT_H
#define POINT_H

/// Simple templatized struct with x and y coordinates.
template<class T = int>
struct Point {
  Point(const T & x = 0, const T & y = 0): x(x), y(y) {}
  T x;
  T y;

  Point operator+=(const Point<T> & o) { x+= o.x; y += o.y; return *this; }
  Point operator-=(const Point<T> & o) { x-= o.x; y -= o.y; return *this; }
};

template<class T>
inline bool operator==(const Point<T> & r, const Point<T> & l) {
  return r.x == l.x && r.y == l.y;
}

template<class T>
inline bool operator<(const Point<T> & r, const Point<T> & l) {
  if (r.x < l.x) return true;
  if (r.x > l.x) return false;
  if (r.y < l.y) return true;
  return false;
}

template<class T>
inline Point<T> operator+(const Point<T> & r, const Point<T> & l) {
  return Point<T>(r)+=l;
}

template<class T>
inline Point<T> operator-(const Point<T> & r, const Point<T> & l) {
  return Point<T>(r)-=l;
}

#endif
