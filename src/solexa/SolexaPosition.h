/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/** Class to store solexa metainformation: lane, tile, x, y.
*/

#include<iostream>
#include "math/Point.h"
#include "feudal/BinaryStreamTraits.h"

struct SolexaPosition {
  int lane;
  int tile;
  int x;
  int y;
  SolexaPosition(int lane=-1, int tile=-1, int x=-1, int y=-1):
    lane(lane), tile(tile), x(x), y(y) {}

  operator Point<int>() const {
    return Point<int>(x,y);
  }

  friend istream & operator>>(istream & in, SolexaPosition & s) {
    in >> s.lane >> s.tile >> s.x >> s.y;
    return in;
  }
  friend ostream & operator<<(ostream & out, const SolexaPosition & s) {
    out << s.lane << " " << s.tile << " " << s.x << " " << s.y << " ";
    return out;
  }
  friend bool operator==(const SolexaPosition & l, const SolexaPosition & r) {
    return l.lane == r.lane && l.tile == r.tile && l.x == r.x && l.y == r.y;
  }
  friend bool operator!=(const SolexaPosition & l, const SolexaPosition & r) {
    return !(l == r);
  }
  friend bool operator<(const SolexaPosition & l, const SolexaPosition & r) {
   if (l.lane < r.lane) return true;
   if (l.lane > r.lane) return false;
   if (l.tile < r.tile) return true;
   if (l.tile > r.tile) return false;
   if (l.x < r.x) return true;
   if (l.x > r.x) return false;
   if (l.y < r.y) return true;
   return false;
  }

};
TRIVIALLY_SERIALIZABLE(SolexaPosition);
