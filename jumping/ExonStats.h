/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef EXON_STATS_H
#define EXON_STATS_H

/** Store statistics for an exon.

    \class ExonStats

  The stats stored include accession name, chromosome, start,
  length, GC%, Tm, and distance to nearest exon on the same gene.

*/
#include "String.h"
#include "Vec.h"
#include <istream>
#include <ostream>
#include <climits>
#include <functional>

struct ExonStats {
  ExonStats(const String & cdna = String(), char strand = 0, 
	    int exonIndex = -1, int chromosome = -1,
	    int start = -1, int length = -1, 
	    float gc = -1, float tm = -1, int closestExon = INT_MAX):
    cdna (cdna),
    strand (strand),
    exonIndex (exonIndex),
    chromosome (chromosome),
    start (start),
    length (length),
    gc (gc),
    tm (tm),
    closestExon (closestExon)
  {}

  String cdna;
  char strand;
  int exonIndex;
  int chromosome;
  int start;
  int length;
  float gc;
  float tm;
  int closestExon;

  //Note: these methods do not need to be friends right now, but
  //this helps the doc and the instance variables might be private later.

  friend std::ostream & operator<<(std::ostream & os, const ExonStats & es);
  friend std::istream & operator>>(std::istream & is, ExonStats & es);
  ///Order by chromosome and position on chromosome
  friend bool operator<(const ExonStats & l, const ExonStats & r);
};


///Write to a text file.
inline std::ostream & operator<<(std::ostream & os, const ExonStats & es) {
  os << es.cdna << "\t" << es.exonIndex << "\t" 
     << es.chromosome << "\t" << es.start << "\t" 
     << es.length << "\t" << es.gc << "\t" 
     << es.tm << "\t" << es.closestExon;
  return os;
}

///Read from a text file.
inline std::istream & operator>>(std::istream & is, ExonStats & es) {
  is >> es.cdna >> es.exonIndex 
     >> es.chromosome >> es.start 
     >> es.length >> es.gc 
     >> es.tm >> es.closestExon;
  return is;
}

///Order by location in genome.
inline bool operator<(const ExonStats & l, const ExonStats & r) {
  if (l.chromosome < r.chromosome) return true;
  if (l.chromosome > r.chromosome) return false;
  if (l.start < r.start) return true;
  return false;
}

///Order by accession name and exon number.
struct SortByCdna: public std::binary_function<ExonStats,ExonStats,bool> {
  bool operator()(const ExonStats & l, const ExonStats & r) {
    if (l.cdna < r.cdna) return true;
    if (l.cdna > r.cdna) return false;
    if (l.start < r.start) return true;
    return false;
  }
};

///Load a set of exon information from a browser line.
void ReadFromBrowserLine(const String & line, vec<ExonStats> & exons);

#endif // EXON_STATS_H
