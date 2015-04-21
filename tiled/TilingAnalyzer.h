// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef TILING_ANALYZER_H
#define TILING_ANALYZER_H

#include "PackAlign.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "Vec.h"
#include "tiled/CPos.h"
#include "tiled/PaddedSeq.h"
#include "tiled/Tiling.h"



/*
 * class tiling_analyzer
 *
 * It performs basic operation onto a given tiling. To use: instantiate a
 * tiling_analyzer and SetPointers. Then each time SetTiling is called
 * data structures are cleaned up and refilled with the fresh data. You
 * can then use the analyzer to see bases and quals at any given point.
 */
class tiling_analyzer {
  
public:
  
  tiling_analyzer( );
  
  void SetPointers( const vecbasevector *r_bases,
		    const vecqualvector *r_quals = 0,
		    const vecbasevector *c_bases = 0,
		    const vecqualvector *c_quals = 0,
		    const vecbasevector *k_bases = 0,
		    const vecqualvector *k_quals = 0 );
  
  void SetTiling( const tiling *the_tiling );
  
  void PrintKnown( int ii, char &base, int &qual ) const;
  
  void PrintContig( int ii, char &base, int &qual ) const;

  void PrintRead( int read, int ii, char &base, int &qual ) const;
  
  int QualKnown( int ii ) const;

  int QualContig( int ii ) const;
  
  int QualRead( int read, int ii ) const;

  char QualToChar( int n_qual ) const;

  int StartOnContig( int read ) const;

  int StopOnContig( int read ) const;
  
  int WinBegin( ) const;

  int WinEnd( ) const;

  pair<int, int> CWin( ) const;

  pair<int, int> RWin( int pos ) const;

  const tiling &Tiling( ) const;
  
  int PaddedLengthKnown( ) const;
  
  int PaddedLengthContig( ) const;

  int PaddedLengthRead( int pos ) const;

  void TamedRectangles( int post,
			int width,
			const vec< int > &read_pos,
			vec< vec<char> > &rbases,
			vec< vec<int> > &rquals ) const;
  
  void Rectangles( int post,
		   int width,
		   vec< int > &read_pos,
		   vec< vec<char> > &rbases,
		   vec< vec<int> > &rquals ) const;
  
  void ContigToReadPos( int pos, vec<CPos> &rpos ) const;
  
  
private:
  
  void BuildStringsKnown( ) const;
  
  void BuildStringsContig( ) const;
  
  void BuildStringsRead( int pos ) const;
  
  int OnContig( int coordinate ) const;

  
private:
  
  // bases and quals (may be null)
  const vecbasevector *k_bases_;
  const vecbasevector *c_bases_;
  const vecbasevector *r_bases_;
  const vecqualvector *k_quals_;
  const vecqualvector *c_quals_;
  const vecqualvector *r_quals_;
  
  // the tiling
  const tiling *the_tiling_;
  
  // bases and quals as printable objects (filled on demand)
  mutable vec<char> str_k_bases_;
  mutable vec<int> str_k_quals_;
  mutable vec<char> str_c_bases_;
  mutable vec<int> str_c_quals_;
  mutable vec< vec<char> > str_r_bases_;
  mutable vec< vec<int> > str_r_quals_;

  // windows on master of aligned sequences (full, not just aligned parts)
  pair<int, int> c_win_;
  vec< pair<int, int> > r_win_;
  
};



#endif
