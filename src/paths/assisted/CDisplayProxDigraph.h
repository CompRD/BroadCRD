///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__ASSISTED__C_DISPLAY_PROX_DIGRAPH__H
#define PATHS__ASSISTED__C_DISPLAY_PROX_DIGRAPH__H

#include "Basevector.h"
#include "graph/Digraph.h"
#include "paths/assisted/CProx.h"
#include "paths/assisted/Proxgraph.h"
#include "util/RunCommand.h"

/**
 * CDisplayProxDigraph
 *
 * Class to display a proxgraph (interactive, see also
 * DisplayHyperKmerPath). Notice that the vertices are oriented
 * contigs (0[+], 0[-], 1[+], 1[-], etc.)
 */
class CDisplayProxDigraph {

public:
  
  CDisplayProxDigraph( const proxgraph &graph, const vecbvec &contigs );
  
  // Entry point (wait for input from cin).
  void GoDisplay( ) const;
  
  
private:
  
  void Setup( );
  
  void DefaultArgs( );

  void PrintArgs( ) const;

  bool FetchCommand( istream &in ) const;
  
  bool ExecCommand( const String &cmnd ) const;

  void AllPredecessors( vec<int> &data, int &level, int rad, int vtx ) const;

  void AllSuccessors( vec<int> &data, int &level, int rad, int vtx ) const;
  
  void ContigLengthFilter( vec<int> &edges ) const;

  void SmartRefFilter( vec<int> &edges ) const;

  void BasicStats( ) const;

  void DisplayRadius( int vtx, int rad ) const;

  void DisplayComponent( int *c_id = 0, int *v_id = 0 ) const;

  void DisplayRegion( int v1, int v2 ) const;

  void DotAndPlot( const vec<int> &vertices ) const;
  
  
private:
  
  // Constant input data (graph and contigs).
  const digraphE<CProx> &graph_;
  const vecbvec &contigs_;
  
  // Arguments (see documentation in DefaultArgs( ), in the .cc file).
  mutable bool SMART_REF_FILTER_;
  mutable int MIN_CLEN_;
  
  // Standard maps from edges to left and right vertices.
  vec<int> to_left_;
  vec<int> to_right_;
  
};

#endif

