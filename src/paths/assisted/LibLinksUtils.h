/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__ASSISTED__LIB_LINKS_UTILS__H
#define PATHS__ASSISTED__LIB_LINKS_UTILS__H

#include "PairsManager.h"
#include "graph/Digraph.h"
#include "lookup/LookAlign.h"
#include "paths/Alignlet.h"
#include "paths/assisted/CLibLinks.h"
#include "paths/assisted/CProxInfo.h"
#include "paths/assisted/CProx.h"
#include "paths/assisted/Proxgraph.h"
#include "paths/reporting/CLinkBundle.h"

/**
 * FillLibLinks
 *
 * Parse the given aligns and collect all pairs of linked oriented
 * contigs (sorted). Remark: it will also build the index first_.
 */
void FillLibLinks( const int MAX_CN,
		   const PairsManager &pairs,
		   const vec<bool> &isrep,
		   const vec<int> &cns,
		   const vec<alignlet> &aligns,
		   const vec<int> &idx,
		   vec<CLibLinks> &lib_links );

/**
 * FillRefLinks
 *
 * Parse the given aligns of contigs onto a ref, and collect all pairs
 * of linked oriented contigs. WARNING: it assumes hits is sorted by
 * start of contig on reference, but it does not check this is the
 * case. If you pass an unsorted hits, you will get meaningless
 * output. Output will be sorted, and the index first_ will be built.
 *
 * Only uniquely aligned contigs with copy number <= MAX_CN will be
 * taken into account.
 */
void FillRefLinks( const int MIN_REF_GAP,
		   const int MAX_REF_GAP,
		   const int MAX_CN,
		   const vec<look_align> &hits,
		   const vec<int> &cns,
		   CLibLinks &ref_links );

/**
 * LibLinksToDigraphE
 *
 * Genererate a proxgraph from sets of CLibLinks. One or more of
 * jumps, Jumps, and ref may be null, but not all of them. The graph
 * will contain 2 * n_contigs vertices (fw and rc, as 0[+], 0[-],
 * 1[+], 1[-], etc.). 
 */
void LibLinksToDigraphE( const int n_contigs,
			 const vec<CLibLinks> *jlinks,
			 const vec<CLibLinks> *Jlinks,
			 const CLibLinks *rlinks,
			 const CProxInfo *info,
			 proxgraph &graphE );



void AddBundlesToDigraphE( const digraphE<CLinkBundle> &bundles,
			   const CProxInfo *info,
			   proxgraph &graphE );
#endif
