/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "jumping/EdgeListGraph.h"

std::ostream & operator<< (std::ostream &out, EdgeListGraph & G)
{
  EdgeListGraph::iterator it=G.begin(), last=G.end();
  for ( ; it!=last; ++it) {
    out << it->first << ": ";
    EdgeListGraph::Set::iterator sit=it->second.begin(),
      slast=it->second.end();
    for ( ; sit!=slast; ++sit) {
      out << *sit << " ";
    }
    out << "\n";
  }
  return out;
}
