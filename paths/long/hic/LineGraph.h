///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Aug 8, 2014 - <crdhelp@broadinstitute.org>
//

#ifndef LINEGRAPH_H_
#define LINEGRAPH_H_

#include "paths/long/hic/HiCDefs.h"

struct EdgeEvidence {
    double ratio;
    size_t count;
    size_t sumLen;
    size_t alt_config;
};
TRIVIALLY_SERIALIZABLE(EdgeEvidence);

extern template class digraphVE<TaggedLineId,EdgeEvidence>;
typedef digraphVE<TaggedLineId,EdgeEvidence> LineGraph;

void dumpLineGraph( LineGraph const& graph, String const& filename, int cseed = -1, Bool singletons = False,
        Bool edgeLabels = True);

#endif /* LINEGRAPH_H_ */
