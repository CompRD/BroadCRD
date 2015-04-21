// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology
//

#ifndef IDENTIFY_DIFFERENCES_H
#define IDENTIFY_DIFFERENCES_H

#include "tiled/CKDiff.h"
#include "tiled/ColumnConsensus.h"
#include "tiled/Tiling.h"
#include "tiled/TilingAnalyzer.h"
#include "tiled/TilingPrinter.h"
#include "tiled/VecTiling.h"



/*
 * IdentifyDifferences
 *
 * Tag as bad all bases for which known and consensus do not match (indels
 * included). Returns basic statistics for the differences found.
 */
ck_diff IdentifyDifferences( int radius,
			     int tile_id,
			     tiling_analyzer &analyzer,
			     ostream &out );



#endif
