// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 

#ifndef LOAD_ON_HUMAN_ALIGNS_H
#define LOAD_ON_HUMAN_ALIGNS_H

#include "Vec.h"
#include "lookup/LookAlign.h"



/*
 * LoadOnHumanAligns
 *
 * Loads on-human alignments on a chromosome based level, and sorts
 * them. Returns the number of loaded alignments (which may be
 * smaller than the number of alignments in the original file).
 *
 * hits: output vec of hits (the output)
 * hits_file: input file with alignments
 * individuals: use only reads coming from individual 0 (may be null)
 * input_log: save ids of all reads found in hits (may be null)
 * excluded_log: save ids of excluded reads (may be null)
 */
int LoadOnHumanAligns( vec<look_align_plus> &hits,
		       const String &hits_file,
		       const vec<int> *individual = 0,
		       ostream *input_log = 0,
		       ostream *excluded_log = 0 );



#endif
