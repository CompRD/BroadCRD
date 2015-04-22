// Copyright (c) 2003 Broad Institute/Massachusetts Institute of Technology

#ifndef GENERATE_Q20_PLUSES_H
#define GENERATE_Q20_PLUSES_H

#include "math/Functions.h"
#include "Qualvector.h"
#include "system/System.h"
#include "Vec.h"



/*
 * GenerateQ20Pluses
 *
 * It generates a vec of q20pluses, where q20pluses[ii] is the number of
 * bases in quals[ii] with quality 20 or better. If a q20pluses file is
 * not found, then GenerateQ20Pluses will generate it by loading quals
 * in chunks (in this case it will save the newly generater q20pluses
 * in a file parallel to quals_file).
 *
 * q20pluses: the output
 * quals_file: qualb file
 * log: if not null, then it will send to log some progress output
 */
void GenerateQ20Pluses( vec<int> &q20pluses,
			const String &quals_file,
			ostream *log = 0 );



#endif
