// Copyright (c) 2003 Broad Institute/Massachusetts Institute of Technology
//

#ifndef PARSE_FASTA_H
#define PARSE_FASTA_H

#include "FastIfstream.h"
#include "String.h"
#include "system/System.h"
#include "Vec.h"



/*
 * ParseFasta
 *
 * It selects and sends to out some of the entries from a fasta file (bases
 * or quality scores). It saves the number of found objects and the number
 * of parsed objects.
 */
void ParseFasta( const String &in_file,
		 const vec<String> &names,
		 int &n_found,
		 int &n_parsed,
		 ostream &out );



#endif
