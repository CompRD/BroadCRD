// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#ifndef LOOK_ALIGN_INDEX_H
#define LOOK_ALIGN_INDEX_H

#include "feudal/BinaryStreamTraits.h"

// The following class is used to index a file of QueryLookupTable aligns, ordered
// by (target_id, start_on_target).  The stop_on_target entry is the maximum 
// stop_on_target for all aligns in the file, up to but not including the next
// index.

class QLT_index {

     public:

     int target_id;
     int start_on_target;
     int stop_on_target;
     int alignment_id;
     longlong start_on_hits_file;

     QLT_index( ) { }

     QLT_index( int target_id_arg,
		int start_on_target_arg,
		int stop_on_target_arg,
		int alignment_id_arg,
		longlong start_on_hits_file_arg )
       : target_id(target_id_arg),
	 start_on_target(start_on_target_arg),
	 stop_on_target(stop_on_target_arg),
	 alignment_id(alignment_id_arg),
	 start_on_hits_file(start_on_hits_file_arg)
     { }

};
TRIVIALLY_SERIALIZABLE(QLT_index);

// PullAligns( dir, target, start, stop, buf ): given a directory dir as described
// in MakeLookAlignsIndex, pull all alignments of query sequences, which
// MIGHT align to target between start and stop, placing output in buf.

void PullAligns( const String& dir, int target, int start, int stop, 
     vec<char>& buf, vec<unsigned char>& mult, vec<QLT_index>* indexp = 0 );

#endif
