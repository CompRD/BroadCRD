// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

/** SubdirCopy.cc
    \file SubdirCopy.cc

    Copy the following files from the given subdir to the given outdir:

    supercontig structure:
      mergedcontigs.answer.gz
      mergedcontigs.superb
      mergedcontigs.superb_index (optional)
      mergedcontigs.summary (optional)

    contig data:
      mergedcontigs.fastb
      mergedcontigs.qualb

    read location data:
      mergedcontigs_orig.locs
      mergedcontigs_orig.locs_index (optional)
      mergedcontigs_orig.locs_indexr (optional)
      mergedcontigs.pma (optional)
      mergedcontigs.locs (optional)


    Needed files that do not exist will cause the function to assert
    before any files are copied.  Optional files that do not exist
    will be noted to cout, and execution will continue.
**/

#include "String.h"

void SubdirCopy( const String& fullSubdirPath, const String& fullOutdirPath );

void SubdirCopy( const String& fullRundirPath, const String& subdir, const String& outdir );
