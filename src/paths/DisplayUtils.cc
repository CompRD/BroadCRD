/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   C file: DisplayUtils.cc

   Utilities to simplify the display of assembly information
   (both the final assembly and various intermediate stages)
   as HTML.
   
   @file
*/

#include <sstream>
#include <algorithm>
#include "paths/DisplayUtils.h"
#include "paths/DisplayUtilsTemplate.h"
#include "paths/simulation/GenomePlacement.h"
#include "lookup/LookAlign.h"

DEFINE_DISPLAY_UTILS(genome_placement);
DEFINE_DISPLAY_UTILS(look_align);
DEFINE_DISPLAY_UTILS(GaplessAlign);

