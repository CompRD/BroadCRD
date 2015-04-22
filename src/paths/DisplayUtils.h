/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   Header file: DisplayUtils.h

   Utilities to simplify the display of assembly information
   (both the final assembly and various intermediate stages)
   as HTML.
   
   @file
*/

#ifndef __INCLUDE_paths_DisplayUtils_h
#define __INCLUDE_paths_DisplayUtils_h

#include "Vec.h"
#include "CommonSemanticTypes.h"
#include "STLExtensions.h"
#include "PackAlign.h"
#include "system/HTMLUtils.h"

typedef triple< align_id_t, align_id_t, html_t > placement_link_t;

/**
   FuncDecl: ShowGenomePlacements

   Given a set of possibly overlapping ranges on the reference,
   construct an HTML representation of how these ranges align
   to the reference (all ranges must align to the same <genome part).
   The ranges can be: reads, <unipaths>,
   <trusted paths>, or any other ranges.   For each range
   you just need to say where it aligns to the reference,
   and how to represent that range in HTML.
   This function will then construct an HTML table showing
   all the ranges.  The size of the table is proportional not
   to the lengths of the ranges, but to the number of overlap
   boundaries between the ranges.

   Each range is represented by an object of a class that models
   the BasicAlign concept defined in PackAlign.h .
   For each range, in parallel arrays, the caller needs to provide
   the HTML that will be put in the table cell(s) representing
   the range (for example, unipath number or trusted path name).
   The HTML contents for each range is given in placement_td_contents,
   and the HTML attributes are given in placement_td_attrs.

   *Important*: note, once more, that all ranges must align to the
   same <genome part> (i.e. to the same chromosome or more generally
   to the same record of the fastb file).   Thus, if the ranges are are
   represented by look_align's, the TargetId() method of each range
   must return the same value.

   The method returns an HTML table representing the ranges.
   You have to add valid HTML headers and footers to create a valid
   HTML page; see HTMLHead() and HTMLTail() in system/HTMLUtils.h .

   If you have several <genome parts>, right now you have to split
   your aligns into several vectors and call ShowGenomePlacements()
   separately for each.   For each genome part you then get
   an HTML table representing aligns to that part.
   A wrapper should be added which does the splitting and returns
   several HTML tables.

   Parameters:

      placements - for each of the overlapping ranges, its placement
         on the reference.  Each placement is represented by
	 an object of a class that models the BasicAlign concept
	 defined in PackAlign.h

      placement_td_contents - for each range, the HTML content
         representing the range.  this will be put between
	 <TD> and </TD> tages in the table cell(s) representing
	 the range.

      placement_td_attrs - for each range, the (optional) HTML
         attributes representing the range.  these will be added
	 as attributes to the <TD ... > element.  you can use this
	 for example to color your ranges into random colors.

      GPART - the <genome part> number for the genome part to which all
         the ranges in 'placements' align; this is displayed to the left
	 of the table.  if negative, interpreted as the rc of
	 -GPART-1.  this parameter only affects the display of the genome
	 part number to the left of the table.  however, all aligns in
	 'placements' should really be aligns to this genome part.


       genomePartSize - the size of the genome part, to which all
         the ranges align.  this is so we can show the correct
	 right margin (how many bases are to the right of the rightmost
	 range).

       placement_links - optional, used to show links between the ranges
         (for example scaffold links between contigs).
	 Each links is represented as
	 triple< align_id_t, align_id_t, html_t >
	 where the two align_id's give the positions of the two linked
	 ranges in 'placements', and the html is the HTML representation
	 of the link (for example giving the names of the two linked contigs
	 and information about link strength.)
      
*/
template <class GENOME_PLACEMENT>
html_t ShowGenomePlacements( const vec< GENOME_PLACEMENT >& placements,
			     const vec< html_t >& placement_td_contents,
			     const vec< html_attrs_t >& placement_td_attrs,
			     genome_part_id_t GPART, nbases_t genomePartSize,
			     const vec< placement_link_t > *placement_links = NULL
			     );


#endif
// #ifndef  __INCLUDE_paths_DisplayUtils_h

