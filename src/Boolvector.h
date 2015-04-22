/*
  The Broad Institute
  SOFTWARE COPYRIGHT NOTICE AGREEMENT
  This software and its documentation are copyright 2008 by the
  Broad Institute/Massachusetts Institute of Technology. All rights are
  reserved.

  This software is supplied without any warranty or guaranteed support
  whatsoever. Neither the Broad Institute nor MIT can be responsible for its
  use, misuse, or functionality.
  
* FILE Boolvector.h
*/

#ifndef _BOOLVECTOR_H
#define _BOOLVECTOR_H

#include "feudal/TrackingAllocator.h"
#include <vector>

// Eventually this should replace bitvector, once all the bitvector functionality
// has been added.  It is significantly faster than bitvector, especially with
// respect to resizing.
class boolvector: public std::vector<bool,typename DefaultAllocator<bool>::type>
{
  public:
    int isize() const    { return size();     }
};

#endif /* _BOOLVECTOR_H */

/******************************************************************/
/**************************[END OF Boolvector.h]**********************/
/******************************************************************/
/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
