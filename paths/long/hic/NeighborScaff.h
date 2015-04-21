///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Oct 8, 2014 - <crdhelp@broadinstitute.org>
//

#ifndef NEIGHBORSCAFF_H_
#define NEIGHBORSCAFF_H_


void neighborScaff( int scaff_seed,
        ScaffoldPile const& scaffPile,
        HiCExperiment const& hic,
        int DEPTH,
        int MINLINESIZE,
        vec<int>& testStarts,
        size_t head = 0,
        Bool verb = False);


#endif /* NEIGHBORSCAFF_H_ */
