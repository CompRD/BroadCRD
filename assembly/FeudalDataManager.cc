///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "assembly/FeudalDataManagerTemplate.h"

template class FeudalDataManager<vecbasevector,basevector>;
template class FeudalDataManager<vecqvec,qvec>;
template class FeudalDataManager<veccompseq,CompressedSequence>;
template class FeudalDataManager<vecString,String>;

