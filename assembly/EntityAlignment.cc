// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology
// 

#include "assembly/EntityAlignment.h"

template class BoundAlignment<Read,Read>;
template class BoundAlignment<Contig,Read>;
template class BoundAlignment<Read,Contig>;
template class BoundAlignment<Contig,Contig>;
