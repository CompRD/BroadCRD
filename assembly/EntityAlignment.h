// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology
// 

// Instantiations of BoundAlignment for Reads and Contigs.

#ifndef ASSEMBLY_ENTITYALIGNMENT
#define ASSEMBLY_ENTITYALIGNMENT

#include "assembly/BoundAlignment.h"

#include "assembly/AssemblyContig.h"
#include "assembly/AssemblyRead.h"

typedef BoundAlignment<Read,Read>     ReadReadAlignment;
typedef BoundAlignment<Contig,Read>   ContigReadAlignment;
typedef BoundAlignment<Read,Contig>   ReadContigAlignment;
typedef BoundAlignment<Contig,Contig> ContigContigAlignment;
extern template class BoundAlignment<Read,Read>;
extern template class BoundAlignment<Contig,Read>;
extern template class BoundAlignment<Read,Contig>;
extern template class BoundAlignment<Contig,Contig>;

#endif
