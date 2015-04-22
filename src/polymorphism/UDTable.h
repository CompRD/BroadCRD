/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// UeberDrueber.  Find polymorphisms in test strain REFA vs finished strain REFB.
#include "MainTools.h"
#include "Vec.h"
#include "Map.h"
#include "Basevector.h"
#include "String.h"

class UDTable {
public:
   // Replace below by hash map once we know this works.
   vec< StdMap<basevector, vec<int> > > theData;
   // contig   string     list of pos
   void Read(String filename);
   void Write(String filename);
   void GenerateTable(vecbasevector &ref);
   int blockSize;
   // The initialization # of the below is purely arbitrary.  
   // Smaller will find you doing longer Smith-Waterman at the end, but 
   // less large-scale aligning. 
   UDTable() { blockSize = 20; }
   UDTable(int bl) { blockSize = bl; }
};

