///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Sep 27, 2012
//
// This is a class-ification of ReadSimulatorSimple.cc that works with Framework.h
//

#ifndef READSIMULATORSIMPLEGENERATOR_H_
#define READSIMULATORSIMPLEGENERATOR_H_

#include "simulation/Framework.h"
#include "String.h"
#include "bam/CachedBAMFile.h"
#include <vector>


////////////////////////////////////////////////////////////////////
// BamReadGenerator - extract reads from a .BAM file by wrapping
// ExtractFromBAM.pl
////////////////////////////////////////////////////////////////////

class BamReadGenerator : public Generator
{
public:
  // generate a read set from a bam file.  if you want to, you can specify
  // regions of interest as a string in the format that "samtools view" expects,
  // e.g., "chr1:1-2000000".
  explicit BamReadGenerator( CachedBAMFile const& bam,
                              String const& regions = "" );

  // generate a read set from some bam files.  region arg is as above.
  explicit BamReadGenerator( std::vector<CachedBAMFile> const& bams,
                              String const& regions = "" );

  // compiler-supplied copying and destructor is OK

  // output getters

  BaseVecVec const& getReads() const { return _bvv; }
  QualVecVec const& getQuals() const { return _qvv; }
  RefLocusVec const& getReadLocs() const { return _refv; }
  PairsManager const& getPairsManager() const { return _pm; }

private:
  void init();
  void processBAM( CachedBAMFile const& bam, String const& regions );

  BaseVecVec _bvv;
  QualVecVec _qvv;
  RefLocusVec _refv;
  PairsManager _pm;
};


#endif /* READSIMULATORSIMPLEGENERATOR_H_ */
