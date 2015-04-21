///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Sep 2012 
//



#ifndef FRAMEWORK_H_
#define FRAMEWORK_H_

#include "Basevector.h"
#include "Qualvector.h"
#include "PairsManager.h"
#include "ReadError.h"
#include "RefLocus.h"



////////////////////////////////////////////////////////////////////
// Base classes
////////////////////////////////////////////////////////////////////

class Generator {
public:
  virtual ~Generator() {};
  virtual const BaseVecVec& getReads() const = 0;
  virtual const QualVecVec& getQuals() const = 0;
  virtual const RefLocusVec& getReadLocs() const = 0;
  virtual const PairsManager& getPairsManager() const = 0;
};


class Corrector {
  virtual ~Corrector() {};
};

class Analyzer {
  virtual ~Analyzer() {};
};


////////////////////////////////////////////////////////////////////
// Generators not specified elsewhere...
////////////////////////////////////////////////////////////////////

class RealDataGenerator : public Generator {
public:
  RealDataGenerator( const String& bamfile /*, specification of a genome region */);
private:
  RealDataGenerator();	 // not implemented
};


class SavedDataGenerator : public Generator {
  SavedDataGenerator( const String& fastbfile, const String& qualbfile, const String& pairsfile /* , some metadata */);
private:
  SavedDataGenerator();  // not implemented
};


#endif /* FRAMEWORK_H_ */
