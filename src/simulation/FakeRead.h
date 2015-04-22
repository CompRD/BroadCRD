// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology


#ifndef _FAKEREAD_H_
#define _FAKEREAD_H_


#include "Vec.h"
#include "String.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "Qualvector.h"

class CFakeRead
{
public:
  CFakeRead() {
    m_size = 0;
    m_bPaired = 0;
  }
  ~CFakeRead() {}

  void Set(const basevector & bases,
	   const qualvector & quals,
	   int from, int to,
	   bool r,
	   int size,
	   const String & lib,
	   const String & plate,
	   const String & name,
	   const String & tmp,
	   const String & well = "",
	   bool bPaired = true);
 
  void Print( ostream& xmls, ostream& bases, ostream& quals ) const;

private:
  basevector m_bases;
  qualvector m_quals;
  String m_fr;
  String m_lib;
  String m_plate;
  String m_name;
  String m_template;
  String m_well;
  int m_size;
  bool m_bPaired;
};

void PrintFakeReads( vec<CFakeRead> const& fakeReads, String const& xmlFile,
                            String const& baseFile, String const& qualFile );

#endif


