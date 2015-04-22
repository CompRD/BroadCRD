// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology


#ifndef AMINO_ACID_CONVERTER_H_
#define AMINO_ACID_CONVERTER_H_

// Simple function that converts kmers (must be a multiple of 3!)
// into amino acid codons, e.g. ATC will be translated into ATT
// since they both code ile.
// NOTE: there is only one table and it's for mammalians!

#include "Basevector.h"
#include "Vec.h"


struct AminoAcid
{
  char m_b1;
  char m_b2;
  char m_b3;
};


class AminoAcidTable
{
 public:
  AminoAcidTable() {
    m_table.resize(64);

    Set('A', 'T', 'T',
	'A', 'T', 'T');
    Set('A', 'T', 'C',
	'A', 'T', 'T');
    Set('A', 'T', 'A',
	'A', 'T', 'T');
    Set('C', 'T', 'T',
	'C', 'T', 'T');
    Set('C', 'T', 'C',
	'C', 'T', 'T');
    Set('C', 'T', 'A',
	'C', 'T', 'T');
    Set('C', 'T', 'G',
	'C', 'T', 'T');
    Set('T', 'T', 'A',
	'C', 'T', 'T');
    Set('T', 'T', 'G',
	'C', 'T', 'T');
    Set('G', 'T', 'T',
	'G', 'T', 'T');
    Set('G', 'T', 'C',
	'G', 'T', 'T');
    Set('G', 'T', 'A',
	'G', 'T', 'T');
    Set('G', 'T', 'G',
	'G', 'T', 'T');
    Set('T', 'T', 'T',
	'T', 'T', 'T');
    Set('T', 'T', 'C',
	'T', 'T', 'T');
    Set('A', 'T', 'G',
	'A', 'T', 'G');
    Set('T', 'G', 'T',
	'A', 'T', 'G');
    Set('T', 'G', 'C',
	'A', 'T', 'G');
    Set('G', 'C', 'T',
	'G', 'C', 'T');
    Set('G', 'C', 'C',
	'G', 'C', 'T');
    Set('G', 'C', 'A',
	'G', 'C', 'T');
    Set('G', 'C', 'G',
	'G', 'C', 'T');
    Set('G', 'G', 'T',
	'G', 'G', 'T');
    Set('G', 'G', 'C',
	'G', 'G', 'T');
    Set('G', 'G', 'A',
	'G', 'G', 'T');
    Set('G', 'G', 'G',
	'G', 'G', 'T');
    Set('C', 'C', 'T',
	'C', 'C', 'T');
    Set('C', 'C', 'C',
	'C', 'C', 'T');
    Set('C', 'C', 'A',
	'C', 'C', 'T');
    Set('C', 'C', 'G',
	'C', 'C', 'T');
    Set('A', 'C', 'T',
	'A', 'C', 'T');
    Set('A', 'C', 'C',
	'A', 'C', 'T');
    Set('A', 'C', 'A',
	'A', 'C', 'T');
    Set('A', 'C', 'G',
	'A', 'C', 'T');
    Set('T', 'C', 'T',
	'T', 'C', 'T');
    Set('T', 'C', 'C',
	'T', 'C', 'T');
    Set('T', 'C', 'A',
	'T', 'C', 'T');
    Set('T', 'C', 'G',
	'T', 'C', 'T');
    Set('A', 'G', 'T',
	'T', 'C', 'T');
    Set('A', 'G', 'C',
	'T', 'C', 'T');
    Set('T', 'A', 'T',
	'T', 'A', 'T');
    Set('T', 'A', 'C',
	'T', 'A', 'T');
    Set('T', 'G', 'G',
	'T', 'G', 'G');
    Set('C', 'A', 'A',
	'C', 'A', 'A');
    Set('C', 'A', 'G',
	'C', 'A', 'A');
    Set('A', 'A', 'T',
	'A', 'A', 'T');
    Set('A', 'A', 'C',
	'A', 'A', 'T');
    Set('C', 'A', 'T',
	'C', 'A', 'T');
    Set('C', 'A', 'C',
	'C', 'A', 'T');
    Set('G', 'A', 'A',
	'G', 'A', 'A');
    Set('G', 'A', 'G',
	'G', 'A', 'A');
    Set('G', 'A', 'T',
	'G', 'A', 'T');
    Set('G', 'A', 'C',
	'G', 'A', 'T');
    Set('A', 'A', 'A',
	'A', 'A', 'A');
    Set('A', 'A', 'G',
	'A', 'A', 'A');
    Set('C', 'G', 'T',
	'C', 'G', 'T');
    Set('C', 'G', 'C',
	'C', 'G', 'T');
    Set('C', 'G', 'A',
	'C', 'G', 'T');
    Set('C', 'G', 'G',
	'C', 'G', 'T');
    Set('A', 'G', 'A',
	'C', 'G', 'T');
    Set('A', 'G', 'G',
	'C', 'G', 'T');
    Set('T', 'A', 'A',
	'T', 'A', 'A');
    Set('T', 'A', 'G',
	'T', 'A', 'A');
    Set('T', 'G', 'A',
	'T', 'A', 'A');
    
  }
  
  const AminoAcid & Get(char b1, char b2, char b3) {
    int index = b1 + 4 * b2 + 16 * b3;
    return m_table[index];
  }

 private:
  void Set(char b1, char b2, char b3,
	   char c1, char c2, char c3) {
    int index = as_char(b1) + 4 * as_char(b2) + 16 * as_char(b3);
    AminoAcid & a = m_table[index];
    a.m_b1 = as_char(c1);
    a.m_b2 = as_char(c2);
    a.m_b3 = as_char(c3);
  }

  vec<AminoAcid> m_table;
};



inline bool ReduceToAminoAcids(basevector & b) 
{
  static AminoAcidTable table;

  for (unsigned int i=0; i<b.size(); i+= 3) {
    const AminoAcid & a = table.Get(b[i], b[i+1], b[i+2]);
    b.Set(i,   a.m_b1);
    b.Set(i+1, a.m_b2);
    b.Set(i+2, a.m_b3);
  }
  return true;

}



#endif //AMINO_ACID_CONVERTER_H_

