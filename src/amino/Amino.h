/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef AMINO_AMINO_H
#define AMINO_AMINO_H

class AminoAcid
{
public:
    /// The 1-letter name.
    char getSymbol() const { return mSymbol; }

    /// The 3-letter name.
    char const* getAbbreviation() const { return mAbbreviation; }

    /// The full name.
    char const* getName() const { return mName; }

    unsigned getIndex() const { return this-gTable; }
    static AminoAcid const& forIndex( unsigned idx ) { return gTable[idx]; }

    /// Iterators for the entire table of Amino Acids.
    static AminoAcid const* begin()
    { return gTable; }
    static AminoAcid const* end()
    { return gTable+sizeof(gTable)/sizeof(gTable[0]); }

    /// Get the AminoAcid associated with a codon.
    /// Codons are 6-bit numbers with the top two bits encoding the leading
    /// base (0=A,1=C,2=G,3=T/U, as usual), the middle two encoding the 2nd
    /// base, and the lowest two bits encoding the 3rd base.
    static AminoAcid const& forCodon( unsigned codon )
    { return gTable[gCodonToAA[codon]]; }

    /// Get the AminoAcid associated with the next three bases.
    template <class Itr> // Itr is an iterator that returns base codes.
    static AminoAcid const& forBases( Itr itr )
    { unsigned codon = *itr;
      codon = (codon << 2) | *++itr;
      codon = (codon << 2) | *++itr;
      return forCodon(codon); }

    /// Get the AminoAcid for a symbol (1-letter abbreviation).
    static AminoAcid const& forSymbol( char symbol )
    { return gTable[gSymbolToAA[symbol&0x1f]]; }

private:
    AminoAcid( char symbol, char const* abbreviation, char const* name )
    : mSymbol(symbol), mAbbreviation(abbreviation), mName(name) {}
    AminoAcid( AminoAcid const& ) = delete;

    char mSymbol;
    char const* mAbbreviation;
    char const* mName;
    static AminoAcid gTable[24];
    static unsigned gCodonToAA[64];
    static unsigned gSymbolToAA[32];
};

#endif
