/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef SNP_H
#define SNP_H
#include "feudal/BinaryStreamTraits.h"

class snp {

     public:

     snp( ) { }

     snp( const int tig, const int pos, const char alt_allele )
          : tig_(tig), pos_(pos), alt_allele_(alt_allele) { }

     int Tig( ) const { return tig_; }
     int Pos( ) const { return pos_; }
     char AltAllele( ) const { return alt_allele_; }

     friend Bool operator==( const snp& s1, const snp& s2 )
     {    return s1.Tig( ) == s2.Tig( ) && s1.Pos( ) == s2.Pos( )
               && s1.AltAllele( ) == s2.AltAllele( );    }

     friend Bool operator<( const snp& s1, const snp& s2 )
     {    if ( s1.Tig( ) < s2.Tig( ) ) return True;
          if ( s1.Tig( ) > s2.Tig( ) ) return False;
          if ( s1.Pos( ) < s2.Pos( ) ) return True;
          if ( s1.Pos( ) > s2.Pos( ) ) return False;
          if ( s1.AltAllele( ) < s2.AltAllele( ) ) return True;
          return False;    }

     private:

     int tig_;          // zero-based index of fasta record in reference
     int pos_;          // zero-based index of base on fasta record
     char alt_allele_;  // alternate allele = 0, 1, 2, or 3 (for A, C, G, T)

};
TRIVIALLY_SERIALIZABLE(snp);

#endif
