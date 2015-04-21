/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// This defines the guts of TruePoly.  See TruePoly.cc for what it does.

#ifndef TRUE_POLY_CORE_H
#define TRUE_POLY_CORE_H

#include "Basevector.h"
#include "CoreTools.h"
#include "FastIfstream.h"
#include "Bitvector.h"
#include "feudal/BinaryStream.h"
#include <limits.h>

class poly {

     public:

     int tig;
     int ref_start;
     String left, sample, ref, right;
     int orig_ref_start;
     int orig_sample_start;
     int sample_tig;
     int orig_sample_len;
     int orig_ref_len;

     poly( ) { }
     poly( const int tig, const int ref_start, const String& left,
          const String& sample, const String& ref, const String& right )
          : tig(tig), ref_start(ref_start), left(left), sample(sample),
          ref(ref), right(right) { }
     poly( const int tig, const int ref_start, const String& left,
          const String& sample, const String& ref, const String& right,
          const int orig_ref_start, const int orig_sample_start, 
          const int sample_tig, const int orig_sample_len, const int orig_ref_len )
          : tig(tig), ref_start(ref_start), left(left), sample(sample),
          ref(ref), right(right), orig_ref_start(orig_ref_start),
          orig_sample_start(orig_sample_start), sample_tig(sample_tig),
          orig_sample_len(orig_sample_len), orig_ref_len(orig_ref_len) { }

     friend Bool operator<( const poly& p1, const poly& p2 )
     {    if ( p1.tig < p2.tig ) return True;
          if ( p1.tig > p2.tig ) return False;
          if ( p1.ref_start < p2.ref_start ) return True;
          return False;    }

     friend ostream& operator<<( ostream& out, const poly& p )
     {    return out << p.tig << " " << p.ref_start << " left=" << p.left 
               << " sample=" << p.sample << " ref=" << p.ref << " right=" << p.right 
               << "\n";    }

     int Class( ) const
     {    if ( sample.size( ) == 1 && ref.size( ) == 1 ) return 1;
          if ( sample.size( ) == 0 ) return 2;
          return 3;    }

     void writeBinary( BinaryWriter& bw ) const
     {
         bw.write(tig);
         bw.write(ref_start);
         bw.write(orig_ref_start);
         bw.write(orig_sample_start);
         bw.write(sample_tig);
         bw.write(orig_sample_len);
         bw.write(orig_ref_len);
     }

     void readBinary( BinaryReader& br )
     {
         br.read(&tig);
         br.read(&ref_start);
         br.read(&orig_ref_start);
         br.read(&orig_sample_start);
         br.read(&sample_tig);
         br.read(&orig_sample_len);
         br.read(&orig_ref_len);
         left.clear();
         sample.clear();
         ref.clear();
         right.clear();
     }

     static size_t externalSizeof()
     { return sizeof(tig)+sizeof(ref_start)+sizeof(orig_ref_start) +
             sizeof(orig_sample_start) + sizeof(sample_tig) +
             sizeof(orig_sample_len); }
};
SELF_SERIALIZABLE(poly);

void SummarizeClasses( ostream& out, const vec<poly>& polys );

void DumpHeader( ostream& out, const String& REFB );

void TruePolyCore( vec<poly>& polys, vecbasevector refA, 
     const vecbasevector& refB, const String& TRUSTEDA = "", 
     const String& TRUSTEDB = "", const int K = 20, 
     const String& CLASSES_TO_PRINT = "{1,2,3}", 
     vec< vec<poly> > *compound_polys = NULL, vecbitvector *untrusted_bases = NULL,
          vecbitvector* nonpolymorphic_bases = NULL );

void TruePolyCore( vec<poly>& polys, const String& REFA, const String& REFB,
     const String& TRUSTEDA = "", const String& TRUSTEDB = "", const int K = 20, 
     const String& CLASSES_TO_PRINT = "{1,2,3}",
     vec< vec<poly> > *compound_polys = NULL, vecbitvector *untrusted_bases = NULL,
          vecbitvector* nonpolymorphic_bases = NULL );

inline void TruePolyCore( ostream& out, const String& REFA, const String& REFB,
     const String& TRUSTEDA = "", const String& TRUSTEDB = "", const int K = 20, 
     const String& CLASSES_TO_PRINT = "{1,2,3}", 
     const String& binary_outfile = "",
     const String& COMPOUND_FILE = "", vecbitvector* untrusted_bases = NULL,
          vecbitvector* nonpolymorphic_bases = NULL )
{    DumpHeader( out, REFB );
     vec<poly> polys;
     vec< vec<poly> > compound_polys;
     if (COMPOUND_FILE == "") 
        TruePolyCore( polys, REFA, REFB, TRUSTEDA, TRUSTEDB, K, CLASSES_TO_PRINT, 
             NULL, NULL, nonpolymorphic_bases );
     else {
        vecbitvector untrusted_bases;
        TruePolyCore( polys, REFA, REFB, TRUSTEDA, TRUSTEDB, K, CLASSES_TO_PRINT, 
             &compound_polys, &untrusted_bases, nonpolymorphic_bases );
        String untrusted_filename = "";
        if ( COMPOUND_FILE.Contains(".compound") )
            untrusted_filename = COMPOUND_FILE.Before(".compound");
        else 
            untrusted_filename = COMPOUND_FILE;
        untrusted_bases.WriteAll(untrusted_filename + ".untrusted_bases.vecbitvector");
     }
     for ( int i = 0; i < polys.isize( ); i++ )
          out << polys[i];
     SummarizeClasses( out, polys );
     if ( binary_outfile != "" )
         BinaryWriter::writeFile( binary_outfile, polys );

     if (COMPOUND_FILE != "") {
        Ofstream(cpout, COMPOUND_FILE);
        for (unsigned int ii = 0; ii < compound_polys.size( ); ++ii) {
            vec<poly>& temp_compound = compound_polys[ii];
            cpout << "Group " << ii << " : " << temp_compound[0];
            for (unsigned int jj = 1; jj < temp_compound.size( ); ++jj) {
                cpout << temp_compound[jj];
            }
            cpout << endl << endl;
        }
        cpout.close( );
     }
}



#endif
