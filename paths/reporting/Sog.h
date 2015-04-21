///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// A sog is a "scaffold on genome".  It represents a scaffold, in alignment to a 
// reference sequence via a set of perfect alignments.

#ifndef SOG_H
#define SOG_H

#include "Basevector.h"
#include "Bitvector.h"
#include "CoreTools.h"
#include "math/HoInterval.h"
#include "paths/reporting/PerfAlign.h"

class sog {

     public:

     sog( ) { }
     sog( const basevector& s, const bitvector& g, const vec<perf_align>& aligns );

     void SetGenome( vecbasevector* genome, vecbitvector* genome_amb ) 
     {    genome_ = genome;
          genome_amb_ = genome_amb;    }

     char GenomeBase( const int t, const int p ) const
     {    if ( (*genome_amb_)[t][p] ) return 'N';
          else return as_base( (*genome_)[t][p] );    }

     int AlignsCount( ) const { return aligns_.size( ); }
     const perf_align& Align( const int i ) const { return aligns_[i]; }
     perf_align& Align( const int i ) { return aligns_[i]; }
     vec<perf_align>& Aligns( ) { return aligns_; }

     const vec<char>& Scaff( ) const { return s_; }
     vec<char>& Scaff( ) { return s_; }

     char Base( const int i ) const { return s_[i]; }
     char CompBase( const int i ) const { return Comp( Base(i) ); }
     void SetBase( const int i, const char c ) { s_[i] = c; }
     int BaseCount( ) const { return s_.size( ); }

     Bool Gap( const int i ) const { return s_[i] == 'N'; }

     void Print( const int id ) const;

     // PrintFasta.  Print a range of bases, folding lines.  If title present,
     // print fasta header line.

     void PrintFasta( ostream& out, const int start, const int stop, 
          const String title = "" ) const;

     vec<ho_interval> Tigs( ) const;

     void Reverse( );

     void Validate( ) const;

     // Replace.  Replace a given range of bases in the scaffold by an alternate
     // sequence, and adjust the coordinates of the alignments accordingly.  We
     // assume that two adjacent alignments are to be replaced by one.

     void Replace( const int start, const int stop, const vec<char>& replacement,
          const int first_align_index, const perf_align& new_align );

     Bool ReplaceOK( const int start, const int stop, const int first_align_index );

     private:

     vec<char> s_;                      // scaffold
     vec<perf_align> aligns_;           // alignments of scaffolds to reference
     static vecbasevector* genome_;     // pointer to the genome
     static vecbitvector* genome_amb_;  // pointer to the genome

};

#endif
