// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#ifndef GENOME_COORD_MAP_H
#define GENOME_COORD_MAP_H

// A genome_coordinate_map defines a mapping from a reference genome R to a draft
// genome X, which lives as an Arachne assembly.

// To create a genome_coordinate_map, you need an Arachne-style AGP for X
// (created by CreateStandardOutputs), and alignments from R to X, in extended
// chain format.  Then use "Build" to build the genome_coordinate_map, and "Write"
// to write it to a file.  Then you can read it back in with "Read".

// The Lookup function may be used to translate R coordinates to X.  Let chr
// be an object in R's agp (e.g. a chromosome).  Let [pos1, pos2] be a closed
// interval on chr (using zero-based coordinates).  Then
//
//      Lookup( /* inputs */    chr, pos1, pos2,
//              /* outputs: */  OK, indent1, indent2, tig1, tig2, 
//                              tig1_pos, tig2_pos, tig1_rc, tig2_rc )
//
// attempts to map the interval to X.  It does so by moving forward from pos1 and
// backward from pos2 until it finds aligning bases (as specified in the original
// alignments from R to X).  The amount moved is specified by indent1 and indent2.
// The corresponding bases are on Arachne X contigs tig1 and tig2 (respectively),
// at positions tig1_pos and tig2_pos (respectively).  Orientations of the 
// correspondence are recorded via tig1_rc and tig2_rc (but the positions always
// refer to positions on the contig itself).
//
// If pos1 > pos2 or [pos1, pos2] lies in a gap or is off the end of the mapped
// region of chr, then Lookup fails, and indicates this by returning OK = False.
// Lookup will also fail if indent1 + indent2 > pos2 - pos1.  However, Lookup does
// not require that tig1 and tig2 lie on the same supercontig or satisfy any other
// consistency conditions.

#include "CoreTools.h"
#include "math/HoInterval.h"

// A balign is a gap-free alignment from bases range1 on sequence id1 to
// bases range2 on sequence id2.  It is reverse-complemented if rc != 0.

class balign {

     public:

     int id1, id2;
     ho_interval range1, range2;
     int rc;

     balign( ) { }
     balign( int id1_arg,
	     int id2_arg,
	     ho_interval range1_arg,
	     ho_interval range2_arg,
	     int rc_arg )
       : id1(id1_arg),
	 id2(id2_arg),
	 range1(range1_arg),
	 range2(range2_arg),
	 rc(rc_arg)
     { }

     int Start1( ) { return range1.Start( ); }
     int Stop1( ) { return range1.Stop( ); }
     int Start2( ) { return range2.Start( ); }
     int Stop2( ) { return range2.Stop( ); }

     friend Bool operator<( const balign& a1, const balign& a2 )
     {    if ( a1.id1 < a2.id1 ) return True;
          if ( a1.id1 > a2.id1 ) return False;
          return a1.range1.Start( ) < a2.range1.Start( );    }

};

class genome_coordinate_map {

     public:

     vec<balign> RtoX;
     vec<String> Rchr;

     genome_coordinate_map( ) { }

     // Valid: check to see if genome_coordinate_map has overlaps.

     Bool Valid( Bool verbose );

     void Build( const String& X_agp, const String& R_to_X );

     void Read( const String& fn );

     void Write( const String& fn );

     void Print( ostream& out, int i );

     // Lookup:
     // [pos1,pos2] is a closed interval on chr

     void Lookup( 
          /* inputs: */ const String& chr, int pos1, int pos2,
          /* outputs: */ Bool& OK, int& indent1, int& indent2,
                         int& tig1, int& tig2, int& tig1_pos, int& tig2_pos,
                         Bool& tig1_rc, Bool& tig2_rc );

     private:

     vec<int> index_;

};

#endif
