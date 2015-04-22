/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef NOBBIT
#define NOBBIT

#include <fstream>
#include "PackAlign.h"
#include "system/System.h"

// Comments for class nobbit:
//
// On sequence id1, the match goes from pos1 to Pos1;
// on sequence id2, the match goes from pos2 to Pos2.

class nobbit {    

     public:

     align a;
     int RC;
     int length1, length2;
     int id1, id2;
     int pos1, pos2;
     float score;
     int Pos1, Pos2;

     nobbit( const align& a_arg,
	     int RC_arg,
	     int length1_arg,
	     int length2_arg,
	     int id1_arg,
	     int id2_arg,
	     int pos1_arg,
	     int pos2_arg,
	     float score_arg,
	     int Pos1_arg,
	     int Pos2_arg )
       : a(a_arg), 
	 RC(RC_arg),
	 length1(length1_arg),
	 length2(length2_arg),
	 id1(id1_arg),
	 id2(id2_arg),
	 pos1(pos1_arg),
	 pos2(pos2_arg),
	 score(score_arg),
	 Pos1(Pos1_arg),
	 Pos2(Pos2_arg)
     { }

     nobbit( ) { }

     void Kill() { id1 = -1; }

     Bool IsDead() const { return ( id1 == -1 ); }

     friend ostream& operator<<( ostream& out, const nobbit& n );
     friend istream& operator>>( istream& in, nobbit& n );

     // Print in a human readable form.
     void Print( ostream& out );

};

/**
 * order_nobbit_Contig_Known
 * Order nobbits by contig_id, known_id, begin_on_known
 */
struct order_nobbit_Contig_Known :
  public binary_function<const nobbit&, const nobbit&, bool>
{
  bool operator( ) ( const nobbit &left, const nobbit &right ) {
    if ( left.id1 == right.id1 ) {
      if ( left.id2 == right.id2 ) return ( left.pos2 < right.pos2 );
      else return ( left.id2 < right.id2 );
    }
    return ( left.id1 < right.id1 );
  }
};

/**
 * order_nobbit_Known_Arachne
 * Order nobbits by known contig id, Arachne contig id.
 */
struct order_nobbit_Known_Arachne :
  public binary_function<const nobbit&, const nobbit&, bool>
{
  bool operator() ( const nobbit &left, const nobbit &right ) {
    if ( left.id2 == right.id2 )
      return ( left.id1 < right.id1 );

    return ( left.id2 < right.id2 );
  }
};

/*
 * order_nobbit_Known_StartOnKnown
 * Order nobbits by known contig id, start position on known contig.
 */
struct order_nobbit_Known_StartOnKnown :
  public binary_function<const nobbit&, const nobbit&, bool>
{
  bool operator() ( const nobbit &left, const nobbit &right ) {
    if ( left.id2 == right.id2 )
      return ( left.pos2 < right.pos2 );
    
    return ( left.id2 < right.id2 );
  }
};

#endif
