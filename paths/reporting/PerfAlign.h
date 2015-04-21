///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PERF_ALIGN_H
#define PERF_ALIGN_H

// A perf_align tracks a perfect alignment between two sequences.  The coordinates
// are always forward on both.  The Rc flag is set if the correspondence if the
// two sequences are reverse complemented with respect to each other.  Note that it
// is not neccesary to specify which is reverse complemented.

class perf_align {

     public:

     perf_align( ) { }

     perf_align( const int id1, const int id2, const int pos1, const int Pos1,
          const int pos2, const int Pos2, const Bool rc ) : id1_(id1), id2_(id2),
          pos1_(pos1), Pos1_(Pos1), pos2_(pos2), Pos2_(Pos2), rc_(rc) { }

     int Id1( ) const { return id1_; }
     int Id2( ) const { return id2_; }
     int pos1( ) const { return pos1_; }
     int Pos1( ) const { return Pos1_; }
     int pos2( ) const { return pos2_; }
     int Pos2( ) const { return Pos2_; }
     Bool Rc( ) const { return rc_; }
     Bool Fw( ) const { return !rc_; }

     int& pos1( ) { return pos1_; }
     int& Pos1( ) { return Pos1_; }
     int& pos2( ) { return pos2_; }
     int& Pos2( ) { return Pos2_; }
     Bool& Rc( ) { return rc_; }

     ho_interval Extent1( ) const { return ho_interval( pos1_, Pos1_ ); }
     ho_interval Extent2( ) const { return ho_interval( pos2_, Pos2_ ); }

     int Len( ) const { return Pos1_ - pos1_; }

     int Offset( ) const { return pos1_ - pos2_; }

     friend int compare( perf_align const& a1, perf_align const& a2 )
     { int result = compare(a1.id1_,a2.id1_);
       if ( !result ) result = compare(a1.id2_,a2.id2_);
       if ( !result ) result = compare(a1.rc_,a2.rc_);
       if ( !result ) result = compare(a1.pos1_,a2.pos1_);
       //if ( !result ) result = compare(a1.pos2_,a2.pos2_);
       return result; }
     friend Bool operator<( const perf_align& a1, const perf_align& a2 )
     {    if ( a1.Id1( ) < a2.Id1( ) ) return True;
          if ( a1.Id1( ) > a2.Id1( ) ) return False;
          if ( a1.Id2( ) < a2.Id2( ) ) return True;
          if ( a1.Id2( ) > a2.Id2( ) ) return False;
          if ( a1.Rc( ) < a2.Rc( ) ) return True;
          if ( a1.Rc( ) > a2.Rc( ) ) return False;
          if ( a1.pos1( ) < a2.pos1( ) ) return True;
          // if ( a1.pos2( ) < a2.pos2( ) ) return True;
          return False;    }

     friend ostream& operator<<( ostream& out, const perf_align& a )
     {    return out << a.Id1( ) << ":" << a.pos1( ) << "-" << a.Pos1( )
               << " --> " << a.Id2( ) << ":" << a.pos2( ) << "-" << a.Pos2( )
               << " (" << ( a.Fw( ) ? "+" : "-" ) << ")";    }

     private:
     
     int id1_, id2_;
     int pos1_, Pos1_;
     int pos2_, Pos2_;
     Bool rc_;

};

inline char Comp( const char c )
{    if ( c == 'A' ) return 'T';
     if ( c == 'C' ) return 'G';
     if ( c == 'G' ) return 'C';
     if ( c == 'T' ) return 'A';
     if ( c == 'N' ) return 'N';
     ForceAssert( 0 == 1 );
     return 0;    }

inline Bool Comparable( const perf_align& a1, const perf_align& a2 )
{    return a1.Id1( ) == a2.Id1( ) && a1.Id2( ) == a2.Id2( )
          && a1.Rc( ) == a2.Rc( );    }

#endif
