// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology
// 

#include "Alignment.h"
#include "Basevector.h"
#include "String.h"
#include "Vec.h"
#include "PrintAlignment.h"
#include "lookup/LookAlign.h"
#include "lookup/PrintLookAlign.h"

void PrintLookAlignPlus( ostream &out,
			 const vecbasevector &target_bases,
			 const vecbasevector &query_bases,
			 const look_align_plus &look_alplus,
			 bool abbreviate )
{
  int t_id = look_alplus.target_id;
  int q_id = look_alplus.query_id;
  bool RC = look_alplus.rc1;
  const align &al = look_alplus.a;
  const basevector &t_bases = target_bases[t_id];
  basevector qbases_rc;
  if ( RC ) {
    qbases_rc = query_bases[q_id];
    qbases_rc.ReverseComplement( );
  }
  const basevector &q_bases = RC ? qbases_rc : query_bases[q_id];

  int err = ActualErrors( q_bases, t_bases, al, 1, 1 );
  int al_len = al.Pos1( ) - al.pos1( );
  float err_rate = ( al_len < 1 ) ? 1.0 : float( err ) / float( al_len );

  int len1 = q_bases.size( );
  int len2 = t_bases.size( );
  int tail = Min( len1 - al.Pos1( ), len2 - al.Pos2( ) );
  int head = Min( al.pos1( ), al.pos2( ) );
  String str_rc = RC ? "rc" : "";
  String str_head = ( head > 0 ) ? " head=" + ToString( head ) : "";
  String str_tail = ( tail > 0 ) ? " tail=" + ToString( tail ) : "";
  
  out << "ALIGN q" << q_id << str_rc
      << " [" << al.pos1( )
      << ", " << al.Pos1( )
      << ")_" << len1
      << " vs t" << t_id
      << " [" << al.pos2( )
      << ", " << al.Pos2( )
      << ")_" << len2
      << " err_rate=" << ToString( err_rate, 2 )
      << "%" << str_head << str_tail
      << "\n";

  PrintVisualAlignment( abbreviate, out, q_bases, t_bases, al );
  
}

