/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef SOMATIC_CALL_TOOLS_H
#define SOMATIC_CALL_TOOLS_H

#include "Basevector.h"
#include "CoreTools.h"
#include "FastIfstream.h"

Bool SomaticMutation( 

     // inputs:

     const int chr, const int pos, const vecbasevector& ref,
     const vec<unsigned char>* TQ, const vec<unsigned char>* NQ, 
     const vec<unsigned char>* TQ_rc, const vec<unsigned char>* NQ_rc,
     const unsigned char refbase, 

     // heuristic constants:

     const int min_mutant_sum_pretest,
     const int min_mutant_sum,
     const double tumor_threshold,
     const double normal_threshold,

     // input and output (if != 4, assume altbase is given):

     unsigned char& altbase,

     // outputs (apart from return value):

     String& info

     );

inline void WriteStrings( const String& fn, const vec<String>& S )
{    Ofstream( out, fn );
     for ( int i = 0; i < S.isize( ); i++ )
          out << S[i] << "\n";    }

inline void ReadStrings( const String& fn, vec<String>& S )
{    String line;
     Ifstream( in, fn );
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          S.push_back(line);    }    }

#endif
