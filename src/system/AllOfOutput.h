///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * AllOfOutput.h
 *
 *  Created on: May 20, 2013
 *      Author: tsharpe
 */

#ifndef SYSTEM_ALLOFOUTPUT_H_
#define SYSTEM_ALLOFOUTPUT_H_

#include "String.h"
#include "Vec.h"

vec<String> AllOfOutput(String const& command);

inline String AllOfOutput1(String const& command)
{    String all;
     vec<String> line = AllOfOutput(command);
     for ( size_t i = 0; i < line.size( ); i++ )
          all += line[i] + "\n";
     return all;    }



#endif /* SYSTEM_ALLOFOUTPUT_H_ */
