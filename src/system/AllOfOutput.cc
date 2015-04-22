///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * AllOfOutput.cc
 *
 *  Created on: May 20, 2013
 *      Author: tsharpe
 */

#include "system/AllOfOutput.h"

vec<String> AllOfOutput(String const& command) {
    procbuf pbuf(command.c_str(),std::ios_base::in);
    istream is(&pbuf);
    vec<String> lines;
    String line;
    while ( getline(is,line) )
        lines.push_back(line);
    return lines;
}
