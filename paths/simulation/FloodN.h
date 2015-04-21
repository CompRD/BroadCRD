///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
#include "MainTools.h"

#include <algorithm>
#include <vector>
#include <fstream>
#include <string>
#include <unordered_map>
*/

#include<string>
#include<unordered_map>
#include<unordered_map>
#include<vector>

#include <sys/types.h>

#include "paths/simulation/Regions.h"

#ifndef FLOODN_H
#define FLOODN_H

namespace FloodN{

//line-by-line filtering of the FASTA file according to the ranges specified by "records"
//and then write the result to a file
void flooded_FASTA(const Regions::region_records& records, const std::string& sOutName ,const std::string& sInName);


};

#endif

