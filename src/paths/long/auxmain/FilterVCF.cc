
///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * FilterVCF.cc
 *
 *  Created on: Sep 19, 2013
 *      Author: blau
 */

#include "paths/long/VariantPostProcess.h"
#include "MainTools.h"


int main(int argc, char *argv[]){
    RunTime();

    BeginCommandArguments;
    CommandArgument_String_Doc(IN_VCF, "Name of input vcf");
    CommandArgument_String_Doc(OUT_VCF, "Name of output vcf");
    EndCommandArguments;
    if(filter_vcf(IN_VCF,OUT_VCF)){
        std::cout << "There has been an error." << std::endl;
        return 1;
    }
    return 0;
};
