///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * AlignToSHBV.cc
 *
 *  Created on: Sep 3, 2013
 *      Author: blau
 */
// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
#include <omp.h>
#include <iomanip>
#include "MainTools.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "paths/long/EvalByReads.h"
#include "paths/long/LocTraceEdges.h"
#include <queue>


int main(int argc, char *argv[])
{
    RunTime( );

    BeginCommandArguments;
    CommandDoc("Performs simple gap-free alignment for each sequence in QUERY_HEAD.{fastb,qualb} against the supported hypberbase vector graph stored in SHBV. Report the placement(s) with the minimal sum of q-score at mismatching bases.");
    CommandArgument_String_Doc(QUERY_HEAD, "specify the input sequence files QUERY_HEAD.fastb and QUERY_HEAD.qualb");
    CommandArgument_String_OrDefault_Doc(SHBV, "", "binary file storing the SupportedHyperBaseVector");
    CommandArgument_String_OrDefault_Doc(HBV, "", "binary file storing the HyperBaseVector");
    CommandArgument_Int_OrDefault_Doc(L, 12, "flanking region used for hashing");
    CommandArgument_Double_OrDefault_Doc(DELTA_QSUM, 0.0, "also report placements with a q-sum within DELTA_QSUM of best placements");
    EndCommandArguments;

    HyperBasevector hbv;
    if ( SHBV != "" ) {
	SupportedHyperBasevector shbv;
	BinaryReader::readFile( SHBV, &shbv );
	hbv = shbv;
    } else if ( HBV != "" ) {
	BinaryReader::readFile( HBV, &hbv );
    } else
	FatalErr("must specify either HBV or SHBV");
    vecbasevector bases( QUERY_HEAD+".fastb");
    vecqualvector quals( QUERY_HEAD+".qualb");

    ForceAssert( L > 0);
    LocTraceEdges(hbv,bases,quals,DELTA_QSUM,L,true);
}
