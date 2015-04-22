///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Be sure to run PreGlumper first.

// MakeDepend: dependency LongProto
// MakeDepend: dependency ReadQGrapher
// MakeDepend: dependency PopRandom
// MakeDepend: dependency FilterIn
// MakeDepend: dependency Morgle
// MakeDepend: dependency FosClean

#include "MainTools.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PERSON);
     CommandArgument_String_OrDefault(REGION, "");
     CommandArgument_String_OrDefault(N, "-1");
     CommandArgument_String_OrDefault(FOS, "");
     CommandArgument_String(INSTANCE);
     CommandArgument_String(EVAL);
     EndCommandArguments;

     SystemSucceed( "PopRandom PERSON=\"" + PERSON + "\" N=" + N 
          + " REGION=" + REGION + " FOS=" + FOS + " INSTANCE=" + INSTANCE );

     SystemSucceed( "FilterIn INSTANCE=" + INSTANCE );

     SystemSucceed( "LongProto HAVE_PICARD=True "
          "SAMPLE=unknown READS=/wga/scr4/jaffe/fos_filter/sub." + INSTANCE 
          + ".bam TMP=/wga/scr4/jaffe/Glumper/tmp.xxx." + INSTANCE 
          + " OUT_INT_HEAD= ~/crd/sss." + INSTANCE
          + " HEURISTICS=\"{POP_BUBBLES=False,K2_FORCE=100,DELETE_SATELLITE=True,"
          "SATELLITE_TARGETS={alpha,two,ebv},MAX_ALPHA_SCORE=25,"
          "REQUIRE_EDGE_MATCH=/wga/scr4/jaffe/fos_filter/tmp.fos/"
          "frag_reads_orig.fastb}\"" );

     SystemSucceed( "FosClean INSTANCE=" + INSTANCE + " N=\"" + EVAL + "\"" );

     SystemSucceed( "Morgle INSTANCE=" + INSTANCE + " N=\"" + EVAL + "\"" );    }
