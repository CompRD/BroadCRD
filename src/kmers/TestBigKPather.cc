///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "kmers/BigKPather.h"

int main( int argc, char* argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Doc(FASTB,"reads to path");
     CommandArgument_String(OUT_HEAD);
     CommandArgument_UnsignedInt(K);
     CommandArgument_UnsignedInt(COVERAGE);
     EndCommandArguments;

     vecbasevector edges(FASTB);
     HyperBasevector hbv;
     ReadPathVec rpv;
     buildBigKHBVFromReads(K,edges,COVERAGE,&hbv,&rpv);

     BinaryWriter::writeFile(OUT_HEAD+".hbv", hbv);
     Ofstream(out_edge, OUT_HEAD+".fasta");
     for ( int i = 0; i < hbv.EdgeObjectCount(); ++i ) {
          hbv.EdgeObject(i).PrintCol(out_edge,ToString(i),80);
     }
     Ofstream(out_paths, OUT_HEAD+".paths");
     for ( size_t i = 0; i < rpv.size(); ++i )
          out_paths << i << ": " << rpv[i] << endl;

     return 0;
}
