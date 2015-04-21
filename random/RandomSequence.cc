/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/// Create a fasta file of random sequence with known ACGT content.
///
/// \file RandomSequence.cc
///
/// There is no T parameter because it = 1- (A+C+G).

#include "MainTools.h"
#include "Basevector.h"

int main( int argc, char *argv[] ) 
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_Abbr(OUT_PREFIX,O);
     CommandArgument_Double_OrDefault(A,0.25);
     CommandArgument_Double_OrDefault(C,0.25);
     CommandArgument_Double_OrDefault(G,0.25);
     CommandArgument_UnsignedInt_OrDefault(LENGTH,1000);
     CommandArgument_UnsignedInt_OrDefault(N,100);
     EndCommandArguments;

     srand48(time(NULL));

     basevector b(LENGTH);
     ForceAssert(A+C+G < 1.0);

     Ofstream(fasta, OUT_PREFIX + ".fasta");
     C = A + C;
     G = G + C;
     for (unsigned int i=0; i != N; ++i) {
       for (unsigned int j=0; j != LENGTH; ++j) {
	 double d = drand48();
	 if (d < A) b.Set(j,0);
	 else if (d < C) b.Set(j,1);
	 else if (d < G) b.Set(j,2);
	 else b.Set(j,3);
       }
       b.Print(fasta, i);
     }
     return 0;
}
	 
       
     
 

