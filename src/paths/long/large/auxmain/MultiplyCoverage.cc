///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Multiply coverage for a given sample by a given constant.  If you screw it up
// sufficiently badly, you can go back and recompute coverage from scratch.
// Keeps a record of edits in a.covs.edits.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include "MainTools.h"
#include "paths/long/large/Lines.h"

int main( int argc, char *argv[] ) 
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault_Doc(DIR, ".", "directory");
     CommandArgument_Int_Doc(S, "sample number, 0 or higher");
     CommandArgument_Double_Doc(M, "multiplier");
     EndCommandArguments;

     cout << Date( ) << ": loading coverage" << endl;
     vec<vec<covcount>> covs;
     BinaryReader::readFile( DIR + "/a.covs", &covs );
     ForceAssertLt( S, covs.isize( ) );
     cout << Date( ) << ": multiplying" << endl;
     #pragma omp parallel for
     for ( int e = 0; e < covs[S].isize( ); e++ )
          covs[S][e].Set( covs[S][e].Cov( ) * M );
     cout << Date( ) << ": writing" << endl;
     BinaryWriter::writeFile( DIR + "/a.covs", covs );
     Echo( Date( ) + ": coverage for sample " + ToString(S) + " multiplied by "
          + ToString(M), DIR + "/a.covs.edits" );
     cout << Date( ) << ": done" << endl;
     Scram(0);    }
