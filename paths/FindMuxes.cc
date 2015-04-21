/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// Finds the minimal extensions (muxes) to the left for all reads.

// An extension is a read that matches the given read perfectly and
// either extends beyond the given read to the left or is coterminal
// with (i.e., extends to the same point as) the given read but is
// shorter.

// The minimal extensions of a given read are the extensions of that
// read that are not an extension of any other extension of that read.

// It is possible for a read to have multiple minimal extensions if
// the genome diverges just to the left of the read.

#include "MainTools.h"
#include "paths/KmerPathDatabase.h"
#include "paths/OrientedKmerPathId.h"
#include "paths/MuxGraph.h"
#include "STLExtensions.h"
#include "TaskTimer.h"
#include "Histogram.h"
#include "feudal/BinaryStream.h"
#include "paths/MuxFinder.h"

int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
    CommandArgument_String(PRE);
    CommandArgument_String(DATA);
    CommandArgument_String(RUN);
    CommandArgument_UnsignedInt_OrDefault( K, 48 );
    CommandArgument_String_OrDefault( PATHSHQ, "pathshq" );
    CommandArgument_Bool_OrDefault( printHiMuxCount, False );
    CommandArgument_UnsignedInt_OrDefault( minNumKmers, 10 );
    CommandArgument_UnsignedInt_OrDefault( start, 0 );
    CommandArgument_UnsignedInt_OrDefault( stop, 0 );
    CommandArgument_UnsignedInt_OrDefault( chunkSize, 100000 );
    CommandArgument_Bool_OrDefault( WRITE, True );
    CommandArgument_Bool_OrDefault( CONFIRM_WRITE, False );
  EndCommandArguments;

  String runDir = PRE + "/" + DATA + "/" + RUN;

  cout << Date() << ": Loading paths..." << endl;

  vecKmerPath pathsFw, pathsRc;
  pathsFw.ReadAll( runDir + "/reads." + PATHSHQ + ".k" + ToString(K) );
  pathsRc.ReadAll( runDir + "/reads." + PATHSHQ + "_rc.k" + ToString(K) );

  cout << Date() << ": Done." << endl;

  cout << Date() << ": Loading paths database..." << endl;

  vec<tagged_rpint> rawPathsDB;
  BinaryReader::readFile( runDir + "/reads." + PATHSHQ + "db.k" + ToString(K),
                          &rawPathsDB );

  cout << Date() << ": Done." << endl;

  KmerPathDatabase pathsDb( &rawPathsDB );

  // Find minimal extensions.

  String logfile = runDir + "/FindMuxes." + PATHSHQ + ".log";
  ofstream logstrm( logfile.c_str() );

  MuxFinder theMuxFinder( pathsFw, pathsRc, pathsDb, &logstrm );
  theMuxFinder.SetVerbose( true );
  theMuxFinder.StudyReads();
  theMuxFinder.SetMinNumKmers( minNumKmers );
  theMuxFinder.SetPrintHiMuxCount( printHiMuxCount );

  MuxGraph theMuxGraph( pathsFw.size() );

  longlong numMuxes = theMuxFinder.FindMuxes( theMuxGraph, 
					      chunkSize,
					      start, 
					      ( stop>0 ? stop : pathsDb.size()-1 ) );
  
  cout << Date() << ": Found " << numMuxes << " minimal extensions." << endl;

  if ( WRITE )
  {
    const String& muxGraphFile( runDir + "/reads." + PATHSHQ + "_muxgraph.k" + ToString(K) );
    theMuxGraph.Write( muxGraphFile );

    if ( CONFIRM_WRITE )
    {
      MuxGraph graphOnDisk;
      graphOnDisk.Read( muxGraphFile );
      graphOnDisk.VerifySameAs( theMuxGraph );
    }
  }

  return 0;
}

