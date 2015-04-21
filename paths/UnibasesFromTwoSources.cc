/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// UnibasesFromTwoSources.

#include "Basevector.h"
#include "MainTools.h"
#include "feudal/BinaryStream.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/Unipath.h"

int main( int argc, char *argv[] ) 
{              
     RunTime( );
     
     BeginCommandArguments;
     CommandArgument_String(IN_HEAD1);
     CommandArgument_String(IN_HEAD2);
     CommandArgument_String(OUT_HEAD);
     CommandArgument_Int(K);
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
       "Number of threads to use (use all available processors if set to 0)");
     EndCommandArguments;

     // Thread control
   
     NUM_THREADS = configNumThreads(NUM_THREADS);

     // Load in the corrected reads
     vecbasevector bases( IN_HEAD1 + ".fastb" );
     bases.ReadAll( IN_HEAD2 + ".fastb" , true);

     // Find out which ones to keep
     vec<Bool> keep;
     BinaryReader::readFile( IN_HEAD1 + ".keep", &keep );
     keep.appendFromBinaryFile(IN_HEAD2 + ".keep");

     // Get rid of the bad reads
     vec<Bool> nkeep( bases.size( ) );
     for ( size_t i = 0; i < bases.size( ); i++ )
          nkeep[i] = !keep[i];
     bases.EraseIf(nkeep);

     // Path and Unipath
     vecKmerPath paths, pathsrc, unipaths;
     vec<tagged_rpint> pathsdb, unipathsdb;
     ReadsToPathsCoreY( bases, K, paths, pathsrc, pathsdb, 
                        OUT_HEAD + "/UnibasesFromTwoSources.bases", NUM_THREADS );
     Unipath( paths, pathsrc, pathsdb, unipaths, unipathsdb );
     KmerBaseBroker kbb( K, paths, pathsrc, pathsdb, bases );

     // Write out the unibase
     vecbasevector unibases;
     for ( size_t i = 0; i < unipaths.size( ); i++ )
          unibases.push_back_reserve( kbb.Seq( unipaths[i] ) );
     unibases.WriteAll( OUT_HEAD + ".unibases.k" + ToString(K) );    }
