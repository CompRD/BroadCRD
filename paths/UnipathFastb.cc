/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/// Creates a unibases fastb file from a set of unipaths.
/// The name of the unipath file can be overridden with UNIPATHS_IN, the 
/// default name is reads.unipaths.kxx
/// The output fastb filename can be set in UNIBASES_OUT, the default name
/// is reads.unipaths.kxx.fastb


#include "Basevector.h"
#include "MainTools.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"
#include "paths/Unipath.h"


int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(READS, "reads");
     CommandArgument_String_OrDefault(UNIPATHS_IN, "reads.unipaths");
     CommandArgument_String_OrDefault(UNIBASES_OUT, "reads.unibases");
     CommandArgument_Int(K);
     EndCommandArguments;

     // Define directories.

     String data_dir = PRE + "/" + DATA;
     String run_dir = PRE + "/" + DATA + "/" + RUN;

     // Read in data.
     cout << Date( ) << ": loading data\n";
     String KS = ToString(K);
     vecbasevector bases( run_dir + "/" + READS + ".fastb" );
     vecKmerPath paths( run_dir + "/" + READS + ".paths.k" + KS );
     vecKmerPath paths_rc( run_dir + "/" + READS + ".paths_rc.k" + KS );
     BREAD2( run_dir + "/" + READS + ".pathsdb.k" + KS, vec<tagged_rpint>, pathsdb );
     vecKmerPath unipaths( run_dir + "/" + UNIPATHS_IN  + ".k" + KS );

     cout << Date( ) << ": making unibases\n";
     vecbasevector unibases;
     {    KmerBaseBroker kbb;
          kbb.Initialize( K, bases, paths, paths_rc, pathsdb );
          for ( size_t i = 0; i < unipaths.size( ); i++ )
               unibases.push_back( kbb.Seq( unipaths[i] ) );    }

     cout << Date( ) << ": writing unibases fastb\n";
     unibases.WriteAll(run_dir + "/" + UNIBASES_OUT + ".k" + KS + ".fastb");
 
}
