/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/// UnipathFasta.  Print fasta file of unipaths.

#include "Basevector.h"
#include "MainTools.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN); 
     CommandArgument_Int(K);
     EndCommandArguments;

     String predata = PRE + "/" + DATA;
     String run_dir = PRE + "/" + DATA + "/" + RUN;

     String KS = ToString(K);
     vecKmerPath unipaths( run_dir + "/reads.unipaths.k" + KS );
     vecKmerPath paths( run_dir + "/reads.paths.k" + KS );
     vecKmerPath paths_rc( run_dir + "/reads.paths_rc.k" + KS );
     BREAD2( run_dir + "/reads.pathsdb.k" + KS, vec<tagged_rpint>, pathsdb );
     KmerBaseBroker kbb( run_dir, K, paths, paths_rc, pathsdb );
     for ( size_t i = 0; i < unipaths.size( ); i++ )
          kbb.Seq( unipaths[i] ).Print( cout, "unipath_" + ToString(i) );    }
