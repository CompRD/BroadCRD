///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "TokenizeString.h"
#include "math/Functions.h"
#include "paths/long/large/ReadNameLookup.h"

int main( int argc, char *argv[] ) 
{
     RunTime( );
     
     BeginCommandArguments;
     CommandArgument_Int_Doc(ID, "local read id");
     CommandArgument_String_OrDefault(SAMPLE, "human");
     EndCommandArguments;

     vecString names;
     names.ReadOne( 
          // "tmp.xxx/frag_reads_orig.names", 
          "/wga/dev/jaffe/BroadCRD/tmp.yyy/frag_reads_orig.names", 
          ID );
     cout << Date( ) << ": loading look" << endl;
     readname_lookup look;
     String global;
     if ( SAMPLE == "human" ) global = "/wga/scr4/wg_projects/H.sapien/NA12878";
     else if ( SAMPLE == "rhody" ) global = "/wga/scr4/jaffe/GapToy/bugs/rhody/data";
     else 
     {    cout << "Unknown sample." << endl;
          Scram(1);    }
     BinaryReader::readFile( global + "/frag_reads_orig.names.idx", &look );
     cout << Date( ) << ": done" << endl;
     int64_t id = look.GetReadId( names[0] );
     cout << id << endl;    }
