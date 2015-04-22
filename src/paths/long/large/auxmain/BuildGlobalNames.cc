///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Build a read names index.  Hardcoded for now.

#include "MainTools.h"
#include "paths/long/large/ReadNameLookup.h"

int main( )
{
     RunTime( );

     // String head = "/wga/scr4/wg_projects/H.sapien/NA12878/frag_reads_orig";

     String head = "/wga/dev/jaffe/BroadCRD/tmp.xxx/frag_reads_orig";

     cout << Date( ) << ": loading names" << endl;
     vecString names( head + ".names" );
     readname_lookup look(names);
     cout << Date( ) << ": writing" << endl;
     BinaryWriter::writeFile( head + ".names.idx", look );
     cout << Date( ) << ": done" << endl;    }
