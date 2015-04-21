//////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Imports PacBio reads (the bare minimum)
//
// RUNS notation e.g.  10,12,13-16,19,24

#include "Basevector.h"
#include "FetchReads.h"
#include "MainTools.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PACBIO_DATA_DIR);
     CommandArgument_String_OrDefault(HEAD_OUT, "");
     CommandArgument_String(RUNS);
     EndCommandArguments;

     // Import PacBio Reads

     cout << "Importing PacBio reads from runs:" << endl;

     vecbasevector reads;
     vec<int> runs;

     Bool done(False);
     while( !done )
     {    String first;
          if ( RUNS.Contains( "," ) ) 
          {    first = RUNS.Before( "," );
               RUNS = RUNS.After( "," );    }
          else 
          {    first = RUNS;
               done = True;    }
          if ( !first.Contains( "-" ) ) runs.push_back( first.Int( ) );
          else
          {    int n1 = first.Before( "-" ).Int( ), n2 = first.After( "-" ).Int( );
               for ( int j = n1; j <= n2; j++ )
                    runs.push_back(j);    }    }

     for ( int i = 0; i < runs.isize( ); i++ ) {
       cout << runs[i] << endl;
       vecbasevector x;
       String prefix = "0" + ToString( runs[i] );
       prefix.resize(3);
       String dir = PACBIO_DATA_DIR + "/" + prefix
            + "/0" + ToString( runs[i] ) + "/data";
       String fn;
       if ( IsRegularFile( dir + "/filtered_subreads.fa" ) )
            fn = dir + "/filtered_subreads.fa";
       else fn = dir + "/filtered_subreads.fasta";

       FetchReads( x, 0, fn );
       reads.Append(x);    
     }

     reads.WriteAll(HEAD_OUT + ".fastb");
     
     cout << "Done" << endl;
}
