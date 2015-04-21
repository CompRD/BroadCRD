/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "Bitvector.h"
#include "FeudalMimic.h"
#include "MainTools.h"

#include "FastIfstream.h"
#include "String.h"
#include "TokenizeString.h"


int main( int argc, char *argv[] )
{
    RunTime( );

    BeginCommandArguments;
    CommandArgument_String_Doc(GENOME, "GENOME.fastb will loaded to size the vecbitvector");
    CommandArgument_String_Doc(DBSNP, "file to use for marking.  column 1 should be contig, column 2 should be 1-based position");
    CommandArgument_String(OUT);
    EndCommandArguments;

    vecbasevector genome( GENOME + ".fastb" );
    vecbitvector dbsnp;
    Mimic( genome, dbsnp);
   
    fast_ifstream snp_in(DBSNP);
    String line;
    unsigned int line_number = 0;
    
    while ( 1 ) {
        line_number++;
        getline( snp_in, line);
        if ( snp_in.fail() ) break;
        vec<String> tokens;
        Tokenize(line, tokens);

        unsigned int contig = tokens[0].Int();
        unsigned int pos = tokens[1].Int() - 1;

        dbsnp[contig].Set(pos, True);
    }

    dbsnp.WriteAll(OUT);
    
    cout << "done. " << endl;
    cout << "processed " << line_number << " lines " << endl;
}
