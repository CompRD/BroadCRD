/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// Turns a vecbitvector file into a human-readable text...


#include <ctype.h>

#include "FastIfstream.h"
#include "system/ParsedArgs.h"
#include "system/RunTime.h"
#include "String.h"
#include "system/System.h"
#include "Bitvector.h"

int main( int argc, char *argv[] )
{
     RunTime();

     BeginCommandArguments;
     CommandArgument_StringSet_Doc(IN, "We can bring in multiple vecbitvectors and merge, but the must be all the same size.");
     CommandArgument_String_OrDefault_Doc(MERGEOP, "OR", "In merging the vectors, what kind of operation does one do.  Currently have OR (default), AND, and XOR.");
     CommandArgument_String_OrDefault_Doc(OUT, "", "By default goes to cout, but can specify a file output instead.");
     // Stupid to do PRINT = 2, since that prints all coordinates 
     // in all contigs.
     CommandArgument_Int_OrDefault_Doc(PRINT, 1, "Prints contig and position.  0 prints the 0-bit locations, 1 the 1-bit location." );
     CommandArgument_Bool_OrDefault_Doc(ZERO, True, "If true, prints 0-based coords, else 1-based.");
     CommandArgument_Bool_OrDefault_Doc(PRINTASIS, False, "If true, prints the vecbitvector with a line of zeros and ones per contig.  Otherwise print coordinates.");
     CommandArgument_Bool_OrDefault_Doc(JUST_STATS, False, "If true, just print stats (% of positions with a 1)\n");
     CommandArgument_Int_OrDefault_Doc(COLBREAK, 0, "If > 0 and if PRINTASIS is True, then insert new lines at this many columns.");
     EndCommandArguments;

     Assert(PRINT == 0 || PRINT == 1);

     AssertGt(IN.size( ),0u);

     vecbitvector theData(IN[0]);

     if (JUST_STATS)
     {
        long num_1s = 0;
        long total = 0;
        for (size_t i = 0; i < theData.size(); i++)
        {
            for (unsigned int j = 0; j < theData[i].size(); j++)
            {
                if (theData[i][j] == 1) { num_1s += 1; }
                total += 1; 
            }
        }

        printf("Total positions : %ld\n", total);
        printf("Number of 1's   : %ld (%0.02f%%)\n", num_1s, 100*((float)num_1s/(float)total));
        fflush(stdout);
        exit(0);
     }

     int optype;
     if (MERGEOP == "OR") optype = 0;  
     else if (MERGEOP == "AND") optype = 1;
     else if (MERGEOP == "XOR") optype = 2;
     else {
        cout << MERGEOP << " is an unknown MERGEOP.\n";
        cout << "Currently choices are OR, AND, and XOR.\n";
        exit(1);
     }

     for (unsigned int i = 1; i < IN.size( ); ++i) {
        vecbitvector secondData(IN[i]);
        Assert(SameSizes(theData, secondData));
        switch(optype) {
           case 0:
               theData |= secondData;
               break;
           case 1:
               theData &= secondData;
               break;
           case 2:
               theData &= secondData;
               break;
        }
     }

     ofstream fout; 
     if (OUT != "")
        OpenOfstream( fout, "fout", OUT ); 

     unsigned posPrint = 0;
     for (size_t contig = 0; contig < theData.size( ); ++contig) {
        for (unsigned int pos = 0; pos < theData[contig].size( ); ++pos) {
           if (!PRINTASIS) {
              if (ZERO) posPrint = pos; else posPrint = pos+1;
              if ( (theData[contig][pos] && PRINT == 1) ||
                   (!(theData[contig][pos]) && PRINT == 0) )
                if (OUT != "")
                   fout << contig << "    " << posPrint << endl;
                else
                   cout << contig << "    " << posPrint << endl;
           } else {
                int result = 0;
                if (theData[contig][pos]) result = 1;
                if (OUT != "") {
                   fout << result;
                   if ((COLBREAK > 0) && ((pos+1) % COLBREAK == 0)) fout << endl;
                } else {
                   cout << result;
                   if ((COLBREAK > 0) && ((pos+1) % COLBREAK == 0)) cout << endl;
                }
           }
        }
        if (PRINTASIS) {
           if (OUT != "")
              fout << endl;
           else
              cout << endl;
        }
     }
     if (OUT != "") fout.close( );

}
