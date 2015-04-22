///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Run this three times with PERSON = child, mom, and pop.

// MakeDepend: dependency SearchHuman

#include "MainTools.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PERSON);
     EndCommandArguments;

     // The following to be run once:

     /*
     SystemSucceed( "LongProto SAMPLE=unknown READS="
          "/wga/scr4/human_data/NA12878_A2925/Solexa-127359_aligned.sorted.bam,"
          "/wga/scr4/human_data/NA12878_A2925/Solexa-127365_aligned.sorted.bam "
          "TMP=/wga/scr4/jaffe/fos_filter/tmp.fos EXIT=LOAD" );
     */

     String dir = "/wga/scr4/jaffe/fos_filter";

     cout << Date( ) << ": 0. finding reads" << endl;
     SystemSucceed( "SearchHuman MATCH="
          "/wga/scr4/jaffe/fos_filter/tmp.fos/frag_reads_orig.fastb BAMS=free "
          "GET_PARTNERS=True PERSON=" + PERSON + " PRINT_HEADERS=True NH=True "
          + "> " + dir + "/" + PERSON + ".sam" );

     cout << Date( ) << ": 1. getting headers" << endl;
     SystemSucceed( "cd " + dir + "; cat " + PERSON + ".sam | grep '^@' > "
          + PERSON + ".sam.header" );

     cout << Date( ) << ": 2. getting main" << endl;
     SystemSucceed( "cd " + dir + "; cat " + PERSON + ".sam | grep -v '^@' > "
          + PERSON + ".sam.main" );

     cout << Date( ) << ": 3. combining" << endl;
     SystemSucceed( "cd " + dir + "; cat " + PERSON + ".sam.header " + PERSON
          + ".sam.main > " + PERSON + ".sam.fixed" );

     cout << Date( ) << ": 4. bamming" << endl;
     SystemSucceed( "cd " + dir + "; samtools view -b -S " + PERSON + ".sam.fixed > "
          + PERSON + ".bam" );

     cout << Date( ) << ": 5. sorting" << endl;
     SystemSucceed( "cd " + dir + "; samtools sort " + PERSON + ".bam "
          + PERSON + ".sorted" );

     cout << Date( ) << ": 6. indexing" << endl;
     SystemSucceed( "cd " + dir + "; samtools index " + PERSON + ".sorted.bam" );

     cout << Date( ) << ": 7. done" << endl;    }
