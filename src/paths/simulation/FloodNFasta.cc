///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#include "MainTools.h"

#include "paths/simulation/FloodN.h"
#include "paths/simulation/Regions.h"

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;

     CommandArgument_String_Doc(OUTPUT_FILE, "output file name");
     CommandArgument_String_Doc(REGION_FILE, "*.regions file name");
     CommandArgument_String_Doc(FASTA_FILE, "*.fasta file name");
     CommandArgument_Int_OrDefault_Doc(PAD_LENGTH,0,"length of padding region");
     CommandArgument_Int_OrDefault_Doc(SKIPPED_REGION_COLUMNS,1,"length of padding region");

     std::cout << "Trimming:       " << FASTA_FILE << std::endl;
     std::cout << "According to:   " << REGION_FILE << std::endl;
     std::cout << "in which first  " << SKIPPED_REGION_COLUMNS << " columns are omitted"<<std::endl;
     std::cout << "With pad width: " << PAD_LENGTH << std::endl;
     std::cout << "Writing to:     " << OUTPUT_FILE << std::endl;

     EndCommandArguments;

     Regions::region_records records(REGION_FILE
                                   ,PAD_LENGTH
                                   ,SKIPPED_REGION_COLUMNS);
                                   // this is to skip a number of columns before parsing the region file


     // before Feb4
     // flood all the [start,end) range of chromosome with 'N'

     // after Feb4
     // flood outside the [start,end) range of chromosome with 'N'
     FloodN::flooded_FASTA(records,OUTPUT_FILE,FASTA_FILE);


//     records.Print();

};
