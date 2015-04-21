
/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// Given the probabilities of an error in base 1,2,...,n of a read,
// how many reads do we expect will have 0,1,2,... errors?
// Creates a table of probabilities for each base/error combination
// sorted by probability, most probably first.
//
// Command-line arguments:
//
//   FILE:      Name of a file with percent error for each base,
//   NUM_BASES: How long are the reads?
//              If 0, read length = amount of data from FILE.
//              If nonzero, we'll ignore the tail of the file.
//   TABLE_OUT: Filename for table of probabilities
//   MIN_ERRORS: Minimum number of errors per read
//               Default: 1
//   MAX_ERRORS: Maximum number of errors per read, or...
//               Default: 3
//   MAX_ENTRIES: Maximum number of entries in table
//                Default: 0, no limit



#include "MainTools.h"
#include "system/System.h"
#include "math/Functions.h"
#include "math/Combinatorics.h"
#include "paths/BaseErrorProb.h"


int main( int argc, char *argv[] )
{
  BeginCommandArguments;

  CommandArgument_String(FILE);
  // values in the file are given as percents.

  CommandArgument_UnsignedInt_OrDefault(NUM_BASES, 0);
  // 0 = use the whole file; nonzero = ignore later values.

  CommandArgument_UnsignedInt_OrDefault(MIN_ERRORS, 1);
  // Minimum number of errors per read

  CommandArgument_UnsignedInt_OrDefault(MAX_ERRORS, 2);
  // Maximum number of errors per read

  CommandArgument_UnsignedInt_OrDefault(MAX_ENTRIES, 0);
  // Maximum number of entries in table

  CommandArgument_Bool_OrDefault(PEEK_AHEAD, False);
  // Add entry to table for most prob. MAX_ERRORS+1 error

  CommandArgument_String_OrDefault(TABLE_OUT, "stdout");
  // Write table to text file

  EndCommandArguments;


  // Import per base error probability distribution
  BaseErrorProbProfile prob(FILE);

  
  // Generate error probability table
  BaseErrorTable error_table;
  error_table = prob.getErrorTable(NUM_BASES, MIN_ERRORS, MAX_ERRORS, MAX_ENTRIES,
				   PEEK_AHEAD);

  // Write table to stdout or file
  ofstream file_stream_out;
  if( TABLE_OUT != "stdout" )
    OpenOfstream( file_stream_out, "file_stream_out" , TABLE_OUT );
  ostream& out = (TABLE_OUT=="stdout" ? cout : file_stream_out);
  
  out << error_table << endl;

}
