/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

///
/// \file FeudalAppend.cc
///

const char* DOC =
"Appends one feudal file to another and writes the result to disk.  "
"The output file can be the same as either input file."
"All work is done on disk, so this will work in low memory and is"
"reasonably fast."
"Minimal checking is done to see that IN1 and IN2 are the same type of file.";

#include "MainTools.h"
#include "feudal/FeudalControlBlock.h"
#include "feudal/FeudalTools.h"

int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandDoc(DOC);
  CommandArgument_String_Doc(IN1, "The contents of this file will be first.");
  CommandArgument_String_Doc(IN2, "The contents of this file will be second.");
  CommandArgument_String_Doc(OUT, "The result will be written to this file.");
  EndCommandArguments;

  FeudalControlBlock fcb1(IN1.c_str(),true);
  FeudalControlBlock fcb2(IN2.c_str(),true);

  if ( fcb1.getSizeofA() != fcb2.getSizeofA() )
      FatalErr("IN1 has a sub-element size of " << fcb1.getSizeofA() <<
               " but IN2 has a sub-element size of " << fcb2.getSizeofA());
  if ( fcb1.getSizeofX() != fcb2.getSizeofX() )
      FatalErr("IN1 has an element size of " << fcb1.getSizeofX() <<
               " but IN2 has a element size of " << fcb2.getSizeofX());
  if ( fcb1.getSizeofFixed() != fcb2.getSizeofFixed() )
      FatalErr("IN1 has an fixed-data size of " << fcb1.getSizeofFixed() <<
               " but IN2 has a fixed-data size of " << fcb2.getSizeofFixed());

  String tmpFile(OUT+".tmp");
  MergeMastervecs(IN1,IN2,tmpFile);
  Remove( OUT );
  Mv( tmpFile, OUT );

  SystemSucceed("FeudalSize F="+OUT);

  return 0;
}
