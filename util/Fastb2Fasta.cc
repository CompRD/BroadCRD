///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

const char *DOC =
  "Converts a fastb file into fasta format.";

#include "Basevector.h"
#include "MainTools.h"

int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandDoc(DOC);
  CommandArgument_String_Doc(IN,
    "Fastb file to convert.");
  CommandArgument_String_OrDefault_Doc(OUT, "", 
    "Output fasta filename. Default is to use IN filename, replacing extension "
    "with 'fasta'. (To output to standard out use '-')");
  CommandArgument_Bool_OrDefault_Doc(ONE_LINE_PER_SEQUENCE, False,
    "Output fasta sequence with no line breaks and no record headers");
  CommandArgument_Int_OrDefault_Doc(BREAKCOL, 80,
    "Break lines that exceed this length.");
  CommandArgument_UnsignedInt_OrDefault_Doc(READ_COUNT, 0,
    "Only convert the first READ_COUNT basevectors.");
  EndCommandArguments;

  if (OUT.empty()) 
    OUT = IN.SafeBefore(".fastb") + ".fasta";

  // Setup virtual master vec to stream in the fastb
  typedef VirtualMasterVec<basevector> VmvBv_t;
  VmvBv_t reads( IN.c_str() );

  const bool to_stdout = (OUT == "-");

  size_t count = 0;
  ofstream ofs;
  if (!to_stdout) ofs.open(OUT.c_str());
  ostream & os = (to_stdout) ? cout : ofs;

  for (VmvBv_t::Itr reads_itr = reads.begin(); reads_itr != reads.end() ; ++reads_itr) {
    if (ONE_LINE_PER_SEQUENCE) {
      os << reads_itr->ToString() << endl;
      count++;
    } 
    else {
      reads_itr->PrintCol(os, ToString(count++), BREAKCOL);
    }
    if (count == READ_COUNT)
      break;
  }

  if (!to_stdout) ofs.close();
}
