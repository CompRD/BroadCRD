/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// SelectBases: select some bases from a CONTIG in a fasta file INPUT.
//
// START - mandatory
//
// Exactly one of STOP or LEN must be supplied.
//
// [START,STOP) is half-open interval.

#include <ctype.h>
#include <fstream>

#include "MainTools.h"

namespace {
// largely via: http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
std::istream& safeGetline(std::istream& is, std::string& t) {
     t.clear();

     // The characters in the stream are read one-by-one using a std::streambuf.
     // That is faster than reading them one-by-one using the std::istream.
     // Code that uses streambuf this way must be guarded by a sentry object.
     // The sentry object performs various tasks,
     // such as thread synchronization and updating the stream state.

     std::istream::sentry se(is, true);
     std::streambuf* sb = is.rdbuf();

     for (;;) {
	  int c = sb->sbumpc();
	  switch (c) {
	  case '\n':
	       return is;
	       break; //NOTREACHED
	  case '\r':
	       if (sb->sgetc() == '\n')
		    sb->sbumpc();
	       return is;
	       break; // NOTREACHED
	  case EOF:
	       // Also handle the case when the last line has no line ending
	       if (t.empty())
		    is.setstate(std::ios::eofbit);
	       return is;
	       break; // NOTREACHED
	  default:
	       t += static_cast<char>(c);
	       break;
	  }
     }

     return is; // NOTREACHED
}
};

int main(int argc, char *argv[]) {
     RunTime();
     BeginCommandArguments;
     CommandArgument_String_Doc(PRE, "Source fasta file is PRE/INPUT");
     CommandArgument_String(INPUT);
     CommandArgument_String(OUTPUT);
     CommandArgument_UnsignedInt_OrDefault_Doc(CONTIG, 0,
	       "Contig (record in fasta file) to extract bases from (0-based)");
     CommandArgument_String_OrDefault_Doc(CONTIG_NAME, "",
	       "Contig name (record in fasta file) to extract bases from "
	       "(overrides CONTIG)");
     CommandArgument_UnsignedInt_Doc(START,
	       "Extract bases [START,STOP) (0-based) from CONTIG");
     CommandArgument_Int_OrDefault(STOP, -1);
     CommandArgument_Int_OrDefault_Doc(LEN, -1,
	       "Extract LEN bases starting at START (0-based) from CONTIG");
     CommandArgument_Bool_OrDefault_Doc(UPCASE, False,
	       "Optionally set bases to upper case.");
     EndCommandArguments;

     ForceAssert( STOP >= 0 ^ LEN >= 0);
     if (LEN >= 0)
	  STOP = START + LEN;

     --STOP;

     std::ofstream os(OUTPUT.c_str());

     int contig_count = -1, base_count = 0;
     Bool at_contig = False;
     Ifstream( in, PRE + "/" + INPUT);

     // then pull the bases
     const char c_start = '>';
     if (CONTIG_NAME != "")
	  CONTIG_NAME += c_start;

     while (1) {
	  std::string s;

	  if (!safeGetline(in, s))
	       FatalErr("Input failure: at end of source fasta file, did not get requested bases");

	  if (s[0] == c_start) {

	       if (at_contig)
		    FatalErr("Input failure: bases selected stretched past a contig");

	       ++contig_count;

	       if (!CONTIG_NAME.empty()) {

		    if (s == CONTIG_NAME) {
			 os << '>' << s << ", bases " << START << " -- " << STOP
				   << "\n";
			 at_contig = True;
		    }

	       } else {

		    if (contig_count == (int) CONTIG) {
			 os << s << ", bases " << START << " -- " << STOP
				   << "\n";
			 at_contig = True;
		    } else if (contig_count > (int) CONTIG)
			 FatalErr("Input error: requested contig outside of the contig id range");
		    }

	       continue;
	  }
	  if (!at_contig) {
	       // if the line started with '>' but was not the right contig,
	       // then 'continue' was executed in the 'if' above; we would not get here.
	       // we are here only if we are currently at sequence line (not '>') AND
	       // the contig is wrong
	       continue;
	  }

	  for (size_t i = 0; i < s.size(); ++i) {
	       const char& c = s[i];
	       if (isalpha(c)) {
		    if ((int) START <= base_count && base_count <= (int) STOP) {
			 char uc = UPCASE ? toupper(c) : c;
			 os << uc;

			 if (base_count == (int) STOP)
			      break; // DONE! quit the loop and (later) the program!
			 if ((base_count + 1 - START) % 80 == 0)
			      os << "\n";
		    }
		    ++base_count;
	       }
	  }
	  if (base_count == (int) STOP)
	       break;
     }
     os << "\n";
     os.close();
     cout << Date( ) << ": SelectBases done" << endl;
}
