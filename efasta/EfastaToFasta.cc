///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// EfastaToFasta: convert a file of nucleotides from efasta to fasta format.
// Ambiguities that do not translate to ambiguous base codes are lost: we just
// make the first choice.  This code also turns Ns into ns.

#include "FastIfstream.h"
#include "MainTools.h"
#include "TokenizeString.h"
#include "efasta/EfastaTools.h"

#define Err(message)                                      \
{    cout << message << endl << "\nInvalid.\n" << endl;   \
     return 1;    }

int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String(HEAD);
  CommandArgument_String_OrDefault(SPLIT_DIR, "");
  CommandArgument_Bool_OrDefault_Doc(IUPAC, True,
      "Keep ambiguities that can be expressed by IUPAC code" );
  CommandArgument_Bool_OrDefault_Doc(EXPAND, False, "instead expand efasta" );
  CommandArgument_Bool_OrDefault_Doc(ALLOW_EMPTY, False, "allow empty records" );
  EndCommandArguments;

  ofstream out;
  ofstream tbl;

  if (SPLIT_DIR == "") {
    cout << "Creating '" << HEAD << ".{fasta,tbl}'." << endl;
    OpenOfstream(out, "out", HEAD + ".fasta");
    OpenOfstream(tbl, "tbl", HEAD + ".tbl");
  }
  else {
    Mkdir777(SPLIT_DIR);
    cout << "Creating many '" << SPLIT_DIR << "/<name>.{fsa,tbl}'." << endl;
  }

  fast_ifstream in(HEAD + ".efasta");
  String line;
  while(1) {
    getline(in, line);
    if (in.fail()) {
      if (line.size() == 0) Err("Illegal empty file.");
      break;
    }
    if (line.size() == 0) continue;
    if (line[0] != '>') {
      Err("See line = '" << line << "', which was expected "
           << "to start with >.");
    }
    String name = "";
    for (unsigned i = 1; i != line.size(); i++) 
      name.push_back(line[i]);
    
    if (SPLIT_DIR != "") {
      String head = SPLIT_DIR + "/" + name;
      OpenOfstream(out, "out_i", head + ".fsa");
      OpenOfstream(tbl, "tbl_i", head + ".tbl");
    }

    String iline = line;

    vec<String> lines;
    Bool eof = False;
    while (1) {
      char c;
      in.peek(c);
      if (in.fail()) { eof = True; break; }
      if (c == '>') break;
      getline(in, line);
      lines.push_back(line);
    }
    if (lines.empty()) 
    {    if (ALLOW_EMPTY)
         {    out << iline << "\n";
              if (eof) break;
              continue;    }
         Err("Illegal record of empty length.");    }

    String all;
    int64_t all_size = 0;
    for (size_t i = 0; i < lines.size(); i++)
      all_size += lines[i].size();
    all.reserve(all_size);
    for (size_t i = 0; i < lines.size(); i++)
      all.append(lines[i]);
    ValidateEfastaRecord(lines);

    if (EXPAND)
    {    vec<basevector> bases;
         efasta(all).ExpandTo(bases);
         for ( int j = 0; j < bases.isize( ); j++ )
         {    out << iline << "[" << j << "]\n";
              const basevector& b = bases[j];
              for ( int l = 0; l < b.isize( ); l++ )
              {    if ( l > 0 && l % 80 == 0 ) out << "\n";
                   out << as_base( b[l] );    }
              out << "\n";    }
         continue;    }

    out << iline << "\n";

    fastavector v;
    vec<Ambiguity> va;
    if ( IUPAC )
      efasta(all).FlattenTo(v, va);
    else
      efasta(all).FlattenTo(v, va, False );

    for (size_t i = 0; i < v.size(); i++) {
      if (i > 0 && i % 80 == 0) out << "\n";
      out << v[i];
    }
    out << "\n";
    

    if (va.size()) 
      tbl << ">Feature    " << name << "\n";
    for (size_t i = 0; i < va.size(); i++)
      tbl << va[i].to_annotation() << endl;


    if (SPLIT_DIR != "") {
      out.close();
      tbl.close();
    }

    if (eof) break;
  }


  if (SPLIT_DIR == "") {
    out.close();
    tbl.close();
  }

}
