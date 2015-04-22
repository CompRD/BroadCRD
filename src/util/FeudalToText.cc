/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
///  Transform a feudal file to text format.  Currently implemented
///  for Basevector and Qualvector but the generalization for any file
///  type that has a standard extension and a << operator should be
///  obvious.
///
///  Note that no line breaking is done within a sequence, so the
///  lines can get very long.
///
///  - IN: the name of the file to read.  Results to stdout.
///  - FIELD_DELIM (or F): field delimiter; default is type-dependent.
///  - LINE_DELIM (or L): line delimiter; default is \n.
///  - TYPE (or T): file type (supported: fastb, qualb, intensities); default to check file extension.
///  - START (or S): first line to print (default 0); if negative, count back from end (so -1 is last line).
///  - END (or E): last line to print (default -1 which means end of file).
///
/// \file FeudalToText.cc
///

#include "MainTools.h"
#include "Basevector.h"
#include "Qualvector.h"

template <class FEUDALTYPE>
void ReadFeudal(FEUDALTYPE &seqs, const String &IN, long start, long end)
{
  if ( end < 0 )
  {
      FeudalControlBlock fcb(IN.c_str());
      end = fcb.getNElements();
      if ( start < 0 )
          start += end; // so start==-1 end==-1 indicates printing just the last element
  }
  seqs.ReadRange(IN, start, end);
}

template <class FEUDALTYPE, class T>
void PrintFeudal(const String &IN, const String &fieldDelim,
		 const String &lineDelim, long start, long end)
{
  typedef typename FEUDALTYPE::iterator VBVITR;

  FEUDALTYPE seqs;
  ReadFeudal(seqs, IN, start, end);

  VBVITR vbvend(seqs.end());
  for ( VBVITR vbvitr(seqs.begin()); vbvitr != vbvend; ++vbvitr )
  {
    copy(vbvitr->begin(),vbvitr->end(),
         ostream_iterator<T>(cout,fieldDelim.c_str()));
    cout << lineDelim;
  }
}

// Need to implement this specially because basevector doesn't
// implement begin/end and because we need to run as_base on the
// entries.
template < >
void PrintFeudal<vecbasevector, unsigned int>(const String &IN,
                                const String &fieldDelim,
				const String &lineDelim, long start,
				long end)
{
  typedef vecbasevector::const_iterator VBVITR;
  typedef basevector::const_iterator BVITR;

  vecbasevector seqs;
  ReadFeudal(seqs, IN, start, end);

  VBVITR vbvend(seqs.cend());
  for ( VBVITR vbvitr(seqs.cbegin()); vbvitr != vbvend; ++vbvitr )
  {
      BVITR bvend(vbvitr->End());
      for ( BVITR bvitr(vbvitr->Begin()); bvitr != bvend; ++bvitr )
          cout << as_base(*bvitr) << fieldDelim;
    cout << lineDelim;
  }
}

bool tryType(const String IN, String &TYPE, String &fieldDelim,
             const String type, const String newFieldDelim)
{
  if (IN.Contains(type)) {
    TYPE=type;
    fieldDelim = newFieldDelim;
    return true;
  }
  return false;
}

int main(int argc, char *argv[])
{
  RunTime();

  BeginCommandArguments;
  CommandArgument_String(IN);
  CommandArgument_String_Abbr_OrDefault(FIELD_DELIM, F, "DEFAULT");
  CommandArgument_String_Abbr_OrDefault(LINE_DELIM, L, "DEFAULT");
  CommandArgument_String_Abbr_OrDefault(TYPE, T, "");
  CommandArgument_LongLong_Abbr_OrDefault(START, S, 0);
  CommandArgument_LongLong_Abbr_OrDefault(END, E, -1);

  EndCommandArguments;

  // Here are the default delimiters
  String fieldDelim = "\t", lineDelim = "\n";

  if (TYPE=="") { // Attempt type deduction, override default delimiters
    tryType(IN, TYPE, fieldDelim, "fastb", "")
      || tryType(IN, TYPE, fieldDelim, "qualb", " ");
  }

  // If delims given at command line, they win
  if ("DEFAULT"!=FIELD_DELIM)
    fieldDelim = FIELD_DELIM;
  if ("DEFAULT"!=LINE_DELIM)
    lineDelim = LINE_DELIM;

  // Now print it!
  if (TYPE=="fastb")
    PrintFeudal<vecbasevector, unsigned int>(IN,fieldDelim,lineDelim,START,END);
  else if (TYPE=="qualb")
    PrintFeudal<vecqualvector, unsigned int>(IN,fieldDelim,lineDelim,START,END);
  else
    cout << "Unrecognized TYPE " << TYPE << endl;

  return 0;
}
