/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


#include "jumping/ExonStats.h"
#include "ParseSet.h"

void cleanup(char * & str) {
  if (str[0] == '"') ++str;
  if (str[strlen(str)-1] == '"') str[strlen(str)-1] = 0;
  if (str[strlen(str)-1] == ',') str[strlen(str)-1] = 0;
}
  
void ReadFromBrowserLine(const String & line, vec<ExonStats> & exons) {
  ExonStats es;
  char name[20], strchrom[3];
  int txStart, txEnd, cdsStart, cdsEnd, exonCount;
  char strExonStarts[2000], strExonEnds[2000];
  int count = sscanf(line.c_str(), 
		     "%s\tchr%s\t%c%d\t%d\t%d\t%d\t%d\t%s\t%s\t",
		     name, strchrom, &es.strand, &txStart, &txEnd, 
		     &cdsStart, &cdsEnd, &exonCount, strExonStarts,
		     strExonEnds);
  if (10 != count) {
    cerr << "Bad line: " << line;
    return;
  }
  if (strcmp(strchrom,"X") == 0) {
    es.chromosome = 23;
  } else if (strcmp(strchrom,"Y") == 0) {
    es.chromosome = 24;
  } else if (strcmp(strchrom,"M") == 0) {
    es.chromosome = 0;
  } else es.chromosome = atoi(strchrom);

  es.cdna = name;

  vec<int> exonStarts, exonEnds;
  //remove " from strings if they are there, and parse into vectors.
  char * tempcstr = strExonStarts; //needed so we can modify the array start;
  cleanup(tempcstr);
  ParseIntSet(String("{") + tempcstr + "}", exonStarts, true); 

  tempcstr = strExonEnds; 
  cleanup(tempcstr);
  ParseIntSet(String("{") + tempcstr + "}", exonEnds, true); 

  for (int i=0; i != exonCount; ++i) {
    es.exonIndex = i;
    es.start = exonStarts[i];
    es.length = exonEnds[i] - exonStarts[i];
 
    es.closestExon = INT_MAX;
    if (i!= 0) es.closestExon = es.start - exonEnds[i-1];
    if (i+1 < exonCount) {
      es.closestExon = min(es.closestExon, exonStarts[i+1] - exonEnds[i]);
    }    
    exons.push_back(es);
  }
}
