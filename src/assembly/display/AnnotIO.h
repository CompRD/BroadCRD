// Copyright (c) 2004 Broad Institute/Massachusetts Institute of Technology

#ifndef _ANNOTIO_H_
#define _ANNOTIO_H_

#include "String.h"
#include "feudal/BinaryStream.h"

// this is a goofy interface: you have to call the read functions, or the write
// functions, in the order given or all hell will break loose.


bool ReadAnnotHeader( BinaryReader& br )
{ int header = 0; br.read(&header); return header == 0x4F4E4E41; }

int ReadAnnotLength( BinaryReader& br )
{ int len; br.read(&len); return len; }

bool ReadOneAnnot( BinaryReader& br, int & cID, int & begin, int & end, String& name)
{ br.read(&cID); br.read(&begin); br.read(&end); br.read(&name); return true; }



void WriteAnnotHeader( BinaryWriter& bw )
{ int header = 0x4F4E4E41; bw.write(header);
  int len = 0; bw.write(len); }

void WriteOneAnnot( BinaryWriter& bw, int cID, int begin, int end, String const& name)
{ bw.write(cID); bw.write(begin); bw.write(end); bw.write(name); }

void WriteAnnotLength( BinaryWriter& bw, int len )
{ bw.seek(sizeof(int)); bw.write(len); }


#endif 
