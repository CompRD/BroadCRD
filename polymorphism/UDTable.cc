/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// UeberDrueber.  Find polymorphisms in test strain REFA vs finished strain REFB.
#include "polymorphism/UDTable.h"

#define DOT_PRINT 10000

// vec< map<basevector, vec<unsigned int> > > theData;

void UDTable::Read(String filename) {
   Ifstream(fin, filename);
   unsigned int numcontig;
   int temp_blockSize;
   fin >> temp_blockSize;
   if (temp_blockSize != blockSize) {
      cerr << "ERROR!  Block size in file is not the same as expected block size.  Either delete and regenerate the table or change the expected block size.\n";
      exit(1);
   }
   fin >> numcontig;
   theData.clear_and_resize(numcontig);
   for (unsigned int contig = 0; contig < numcontig; ++contig) {
      int numData;
      fin >> numData;
      for (int miter = 0; miter < numData; ++miter) {
         String theseq;
         fin >> theseq;
         basevector thebv;
         thebv.SetFromString(theseq);
         unsigned int numpos;
         fin >> numpos;
         vec<int> temppos;
         for (unsigned int i = 0; i < numpos; ++i) {
            unsigned int thispos;
            fin >> thispos;
            temppos.push_back(thispos);
         }
         theData[contig][thebv] = temppos;
      }
   }
   fin.close( );
}

void UDTable::Write(String filename) {
   Ofstream(fout, filename);
   fout << blockSize << endl;
   fout << theData.size( ) << endl;
   for (unsigned int contig = 0; contig < theData.size( ); ++contig) {
      fout << theData[contig].size( ) << endl;
      for (map< basevector, vec<int> >::iterator miter = theData[contig].begin();
           miter != theData[contig].end(); ++miter) {
         String theseq = ToString(miter->first);
         fout << theseq << endl;
         vec<int>& temppos = miter->second;
         fout << temppos.size( ) << endl;
         for (unsigned int i = 0; i < temppos.size( ); ++i)
            fout << temppos[i] << "  ";
         fout << endl;
      }
   }
   fout.close( );
}

// Later, generate table contig by contig and read/write by contig for memory
// useage
void UDTable::GenerateTable(vecbasevector& ref) {
    theData.resize(ref.size());
    basevector temp;
    for (size_t contig = 0; contig < ref.size( ); contig++) {
       for (unsigned int pos = 0; pos < ref[contig].size( ) - blockSize; pos++) {
          temp.SetToSubOf(ref[contig], pos, blockSize);
          theData[contig][temp].push_back(pos);
       }
    }
}

