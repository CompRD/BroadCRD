///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Convert fasta bases to amino acids.

#include "Basevector.h"
#include "FetchReads.h"
#include "MainTools.h"
#include "amino/Amino.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(IN);
     CommandArgument_Int_Doc(FRAME, "between 1 and 6");
     CommandArgument_String(OUT);
     CommandArgument_UnsignedInt_OrDefault(COLS,240);
     EndCommandArguments;

     if ( FRAME < 1 || FRAME > 6 )
         FatalErr("FRAME must be between 1 and 6");
     bool rc = ( FRAME >= 4 );
     if ( rc )
         FRAME -= 3;
     FRAME -= 1;

     Ofstream( out, OUT );
     vecbasevector bases;
     vecString names;
     FetchReads( bases, names, IN );
     for ( int i = 0; i < (int) bases.size( ); i++ )
     {    out << ">" << names[i] << "\n";
          basevector const& b = bases[i];
          unsigned col = 0;
          if ( rc )
              for ( auto itr=b.rcbegin(FRAME),end=b.rcend(); itr+2<end; itr+=3 )
              {
                  out << AminoAcid::forBases(itr).getSymbol();
                  if ( !(++col % COLS) ) out << '\n';
              }
          else
              for ( auto itr=b.begin(FRAME),end=b.end(); itr+2<end; itr+=3 )
              {
                  out << AminoAcid::forBases(itr).getSymbol();
                  if ( !(++col % COLS) ) out << '\n';
              }
          if ( col % COLS )
              out << '\n';
     }
}
