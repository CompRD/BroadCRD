// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research
// 


// Fasta: convert a file reads.fastb into a human-readable reads.fasta.
// If run with CLEAN=True, generate a blastable file (no blank lines,
// no empty reads).  If NAMES=True, use original read names.
// If MAXREADS > 0, print no more than MAXREADS reads.

#include "CompressedSequence.h"
#include "math/Functions.h"
#include "system/ParsedArgs.h"
#include "String.h"
#include "system/System.h"
#include "Vec.h"

int main( int argc, char *argv[] )
{
     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(HEAD, "reads_orig");
     CommandArgument_Bool_OrDefault(CLEAN, False);
     CommandArgument_Bool_OrDefault(NAMES, False);
     CommandArgument_UnsignedInt_OrDefault(MAXREADS, 0);
     EndCommandArguments;

     String run_dir = PRE + "/" + DATA + "/" + RUN;

     vecString ids;
     if (NAMES) ids.ReadAll( run_dir + "/reads.ids" );

     veccompseq EE;
     EE.ReadAll( run_dir + "/" + HEAD + ".fastn" );

     Ofstream( out, run_dir + "/" + HEAD + ".fasta" );
     if ( MAXREADS == 0 ) 
       MAXREADS = EE.size( );

     vec<char> bases;
     for ( int id = 0; id < Min( (int) MAXREADS, (int) EE.size( ) ); id++ )
     { 
       if ( CLEAN && EE[id].size( ) == 0 ) 
	 continue;
       if ( !NAMES ) 
	 out << ">read_" << id+1 << "\n";
       else
	 out << ">" << ids[id] << "\n";
       bases = EE[id].asVecChar();
       for ( unsigned int j = 0; j < bases.size(); j++ )
       {
	 if ( j > 0 && j % 80 == 0 ) 
	   out << "\n";
	 out << bases[j];
       }
       out << ( CLEAN ? "\n" : "\n\n" );
     }
}
