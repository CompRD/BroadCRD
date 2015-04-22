///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// FastaToFastq.  Convert a fasta file into a fastq file.  Assign
// constant quality scores.

#include "Fastavector.h"
#include "FastIfstream.h"
#include "MainTools.h"

int main( int argc, char *argv[] )
{
    String emptyString;
    RunTime();
    BeginCommandArguments;
    CommandArgument_String_Doc(FASTA, "name of fasta file, must end in .fasta");
    CommandArgument_String_OrDefault_Doc(FASTQ, emptyString,
            "name of fastq file (defaults to basename of fasta filename");
    CommandArgument_Int_OrDefault_Doc(Q, 40, "quality score to assign");
    EndCommandArguments;

    vec<fastavector> bases;
    vec<String> names;
    LoadFromFastaFile(FASTA, bases, names);

    ofstream out;
    if ( FASTQ == emptyString )
        FASTQ = FASTA.ReplaceExtension(".fasta",".fastq");
    OpenOfstream(out, "FASTQ", FASTQ);

    char qChar = Q + 33;
    String qString;
    for ( size_t i = 0; i < bases.size(); i++ )
    {
        out << "@" << names[i].SafeBefore(" ") << "\n";
        fastavector const& fv = bases[i];
        typedef fastavector::const_iterator Itr;
        for ( Itr itr(fv.begin()), end(fv.end()); itr != end; ++itr )
            out << *itr;
        out << "\n+\n";
        qString.resize(fv.size(),qChar);
        out << qString;
        out << '\n';
    }
}
