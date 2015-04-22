/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file CreateRefDict.cc
 * \author tsharpe
 * \date Dec 11, 2008
 *
 * \brief Reads a FASTA file, and writes a .dict file.
 */
#include "MainTools.h"
#include "util/RefDesc.h"

int main( int argc, char** argv )
{
    RunTime();

    BeginCommandArguments;
    CommandArgument_String_Doc(REF,"The reference FASTA file.");
    EndCommandArguments;

    String const FASTA(".fasta");
    String const DICT(".dict");

    String outFile(REF.ReplaceExtension(FASTA,DICT));
    if ( IsSomeSortOfFile(outFile) )
    {
        cout << "The dictionary already exists in file " << outFile << endl;
        cout << "Remove it and run this program again if you'd like to regenerate it." << endl;
        exit(1);
    }

    String trueOutFile(outFile);
    String trueRefFile(REF);
    if ( IsSymbolicLink(REF) )
    {
        trueRefFile = ReadSymbolicLink(REF);
        trueOutFile = trueRefFile.ReplaceExtension(FASTA,DICT);
        if ( IsSomeSortOfFile(trueOutFile) )
        {
            Symlink(trueOutFile,outFile);
            cout << "The dictionary already existed in file " << trueOutFile << endl;
            cout << "We've just soft-linked you to it." << endl;
            exit(0);
        }
    }

    vec<RefDesc> dict = RefDesc::readFASTA(trueRefFile);
    if ( !dict.size() )
    {
        cout << "We didn't find any sequence in the reference file " << REF << endl;
        cout << "So we didn't write a dictionary.  Sorry.  Maybe check that reference file?" << endl;
        exit(1);
    }

    if ( !RefDesc::writeDict(trueOutFile,dict) )
    {
        cout << "Writing the file " << trueOutFile << " failed for some reason.  Sorry." << endl;
        cout << "Check to see that the location is writable, that the disk is not full, etc." << endl;
        exit(1);
    }

    if ( outFile != trueOutFile )
    {
        Symlink(trueOutFile,outFile);
    }

    cout << "Dictionary written to " << outFile << endl;
    cout << "You'll want to edit this file." << endl;
    cout << "FASTA comments are not standardized, so this new dictionary" << endl;
    cout << "will almost certainly require some manual clean-up." << endl;
    cout << "You may also want to reorder the descriptors:  their order in this file" << endl;
    cout << "controls the order of sorting in a BAM file." << endl;
    cout << endl;

    return 0;
}
