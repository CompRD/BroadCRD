////////////////////////////////////////////////////////////////////////////
//                  SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//      This software and its documentation are copyright (2009) by the   //
//  Broad Institute.  All rights are reserved.  This software is supplied //
//  without any warranty or guaranteed support whatsoever. The Broad      //
//  Institute is not responsible for its use, misuse, or functionality.   //
////////////////////////////////////////////////////////////////////////////
/*
 * \file UpdateFastn.cc
 * \author tsharpe
 * \date Nov 13, 2009
 *
 * \brief
 */
#include "MainTools.h"
#include "feudal/FeudalFileReader.h"
#include "feudal/IncrementalWriter.h"
#include "CompressedSequence.h"
#include <cstring>

int main( int argc, char** argv )
{
    RunTime( );

    BeginCommandArguments;
    CommandArgument_String_Doc( FASTN, "Name of fastn file to convert." );
    EndCommandArguments;

    FeudalFileReader rdr(FASTN.c_str());
    if ( rdr.getFCB().getSizeofA() == 1 )
    {
        std::cout << "The file " << FASTN << " has already been converted." << std::endl;
        return 0;
    }

    if ( rdr.getFCB().getSizeofA() != 2  )
    {
        FatalErr("The file " << FASTN << " doesn't look like an old compressed sequence feudal file.");
    }

    String newFileName = FASTN + ".new";
    size_t nnn = rdr.getNElements();
    CompressedSequence seq;
    IncrementalWriter<CompressedSequence> wrtr(newFileName.c_str(),nnn);
    for ( size_t idx = 0; idx < nnn; ++idx )
    {
        size_t varDataLen = rdr.getDataLen(idx);
        if ( varDataLen & 1 )
            FatalErr("I'm confused about the var data len of element "<<idx<<" in file "<<FASTN);

        unsigned int nEle;
        memcpy(&nEle,rdr.getFixedData(idx,sizeof(nEle)),sizeof(nEle));
        if ( (nEle+4)/5 != varDataLen/2 )
            FatalErr("I'm confused about the number of bases in element"<<idx<<" in file "<<FASTN);

        seq.clear().reserve(nEle);

        BinaryReader& br = rdr.getData(idx);
        unsigned short val = 0;
        for ( unsigned int iii = 0; iii < nEle; ++iii )
        {
            if ( !(iii % 5) )
	        br.read(&val);

            unsigned char vvv = val & 7;
            if ( vvv > 4 )
                FatalErr("I'm confused about the value of base "<<iii<<" in element "<<idx<<" in file "<<FASTN<<" -- it's "<<vvv);
            if ( vvv == 4 ) vvv = 15;
            else vvv = 1 << vvv;
            seq.push_back(vvv);
            val >>= 3;
        }

        wrtr.add(seq);
    }
    wrtr.close();

    String oldFileName = FASTN + ".old";
    Rename(FASTN,oldFileName);
    Rename(newFileName,FASTN);

    std::cout << "OK." << std::endl;
    return 0;
}
