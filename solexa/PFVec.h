/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file PFVec.h
 * \author tsharpe
 * \date Jan 21, 2009
 *
 * \brief Makes a vector<bool> from a [solexa head].filter.vector file.
 *
 * This class just adds some constructors that read the file.
 */
#ifndef PFVEC_H_
#define PFVEC_H_

#include "String.h"
#include <fstream>
#include <vector>

class PFVec : public vector<bool>
{
public:
    PFVec()
    {}

    PFVec( String const& filterVectorFile )
    {
        ifstream is(filterVectorFile.c_str());
        readStream(is,filterVectorFile.c_str());
    }

    PFVec( string const& filterVectorFile )
    {
        ifstream is(filterVectorFile.c_str());
        readStream(is,filterVectorFile.c_str());
    }

    PFVec( char const* filterVectorFile )
    {
        ifstream is(filterVectorFile);
        readStream(is,filterVectorFile);
    }

private:
    void readStream( ifstream& is, char const* filterVectorFile )
    {
        if ( !is.is_open() || !is.good() )
        { FatalErr("Unable to read file: " << filterVectorFile); }

        reserve(MAX_READS_PER_RUN);

        char buf[8192];
        int lineNo = 0;
        while ( is.good() )
        {
            if ( is.getline(buf,sizeof(buf)).gcount() )
            {
                ++lineNo;
                if ( buf[0] == '0' )
                    push_back(false);
                else if ( buf[0] == '1' )
                    push_back(true);
                else
                { FatalErr("In file " << filterVectorFile << ", unable to interpret line #" << lineNo << ": '" << buf << "'"); }
            }
        }
        if ( !is.eof() )
        { FatalErr("Unable to fully process " << filterVectorFile); }

        is.close();
    }

    static size_t const MAX_READS_PER_RUN = 30000000; // Largest expected number of reads in a Solexa run.
    static size_t const BUF_LEN = 8192; // I always use 8192.  I don't know why.  This buffer only needs to be 2 bytes, actually.
};

#endif /* PFVEC_H_ */
