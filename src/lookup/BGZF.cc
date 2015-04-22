/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file BGZF.cc
 * \author tsharpe
 * \date May 21, 2009
 *
 * \brief Utilities for writing BAM files.
 */
#include "lookup/BGZF.h"
#include <zlib.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>

namespace
{
    void fatalErr( char const* msg )
    {
        std::cerr << "Can't write BAM file: " << msg << std::endl;
        exit(1);
    }
}

unsigned int BGZFBlock::compress( void* data, unsigned int len )
{
    while ( !tryCompress(data,len) )
        ;
    return len;
}

bool BGZFBlock::tryCompress( void* data, unsigned int& len )
{
    z_stream zs;
    zs.zalloc = 0;
    zs.zfree = 0;
    zs.opaque = 0;
    if ( deflateInit2(&zs,Z_DEFAULT_COMPRESSION,Z_DEFLATED,WINDOW_BITS,MEM_LEVEL,Z_DEFAULT_STRATEGY) != Z_OK )
        fatalErr("Can't initialize BAM file's z_stream.");

    zs.next_in = reinterpret_cast<Bytef*>(data);
    zs.avail_in = len;
    zs.next_out = mDataBlock;
    zs.avail_out = sizeof(mDataBlock);
    if ( deflate(&zs,Z_FINISH) != Z_STREAM_END )
    {
        zs.avail_in += 256;
        if ( zs.avail_in >= len )
            fatalErr("Can't deflate even a tiny amount of data for BAM file.");
        len -= zs.avail_in; // back off by the amount we were unable to compress, and a little more
        return false;
    }

    if ( deflateEnd(&zs) != Z_OK )
        fatalErr("Can't finalize BAM file's z_stream.");

    unsigned char* ppp = mDataBlock + zs.total_out;
    unsigned int crc = crc32(crc32(0,0,0),reinterpret_cast<Bytef*>(data),len);
    memcpy(ppp,&crc,sizeof(unsigned int)); // setting mCRC32
    ppp += sizeof(unsigned int);
    memcpy(ppp,&len,sizeof(unsigned int)); // setting mInputSize
    ppp += sizeof(unsigned int);
    unsigned int compressedSize = sizeof(BGZFBlock)-zs.avail_out;
    mBlockSizeLessOne = (ppp - reinterpret_cast<unsigned char*>(this)) - 1U;

    return true;
}

BGZFStreambuf::int_type BGZFStreambuf::overflow( int_type ch )
{
    if ( ch != traits_type::eof() )
    {
        *pptr() = ch;
        pbump(1);
    }

    char* ppp = pptr();
    unsigned int inLen = ppp - mBuf;
    if ( !inLen )
        return 0; // EARLY RETURN!

    setp(mBuf,epptr());

    BGZFBlock output;
    inLen -= output.compress(mBuf,inLen);
    if ( inLen )
    {
        memmove(mBuf,ppp-inLen,inLen);
        pbump(inLen);
    }

    if ( mpSB->sputn(reinterpret_cast<char const*>(&output),output.getBlockSize()) != output.getBlockSize() )
        fatalErr("Can't write to BAM file.");

    return 0;
}

int BGZFStreambuf::sync()
{
    while ( pptr() != pbase() )
        overflow( traits_type::eof() );
    return 0;
}

BAMostream::BAMostream( char const* bamFile, char const* header, RefDict const& dict )
: std::ostream(&mSB), mSB(&mFilebuf)
{
    mFilebuf.open(bamFile,std::ios_base::out|std::ios_base::binary|std::ios_base::trunc);
    init(header,dict);
}

void BAMostream::write( BAMAlignment const& alignment )
{
    std::ostream::write(reinterpret_cast<char const*>(&alignment),alignment.blockSize+sizeof(alignment.blockSize));
}

void BAMostream::close()
{
    mSB.pubsync();
    mFilebuf.pubsync();
    mFilebuf.close();
}

void BAMostream::init( char const* header, RefDict const& dict )
{
    write("BAM\1",4U);
    unsigned int headerLen = strlen(header);
    write(&headerLen,sizeof(headerLen));
    write(header,headerLen);
    unsigned int dictSize = dict.size();
    write(&dictSize,sizeof(dictSize));
    RefDict::const_iterator end(dict.end());
    for ( RefDict::const_iterator itr(dict.begin()); itr != end; ++itr )
    {
        RefDesc const& desc = *itr;
        unsigned int nameLen = desc.getName().size() + 1;
        write(&nameLen,sizeof(nameLen));
        write(desc.getName().c_str(),nameLen);
        unsigned int refLen = desc.getLength();
        write(&refLen,sizeof(refLen));
    }
}
