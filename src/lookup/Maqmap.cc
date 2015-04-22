/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file Maqmap.cc
 * \author tsharpe
 * \date May 22, 2009
 *
 * \brief Utility classes for dealing with the maq aligner's map output.
 */
#include "lookup/Maqmap.h"
#include <iostream>
#include <stdlib.h>
#include <string.h>

namespace
{
    void fatalErr( char const* msg )
    {
        std::cerr << msg << std::endl;
        exit(1);
    }
}

BFQReader::BFQReader( char const* bfqFile )
: mIS(bfqFile,std::ios_base::in|std::ios_base::binary)
{
    memset(mMapRec.seq,0,sizeof(mMapRec.seq));
    mMapRec.se_map_qual_or_indel_len = 0U;
    mMapRec.map_qual_or_indel_pos = 0U;
    mMapRec.mismatch_counts = 0U;
    mMapRec.no_of_errors = 0U;
    mMapRec.count_0_error_matches = 0U;
    mMapRec.count_1_error_matches = 0U;
    mMapRec.flag = MapRec::FLAG_UNMAPPED_READ;
    mMapRec.alt_qual_or_mate_map_qual = 0U;
    mMapRec.seqid = ~0U;
    mMapRec.pos = ~0U;
    mMapRec.dist = 0;
    memset(mMapRec.name,0,sizeof(mMapRec.name));
    if ( !mIS )
    {
        std::cerr << "Unable to open BFQ file: " << bfqFile << std::endl;
        exit(1);
    }
}

MapRec const& BFQReader::next()
{
    unsigned int nameLen;
    if ( !mIS.read(reinterpret_cast<char*>(&nameLen),sizeof(nameLen)) )
        fatalErr("Can't read name length from BFQ file.");

    if ( nameLen > sizeof(mMapRec.name) )
        fatalErr("BFQ name length is too big.");

    if ( !mIS.read(mMapRec.name,nameLen) )
        fatalErr("Can't read BFQ name.");

    if ( strlen(mMapRec.name)+1 != nameLen )
        fatalErr("BFQ name corrupt.");

    unsigned int seqLen;
    if ( !mIS.read(reinterpret_cast<char*>(&seqLen),sizeof(seqLen)) )
        fatalErr("Can't read data length from BFQ file.");

    if ( seqLen > sizeof(mMapRec.seq) )
        fatalErr("BFQ data length is too big.");

    mMapRec.size = seqLen;
    if ( !mIS.read(reinterpret_cast<char*>(mMapRec.seq),seqLen) )
        fatalErr("Can't read BFQ sequence data.");

    return mMapRec;
}

inline unsigned int MaqMapReader::readUInt()
{
    unsigned int val;
    if ( gzread(mGZFile,&val,sizeof(val)) != sizeof(val) )
        fatalErr("Can't read next uint from MAP file.");
    return val;
}

MaqMapReader::MaqMapReader( char const* mapFile, RefDict const& dict )
: mGZFile(gzopen(mapFile,"rb"))
{
    if ( readUInt() != MAQMAP_FORMAT_NEW )
        fatalErr("Wrong MAP file magic number.");

    unsigned int nRefs = readUInt();
    if ( nRefs != dict.size() )
        fatalErr("Number of reference sequences in map file disagrees with no. of references in dictionary.");

    for ( unsigned int iii = 0; iii < nRefs; ++iii )
    {
        unsigned int refNameLen = readUInt();
        std::string const& dictName = dict[iii].getName();
        if ( refNameLen != dictName.size()+1 )
            fatalErr("Reference info inconsistent between MAP file and dictionary (nameLen).");

        char* buf = new char[refNameLen];
        if ( static_cast<unsigned int>(gzread(mGZFile,buf,refNameLen)) != refNameLen )
            fatalErr("Can't read reference name from MAP file.");
        if ( strcmp(buf,dictName.c_str()) )
            fatalErr("Reference info inconsistent between MAP file and dictionary (name).");
        delete [] buf;
    }

    readUInt();
    readUInt();
}

MapRec const& MaqMapReader::next()
{
    if ( gzread(mGZFile,&mMapRec,sizeof(mMapRec)) != sizeof(mMapRec) )
        fatalErr("Unable to read next record from MAP file.");
    if ( mMapRec.size > sizeof(mMapRec.seq) )
        fatalErr("Bogus sequence length in record from MAP file.");
    return mMapRec;
}
