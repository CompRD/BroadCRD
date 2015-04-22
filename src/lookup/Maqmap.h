/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file Maqmap.h
 * \author tsharpe
 * \date May 22, 2009
 *
 * \brief Utility classes for dealing with the maq aligner's map output.
 */
#ifndef LOOKUP_MAQMAP_H_
#define LOOKUP_MAQMAP_H_

#include <algorithm>
#include <fstream>
#include <zlib.h>
#include <string.h>
#include "util/RefDesc.h"

struct MapRec
{
    unsigned char seq[127];
    unsigned char se_map_qual_or_indel_len;
    unsigned char size;
    unsigned char map_qual_or_indel_pos;
    unsigned char mismatch_counts;
    unsigned char no_of_errors;
    unsigned char count_0_error_matches;
    unsigned char count_1_error_matches;
    unsigned char flag;
    unsigned char alt_qual_or_mate_map_qual;
    unsigned int seqid;
    unsigned int pos;
    int dist;
    char name[36];

    static unsigned char const FLAG_FF = 0x01;
    static unsigned char const FLAG_FR = 0x02;
    static unsigned char const FLAG_RF = 0x04;
    static unsigned char const FLAG_RR = 0x08;
    static unsigned char const FLAG_PAIRED = 0x10;
    static unsigned char const FLAG_DIFFCHR = 0x20;
    static unsigned char const FLAG_NOMATCH = 0x40;
    static unsigned char const FLAG_SW = 0x80;
    static unsigned char const FLAG_UNMAPPED_READ = FLAG_NOMATCH|FLAG_SW;
    static unsigned char const FLAG_PROPER_PAIRED = FLAG_PAIRED|FLAG_FR;
    static unsigned char const FLAG_SW_PAIRED = FLAG_SW|FLAG_FR;

    bool operator<( MapRec const& that ) const
    {
        bool result = false;
        int val = strcmp(name,that.name);
        if ( val < 0 )
            result = true;
        else if ( !val )
        {
            val = 0;
            unsigned int nnn = std::min(size,that.size);
            for ( unsigned int iii = 0; !val && iii < nnn; ++iii )
                val = seq[iii]&0xc0 - that.seq[iii]&0xc0;

            if ( val < 0 )
                result = true;
            else if ( !val )
                result = size < that.size;
        }
        return result;
    }
};

/// Reads a BFQ (binary fastq) file, and returns MapRecs.
/// The MapRecs are created to appear as if maq couldn't map the read.
class BFQReader
{
public:
    BFQReader( char const* bfqFile );

    bool hasNext()
    { return mIS.peek() != EOF; }

    MapRec const& next();

    void close()
    { mIS.close(); }

private:
    BFQReader( BFQReader const& ); // unimplemented -- no copying
    BFQReader& operator=( BFQReader const& ); // unimplemented -- no copying

    std::ifstream mIS;
    MapRec mMapRec;
};

/// Reads a maq output file (a map file), and returns MapRecs.
class MaqMapReader
{
public:
    /// Second arg is a dictionary to compare to the map file's header info.
    MaqMapReader( char const* mapFile, RefDict const& dict );

    bool hasNext()
    { bool result = false;
      int ch = gzgetc(mGZFile);
      if ( ch != -1 ) { gzungetc(ch,mGZFile); result = true; }
      return result; }

    MapRec const& next();

    void close()
    { gzclose(mGZFile); }

private:
    MaqMapReader( MaqMapReader const& ); // unimplemented -- no copying
    MaqMapReader& operator=( MaqMapReader const& ); // unimplemented -- no copying

    unsigned int readUInt();

    gzFile mGZFile;
    MapRec mMapRec;

    static unsigned int const MAQMAP_FORMAT_NEW = ~0U;
};

#endif /* MAQMAP_H_ */
