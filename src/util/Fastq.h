///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file Fastq.h
 * \author blau
 * \date Mar 5, 2014
 *
 * \brief Converts a Fastq file to a fastb and qualb.
 */

#ifndef UTIL_FASTQ_H
#define UTIL_FASTQ_H

#include <istream>
#include "String.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "system/ProcBuf.h"

namespace fastq{

// factored out from FastqToFastbQualb.cc at r48873
struct FastqEntry 
{
  char phred_min;
  char phred_max;
  BaseVec sequence;
  QualVec qualities;
  String name;
};

// factored out from FastqToFastbQualb.cc at r48873
class Swizzler
{
public:
    Swizzler( char offset ) : mOffset(offset), mMin('\x7f'), mMax('\x80') {}

    char operator()( char score )
    { if ( score < mMin ) mMin = score;
      if ( score > mMax ) mMax = score;
      return score > mOffset ? score-mOffset : 0; }

    char getMin() const { return mMin; }
    char getMax() const { return mMax; }

private:
    char mOffset;
    char mMin;
    char mMax;
};

// factored out from FastqToFastbQualb.cc at r48873
bool NextFastqEntry( std::istream& is, FastqEntry &entry,
                            vec<size_t> * p_n_amb, const char phred_offset );

// interleave the reads in from fastq1 and fastq2 as pairs, and write to outhead.fastb/qualb/pairs
void Fastq2FastbQualbPair_Interleave( const String& fastq1, const String& fastq2, const char q_offset
                                    , const String& outhead);
// interleave the read 2i and 2i+1 of a fastq file as pairs, and write to outhead.fastb/qualb/pairs
void Fastq2FastbQualbPair( const String& fastq, const char q_offset , const String& outhead);


//mechanism to read line by line
class FastQReader
{
public:
    FastQReader(const String& sFileName, const char q_offset)
        :buffer(),n_amb(4,0),ifs(sFileName.c_str()),phred_offset(q_offset),bRead(true){Step();};

    bool Step(){
        if(bRead){
            bRead=NextFastqEntry(ifs,buffer,&n_amb,phred_offset);
        }
        return bRead;
    }
    bool State()const{return bRead;};
    const FastqEntry& Get()const{return buffer;};


private:
    FastQReader();
    FastQReader(const FastQReader&);
    FastQReader& operator=(const FastQReader&);
    FastqEntry buffer;
    vec<size_t> n_amb;
    std::ifstream ifs;
    const char phred_offset;
    bool bRead;

};


size_t ReadFastq( String const& filename, vecbasevector& bases_out,
        vecqualvector& quals_out, char const phred_offset = 33,
        vec<String>* names_outp = nullptr);

}

#endif

