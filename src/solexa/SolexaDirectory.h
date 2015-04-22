/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file SolexaDirectory.h
 * \author tsharpe
 * \date Jan 22, 2009
 *
 * \brief Navigation of Solexa run folders.
 */
#ifndef SOLEXA_SOLEXADIRECTORY_H_
#define SOLEXA_SOLEXADIRECTORY_H_

#include "system/file/Directory.h"
#include <string>
#include <map>
#include <cctype>

/// This class hands out File objects representing all the stuff
/// one normally finds in a Solexa run folder.
/// Normally, you'd get an instance of this class from the SolexaPipelineDir.
class SolexaDir
{
public:
    enum Aligner { MERLIN, PAIREDREADPARAMSBYREAD, PAIRFINDER };
    enum End { END1, END2 };

    SolexaDir( Directory const& dir, std::string head ) // for single-end runs
    : mDir1(dir), mPrefix(SLASH+head+DOT), mIsPaired(false)
    {}

    SolexaDir( Directory const& dir1, Directory const& dir2, std::string head ) // for paired-end runs
    : mDir1(dir1), mDir2(dir2), mPrefix(SLASH+head+DOT), mIsPaired(true)
    {}

    // compiler-supplied copying and destructor are OK

    bool isPaired() const
    { return mIsPaired; }

    File getQltout( End end = END1, Aligner aligner = MERLIN ) const
    {
        switch ( aligner )
        {
        default:
        case MERLIN:
            return File(dir(end) + QLTOUT_TYPE);
        case PAIREDREADPARAMSBYREAD:
            return File(dir(END1) + pairedReadParamsByReadSuffix(end) + QLTOUT_TYPE);
        case PAIRFINDER:
            return File(dir(END1) + pairFinderSuffix(end) + QLTOUT_TYPE);
        }
    }

    File getRefFasta() const
    { return File(dir(END1)+REFFASTA_TYPE); }

    File getRefDict() const
    { return File(dir(END1)+REFDICT_TYPE); }

    File getMetrics( End end = END1 ) const
    { return File(dir(end)+METRICS_TYPE); }

    File getFastb( End end = END1, Aligner aligner = MERLIN ) const
    {
        switch ( aligner )
        {
        default:
        case MERLIN:
        case PAIRFINDER:
            return File(dir(end) + FASTB_TYPE);
        case PAIREDREADPARAMSBYREAD:
            return File(dir(END1) + pairedReadParamsByReadSuffix(end) + FASTB_TYPE);
        }
    }

    File getQualb( End end = END1, Aligner aligner = MERLIN ) const
    {
        switch ( aligner )
        {
        default:
        case MERLIN:
        case PAIRFINDER:
            return File(dir(end) + QUALB_TYPE);
        case PAIREDREADPARAMSBYREAD:
            return File(dir(END1) + pairedReadParamsByReadSuffix(end) + QUALB_TYPE);
        }
    }

    File getAdaptedQualb( End end = END1 ) const
    { return File(dir(end)+ADAPTED_QUALB_TYPE); }

    File getFilterVector( End end = END1 ) const
    { return File(dir(end)+PASSFILTER_TYPE); }

    /// see if the filename looks like a solexa dir ( i.e., <date>_<flowcell> )
    /// return the flowcell part, if so
    static bool isSolexaDir( File const& file, std::string& flowcell );

    static std::string const FASTB_TYPE;
    static std::string const QUALB_TYPE;
    static std::string const ADAPTED_QUALB_TYPE;
    static std::string const QLTOUT_TYPE;
    static std::string const REFDICT_TYPE;
    static std::string const REFFASTA_TYPE;
    static std::string const METRICS_TYPE;
    static std::string const PASSFILTER_TYPE;
    static std::string const SAM_TYPE;
    static std::string const DICT_TYPE;

private:
    std::string dir( End end ) const;

    std::string const& pairedReadParamsByReadSuffix( End end ) const
    { return PAIREDREADPARAMSBYREAD_SUFFICES[end]; }

    std::string const& pairFinderSuffix( End end ) const
    { return PAIRFINDER_SUFFICES[end]; }

    Directory mDir1;
    Directory mDir2;
    std::string mPrefix;
    bool mIsPaired;

    static std::string const PAIREDREADPARAMSBYREAD_SUFFICES[2];
    static std::string const PAIRFINDER_SUFFICES[2];
    static std::string const SLASH;
    static std::string const DOT;

    friend class SolexaPipelineDir;
};

/// class to find SolexaDirs, given a head (i.e., read group, i.e., "flowcell.lane") string.
/// methods are all static, because the class is a singleton
class SolexaPipelineDir
{
public:
    static SolexaDir getSolexaDir( std::string const& head )
    { if ( !gpInstance ) gpInstance = new SolexaPipelineDir();
      return gpInstance->find(head); }

    static void refresh()
    { delete gpInstance; gpInstance = 0; }

private:
    SolexaPipelineDir(); // private constructor:  singleton

    SolexaDir find( std::string const& head ) const;

    std::multimap<std::string,Directory> mMap; // multimap of flowcell onto Directory

    static char const DEFAULT_LOCATION[];
    static SolexaPipelineDir* gpInstance;
};

#endif /* SOLEXA_SOLEXADIRECTORY_H_ */
