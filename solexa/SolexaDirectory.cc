/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file SolexaDirectory.cc
 * \author tsharpe
 * \date Jan 23, 2009
 *
 * \brief Navigation of Solexa run folders.
 */
#include "solexa/SolexaDirectory.h"
#include "system/file/Directory.h"
#include "system/file/SymLink.h"
#include "system/System.h"
#include <iostream>

using std::string;
using std::multimap;

namespace
{

struct SingletonKiller
{
    ~SingletonKiller() { SolexaPipelineDir::refresh(); }
    static SingletonKiller gSingletonKiller;
};

SingletonKiller SingletonKiller::gSingletonKiller;

}

string const SolexaDir::FASTB_TYPE("fastb");
string const SolexaDir::QUALB_TYPE("qualb");
string const SolexaDir::ADAPTED_QUALB_TYPE("adapted.qualb");
string const SolexaDir::QLTOUT_TYPE("qltout");
string const SolexaDir::REFDICT_TYPE("reference.dict");
string const SolexaDir::REFFASTA_TYPE("reference.fasta");
string const SolexaDir::METRICS_TYPE("metrics");
string const SolexaDir::PASSFILTER_TYPE("filter.vector");
string const SolexaDir::SAM_TYPE("sam");
string const SolexaDir::DICT_TYPE("dict");

string const SolexaDir::PAIREDREADPARAMSBYREAD_SUFFICES[2] = { string(".a"), string(".b") };
string const SolexaDir::PAIRFINDER_SUFFICES[2] = { string(".end1"), string(".end2") };
string const SolexaDir::SLASH("/");
string const SolexaDir::DOT(".");

char const SolexaPipelineDir::DEFAULT_LOCATION[] = "/slxa/pipelineOutput/farm";
SolexaPipelineDir* SolexaPipelineDir::gpInstance;

bool SolexaDir::isSolexaDir( File const& file, string& flowcell )
{
    bool result = false;

    if ( file.isDir() )
    {
        string const& fileName = file.filename();
        if ( fileName.size() == 12 && fileName[6] == '_' )
        {
            result = true;
            for ( uint iii = 0; iii < 6; ++iii )
            {
                if ( !isdigit(fileName[iii]) )
                {
                    result = false;
                    break;
                }
            }

            if ( result )
            {
                flowcell = fileName.substr(7);
            }
        }
    }

    return result;
}

string SolexaDir::dir( End end ) const
{
    if ( end == END2 && !mIsPaired )
    { FatalErr("There is no paired end for flowcell " << mPrefix.substr(1)); }

    return (end == END1 ? mDir1.toString() : mDir2.toString()) + mPrefix;
}

SolexaPipelineDir::SolexaPipelineDir()
{
    char const* dirName = getenv("SOLEXA_PIPELINE_DIR");
    if ( !dirName )
        dirName = DEFAULT_LOCATION;
    Directory dir(dirName);

    string flowCell;
    Directory::const_iterator end = dir.end();
    for ( Directory::const_iterator itr = dir.begin(); itr != end; ++itr )
    {
        File file(*itr);
        if ( SolexaDir::isSolexaDir(file,flowCell) )
        {
            mMap.insert( pair<string const,Directory>(flowCell,Directory(file)) );
        }
    }
}

SolexaDir SolexaPipelineDir::find( string const& head ) const
{
    size_t posn = head.find('.');
    if ( posn == string::npos )
    {
        FatalErr("Can't find flowcell name in " << head);
    }

    string flowcell = head.substr(0,posn);
    typedef multimap<string,Directory>::const_iterator itr_type;
    pair<itr_type,itr_type> range(mMap.equal_range(head.substr(0,posn)));
    if ( range.first == range.second )
    {
        FatalErr("There is no solexa pipeline dir for flowcell " << flowcell);
    }

    Directory const& dir1 = range.first->second;
    ++range.first;
    if ( range.first == range.second )
    {
        return SolexaDir(dir1,head);
    }
    else
    {
        Directory const& dir2 = range.first->second;
        ++range.first;
        if ( range.first != range.second )
        {
            cerr << "There are more than 2 directories for flowcell " << flowcell <<
            ".\nWe're going to pretend it's a single-end run on the first directory, and ignore all the others.";
            return SolexaDir(dir1,head);
        }
        else if ( dir1 < dir2 )
        {
            return SolexaDir(dir1,dir2,head);
        }
        else
        {
            return SolexaDir(dir2,dir1,head);
        }
    }
}
