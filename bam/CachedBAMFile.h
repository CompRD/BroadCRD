///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file CachedBAMFile.h
 * \author tsharpe
 * \date Oct 26, 2012
 *
 * \brief
 */
#ifndef CACHEDBAMFILE_H_
#define CACHEDBAMFILE_H_

#include "String.h"

class CachedBAMFile
{
public:
    // if you know where it is
    CachedBAMFile( String const& filename )
    : mFilename(filename) {}

    // if you know its identifiers
    CachedBAMFile( String const& flowcell, String const& run, int lane,
                    String const& library )
    : mFilename(getPicardDir()+flowcell+'/'+run+'/'+ToString(lane)+'/'+library
                +'/'+flowcell+'.'+ToString(lane)
                +".aligned.duplicates_marked.bam") {}

    String const& getPath() const { return mFilename; }

    static String getPicardDir()
    { return "/wga/scr4/picard/"; }

private:
    String mFilename;
};

#endif /* CACHEDBAMFILE_H_ */
