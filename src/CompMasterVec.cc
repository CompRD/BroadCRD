/////////////////////////////////////////////////////////////////////////////
//                  SOFTWARE COPYRIGHT NOTICE AGREEMENT                    //
//      This software and its documentation are copyright (2009) by the    //
//  Broad Institute.  All rights are reserved.  This software is supplied  //
//  without any warranty or guaranteed support whatsoever. The Broad       //
//  Institute cannot be responsible for its use, misuse, or functionality. //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file CompMasterVec.cc
 * \author palvarez
 * \date Nov. 9, 2004
 */
#include "CompMasterVec.h"
#include "feudal/FeudalControlBlock.h"

bool IsGoodCompmastervecFile( const String& filename )
{
    unsigned long fileSize;
    FeudalControlBlock fcb(filename.c_str(),false,&fileSize);
    return fcb.isValid(filename.c_str(),fileSize,false) && fcb.isCompressed();
}

///Return number of unique objects in file
size_t CompmastervecFileUniqueObjectCount( const String& filename )
{
    unsigned long fileSize;
    FeudalControlBlock fcb(filename.c_str(),false,&fileSize);
    Assert(fcb.isValid(filename.c_str(),fileSize,false) && fcb.isCompressed());
    return fcb.getNElements();
}

///Return total size of unique dynamic data in the file
size_t CompmastervecFileRawCount( const String& filename )
{
    unsigned long fileSize;
    FeudalControlBlock fcb(filename.c_str(),false,&fileSize);
    Assert(fcb.isValid(filename.c_str(),fileSize,false) && fcb.isCompressed());
    Assert(fcb.getSizeofA());
    return fcb.getVarDataLen()/fcb.getSizeofA();
}
