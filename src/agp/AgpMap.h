// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

#ifndef AGP_MAP
#define AGP_MAP

#include "String.h"

#include <set>
#include <map>

class AgpMapContig
{
  public:
    AgpMapContig()
        : mContigId( -1 ), 
          mStartCoord( -1 ), 
          mStopCoord( -1 ),
          mIsRc( false ) {}
    AgpMapContig( const int contigId, 
                  const int startCoord, const int stopCoord, const bool isRc ) 
        : mContigId( contigId ),
          mStartCoord( startCoord ),
          mStopCoord( stopCoord ),
          mIsRc( isRc ) {}
    
    bool IsValid()       const { return mContigId == -1; }

    int  GetContigId()   const { return mContigId; }
    int  GetStartCoord() const { return mStartCoord; }
    int  GetStopCoord()  const { return mStopCoord; }
    bool IsRc()          const { return mIsRc; }

    bool IsAfter ( const int coord ) const { return mStartCoord > coord; }
    bool IsBefore( const int coord ) const { return mStopCoord < coord; }
    bool Contains( const int coord ) const { return ! IsAfter(coord) && ! IsBefore(coord); }
        

    bool operator < ( const AgpMapContig &other ) const
    {
        return mStartCoord < other.mStartCoord; 
    }
    
  private:
    int mContigId;
    int mStartCoord;
    int mStopCoord;
    bool mIsRc;
};


class AgpMapBuilder
{
  public:
    AgpMapBuilder() {}

    void BuildFrom( const String &strAgpFile, 
                    std::map< String, std::set<AgpMapContig> > &mapChromToAgpMap );
};

#endif
