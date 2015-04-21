///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ASSEMBLY_IDMANAGER
#define ASSEMBLY_IDMANAGER

class IdManager
{
  public:
    IdManager( );

    int GetMaxId( ) const;

    int GetNewId( );

    void UseId( const int id );

  private:
    int mMaxId;
};

inline
IdManager::IdManager( )
    : mMaxId( -1 )
{
}

inline
int
IdManager::GetMaxId( ) const
{
    return mMaxId;
}

inline
int
IdManager::GetNewId( )
{
    return ++mMaxId;
}

inline
void
IdManager::UseId( const int id )
{
    mMaxId = std::max( mMaxId, id );
}

#endif
