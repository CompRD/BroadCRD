/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////
/*
 * \file Qltout2SAMUtils.h
 * \author tsharpe
 * \date Jan 6, 2009
 *
 * \brief Some utility classes to facilitate SAM creation.
 */
#ifndef QLTOUT2SAMUTILS_H_
#define QLTOUT2SAMUTILS_H_

#include "Basevector.h"
#include "Qualvector.h"
#include "system/ErrNo.h"
#include "system/file/File.h"
#include <algorithm>
#include <sstream>
#include <unistd.h>
#include <signal.h>

// some of the class names defined here (Block,SubSequence,Alignment) are common as dirt, so they're in a namespace
namespace qltout2sam
{
using std::string;

extern string NULL_STR;

template<class T> T istreamGet( istream& is, T val ) { is >> val; return val; }

/// An alignment block, as it appears in a qltout file:  gapSize, matchSize, mismatchCount
class Block
{
public:
    Block() : mGapSize(0), mMatchSize(0), mMismatchCount(0)
    {}

    Block( istream& is ) : mGapSize(istreamGet(is,mGapSize)), mMatchSize(istreamGet(is,mMatchSize)), mMismatchCount(istreamGet(is,mMismatchCount))
    {}

    // member-wise copying ok, no destructor needed

    int getGapSize() const
    { return mGapSize; }

    uint getMatchSize() const
    { return mMatchSize; }

    uint getMismatchCount() const
    { return mMismatchCount; }

    uint getEditDistance() const
    { return mMismatchCount + max(mGapSize,-mGapSize); }

    uint readLen() const
    { return mMatchSize + max(-mGapSize,0); }

    uint refLen() const
    { return mMatchSize + max(mGapSize,0); }

    /// make up a score for the block, based on the author's idea of a reasonable scoring matrix
    int getScore() const
    { return (mGapSize ? GAP_OPEN_SCORE + GAP_EXTEND_SCORE*(mGapSize-1) : 0) + MATCH_SCORE*(mMatchSize-mMismatchCount) + MISMATCH_SCORE*mMismatchCount; }

private:
    int mGapSize;
    uint mMatchSize;
    uint mMismatchCount;
    friend istream& operator>>( istream&, Block& );

    // scores
    static int const MATCH_SCORE = 4;
    static int const MISMATCH_SCORE = -6;
    static int const GAP_EXTEND_SCORE = -3;
    static int const GAP_OPEN_SCORE = -13;
};

inline istream& operator>>( istream& is, Block& block )
{
    is >> block.mGapSize;
    is >> block.mMatchSize;
    is >> block.mMismatchCount;
    return is;
}

/// a chunk of a larger sequence of bases, as it appears in a qltout file: id, start, end, len
class SubSequence
{
public:
    SubSequence() : mID(0), mStart(0), mEnd(0), mLen(0)
    {}

    SubSequence( istream& is ) : mID(istreamGet(is,mID)), mStart(istreamGet(is,mStart)), mEnd(istreamGet(is,mEnd)), mLen(istreamGet(is,mLen))
    {}

    // member-wise copying ok, no destructor needed

    uint getID() const
    { return mID; }

    uint getStart() const
    { return mStart; }

    uint getEnd() const
    { return mEnd; }

    uint getLen() const
    { return mLen; }

    bool isValid( string& reason, uint maxID ) const
    {
        bool result = false;
        if ( mID >= maxID )
        {
            reason = "ID exceeds max";
        }
        else if ( mStart >= mEnd )
        {
            reason = "start >= end";
        }
        else if ( mEnd > mLen )
        {
            reason = "end > len";
        }
        else
        {
            result = true;
        }
        return result;
    }

    uint alignmentLen() const
    { return mEnd - mStart; }

private:
    uint mID;
    uint mStart;
    uint mEnd;
    uint mLen;
    friend istream& operator>>( istream&, SubSequence& ss );
};

inline istream& operator>>( istream& is, SubSequence& ss )
{
    is >> ss.mID;
    is >> ss.mStart;
    is >> ss.mEnd;
    is >> ss.mLen;
    return is;
}

/// two subsequences, a set of blocks, and a flag indicating whether the subsequences are on the same or opposite strand
class Alignment
{
public:
    Alignment() : mSameStrand(false), mNBlocks(0), mBlocks(0)
    {}

    Alignment( istream& is ) : mRead(is), mSameStrand(!istreamGet(is,mSameStrand)), mRef(is), mNBlocks(istreamGet(is,mNBlocks)), mBlocks(0)
    {
        if ( is.fail() )
        {
            mNBlocks = 0;
        }
        Block* pBlock = allocBlocks();
        Block* pEnd = end();
        while ( pBlock != pEnd )
        {
            is >> *pBlock++;
        }
    }

    Alignment( Alignment const& alignment ) : mSameStrand(false), mNBlocks(0), mBlocks(0)
    { *this = alignment; }

    ~Alignment()
    { freeBlocks(); }

    Alignment& operator=( Alignment const& that )
    {
        if ( this != &that )
        {
            mRead = that.mRead;
            mSameStrand = that.mSameStrand;
            mRef = that.mRef;
            freeBlocks();
            mNBlocks = that.mNBlocks;
            allocBlocks();
            std::copy( that.begin(), that.end(), begin() );
        }
        return *this;
    }

    SubSequence& getRead()
    { return mRead; }

    SubSequence const& getRead() const
    { return mRead; }

    bool isSameStrand() const
    { return mSameStrand; }

    SubSequence& getRef()
    { return mRef; }

    SubSequence const& getRef() const
    { return mRef; }

    uint getNBlocks() const
    { return mNBlocks; }

    Block& getBlock( int idx )
    { return begin()[idx]; }

    Block const& getBlock( int idx ) const
    { return begin()[idx]; }

    Block* begin()
    { return canInline() ? reinterpret_cast<Block*>(&mBlocks) : mBlocks; }

    Block const* begin() const
    { return canInline() ? reinterpret_cast<Block const*>(&mBlocks) : mBlocks; }

    Block* end()
    { return begin() + mNBlocks; }

    Block const* end() const
    { return begin() + mNBlocks; }

    string getCigar() const
    {
        ostringstream oss;
        if ( mRead.getStart() )
        {
            oss << mRead.getStart() << 'S';
        }
        Block const* pEnd = end();
        Block const* pBlock = begin();
        while ( pBlock != pEnd )
        {
            Block const& block = *pBlock++;
            int gapSize = block.getGapSize();
            if ( gapSize > 0 )
            {
                oss << gapSize << 'D';
            }
            else if ( gapSize < 0 )
            {
                oss << -gapSize << 'I';
            }
            oss << block.getMatchSize() << 'M';
        }
        uint skip = mRead.getLen() - mRead.getEnd();
        if ( skip )
        {
            oss << skip << 'S';
        }
        return oss.str();
    }

    uint getEditDistance() const
    {
        return blockSum(&Block::getEditDistance);
    }

    bool isValid( string& reason, uint maxReadID, uint maxRefID ) const
    {
        bool result = false;
        if ( !mRead.isValid(reason,maxReadID) )
        {
            reason = "read " + reason;
        }
        else if ( !mRef.isValid(reason,maxRefID) )
        {
            reason = "ref " + reason;
        }
        else if ( readBlocksLen() != mRead.alignmentLen() )
        {
            reason = "length of blocks disagrees with read start and stop";
        }
        else if ( refBlocksLen() != mRef.alignmentLen() )
        {
            reason = "length of blocks disagrees with ref start and stop";
        }
        // note:  sometimes all the mismatches are recorded in a single, arbitrary block.
        // so the only acccuracy criterion we can impose is that the totals across all
        // blocks stand in the right relation
        else if ( mismatchCount() > matchSize() )
        {
            reason = "mismatch count exceeds match count";
        }
        else
        {
            result = true;
        }
        return result;
    }

    int getScore() const
    {
        int result = 0;
        Block const* pBlock = begin();
        Block const* pEnd = end();
        while ( pBlock != pEnd )
        {
            result += pBlock++->getScore();
        }
        return result;
    }

private:
    bool canInline() const
    { return mNBlocks <= sizeof(mBlocks)/sizeof(Block); }

    Block* allocBlocks()
    {
        Block* result = reinterpret_cast<Block*>(&mBlocks);
        if ( !canInline() )
            result = mBlocks = new Block[mNBlocks];
        return result;
    }

    void freeBlocks()
    { if ( !canInline() ) delete [] mBlocks; }

    uint blockSum( uint (Block::*func)() const ) const
    {
        uint result = 0;
        Block const* pBlock = begin();
        Block const* pEnd = end();
        while ( pBlock != pEnd )
        {
            result += (pBlock++->*func)();
        }
        return result;
    }

    uint readBlocksLen() const
    {
        return blockSum(&Block::readLen);
    }

    uint refBlocksLen() const
    {
        return blockSum(&Block::refLen);
    }

    uint mismatchCount() const
    {
        return blockSum(&Block::getMismatchCount);
    }

    uint matchSize() const
    {
        return blockSum(&Block::getMatchSize);
    }

    // even though it's not compact, these members have to be declared in this order
    // so that initialization from a qltout QUERY line will happen correctly
    SubSequence mRead;
    bool mSameStrand;
    SubSequence mRef;
    uint mNBlocks;
    Block* mBlocks;

    friend istream& operator>>( istream&, Alignment& );
};

inline istream& operator>>( istream& is, Alignment& alignment )
{
    alignment.freeBlocks();
    new (&alignment) Alignment(is);
    return is;
}

/// reads the next alignment from the stream, and stashes it in the alignment reference
/// returns true if successful
bool nextAlignment( istream& is, Alignment& alignment, uint nReads, uint nRefs );

/// a read:  id, seqs, and quals.
/// this base class represents an unmapped read, so there is no alignment data.
class Read
{
public:
    Read( uint readID, basevector const& bvec, qvec const& qvec, bool isPF )
        : mReadID(readID), mBVec(bvec), mQVec(qvec), mFlag(NOT_MAPPED|(isPF?0:NOT_PF)), mpMate(0)
    {}

    // member-wise copying ok, no destructor needed

    void setMate( Read const* pMate, bool isFirst )
    { if ( mReadID != pMate->getReadID() )
      { FatalErr("Trying to pair reads with different IDs: " << mReadID << " vs. " << pMate->getReadID() ); }

      mpMate = pMate;

      mFlag |= IS_PAIRED;
      if ( isFirst )
          mFlag |= FIRST_READ;
      else
          mFlag |= SECOND_READ;

      if ( !pMate->isMapped() )
        mFlag |= MATE_NOT_MAPPED;
      else if ( isMapped() && pMate->getRefID() == getRefID() )
        mFlag |= IS_PROPERLY_PAIRED;

      if ( pMate->isReversed() )
        mFlag |= MATE_IS_REVERSED;
    }

    virtual bool isValid( string& reason ) const
    { bool result = true;
      if ( mBVec.size() && mQVec.size() && mBVec.size() != static_cast<uint>(mQVec.size()) )
      { reason = "bases and quals have different lengths";
        result = false; }
      return result;
    }

    virtual bool isMapped() const
    { return false; }

    bool isReversed() const
    { return mFlag & IS_REVERSED; }

    bool isProperlyPaired() const
    { return mFlag & IS_PROPERLY_PAIRED; }

    uint getReadID() const
    { return mReadID; }

    int getFlag() const
    { return mFlag; }

    virtual uint getRefID() const
    { return gUnmappedRefID; }

    virtual uint getRefPosition() const
    { return static_cast<uint>(-1); }

    virtual uint get5PrimeRefPosition() const
    { return static_cast<uint>(-1); }

    virtual uint getMapQ() const
    { return 0; }

    virtual string getCigar() const
    { return "*"; }

    Read const* getMate() const
    { return mpMate; }

    virtual int getISize() const
    { return 0; }

    virtual uint getEditDistance() const
    { return getReadLen(); }

    basevector const& getBVec() const
    { return mBVec; }

    qvec const& getQVec() const
    { return mQVec; }

    uint getReadLen() const
    { return mBVec.size(); }


    static void setUnmappedRefID( uint refID )
    { gUnmappedRefID = refID; }

    // values for the SAM flag field
    static int const IS_PAIRED = 0x0001;
    static int const IS_PROPERLY_PAIRED = 0x0002;
    static int const NOT_MAPPED = 0x0004;
    static int const MATE_NOT_MAPPED = 0x0008;
    static int const IS_REVERSED = 0x0010;
    static int const MATE_IS_REVERSED = 0x0020;
    static int const FIRST_READ = 0x0040;
    static int const SECOND_READ = 0x0080;
    static int const NOT_PRIMARY = 0x0100;
    static int const NOT_PF = 0x0200;

protected:
    Read( uint readID, basevector const& bvec, qvec const& qvec, int flag )
        : mReadID(readID), mBVec(bvec), mQVec(qvec), mFlag(flag), mpMate(0)
    {}

private:
    uint mReadID;
    basevector const& mBVec;
    qvec const& mQVec;
    int mFlag;
    Read const* mpMate;

    static uint gUnmappedRefID;
};

/// transforms an iterator for bases (encoded in a small integer, as we usually do it) into an iterator of characters
template<class Itr> class BaseAsCharIterator
    : public std::iterator<std::input_iterator_tag,char,void,void,char>
{
public:
    BaseAsCharIterator( Itr itr ) : mItr(itr)
    {}

    // member-wise copying ok, no destructor needed

    char operator*()
    { return "ACGT"[*mItr]; }

    bool operator==( BaseAsCharIterator<Itr> const& that )
    { return mItr == that.mItr; }

    bool operator!=( BaseAsCharIterator<Itr> const& that )
    { return mItr != that.mItr; }

    BaseAsCharIterator<Itr>& operator++()
    { ++mItr; return *this; }

    BaseAsCharIterator<Itr> operator++(int)
    { BaseAsCharIterator<Itr> tmp(*this); ++mItr; return tmp; }

private:
    Itr mItr;
};

/// transforms an iterator for quality scores (encoded in a small integer, as we usually do it) into an iterator of characters
/// the character-encoding is specific to SAM files:  it's the quality score + 33.
template<class Itr> class QualAsCharIterator
    : public std::iterator<std::forward_iterator_tag,char,void,void,char>
{
public:
    QualAsCharIterator( Itr itr ) : mItr(itr)
    {}

    // member-wise copying ok, no destructor needed

    char operator*()
    { using std::min; // max legal qual score is 60
      return min(*mItr,static_cast<unsigned char>(60))+33; }

    bool operator==( QualAsCharIterator<Itr> const& that )
    { return mItr == that.mItr; }

    bool operator!=( QualAsCharIterator<Itr> const& that )
    { return mItr != that.mItr; }

    QualAsCharIterator<Itr>& operator++()
    { ++mItr; return *this; }

    QualAsCharIterator<Itr> operator++(int)
    { QualAsCharIterator<Itr> tmp(*this); ++mItr; return tmp; }

private:
    Itr mItr;
};

inline ostream& operator<<( ostream& os, Read const& read )
{
    Read const* pMate = read.getMate();

    //        QNAME                       FLAG
    os << read.getReadID() << '\t' << read.getFlag() << '\t';

    if ( !read.isMapped() && pMate && pMate->isMapped() )
        //        RNAME                       POS
        os << pMate->getRefID() << '\t' << pMate->getRefPosition()+1;
    else
        //        RNAME                       POS
        os << read.getRefID() << '\t' << read.getRefPosition()+1;

    //        MAPQ                        CIGAR
    os << '\t' << read.getMapQ() << '\t' << read.getCigar() << '\t';

    if ( !pMate || !(read.isMapped() || pMate->isMapped()) ) // unpaired or both unmapped
        //    MRNM          MPOS
        os << '*' << '\t' << 0;
    else if ( read.isMapped() && !pMate->isMapped() ) // mate-only unmapped
        //    MRNM                 MPOS
        os << '=' << '\t' << read.getRefPosition()+1;
    else if ( !read.isMapped() && pMate->isMapped() ) // read-only unmapped
        //    MRNM                 MPOS
        os << '=' << '\t' << pMate->getRefPosition()+1;
    else if ( read.getRefID() == pMate->getRefID() ) // both mapped to same ref
        //    MRNM                 MPOS
        os << '=' << '\t' << pMate->getRefPosition()+1;
    else                                             // both mapped, but to different refs
        //        MRNM                            MPOS
        os << pMate->getRefID() << '\t' << pMate->getRefPosition()+1;

    //              ISIZE
    os << '\t' << read.getISize() << '\t';

    std::ostream_iterator<char> osi(os);

    // SEQ
    basevector const& vb = read.getBVec();
    if ( !vb.size() )
    {
        os << '*';
    }
    else if ( read.isReversed() )
    {
        copy( BaseAsCharIterator<basevector::const_reverse_complement_iterator>(vb.RCBegin()),
                BaseAsCharIterator<basevector::const_reverse_complement_iterator>(vb.RCEnd()),
                osi );
    }
    else
    {
        std::copy( BaseAsCharIterator<basevector::const_iterator>(vb.Begin()),
                    BaseAsCharIterator<basevector::const_iterator>(vb.End()),
                    osi );
    }
    os << '\t';

    // QUAL
    qvec const& vq = read.getQVec();
    if ( !vq.size() )
    {
        os << '*';
    }
    else if ( read.isReversed() )
    {
        copy( QualAsCharIterator<qvec::const_reverse_iterator>(vq.rbegin()),
                QualAsCharIterator<qvec::const_reverse_iterator>(vq.rend()),
                osi );
    }
    else
    {
        std::copy( QualAsCharIterator<qvec::const_iterator>(vq.begin()),
                    QualAsCharIterator<qvec::const_iterator>(vq.end()),
                    osi );
    }

    return os;
}

/// a read, supplemented by alignment data
class MappedRead : public Read
{
public:
    MappedRead( basevector const& bvec, qvec const& qvec, Alignment& alignment, int mapq, bool isPrimary, bool isPF )
     : Read(alignment.getRead().getID(),bvec,qvec,
             (alignment.isSameStrand()?0:IS_REVERSED)|(isPrimary?0:NOT_PRIMARY)|(isPF?0:NOT_PF)),
       mAlignment(alignment), mMapQ(mapq)
    {}

    // member-wise copying ok, no destructor needed

    virtual bool isValid( string& reason ) const
    { bool result = Read::isValid(reason);
      if ( result && getReadLen() && getReadLen() != mAlignment.getRead().getLen() )
      { reason = "read length in alignment disagrees with length of bases and quals";
        result = false; }
      return result; }

    virtual bool isMapped() const
    { return true; }

    virtual uint getRefID() const
    { return mAlignment.getRef().getID(); }

    virtual uint getRefPosition() const
    { return mAlignment.getRef().getStart(); }

    virtual uint get5PrimeRefPosition() const
    { return isReversed() ? mAlignment.getRef().getEnd()-1 : mAlignment.getRef().getStart(); }

    virtual uint getMapQ() const
    { return mMapQ; }

    virtual string getCigar() const
    { return mAlignment.getCigar(); }

    virtual int getISize() const
    { return isProperlyPaired() ? (getMate()->get5PrimeRefPosition()-get5PrimeRefPosition()) : 0; }

    virtual uint getEditDistance() const
    { return mAlignment.getEditDistance(); }

    static uint const MAPQ_MAX = 255;

private:
    Alignment const& mAlignment;
    uint mMapQ;
};

/// constructor creates a named pipe, destructor removes it
class NamedPipe
{
public:
    NamedPipe() : mName(synthesizeName())
    { if ( mkfifo(mName.c_str(),0600) == -1 )
      { ErrNo err;
        FatalErr("Can't make named pipe" << err); } }

    ~NamedPipe()
    { unlink(mName.c_str()); }

    string const& getName()
    { return mName; }

private:
    // no copying
    NamedPipe( NamedPipe const& );
    NamedPipe& operator=( NamedPipe const& );

    string synthesizeName()
    {
        stringstream ss;
        ss << "/tmp/qltout2SAM." << getpid() << '.' << getuid() << '.' << ++gNPipes;
        return ss.str();
    }

    string mName;
    static int gNPipes;
};

void linkRefDict( File const& refDict, File const& refFasta );

/// A couple of little goofy utils to track child processes
class Executor
{
public:
    Executor( string command ) : mCommand(command)
    {
        EXEC_ARGS[2] = command.c_str();
        pid_t pid = vfork();
        if ( !pid )
        {
            execv(EXEC_ARGS[0],const_cast<char*const*>(EXEC_ARGS));
            _exit(errno);
        }

        if ( pid == -1 )
        { ErrNo err;
          FatalErr("Unable to vfork for " << command << err); }

        mPid = pid;
    }

    void kill( int sigNo )
    {
        ::kill(mPid,sigNo);
    }

    // if the exit status is associated with our pid, report on it
    // result is whether the exit was abnormal
    bool reportStatus( pid_t pid, int status )
    {
        bool result = false;
        if ( pid == mPid )
        {
            if ( WIFSIGNALED(status) )
            {
                cout << mCommand << " died on signal " << WTERMSIG(status) << endl;
                result = true;
            }
            else if ( WIFEXITED(status) )
            {
                if ( WEXITSTATUS(status) )
                {
                    cout << mCommand << " exited with status " << WEXITSTATUS(status) << endl;
                    result = true;
                }
            }
        }
        return result;
    }

private:
    // no copying
    Executor( Executor const& );
    Executor& operator=( Executor const& );

    string mCommand;
    pid_t mPid;
    static char const* EXEC_ARGS[];
};

void execChild( string command );
void signalCatcher( int sigNo, siginfo_t* info, void* context );
void waitForChildren( string exitMessage );

} // end of namespace qltout2sam

#endif /* QLTOUT2SAMUTILS_H_ */
