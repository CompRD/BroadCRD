// Copyright (c) 2004 Broad Institute / Massachusetts Institute of Technology

#ifndef MINIMAL_ALIGN_H
#define MINIMAL_ALIGN_H

#include <stdlib.h>
#include <math.h>

#include "CoreTools.h"
#include "math/Functions.h"
#include "FastIfstream.h"

struct pos_score
{
    pos_score(const int target_arg,
	      const int beg_arg,
	      const int end_arg,
	      const int rc1_arg,
	      const float score_arg)
      : targetID(target_arg),
	targetBeg(beg_arg),
	targetEnd(end_arg),
	rc(rc1_arg),
	score(score_arg) 
    {}

    pos_score()
      : targetID(0), targetBeg(0), targetEnd(0), rc(0), score(0.0) 
    {}

    friend 
    bool operator < (const pos_score & r1, const pos_score & r2 )
    { 
        if ( r1.targetID < r2.targetID ) return True; 
        if ( r1.targetID > r2.targetID ) return False; 
        if ( r1.targetBeg < r2.targetBeg ) return True; 
        if ( r1.targetBeg > r2.targetBeg ) return False; 
        if ( r1.targetEnd < r2.targetEnd ) return True; 
        if ( r1.targetEnd > r2.targetEnd) return False; 
        return False;
    }

    friend
    istream& operator>>(istream& in, pos_score & r )
    {
        return in >> r.targetID >> r.targetBeg >> r.targetEnd >> r.rc >> r.score;    
    } 

    friend
    ostream& operator<<(ostream& out, pos_score & r )
    {
        return out << r.targetID << " " <<  r.targetBeg << " " 
		   << r.targetEnd << " " << r.rc << " " << r.score << endl;    
    } 

    int targetID;
    int targetBeg;
    int targetEnd;
    int rc;
    float score;
};

class MinimalAlign 
{
  public:

    MinimalAlign(const int queryID, const Bool uniquelyPlaced, const vec<pos_score> & alignment):
        m_aligns(alignment), m_readID(queryID), m_uniquelyPlaced(uniquelyPlaced)
    { }

    MinimalAlign():
        m_readID(-1),
	m_uniquelyPlaced(True)
    { }

    void SetAligns(const vec<pos_score> & alignment)
    { m_aligns = alignment; }
    void SetReadID(const int readID )
    { m_readID = readID; }

    vec<pos_score> Aligns() const { return m_aligns; }
    int readID() const {return m_readID; }

    unsigned int NumAligns() const { return m_aligns.size(); }
    Bool UniquelyPlaced() const { return m_uniquelyPlaced; }

    void sortAligns()
    {
        if ( !m_aligns.empty() )
            sort(m_aligns.begin(), m_aligns.end());  
    }

    friend 
    bool operator < (const MinimalAlign & r1, const MinimalAlign & r2 )
    { 
        if ( r1.m_readID < r2.m_readID ) return True; 
        if ( r1.m_readID > r2.m_readID ) return False; 
        return False;
    }
   
    friend 
    bool operator == ( const MinimalAlign & r1, const MinimalAlign & r2 )
    {
        if ( r1.m_readID != r2.m_readID ) return False;
        if ( r1.m_uniquelyPlaced != r2.m_uniquelyPlaced ) return False;
        if ( r1.m_aligns.size() != r2.m_aligns.size() ) return False;
        for ( unsigned int i=0; i < r1.m_aligns.size(); ++i )
        {
            if ( r1.m_aligns[i].targetID != r2.m_aligns[i].targetID ) return False;
            if ( r1.m_aligns[i].targetBeg != r2.m_aligns[i].targetBeg ) return False;
            if ( r1.m_aligns[i].targetEnd != r2.m_aligns[i].targetEnd ) return False;
            if ( r1.m_aligns[i].rc != r2.m_aligns[i].rc ) return False;
            if ( r1.m_aligns[i].score != r2.m_aligns[i].score ) return False;
        }
        return True;
    }

    friend
    ostream& operator<<(ostream& out, const MinimalAlign& sa)
    {
        BinWrite( out, sa.m_readID );
	BinWrite( out, sa.m_uniquelyPlaced );
        unsigned int num = sa.m_aligns.size();
        BinWrite( out, num );
        for ( unsigned int i=0; i < sa.m_aligns.size(); ++i )
        {
            BinWrite(out, sa.m_aligns[i].targetID);
            BinWrite(out, sa.m_aligns[i].targetBeg);
            BinWrite(out, sa.m_aligns[i].targetEnd);
            BinWrite(out, sa.m_aligns[i].rc);
            BinWrite(out , sa.m_aligns[i].score);
        }

        return out;
    }

    friend
    istream& operator>>(istream& in, MinimalAlign& sa)
    {
        BinRead( in, sa.m_readID );
	BinRead( in, sa.m_uniquelyPlaced );
        unsigned int num;
        BinRead( in, num );
        sa.m_aligns.resize(num);
        for ( unsigned int i=0; i < sa.m_aligns.size(); ++i )
        {
            BinRead(in, sa.m_aligns[i].targetID);
            BinRead(in, sa.m_aligns[i].targetBeg);
            BinRead(in, sa.m_aligns[i].targetEnd);
            BinRead(in, sa.m_aligns[i].rc);
            BinRead(in , sa.m_aligns[i].score);
        }

        return in;
    }  

    void Print(ostream& out)
    {
        for ( unsigned int i=0; i < m_aligns.size(); ++i )
        {
            out << "QUERY " << m_readID << "\t"
                << m_aligns[i].targetID << "\t"
                << m_aligns[i].targetBeg << "\t"
                << m_aligns[i].targetEnd << "\t"
            	<< m_aligns[i].rc << "\t"
                << m_aligns[i].score << endl;
        }
    }

    void Print(ostream& out, const unsigned int i )
    {
        out << "QUERY " << m_readID << "\t"
            << m_aligns[i].targetID << "\t"
            << m_aligns[i].targetBeg << "\t"
            << m_aligns[i].targetEnd << "\t"
            << m_aligns[i].rc << "\t"
            << m_aligns[i].score << endl;
    }

    void
    Read( const String& in )
    {    
        istrstream ins( in.c_str( ) );
        String field;
        ins >> field;
        ForceAssert( field == "QUERY" );
        unsigned int numAligns;
        ins >> m_readID;
        m_aligns.resize(1);
        ins >> m_aligns[0].targetID >> m_aligns[0].targetBeg >> m_aligns[0].targetEnd
            >> m_aligns[0].rc >> m_aligns[0].score;
    }

    // TODO: this shouldn't be public but since to remove
    // items, made public for easy access (less coding).
    vec<pos_score> m_aligns;

    private:
    
    int m_readID;
    Bool m_uniquelyPlaced;
};

class FlatMinimalAlign 
{
  public:

    FlatMinimalAlign(const int readID, const int targetID, const int targetBeg,
                     const int targetEnd, const int rc)
        : m_readID(readID), m_targetID(targetID), m_targetBeg(targetBeg),
          m_targetEnd(targetEnd), m_rc(rc)
    {}
    FlatMinimalAlign()
        : m_readID(0), m_targetID(0), m_targetBeg(0), m_targetEnd(0), m_rc(0)
    {}

    int ReadID()const { return m_readID; }
    int TargetID() const { return m_targetID; }
    int TargetBeg() const { return m_targetBeg; }
    int TargetEnd() const { return m_targetEnd; }
    int RC() const { return m_rc; }

    friend 
    bool operator < (const FlatMinimalAlign & r1, const FlatMinimalAlign & r2 )
    { 
        if ( r1.m_targetID < r2.m_targetID ) return True; 
        if ( r1.m_targetID > r2.m_targetID ) return False; 
        if ( r1.m_targetBeg < r2.m_targetBeg ) return True; 
        if ( r1.m_targetBeg > r2.m_targetBeg ) return False; 
        if ( r1.m_targetEnd < r2.m_targetEnd ) return True; 
        if ( r1.m_targetEnd > r2.m_targetEnd ) return False;
        if ( r1.m_readID < r1.m_readID ) return True;
        if ( r1.m_readID < r1.m_readID ) return False;
        return False;
    }

    friend
    ostream& operator<<(ostream& out, const FlatMinimalAlign& sa)
    {
        return  out << "QUERY " << sa.m_readID << "\t"
            << sa.m_targetID << "\t"
            << sa.m_targetBeg<< "\t"
            << sa.m_targetEnd << "\t"
            << sa.m_rc << endl;
    }

    void Print( ostream& out )
    {
        out << "QUERY " << m_readID << "\t"
            << m_targetID << "\t"
            << m_targetBeg<< "\t"
            << m_targetEnd << "\t"
            << m_rc << endl;
    }
        
    void
    Read( const String& in )
    {    
        istrstream ins( in.c_str( ) );
        String field;
        ins >> field;
        ForceAssert( field == "QUERY" );
        ins >> m_readID >> m_targetID >> m_targetBeg >> m_targetEnd >> m_rc;
    }

  private:

    int m_readID;
    int m_targetID;
    int m_targetBeg;
    int m_targetEnd;
    int m_rc;
};

class AssemblyAligns
{
  public:
    AssemblyAligns(){}

    void Read(const String & fileName);
    vec<FlatMinimalAlign> & GetAligns()
    {return m_minAligns;}

  private:
    vec<FlatMinimalAlign> m_minAligns;
};

#endif
