///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * DiffDiscovar.cc
 *
 *  Created on: May 9, 2013
 *      Author: tsharpe
 */
#include "MainTools.h"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

class Call
{
public:
    Call( std::string const& line );

    Call( Call&& that )
    : mLine(std::move(that.mLine)), mTabs(std::move(that.mTabs)),
      mID(that.mID), mStatus(that.mStatus) {}

    Call( Call const& that )
    { *this = that; }

    Call& operator=( Call const& that )
    { mLine = that.mLine; mTabs = that.mTabs;
      mID = that.mID; mStatus = that.mStatus;
      return *this; }

    std::string const& getLine() const { return mLine; }
    char getStatus() const { return mStatus; }
    void setStatus( char status )
    { mStatus = status; }

private:
    std::string mLine;
    std::vector<std::string::const_iterator> mTabs;
    unsigned mID;
    char mStatus = ' ';

    friend bool locusComparator( Call const& c1, Call const& c2 )
    { return std::lexicographical_compare(c1.mTabs[1]+1,c1.mTabs[2],
                                          c2.mTabs[1]+1,c2.mTabs[2]); }

    friend bool idComparator( Call const& c1, Call const& c2 )
    { return c1.mID < c2.mID; }

    friend bool sameEvent( Call const& c1, Call const& c2 )
    { using std::distance; using std::equal;
      if ( distance(c1.mTabs[2], c1.mTabs[4]) !=
              distance(c2.mTabs[2], c2.mTabs[4]) ) return false;
      if ( !equal(c1.mTabs[2], c1.mTabs[4], c2.mTabs[2]) ) return false;
      if ( c1.mTabs.size() != c2.mTabs.size() ) return false;
      if ( c1.mTabs.size() < 6 ) return true;
      if ( distance(c1.mTabs[5], c1.mLine.end()) !=
              distance(c2.mTabs[5], c2.mLine.end()) ) return false;
      return equal(c1.mTabs[5], c1.mLine.end(), c2.mTabs[5]); }

    friend void swap( Call& c1, Call& c2 )
    { using std::swap;
      swap(c1.mLine,c2.mLine); swap(c1.mTabs,c2.mTabs); swap(c1.mID,c2.mID); }
};

bool locusComparator( Call const& c1, Call const& c2 );
bool idComparator( Call const& c1, Call const& c2 );

Call::Call( std::string const& line )
: mLine(line)
{
    mTabs.reserve(5);
    std::string const& myLine = mLine;
    auto itr = myLine.begin(), end = myLine.end();
    while ( true )
    {
        itr = std::find(itr,end,'\t');
        if ( itr == end )
            break;
        mTabs.push_back(itr);
        ++itr;
    }
    if ( mTabs.size() < 5 )
        FatalErr("Garbled calls file.  Line is:\n" << line);

    std::istringstream iss(std::string(mTabs[0]+1,mTabs[1]));
    iss >> mID;
    if ( !iss )
        FatalErr("Garbled calls file.  Can't extract ID from:\n" << line);
}

void diff( std::vector<Call> const& exp, std::vector<Call> const& act,
            std::vector<Call>* pDiffs )
{
    std::vector<Call>& diffs = *pDiffs;
    auto eItr = exp.begin(), eEnd = exp.end();
    auto aItr = act.begin(), aEnd = act.end();
    while ( eItr != eEnd && aItr != aEnd )
    {
        if ( locusComparator(*eItr,*aItr) )
        {
            diffs.push_back(*eItr);
            diffs.back().setStatus('-');
            ++eItr;
        }
        else if ( locusComparator(*aItr,*eItr) )
        {
            diffs.push_back(*aItr);
            diffs.back().setStatus('+');
            ++aItr;
        }
        else if ( !sameEvent(*eItr,*aItr) )
        {
            diffs.push_back(*eItr);
            diffs.back().setStatus('-');
            diffs.push_back(*aItr);
            diffs.back().setStatus('+');
            ++eItr;
            ++aItr;
        }
        else
        {
            ++eItr;
            ++aItr;
        }
    }
    while ( eItr != eEnd )
    {
        diffs.push_back(*eItr);
        diffs.back().setStatus('-');
        ++eItr;
    }
    while ( aItr != aEnd )
    {
        diffs.push_back(*aItr);
        diffs.back().setStatus('+');
        ++aItr;
    }
}

int main( int argc, char** argv )
{
    RunTime();
    BeginCommandArgumentsNoHeader;
    CommandArgument_String_Doc(EXPECTED, "A calls file from a previous run.");
    CommandArgument_String_Doc(ACTUAL, "A calls file from a new run.");
    EndCommandArguments;

    std::string line;
    std::ifstream exp(EXPECTED.c_str());
    if ( !exp )
        FatalErr("Can't open file: " << EXPECTED);
    std::vector<Call> expectedCalls;
    while ( getline(exp,line) )
    {
        auto beg = line.begin();
        if ( line.size() > 3 && std::equal(beg,beg+4,"VAR\t") )
            expectedCalls.push_back(Call(line));
    }
    if ( !exp.eof() )
        FatalErr("Can't read file: " << EXPECTED);

    std::ifstream act(ACTUAL.c_str());
    if ( !act )
        FatalErr("Can't open file: " << ACTUAL);
    std::vector<Call> actualCalls;
    while ( getline(act,line) )
    {
        auto beg = line.begin();
        if ( line.size() > 3 && std::equal(beg,beg+4,"VAR\t") )
            actualCalls.push_back(Call(line));
    }
    if ( !act.eof() )
        FatalErr("Can't read file: " << ACTUAL);

    std::sort( expectedCalls.begin(), expectedCalls.end(), locusComparator );
    std::sort( actualCalls.begin(), actualCalls.end(), locusComparator );

    std::vector<Call> diffs;
    diff( expectedCalls, actualCalls, &diffs );
    std::sort( diffs.begin(), diffs.end(), idComparator );

    for ( auto itr=diffs.begin(), end=diffs.end(); itr != end; ++ itr )
        std::cout << itr->getStatus() << itr->getLine() << '\n';

    return !diffs.empty();
}
