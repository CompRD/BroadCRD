///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * CombineSAMs.cc
 *
 *  Created on: Jun 26, 2014
 *      Author: tsharpe
 */
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace
{

inline void tabTokenize( std::string const& str, std::vector<std::string>& vec )
{
    vec.clear();
    auto beg = str.begin();
    auto end = str.end();
    while ( beg != end )
    {
        auto itr = std::find(beg,end,'\t');
        vec.emplace_back(beg,itr);
        if ( itr != end ) ++itr;
        beg = itr;
    }
}

inline int toInt( std::string const& str )
{
    size_t pos;
    int result = std::stoi(str,&pos);
    if ( pos != str.size() )
    {
        std::cerr << "Failed to convert int.\n";
        exit(1);
    }
    return result;
}

}

int main( int argc, char** argv )
{
    std::ifstream if1(argv[1]);
    std::ifstream if2(argv[2]);
    std::string rec1;
    std::string rec2;
    std::vector<std::string> tokens1;
    std::vector<std::string> tokens2;
    while ( getline(if1,rec1) )
    {
        if ( !getline(if2,rec2) )
        {
            std::cerr << "unpaired line at end of file 1\n";
            exit(1);
        }
        if ( rec1.empty() )
        {
            std::cerr << "read empty line from file 1\n";
            exit(1);
        }
        if ( rec2.empty() )
        {
            std::cerr << "read empty line from file 2\n";
            exit(1);
        }
        if ( rec1.front() == '@' )
        {
            if ( rec1 != rec2 )
            {
                std::cerr << "header lines don't match\n";
                exit(1);
            }
            std::cout << rec1 << '\n';
            continue;
        }

        tabTokenize(rec1,tokens1);
        tabTokenize(rec2,tokens2);
        if ( tokens1.size() < 11 )
        {
            std::cerr << "not enough tokens in line from file 1\n";
            exit(1);
        }
        if ( tokens2.size() < 11 )
        {
            std::cerr << "not enough tokens in line from file 2\n";
            exit(1);
        }
        if ( tokens1.front() != tokens2.front() )
        {
            std::cerr << "QNAMES out of sync\n";
            exit(1);
        }
        int flag1 = toInt(tokens1[1]);
        if ( flag1 & 0x4 ) // UNMAPPED
            continue;
        int flag2 = toInt(tokens2[1]);
        if ( flag2 & 0x4 ) // UNMAPPED
            continue;
        tokens1[6] = tokens2[2];
        tokens1[7] = tokens2[3];
        tokens2[6] = tokens1[2];
        tokens2[7] = tokens1[3];
        if ( flag1 & 0x10 ) // RC
            flag2 |= 0x20; // Mate is RC
        if ( flag2 & 0x10 )
            flag1 |= 0x20;
        flag1 |= 0x41; // first segment in pair
        flag2 |= 0x81; // last segment in pair
        tokens1[1] = std::to_string(flag1);
        tokens2[1] = std::to_string(flag2);

        std::cout << tokens1[0];
        for ( auto itr=tokens1.begin()+1, end=tokens1.end(); itr != end; ++itr )
            std::cout << '\t' << *itr;
        std::cout << '\n';
        std::cout << tokens2[0];
        for ( auto itr=tokens2.begin()+1, end=tokens2.end(); itr != end; ++itr )
            std::cout << '\t' << *itr;
        std::cout << '\n';
    }
    if ( getline(if2,rec2) )
    {
        std::cerr << "unpaired line at end of file 2\n";
        exit(1);
    }
}
