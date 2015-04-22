///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include "system/Assert.h"
#include "paths/simulation/FloodN.h"

namespace FloodN{

//void region_records::flooded_FASTA(const std::string& sOutName ,const std::string& sInName){
void flooded_FASTA(const Regions::region_records& records, const std::string& sOutName ,const std::string& sInName){
    std::ifstream infile(sInName);
    std::ofstream outfile(sOutName);
    ForceAssert( infile.is_open() );
    ForceAssert( outfile.is_open() );


    std::string sBuffer;     // read buffer
    size_t uNProcessedCHR=0; // number of chromosome (with flooding specification) processed

    if(std::getline(infile,sBuffer)){// score one line first
        outfile << sBuffer << "\n";
    }
    else{
        infile.close();
        outfile.close();
        return;
    }

    while(infile.good()&&sBuffer.size()==0){//fast forward over all blank lines
        std::getline(infile,sBuffer);
        outfile << sBuffer << "\n";
    }

    std::string sOut;     // output buffer
    while(infile.good()){//&&uNProcessedCHR<records.size()){
        // chromosome
        ForceAssert( sBuffer[0]=='>');
        size_t uStart = sBuffer.find_first_not_of(" \t",1);
        size_t uEnd = sBuffer.find_first_of(" \t\n",uStart);
        std::string sCHR(sBuffer,uStart,uEnd-uStart);

        const Regions::region_records::const_iterator itr_chr = records.find(sCHR);
        if(itr_chr==records.end()){
            sOut.assign(sBuffer.size(),'N');
            while(std::getline(infile,sBuffer)){
                if(sBuffer[0]=='>'){
                    outfile << sBuffer << "\n";
                    break;
                }
                else{
                    //avoiding space and tabs
                    sBuffer.erase(remove(sBuffer.begin(),sBuffer.end(),' '),sBuffer.end());
                    sBuffer.erase(remove(sBuffer.begin(),sBuffer.end(),'\t'),sBuffer.end());
                    if(sOut.size()!=sBuffer.size()) sOut.assign(sBuffer.size(),'N');
                    outfile << sOut << "\n";
                }
            }
        }
        else{
            //if there is flooding instruction
            auto itr = itr_chr->second.begin();// go through all flooding range, assume sorted and merged
            auto end = itr_chr->second.end();

            ssize_t llShift = 0;// number of basis processed so far
            while(std::getline(infile,sBuffer)){
                if(sBuffer[0]=='>'){ // end of chromosome record
                    outfile << sBuffer << "\n";
                    itr=itr_chr->second.end(); // this line should be redundant now, but just to safe guard against future bugs
                    break;
                }
                else{// if not end of chromosome record
                    //compact the string, if only regex could work properly...
                    sBuffer.erase(remove(sBuffer.begin(),sBuffer.end(),' '),sBuffer.end());
                    sBuffer.erase(remove(sBuffer.begin(),sBuffer.end(),'\t'),sBuffer.end());

                    sOut.assign(sBuffer.size(),'N');

                    if( itr!=end){// if there are still special range unchecked

                        ssize_t llFloodStart = std::max(ssize_t(0),itr->llStart - llShift);
                        while(itr!=end && llFloodStart<ssize_t(sBuffer.size())){
                            ssize_t llFloodEnd = itr->llEnd - llShift;
//                            std::fill( sBuffer.begin()+llFloodStart
//                                     , sBuffer.begin()+std::min(ssize_t(sBuffer.size())
//                                                               ,llFloodEnd
//                                                               )
//                                     , 'N');
                            std::copy( sBuffer.begin()+llFloodStart
                                     , sBuffer.begin()+std::min(ssize_t(sBuffer.size())
                                                               ,llFloodEnd
                                                               )
                                     , sOut.begin()+llFloodStart);
                            if( llFloodEnd > ssize_t(sBuffer.size()) ){
                                break;
                            }
                            else{
                                ++itr;
                                llFloodStart = std::max(ssize_t(0),itr->llStart - llShift);
                            }
                        }//while(itr!=end && llFloodStart<sBuffer.size())
                        llShift += sBuffer.size();
                    }//if itr!=end

                    outfile << sOut <<"\n";
                }//if not end of chromosome record
            }//while getline
            ++uNProcessedCHR;
        }//if there is flooding
    }//while(itr!=end && llFloodStart<sBuffer.size())


    // copy over the rest of the files
//    while(std::getline(infile,sBuffer)){ outfile << sBuffer << "\n"; }

    infile.close();
    outfile.close();
}


}//namespace FloodN
