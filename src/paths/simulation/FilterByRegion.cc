///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#include <string>
#include "MainTools.h"

#include "paths/simulation/Regions.h"

//quick and dirty fasta reader
struct fasta_reader{
  fasta_reader(const std::string& sFileName)
      :inFile(sFileName)
      ,llBufferStart(-1)
      ,chromo_start(-1)
      ,sBuffer()
      ,sCurrentChromo(){
    ForceAssert( inFile.is_open() );
  };
  ~fasta_reader(){ inFile.close(); }


  // see if sChromo is logged in the fasta record
  // if yes, and sChromo is not sCurrentChromo, read the first line into buffer
  // if not in record, return an error
  int seek(const std::string& sChromo){
    ForceAssert( inFile.is_open() );
    std::cout << "seeking " << sChromo << std::endl;
    if(sChromo==sCurrentChromo) return 0;
    reset();
    while(  std::getline(inFile,sBuffer) ){
      if(sBuffer[0]=='>'){
        size_t uStart = sBuffer.find_first_not_of(" \t",1);
        size_t uEnd    = sBuffer.find_first_of(" \t",1);
        if( sChromo == sBuffer.substr(uStart,uEnd-uStart)){
          sCurrentChromo=sChromo;
          chromo_start = inFile.tellg();
          int info = nextline();
          if(info) return info;
          llBufferStart=0;
          return 0;
        }
      }
    }
    return 1;
  };

  //returns error code, 0 == no error
  int extract(std::string& sOut                              // output
             , const std::string& sChromo                    // chromosome name
             , const ssize_t llStart , const ssize_t llEnd){ // specifying [ llStart,llEnd )
    ForceAssert( inFile.is_open() );
    sOut.clear();

    int info = seek(sChromo);
    if(info) return info;

    // if llStart < llBufferStart, start reading record from the begining of the chromosome record
    if(llBufferStart>llStart){
        inFile.seekg(chromo_start);
        info = nextline();
        if(info) return info;
        llBufferStart=0;
    }

    //skipp all the lines until llStart is in range
    while(llBufferStart+ssize_t(sBuffer.size()) <= llStart){
        info = nextline();
        if(info) return info;
    }
    //read the records
    while(true){
        ssize_t locStart = std::max(ssize_t(0),llStart-llBufferStart);
        ssize_t locEnd   = llEnd   - llBufferStart;

        if(locEnd<=0) return 0;

        sOut.append(sBuffer ,locStart ,locEnd-locStart);

        if( size_t(locEnd) <= sBuffer.size()) return 0;

        info = nextline();
        if(info) return info;
    }
    return 0;
  };

  void reset(){
    std::cout << "reseting file pointer"<<std::endl;
    ForceAssert( inFile.is_open() );
    inFile.seekg(ios_base::beg);
    llBufferStart=-1;
    chromo_start=-1;
    sCurrentChromo.clear();
  };


  std::ifstream inFile;
  ssize_t llBufferStart;
  std::streampos chromo_start;
  std::string sBuffer;
  std::string sCurrentChromo;
private:
  int nextline(){
    llBufferStart += sBuffer.size();
    if(std::getline(inFile,sBuffer)){
      sBuffer.erase(remove(sBuffer.begin(),sBuffer.end(),' '),sBuffer.end());
      sBuffer.erase(remove(sBuffer.begin(),sBuffer.end(),'\t'),sBuffer.end());
      return 0;
    }
    else{
      return 2;
    }
  }
};

void FilterByRegion( const Regions::region_records&records // a list of records
                   , const std::string&sOutName            // output file's name
                   , const std::string&sFastaName          // input FASTA file's name
                   , const std::string&sInName             // input file's name
                   , const size_t nameColumn        =0    // input column for chromosome name
                   , const ssize_t startColumn      =-1    // input column for region start
                   , const ssize_t endColumn        =-1    // input column for region end, <= startColumn for Michael Mode
                   , const ssize_t llInputRangeShift=0     // shift the input's range before comparing
                   , const std::string sDelim=" \t"
                   ){
  //column tags in Michael's file
  enum MichaelENUM{mike_bin,mike_name,mike_chrom
                  ,mike_strand,mike_txStart,mike_txEnd
                  ,mike_cdsStart,mike_cdsEnd
                  ,mike_exonCount,mike_exonStarts,mike_exonEnds
                  ,mike_score,mike_name2
                  ,mike_cdsStartStat,mike_cdsEndStat,mike_exonFrames
                  ,mike_N};

  std::ifstream infile(sInName);
  std::ofstream outfile(sOutName);
  ForceAssert( infile.is_open() );
  ForceAssert( outfile.is_open() );


  fasta_reader fasta(sFastaName);
  std::string sFASTA;


  std::string sBuffer;

  // are we parsing Michael's format or not
  const bool bMichaelMode = !(startColumn>=0 && endColumn>startColumn);

  // number of columns to be parsed
  size_t nColumns =( (bMichaelMode)
                   ? (size_t(mike_N))
                   : (std::max( std::max(ssize_t(nameColumn),startColumn), endColumn)+1)
                   );

  std::vector<std::string> svFields(nColumns);

  //line-by-line
  while(std::getline(infile,sBuffer)){
      if( sBuffer.size()==0 ) continue;

      //skip over empty or commented line
      size_t uStart = sBuffer.find_first_not_of(sDelim,0);
      if( uStart == std::string::npos || sBuffer[uStart]=='#')
        continue;

      //parsing the needed number of columns
      for(size_t col=0;col<nColumns;++col){
          size_t uEnd = sBuffer.find_first_of(sDelim,uStart);
          svFields[col] = sBuffer.substr( uStart , uEnd-uStart);
          uStart = sBuffer.find_first_not_of(sDelim,uEnd+1);
      }

      // branch depending if we are doing michael's style or not
      if(!bMichaelMode){
          // check if chromosome and range is within the record, spit out the line if within range
          const Regions::Region rBuffer = records.find( svFields[nameColumn]
                                                      , atol(svFields[startColumn].c_str())+llInputRangeShift
                                                      , atol(svFields[endColumn].c_str())+llInputRangeShift
                                                      );
          if( rBuffer.isValid() ){
              outfile << sBuffer << "\n";
          }
          std::cerr << "WARNING: FASTA output has not been implemented for non michael mode" << std::endl;
      }
      else if( records.isInRecord(svFields[mike_chrom])){ //if we're in Mike mode
          // exon entries are in the same row
          const size_t nExon = atol( svFields[mike_exonCount].c_str());
          size_t uStartStart = svFields[mike_exonStarts].find_first_not_of(sDelim,0);
          size_t uEndStart   = svFields[mike_exonEnds]  .find_first_not_of(sDelim,0);
          //
          for(size_t nn=0 ; nn<nExon ;++nn){
            const size_t uStartEnd = svFields[mike_exonStarts].find_first_of(",",uStartStart);
            const size_t uEndEnd   = svFields[mike_exonEnds]  .find_first_of(",",uEndStart);

            const size_t uExonStart = atol(svFields[mike_exonStarts].substr(uStartStart,uStartEnd-uStartStart).c_str())+llInputRangeShift;
            const size_t uExonEnd   = atol(svFields[mike_exonEnds]  .substr(uEndStart,  uEndEnd  -uEndStart  ).c_str())+llInputRangeShift;
            const Regions::Region rBuffer = records.find( svFields[mike_chrom] , uExonStart , uExonEnd);

            if( rBuffer.isValid() ){
                outfile << ">"<<rBuffer.sTag<<"."<<svFields[mike_name]<<"."<<nn+1<<"\t\t#"
                        << "Fosmid="<<rBuffer.sTag << " "
                        << "Name="<<svFields[mike_name] << " "
                        << "Exon="<<nn+1 << " "
                        << "Total_Exon="<<nExon << " "
                        << "Chromo="<<svFields[mike_chrom]<< " "
                        << "Start="<<uExonStart << " "
                        << "End="<<uExonEnd
                        << std::endl;

                int info = fasta.extract(sFASTA,svFields[mike_chrom].substr(3,std::string::npos),uExonStart,uExonEnd);
                ForceAssert( info == 0);
                ForceAssert( sFASTA.size()==uExonEnd-uExonStart );
                outfile<<sFASTA<<"\n";
                sFASTA.clear();
            }
            uStartStart = uStartEnd + 1;
            uEndStart   = uEndEnd + 1;
          }
      }//if(!bMichaelMode)
  }

  infile.close();
  outfile.close();
};

void FilterVCFByRegion( const Regions::region_records&records // a list of records
                      , const std::string&sOutName            // output file's name
                      , const std::string&sInName             // input file's name
                   ){
  const std::string sDelim=" \t";

  std::ifstream infile(sInName);
  std::ofstream outfile(sOutName);
  ForceAssert( infile.is_open() );
  ForceAssert( outfile.is_open() );

  std::string sBuffer;

  //line-by-line header
  while(std::getline(infile,sBuffer)){
      if( sBuffer.size()==0 ) continue;
      size_t uStart = sBuffer.find_first_not_of(sDelim,0);
      if( sBuffer[uStart] != '#'){
          break;
      }
      outfile << sBuffer.substr(uStart,std::string::npos) << std::endl;;
  }

  //now the buffer has the first data line
  do{
      //skip over empty or commented line
      size_t uStart = sBuffer.find_first_not_of(sDelim,0);
      if( uStart == std::string::npos || sBuffer[uStart]=='#')
        continue;

      size_t uEnd = sBuffer.find_first_of(sDelim,uStart);
      const std::string sChromo = sBuffer.substr(uStart,uEnd-uStart);


      uStart = sBuffer.find_first_not_of(sDelim,uEnd+1);
      uEnd = sBuffer.find_first_of(sDelim,uStart);
      const size_t uPOS = atol( sBuffer.substr(uStart,uEnd-uStart).c_str() );

      const Regions::Region rBuffer = records.find( sChromo , uPOS , uPOS+1);
      if(rBuffer.isValid()){
          outfile << sBuffer << std::endl;;
      }

  }while(std::getline(infile,sBuffer));

  infile.close();
  outfile.close();
};

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;

     CommandArgument_String_Doc(OUTPUT_FILE, "output file name");
     CommandArgument_String_Doc(INPUT_FILE, "input file name");
     CommandArgument_String_Doc(REGION_FILE, "*.regions file name");
     CommandArgument_String_OrDefault_Doc(FASTA_FILE, "", "*.fasta file name, OMIT to filter VCF by regions");
     CommandArgument_Int_OrDefault_Doc(SKIPPED_REGION_COLUMNS,1,"length of padding region");

     CommandArgument_Int_OrDefault_Doc(NAME_COLUMN,0,"column number of name field (0 base)");
     CommandArgument_Int_OrDefault_Doc(START_COLUMN,-1,"(<0 for Michael mode) column number of the start field of [start,end) (0 base)");
     CommandArgument_Int_OrDefault_Doc(END_COLUMN,-1,"(<=START_COLUMN for Michael mode) column number of the end field [start,end) (0 base)");


     CommandArgument_Int_OrDefault_Doc(INPUT_RANGE_SHIFT,0,"shift of input files range before comparing to region");
     EndCommandArguments;

     if(FASTA_FILE!=""){
         if(START_COLUMN>=0 && END_COLUMN>START_COLUMN){
             std::cout << "Trimming:          " << INPUT_FILE << std::endl;
             std::cout << "with name column:  " << NAME_COLUMN << std::endl;
             std::cout << "with start column: " << START_COLUMN << std::endl;
             std::cout << "with end column:   " << END_COLUMN << std::endl;
         }
         else{
             std::cout << "Michael Mode       " << std::endl;
             std::cout << "Processing:        " << INPUT_FILE << std::endl;
             std::cout << "with FASTA ref:    " << FASTA_FILE << std::endl;
         }
     }
     else{
         std::cout << "Trimming VCF file:          " << INPUT_FILE << std::endl;
     }


     std::cout << "according to:      " << REGION_FILE << std::endl;
     std::cout << "in which first     " << SKIPPED_REGION_COLUMNS << " columns are omitted"<<std::endl;
     std::cout << "Writing to:        " << OUTPUT_FILE << std::endl;




     Regions::region_records records(REGION_FILE
                                   ,0
                                   ,SKIPPED_REGION_COLUMNS);
                                   // this is to skip a number of columns before parsing the region file
//     records.Print();
     if(FASTA_FILE!=""){
         FilterByRegion(records
                       ,OUTPUT_FILE
                       ,FASTA_FILE
                       ,INPUT_FILE);//,NAME_COLUMN,START_COLUMN,END_COLUMN ,INPUT_RANGE_SHIFT);
     }
     else{
         FilterVCFByRegion(records,OUTPUT_FILE,INPUT_FILE);
     }

};
