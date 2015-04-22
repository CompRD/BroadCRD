/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/* 
   Program: CorrectReadsByKmers

   Given a <.fastb> file of reads, and the table T of strong (trusted)  k-mers in
   them (together with their multiplicities, as produced by
   <FindKmerFrequencies>), edit the reads as follows.
   
   Scan each read R, starting at the left, for strong kmers. Mark base positions that
   belong to strong k-mers as strong. If no strong k-mers can be found within the original 
   read then try to uncover them by mutating the bases in a following way. Scan the read from 
   left to right and try 3k possible mutations for each k-mer. If only one such mutation
   produces a strong k-mer then mark the corresponding base positions as strong. 

   If the original or mutated read contains strong k-mers then take the leftmost one 
   (under the assumption that the beginning of the read has a higher quality) and 
   try to extend that kmer by including (in turn) the neighboring positions to the left 
   and right. If the original base at those positions does not form an adjacent strong kmer
   and was not previously marked as strong then mutate it. If only one of the 3 possible 
   mutations forms the strong kmer then accept it and move to the next position. Otherwise, 
   stop and move to the next seed to the right and follow the same procedure.

   If no read trimming is allowed then accept the modified read if it is composed only of 
   the strong k-mers and the number of mutations is lower than MAX_ERRORS as specified 
   on the command line. 

   If trimming is allowed, then accept an antirely strong sub-read with the highest sum of match
   scores (the single base match score is 1 if the base was not mutated and -1 if it was).
   In that case MAX_ERRORS is disregarded and sub-reads may contain many mutations.

   The actual implementation is in <MakeReadCorrections>.

   In group: Edits
*/

#include "MainTools.h"
#include "math/Functions.h"
#include "ParseSet.h"
#include "Qualvector.h"
#include "kmers/KmerShape.h"
#include "kmer_freq/KmerShortMap.h"
#include "kmer_freq/MakeReadCorrections.h"
#include "paths/BaseErrorProb.h"
#include "feudal/BinaryStream.h"




void ExcludeReads(vec<read_id_t>& i_reads,
                  const vec<Bool>& read_to_exclude)
{
  long nexcluded = 0;
  if (!read_to_exclude.empty()) {
    vec<read_id_t> i_reads2;

    for (int i = 0; i < i_reads.isize(); i++){
      if (!read_to_exclude[i_reads[i]])
        i_reads2.push_back(i_reads[i]);
      else 
	nexcluded++;
    }
    i_reads = i_reads2;
  }
  cout << "Number of protected reads: " << nexcluded << endl;
}

/*****************************************************************************
 *
 *   main starts here
 *
 *****************************************************************************/




int main( int argc, char *argv[] )
{
  RunTime( ); 
  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA); 
  CommandArgument_String(RUN); 
  CommandArgument_String(K);
  CommandArgument_String(READS_IN);
  CommandArgument_String_OrDefault_Doc(PROTECTED_IN, "",
        "a binary file indicating which reads should not be error corrected by kept as they are");
  CommandArgument_String_OrDefault(QUALS_IN, "");
  CommandArgument_String_OrDefault(READS_OUT, "");
  CommandArgument_String_OrDefault(STRONG_IN, "");
  CommandArgument_Int_OrDefault(MAX_ERRORS, 4);
  CommandArgument_Bool_OrDefault(TRIM_READS, False);
  CommandArgument_Bool_OrDefault(THOROUGH, False);
  CommandArgument_Bool_OrDefault_Doc(KEEP_PARTIAL, False,
	"keep corrections even though not entire read is composed of strong kmers");
  CommandArgument_Bool_OrDefault(VERBOSE, False);
  EndCommandArguments; 



  /*****************************************************************************
   *
   *        Load data files
   *
   ****************************************************************************/

  if ( KEEP_PARTIAL && TRIM_READS ){
    cout << "ERROR! inconsistent option setting: you can either TRIM_READS or KEEP_PARTIAL" << endl;
    exit(1);
  }

  String run_dir = PRE + "/" + DATA + "/" + RUN;


  cout << Date() << " Reading " + READS_IN + ".fastb" << endl;
  vecbasevector reads( run_dir + "/" + READS_IN + ".fastb");     // Load reads
  cout << "Number of reads: " << reads.size() << endl;


  // Load Probability Distribution Table
  //BaseErrorProbProfile probProfile;


  // Load Quality Score (optional)
  vecqualvector quals;
  if (QUALS_IN != "") {
    //cout << Date() << "Reading " + QUALS_IN << endl;
    //quals =  vecqualvector( run_dir + "/" + QUALS_IN );
    cout << "FYI: quality score information is not currently used in this algorithm\n";
  } 

  if ( TRIM_READS )
    cout << "FYI: when trimming reads is allowed, MAX_ERRORS parameter is disregarded\n";


  // Load strong kmer table
  String strongKmersFileBase = run_dir + "/" + READS_IN + "." + "strong" + ".k";
  if (STRONG_IN != "")
    strongKmersFileBase =  run_dir + "/" + STRONG_IN + ".k";


  KmerShapeId Ks = (KmerShapeId) K;
  cout << Date() << " Reading " + strongKmersFileBase + ToString(Ks) << endl;
  KmerShortMap StrongTable(Ks, strongKmersFileBase + ToString(Ks));



  
  cout <<  Date() << " Data loaded." << endl;


  /*****************************************************************************
   *
   *         setup vector with the indices of the reads to process
   *
   ****************************************************************************/

  vec<Bool> read_is_strong(reads.size(), False);   // read_is_strong initialized with False

  vec<int> i_reads_todo( reads.size() );   
  for (int i = 0; i < i_reads_todo.isize(); i++) // by default all of them are to be done
      i_reads_todo[i] = i;

  if (!PROTECTED_IN.empty()) { //all non-protected reads will be processed
    vec<Bool> read_is_protected;
    
    BinaryReader::readFile(run_dir + "/" + PROTECTED_IN, &read_is_protected);
    ExcludeReads(i_reads_todo, read_is_protected);
    for ( int i = 0; i < read_is_protected.isize(); i++ ){
      if ( read_is_protected[ i ] )
	read_is_strong[ i ] = True; // we assume that protected => strong
    }

  }

  /*****************************************************************************
   *
   *       correct the reads
   *
   ****************************************************************************/

  // Perform Error Correction
  MakeReadCorrections(reads, 
                      i_reads_todo,
                      quals, 
                      StrongTable, 
                      read_is_strong,
                      VERBOSE, 
                      MAX_ERRORS, 
                      TRIM_READS,
                      THOROUGH,
		      KEEP_PARTIAL);
  
  
  /*****************************************************************************
   *
   *       save results
   *
   ****************************************************************************/


  // output reads (corrected and not).
  if ( READS_OUT != "" ) 
    reads.WriteAll( run_dir + "/" + READS_OUT + ".fastb" );


  // output read_is_strong info
  BinaryWriter::writeFile(run_dir + "/" + READS_OUT + ".is_strong",
                              read_is_strong);


  cout << Date() << " CorrectReadsByKmers finished." << endl;
  return 0;
}
