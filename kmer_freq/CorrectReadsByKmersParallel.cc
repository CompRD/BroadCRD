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
#include "system/WorklistN.h"
#include "feudal/BinaryStream.h"




void ExcludeReads(vec<read_id_t>& i_reads,
                  const vec<Bool>& read_to_exclude)
{
  if (!read_to_exclude.empty()) {
    vec<read_id_t> i_reads2;
    
    for (int i = 0; i < i_reads.isize(); i++)
      if (!read_to_exclude[i_reads[i]])
        i_reads2.push_back(i_reads[i]);

    i_reads = i_reads2;
  }

}

/*****************************************************************************
 *
 *   Processor class for handling multi-threading 
 *   For an explanation of this, see system/Worklist.h.
 *
 *****************************************************************************/



class Processor
{
public:
  // operator(): This is the function called in parallel by the Worklist
  // When called, it will correct a specified set of reads
  void operator() ( unsigned int ipe ) {
    
   

    cout << ( "proc ipe = " + ToString(ipe) + " n_reads_todo = " +
              ToString((*mp_i_reads_todo_pe)[ipe].isize()) ) << endl;



    // Perform Error Correction
    MakeReadCorrections(*mp_reads, 
                        (*mp_i_reads_todo_pe)[ipe],
                        *mp_quals, 
                        *mp_StrongTable, 
                        *mp_read_is_strong,
                        m_VERBOSE, 
                        m_MAX_ERRORS, 
                        m_TRIM_READS,
                        m_THOROUGH);
  

  }
  
  static vecbasevector* mp_reads;
  static vecqualvector* mp_quals;
  static KmerShortMap* mp_StrongTable;
  static vec<Bool>* mp_read_is_strong;
  static Bool m_VERBOSE;
  static int m_MAX_ERRORS;
  static Bool m_TRIM_READS;
  static Bool m_THOROUGH;

  static vec< vec<int> >* mp_i_reads_todo_pe;
};
  


/*****************************************************************************
 *
 *   Declarations of static variables in the Processor class.
 *
 *****************************************************************************/
  
vecbasevector* Processor::mp_reads;
vecqualvector* Processor::mp_quals;
KmerShortMap* Processor::mp_StrongTable;
vec<Bool>* Processor::mp_read_is_strong;
Bool Processor::m_VERBOSE;
int Processor::m_MAX_ERRORS;
Bool Processor::m_TRIM_READS;
Bool Processor::m_THOROUGH;

vec< vec<int> >* Processor::mp_i_reads_todo_pe;
 

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
  CommandArgument_String_OrDefault(QUALS_IN, "");
  CommandArgument_String_OrDefault(READS_OUT, "");
  CommandArgument_String_OrDefault(READ_ID,"");
  CommandArgument_String_OrDefault(STRONG_IN, "");
  CommandArgument_Int_OrDefault(MAX_ERRORS, 4);
  CommandArgument_Bool_OrDefault(TRIM_READS, False);
  CommandArgument_Bool_OrDefault(THOROUGH, False);
  CommandArgument_Bool_OrDefault(VERBOSE, False);

  // Parallelization parameter
  CommandArgument_UnsignedInt_OrDefault_Doc( N_PROCS, 1, "Use this many processes in parallel (or 1 if not parallelizing.)" );
  EndCommandArguments; 





  /*****************************************************************************
   *
   *        Load data files
   *
   ****************************************************************************/

  String run_dir = PRE + "/" + DATA + "/" + RUN;


  cout << Date() << " Reading " + READS_IN + ".fastb" << endl;
  vecbasevector reads( run_dir + "/" + READS_IN + ".fastb");     // Load reads



  // Load Probability Distribution Table
  //BaseErrorProbProfile probProfile;

  
  // Load Quality Score (optional)
  vecqualvector quals;
  if (QUALS_IN != "") {
    cout << Date() << "Reading " + QUALS_IN << endl;
    quals =  vecqualvector( run_dir + "/" + QUALS_IN );
  } 
  else 
    cout << "FYI: quality score information is not used in this algorithm\n";

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

  vec<int> i_reads_todo;

  if (!READ_ID.empty()) {  // READ_ID specified

    cout << Date() << " reads to do: " + READ_ID << endl;
    ParseIntSet( READ_ID, i_reads_todo );
  } 
  else {  // READ_ID not specified => all non-protected reads will be processed

    i_reads_todo.resize(reads.size());
    for (int i = 0; i < i_reads_todo.isize(); i++)
      i_reads_todo[i] = i;

    // Check for read protection
    if ( IsRegularFile(run_dir + "/reads.protected.orig") ) {
      vec<Bool> read_is_protected;
      BinaryReader::readFile(run_dir + "/reads.protected.orig" , &read_is_protected);

      ExcludeReads(i_reads_todo, read_is_protected);

      cout << Date() << " reads to do: all but protected: reads.protected.orig" << endl;
    }
    else {
      cout << Date() << " reads to do: all" << endl;
      
    }
  } 


  vec< vec<int> > i_reads_todo_pe(N_PROCS);
  for (unsigned int i = 0, ipe = 0; i < i_reads_todo.size(); i++, ipe++ ) {
    ipe %= N_PROCS;
    i_reads_todo_pe[ipe].push_back(i_reads_todo[i]);
  }

  
  for (unsigned int ipe = 0; ipe < N_PROCS; ipe++ ) {
    cout << "ipe = " << ipe << "n_reads_todo = " << i_reads_todo_pe[ipe].isize() << endl;
  }



  /*****************************************************************************
   *
   *        load data into processor class
   *
   ****************************************************************************/

  // read_is_strong initialized with False
  vec<Bool> read_is_strong(reads.size(), False); 


  // Load data structures into the Processor class
  Processor::mp_reads = &reads;
  Processor::mp_quals = &quals;
  Processor::mp_StrongTable = &StrongTable;
  Processor::mp_read_is_strong = &read_is_strong;
  Processor::m_VERBOSE = VERBOSE;
  Processor::m_MAX_ERRORS = MAX_ERRORS;
  Processor::m_TRIM_READS = TRIM_READS;
  Processor::m_THOROUGH = THOROUGH;
  
  Processor::mp_i_reads_todo_pe = &i_reads_todo_pe;


  {
      parallelFor(0u,N_PROCS,Processor(),N_PROCS);
  }




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


  cout << Date() << " CorrectReadsByKmers finished" << endl;
  return 0;
}
