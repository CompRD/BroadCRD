// // Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

/*
  Program: EvalReadErrorEdits
  
  
  Given a set of reads and corresponding "erro corrected" versions find out how many of the reads were corrected and how many were miscorrected

*/

#include <strstream>
#include <limits>
#include "MainTools.h"
#include "FastIfstream.h"
#include "FeudalMimic.h"
#include "math/Functions.h"
#include "PackAlign.h"
#include "ReadPairing.h"
#include "ReadLocation.h"
#include "CommonSemanticTypes.h"
#include "Bitvector.h"
#include "lookup/ImperfectLookup.h"
#include "lookup/PerfectLookup.h"
#include "lookup/LookAlign.h"

int main( int argc, char *argv[] )
{
  cout << Date() << " Starting! ... " << endl;
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_String(ORIGINAL_READS);
  CommandArgument_String(CORRECTED_READS);
  CommandArgument_String_OrDefault(GENOME_LOOKUP,"genome");
  CommandArgument_Int_OrDefault(K_LOOKUP,12);
  CommandArgument_Int_OrDefault_Doc(ERRORS, 0, 
       "Evaluate allowing this number of ERRORS in alignments");
  EndCommandArguments;

  ForceAssertGe( ERRORS, 0 );

  String data_dir = PRE + "/" + DATA;
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  
  
  vecbasevector original_reads( run_dir + "/" + ORIGINAL_READS + ".fastb" );
  size_t nreads    = original_reads.size();
  vecbasevector crtd_reads( run_dir + "/" + CORRECTED_READS + ".fastb" );
  ForceAssertEq( nreads, crtd_reads.size() );

  //vecbasevector genome( data_dir + "/" + GENOME + ".fastb" );
  vec<look_align> crtd_aligns;
  if ( ERRORS == 0 ){
    cout << Date() << " running PerfectLookup..." << endl;
    PerfectLookup( K_LOOKUP /* seed alignments to the reference on kmers of this size */, 
		   crtd_reads, data_dir + "/"+GENOME_LOOKUP+".lookup", crtd_aligns, FW_OR_RC );
    cout << Date() << " PerfectLookup done. " << endl;
  }
  else{
    vec< int > best_errors;
    cout << Date() << " running ImperfectLookup..." << endl;
    ImperfectLookup( K_LOOKUP /* seed alignments to the reference on kmers of this size */, 
		     crtd_reads, data_dir + "/"+GENOME_LOOKUP+".lookup", crtd_aligns, best_errors, FW_OR_RC, 1000 );
    cout << Date() << " ImperfectLookup done. " << endl;
  }

  long nCrtdGenomic = 0;
  vec<Bool> crtd_genomic( crtd_reads.size(), False);
  for ( size_t ialign = 0; ialign < crtd_aligns.size( ); ialign++ ){    
    const look_align& la = crtd_aligns[ialign];
    int readId = la.query_id;
    int errors = la.mutations;
    if ( crtd_genomic[readId] == True ) continue;
    if ( la.FullLength() && errors <= ERRORS){
      crtd_genomic[readId] = True;
    }
  }
  crtd_aligns.resize( 0 );
  for ( int i = 0; i < crtd_genomic.isize(); i++ )
    if ( crtd_genomic[i] ) nCrtdGenomic++;
  cout << "number of genomic reads among corrected: " << nCrtdGenomic << endl;

  vec<look_align> orig_aligns;
  if ( ERRORS == 0 ){
    cout << Date() << " running PerfectLookup..." << endl;
    PerfectLookup( K_LOOKUP /* seed alignments to the reference on kmers of this size */, 
		   original_reads, data_dir + "/"+GENOME_LOOKUP+".lookup", orig_aligns, FW_OR_RC );
    cout << Date() << " PerfectLookup done. " << endl;
  }
  else{
    vec< int > best_errors;
    cout << Date() << " running ImperfectLookup..." << endl;
    ImperfectLookup( K_LOOKUP /* seed alignments to the reference on kmers of this size */, 
		     original_reads, data_dir + "/"+GENOME_LOOKUP+".lookup", orig_aligns, best_errors, FW_OR_RC, 1000 );
    cout << Date() << " ImperfectLookup done. " << endl;
  }


  long nOrigGenomic = 0;
  vec<Bool> orig_genomic( crtd_reads.size(), False);
  for ( int ialign = 0; ialign < orig_aligns.isize( ); ialign++ ){    
    const look_align& la = orig_aligns[ialign];
    int readId = la.query_id;
    int errors = la.mutations;
    if ( orig_genomic[readId] == True ) continue;
    if ( la.FullLength() && errors <= ERRORS ){
      orig_genomic[readId] = True;
    }
  }
  orig_aligns.resize( 0 );
  for ( int i = 0; i < orig_genomic.isize(); i++ )
    if ( orig_genomic[i] ) nOrigGenomic++;
  cout << "number of genomic reads among original: " << nOrigGenomic << endl;

  long nreadsCrtd        = 0; // will hold number of reads that were changed
  long nreadsWellCrtd    = 0; // will hold number of reads that were corrected well
  long nreadsMissCrtd    = 0; // will hold number of reads not corrected well enough
  long nreadsCrtdUnknown = 0;
  for ( size_t id1 = 0; id1 < nreads; id1++ ){
    if ( original_reads[id1] != crtd_reads[id1] ){ // corrected read
      ++nreadsCrtd;
      if ( crtd_genomic[id1] )
	++nreadsWellCrtd;
      else if ( orig_genomic[id1] )
	++nreadsMissCrtd;
      else 
	++nreadsCrtdUnknown;
    }
  }
  cout << "reads:                                 " << nreads << "\n";
  cout << "corrected reads:                       " << nreadsCrtd << "\n";
  cout << "corrected genomic:                     " << nreadsWellCrtd << "\n";
  cout << "misscorrected (genomic to non-genomic):" << nreadsMissCrtd << "\n";
  cout << "modified but still non-genomic:        " << nreadsCrtdUnknown << "\n";
 

}
