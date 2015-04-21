///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
//
//  Author: Neil Weisenfeld - Sep 2012 
//


#include "MainTools.h"
#include "simulation/Framework.h"
#include "simulation/ReferenceTools.h"
#include "simulation/ConstructionDReadGenerator.h"
#include "Qualvector.h"
#include <memory>

int main(int argc, char* argv[])
{
  RunTime();
  BeginCommandArguments;
  CommandArgument_String_OrDefault_Doc(FASTB_REF, "", "Reference genome for simulation.");
  CommandArgument_String_OrDefault_Doc(FASTA_REF, "", "Reference genome for simulation.");
  CommandArgument_Int_OrDefault_Doc(REF_ID, -1, "If specified, use only this record from reference fastb file.");
  CommandArgument_Int_OrDefault_Doc(REF_START, -1, "If specified, use only bases starting at this position in reference fastb file. For use with REF_ID and REF_STOP.");
  CommandArgument_Int_OrDefault_Doc(REF_STOP, -1, "If specified, use only bases stopping at this position in reference fastb file. For use with REF_ID and REF_START.");
  CommandArgument_UnsignedInt_OrDefault_Doc(G_RANDOM, 0, "Random genome size.");
  CommandArgument_UnsignedInt_OrDefault_Doc(RANDOM_SEED, 1252723, "Random seed.");
  CommandArgument_String_OrDefault_Doc(HEAD_RANDOM_REF_OUT, "", "If generating a random genome, save it to '<HEAD_RANDOM_REF_OUT>.fastb'.");
  CommandArgument_String_Doc(HEAD_IN_ERRGEN, "The input error generator pathname head.");
  CommandArgument_String_OrDefault_Doc(HEAD_OUT_PERFECT, "", "Perfect reads output head.");
  CommandArgument_String_OrDefault_Doc(HEAD_OUT_ERROR, "", "Error reads output head.");
  CommandArgument_String_OrDefault_Doc(HEAD_OUT_LOCS, "", "Perfect reads output head.");
  CommandArgument_UnsignedInt_OrDefault_Doc(LEN, 40000, "Mean read length.");
  CommandArgument_UnsignedInt_OrDefault_Doc(LEN_SIG, 10000, "Read lenght standard deviation.");
  CommandArgument_UnsignedInt_OrDefault_Doc(COVERAGE, 50, "Coverage.");
  CommandArgument_UnsignedInt_OrDefault_Doc(PAIRED_READ_LEN, 0,  "If paired, this is a distribution on read length.");
  CommandArgument_Bool_OrDefault_Doc(FW_ONLY, False, "Whether to generate reads in the FW direction only.");
  CommandArgument_Bool_OrDefault_Doc(PAIRED, False, "Whether to generate paired reads from the fragment.");
  CommandArgument_String_OrDefault_Doc(HEAD_OUT_PAIRS, "", "output head for .pairs file (when PAIRED=true)");
  CommandArgument_UnsignedInt_OrDefault_Doc(VERBOSITY, 1, "Verbosity level.");


  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0,
      "Number of threads to use (use all available processors if set to 0)");
  EndCommandArguments;

  // check that there's work to do
  if (HEAD_OUT_PERFECT == "" && HEAD_OUT_ERROR == "") {
    cout << Tag() << "Nothing to do. (specify HEAD_OUT_PERFECT and/or HEAD_OUT_ERROR)" << endl;
    return 1;
  }

  // type of reference
  BaseVecVec bvv_ref;

  if (FASTB_REF != "" ) {
    cout << Tag() << "reading reference from " << FASTB_REF << endl;
    ReferenceTools::fromFastb(&bvv_ref, FASTB_REF, REF_ID, REF_START, REF_STOP);
  }
  else if ( FASTA_REF != "" ) {
    cout << Tag() << "reading reference from " << FASTA_REF << endl;
    ReferenceTools::fromFasta(&bvv_ref,FASTA_REF, REF_ID, REF_START, REF_STOP);
  }
  else if ( G_RANDOM != 0 ) {
    cout << Tag() << "creating random genome of size " << G_RANDOM << endl;
    ReferenceTools::fromRandom(&bvv_ref, G_RANDOM, RANDOM_SEED);
  }
  cout << Tag() << "done creating/reading reference." << endl;

  if ( HEAD_RANDOM_REF_OUT != "" ) {
    cout << Tag() << "writing reference... " << HEAD_RANDOM_REF_OUT << endl;
    ReferenceTools::toFastb(&bvv_ref, HEAD_RANDOM_REF_OUT);
  }

  ErrorGenerator errgen(HEAD_IN_ERRGEN);

  ConstructionDReadGenerator cddg( bvv_ref, !FW_ONLY, &errgen, NUM_THREADS, 
	  				HEAD_OUT_PERFECT != "", 
					HEAD_OUT_LOCS != "", 
					HEAD_OUT_PAIRS != "" );

  if ( PAIRED ) {
    cout << Tag() << "generating paired reads" << endl;
    cddg.buildPairedReads(PAIRED_READ_LEN,COVERAGE,LEN,LEN_SIG);
  } else {
    cout << Tag() << "generating unpaired reads" << endl;
    cddg.buildUnpairedReads(LEN,COVERAGE);
  }

  if ( HEAD_OUT_PERFECT != "" ) {
    cout << Tag() << "writing perfect reads..." << endl;
    cddg.getPerfect().WriteAll(HEAD_OUT_PERFECT + ".fastb" );
  }

  if ( HEAD_OUT_ERROR != "" ) {
    cout << Tag() << "writing out reads with errors..." << endl;
    cddg.getReads().WriteAll(HEAD_OUT_ERROR + ".fastb" );
    cddg.getQuals().WriteAll(HEAD_OUT_ERROR + ".qualb" );
  }

  if ( HEAD_OUT_LOCS != "" ) {
    const RefLocusVec& locs = cddg.getReadLocs();
    cout << Tag() << "Saving genomic locations." << endl;
    BinaryIteratingWriter<RefLocusVec> writer((HEAD_OUT_LOCS+".locs").c_str());
    size_t nnn = locs.size();
    for ( size_t idx = 0; idx != nnn; ++idx )
      writer.write(locs[idx]);
    writer.close();
  }

  if ( PAIRED && HEAD_OUT_PAIRS != "" ) {
    cout << Tag() << "Saving pairs." << endl;
    cddg.getPairsManager().Write(HEAD_OUT_PAIRS + ".pairs" );
  }


  return 0;

}
