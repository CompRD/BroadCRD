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
#include "simulation/ReadSimulatorSimpleGenerator.h"
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
  CommandArgument_String_OrDefault_Doc(HEAD_IN_PERFECT, "", "Perfect reads input head. If specified, use to perturb.");
  CommandArgument_String_OrDefault_Doc(HEAD_OUT_PERFECT, "", "Perfect reads output head.");
  CommandArgument_String_OrDefault_Doc(HEAD_OUT_ERROR, "", "Error reads output head.");
  CommandArgument_String_OrDefault_Doc(HEAD_OUT_ERROR_LOCS, "", "Error locations in reads output head.");
  CommandArgument_String_OrDefault_Doc(HEAD_OUT_LOCS, "", "Perfect reads output head.");
  CommandArgument_Double_OrDefault_Doc(TILING, 0.0, "Probability of tiled reads being sequenced.");
  CommandArgument_LongLongSet_OrDefault_Doc(CIRCULARS, "", "List of indexes of the circular chromosomes in the reference.");
  CommandArgument_DoubleSet_OrDefault_Doc(COPY_NUMBERS, "", "List of the copy numbers of each chromosome in the reference.");
  CommandArgument_UnsignedInt_OrDefault_Doc(LEN, 40000, "Mean read length.");
  CommandArgument_UnsignedInt_OrDefault_Doc(LEN_SIG, 10000, "Read lenght standard deviation.");
  CommandArgument_UnsignedInt_OrDefault_Doc(LEN_MIN, LEN/100, "Minimum read length.");
  CommandArgument_UnsignedInt_OrDefault_Doc(COVERAGE, 50, "Coverage.");
  CommandArgument_UnsignedInt_OrDefault_Doc(PAIRED_READ_LEN, 0,  "If paired, this is a distribution on read length.");
  CommandArgument_UnsignedInt_OrDefault_Doc(PAIRED_READ_LEN_SIG, 0,  "If paired, this is a distribution on read length.");
  CommandArgument_UnsignedInt_OrDefault_Doc(PAIRED_READ_LEN_MIN, 0,  "If paired, this is a distribution on read length.");
  CommandArgument_Double_OrDefault_Doc(ERR_DEL, 0.0, "Deletion error rate.");
  CommandArgument_Double_OrDefault_Doc(ERR_SUB, 0.0, "Substitution error rate.");
  CommandArgument_Double_OrDefault_Doc(ERR_INS, 0.0, "Insertion error rate.");
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


  if ( HEAD_RANDOM_REF_OUT != "" ) {
    cout << Tag() << "writing reference... " << HEAD_RANDOM_REF_OUT << endl;
    ReferenceTools::toFastb(&bvv_ref, HEAD_RANDOM_REF_OUT);
  }

  vec<bool> circular( bvv_ref.size(), false );

  ReadSimulatorSimpleGenerator rssg2( bvv_ref, circular, NUM_THREADS, PAIRED, true, true, false );

  RandomLength insert_len( LEN, LEN_SIG, LEN_MIN, RANDOM_SEED );

  if ( PAIRED ) {
    cout << Tag() << "generating paired read locs" << endl;

    RandomLength read_len( PAIRED_READ_LEN, PAIRED_READ_LEN_SIG, PAIRED_READ_LEN_MIN, RANDOM_SEED );

    for (size_t ref_id = 0; ref_id < bvv_ref.size(); ++ref_id )
      rssg2.computeReadsSizeAndPositionPaired(ref_id, COVERAGE, RANDOM_SEED, FW_ONLY, TILING, insert_len, read_len);

  } else {
    cout << Tag() << "generating unpaired read locs" << endl;

    for (size_t ref_id = 0; ref_id < bvv_ref.size(); ++ref_id )
      rssg2.computeReadsSizeAndPositionUnpaired(ref_id, COVERAGE, RANDOM_SEED, FW_ONLY, TILING, insert_len );

  }


  ForceAssert( HEAD_OUT_PERFECT != "" || HEAD_OUT_ERROR != "" );

  if ( HEAD_OUT_PERFECT != "" ) {
    cout << Tag() << "building perfect reads" << endl;
    rssg2.buildReadsFromPositions(RANDOM_SEED);
    rssg2.getReads().WriteAll( HEAD_OUT_PERFECT + ".fastb" );
    rssg2.getQuals().WriteAll( HEAD_OUT_PERFECT + ".qualb" );
  }


  if ( HEAD_OUT_ERROR != "" ) {
    cout << Tag() << "building reads with errors" << endl;
    rssg2.buildReadsFromPositionsWithErrors(RANDOM_SEED,ERR_DEL, ERR_INS, ERR_SUB );
    rssg2.getReads().WriteAll( HEAD_OUT_ERROR + ".fastb" );
    rssg2.getQuals().WriteAll( HEAD_OUT_ERROR + ".qualb" );
  }


  if ( HEAD_OUT_LOCS != "" ) {
    const RefLocusVec& locs = rssg2.getReadLocs();
    cout << Tag() << "Saving genomic locations." << endl;
    BinaryIteratingWriter<RefLocusVec> writer((HEAD_OUT_LOCS+".locs").c_str());
    size_t nnn = locs.size();
    for ( size_t idx = 0; idx != nnn; ++idx )
      writer.write(locs[idx]);
    writer.close();
  }

  if ( HEAD_OUT_ERROR_LOCS != "" ) {
    cout << Tag() << "Saving individual errors." << endl;
    const ReadErrorVecVec& err_locs = rssg2.getReadErrors();
    err_locs.WriteAll( HEAD_OUT_ERROR_LOCS + ".err_locs" );
  }


  if ( PAIRED && HEAD_OUT_PAIRS != "" ) {
    cout << Tag() << "Saving pairs." << endl;
    rssg2.getPairsManager().Write(HEAD_OUT_PAIRS + ".pairs" );
  }

  return 0;

}
