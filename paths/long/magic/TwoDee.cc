///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "paths/long/magic/NanoData.h"
#include "IntPairVec.h"
#include "lookup/LookAlign.h"
#include "paths/long/CreateGenome.h"
#include "pairwise_aligners/ClusterAligner.h"
#include "kmers/KmerRecord.h"
#include "paths/RemodelGapTools.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "PrintAlignment.h"
#include "paths/long/ultra/MultipleAligner.h"
#include "paths/long/ultra/FounderAlignment.h"
#include "paths/long/magic/KmerProfileHMM.h"
#include "paths/long/magic/NanoScaler.h"
#include "paths/long/magic/KmerHHAlign.h"
#include "paths/long/magic/NanoRCModel.h"
#include "pairwise_aligners/SmithWatAffine.h"
#include "PackAlign.h"

namespace {

};

int main(int argc, char* argv[]) {
    RunTime();

    BeginCommandArguments;
    CommandArgument_String_Doc(NANO, ".nano dataset");
    CommandArgument_Int_Doc(READ, "read to process");
    CommandArgument_Int_OrDefault_Doc(MINHITS, 5,
            "minimum number of 12-mer matches between reads");
    CommandArgument_Double_OrDefault_Doc(ERRTOL,0.80,
            "dump reads with score > ERRTOL*overlap_length");
    CommandArgument_Int_OrDefault_Doc(MIN_OVERLAP_BASES,300,
            "don't include aligns with overlap < MIN_OVERLAP_BASES")
    CommandArgument_Bool_OrDefault_Doc(VISUAL, False,
            "print SmithWatBandedA visual alignment");
    CommandArgument_Bool_OrDefault_Doc(DEBUG, False, "verbose debugging");

    EndCommandArguments;

    ForceAssertGe(READ,0);

    cout << Date() << ": reading nano file " << NANO << endl;
    NanoData nano;
    vec<align> aligns;
    vecbasevector gang;

    BinaryReader::readFile(NANO, &nano);
    ForceAssertGt((int)nano.GetSize(), READ);
    cout << Date() << ": read " << nano.GetSize() << " reads" << endl;


    // Instantiate profile HMM
    auto& e0_raw = nano.GetEvents(NanoData::TEMP)[READ];
    auto& m0 = nano.GetModels(NanoData::TEMP)[READ];


    auto& e1_raw = nano.GetEvents(NanoData::COMP)[READ];
    // model is already RC'd, or so it seems
    auto m1 = nano.GetModels(NanoData::COMP)[READ];
//    auto m1_fw = nano.GetModels(NanoData::COMP)[READ];
//    auto m1 = NanoRCModel( m1_fw ).Get();

    // dump FASTA
    auto fw_bases = nano.GetBases(NanoData::TEMP)[READ];
    auto rc_bases = nano.GetBases(NanoData::COMP)[READ];
    rc_bases.ReverseComplement();
    auto co_bases = nano.GetBases(NanoData::CONS)[READ];

    cout << Date() << ": read " << READ << ", template sz=" << fw_bases.size()
            << ", comp sz=" << rc_bases.size() << ", cons sz="
            <<  co_bases.size() << endl;

    Ofstream(OUT, "test.fasta");
    fw_bases.Print(OUT, "fw"+ToString(READ));
    rc_bases.Print(OUT, "rc"+ToString(READ));
    co_bases.Print(OUT, "co"+ToString(READ));

    alignment al;
    unsigned int errs;
//    cout << "=== FW->CONS ===========================" << endl;
//    errs = SmithWatAffine(fw_bases, co_bases, al, true, true);
//    cout << "Errors=" << errs << endl;
//    cout << al << endl;
//    PrintVisualAlignment(False, cout, fw_bases, co_bases, al );

//    cout << "=== RC->CONS ===========================" << endl;
//    errs = SmithWatAffine(rc_bases, co_bases, al, true, true);
//    cout << "Errors=" << errs << endl;
//    cout << al << endl;
//    PrintVisualAlignment(False, cout, rc_bases, co_bases, al );

    cout << "=== RC->FW ===========================" << endl;
    errs = SmithWatAffine(rc_bases, fw_bases, al, true, true);
    cout << "Errors=" << errs << endl;
    cout << al << endl;
    PrintVisualAlignment(False, cout, rc_bases, fw_bases, al );

    auto hairpin = nano.GetHairpin()[READ];
    cout << "hairpin:" << endl;
    hairpin.Dump();

    cout << "template:" << endl;
    nano.GetTempEvents()[READ].Dump();

    cout << "complement:" << endl;
    nano.GetCompEvents()[READ].Dump();


//    cout << "fw model:" << endl;
//    m1_fw.Dump();
//
//    cout << "rc model:" << endl;
//    m1.Dump();

    auto e0 = NanoScaler(
            m0("shift")->get_double(),
            m0("scale")->get_double(),
            m0("drift")->get_double() ).Scale(
                    e0_raw["mean"]->get_vec_double(),
                    e0_raw["start"]->get_vec_double());

    auto e1 = NanoScaler(
            m1("shift")->get_double(),
            m1("scale")->get_double(),
            m1("drift")->get_double() ).ReverseAndScale(
                    e1_raw["mean"]->get_vec_double(),
                    e1_raw["start"]->get_vec_double() );

    KmerProfileHMM q( e0,
                      m0["level_mean"]->get_vec_double(),
                      m0["level_stdv"]->get_vec_double() );

    KmerProfileHMM p( e1,
                      m1["level_mean"]->get_vec_double(),
                      m1["level_stdv"]->get_vec_double() );

    KmerHHAlign aligner;
    align a;
    auto score = aligner.Align(q, p, a, DEBUG);

    cout << Date() << ": done" << endl;
}
