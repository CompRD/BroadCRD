////////////////////////////////////////////////////////////////////////////
//                  SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//      This software and its documentation are copyright (2009) by the   //
//  Broad Institute.  All rights are reserved.  This software is supplied //
//  without any warranty or guaranteed support whatsoever. The Broad      //
//  Institute is not responsible for its use, misuse, or functionality.   //
////////////////////////////////////////////////////////////////////////////
/*
 * \file AlignPairsToHyperTest.cc
 * \author tsharpe
 * \date Dec 15, 2009
 *
 * \brief
 */
#include "MainTools.h"
#include "PairsManager.h"
#include "String.h"
#include "feudal/BaseVec.h"
#include "paths/AlignPairsToHyper.h"
#include "paths/HyperBasevector.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"
#include "paths/SeqOnHyper.h"
#include "feudal/BinaryStream.h"

HyperBasevector getHyperBasevector( HyperKmerPath const& h,
                                    String const& wrun_dir, int K )
{
    KmerBaseBroker kbb( wrun_dir, K );
    return HyperBasevector( h, kbb );
}

int main( int argc, char** argv )
{
    RunTime();
    BeginCommandArguments;
    CommandArgument_String_OrDefault(DIR,"/wga/dev/WGAdata/projects/ALLPATHS/N.crassa/13327.5,13350.3,5,6,202GC.1,202EA.1,2,3,5,6,7,8,201FN.1,2,3,5,6,7,8,13174.1,2,3,5,6,7,8.iainm.1apr/iainm.1apr");
    CommandArgument_String_OrDefault(SUBDIR, "test");
    CommandArgument_Int_OrDefault(K, 20);
    CommandArgument_Bool_OrDefault(SHOW_PAIR_ALIGNS, False);
    EndCommandArguments;

    String sub_dir = DIR + "/ASSEMBLIES/" + SUBDIR;
    String wrun_dir = sub_dir + "/run";

    HyperKmerPath h;
    BinaryReader::readFile( sub_dir + "/hyper.initial", &h );

    HyperBasevector hb = getHyperBasevector(h,wrun_dir,K);

    size_t starting_read_id = 79747142UL;
    size_t ending_read_id = 89506932UL;

    vecbvec reads( DIR + "/reads.fastb" );
    PairsManager pairs( DIR + "/reads.pairs" );

    vec<CompressedSeqOnHyper> csaligns;
    vec<vec<int> > csaligns_index;
    vec<IndexPair> pair_places;

    AlignPairsToHyper( h, hb, reads, pairs, sub_dir, csaligns,
                        csaligns_index, pair_places, SHOW_PAIR_ALIGNS );
}
