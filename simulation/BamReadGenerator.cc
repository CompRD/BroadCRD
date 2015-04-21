///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file BamReadGenerator.cc
 * \author tsharpe
 * \date Oct 26, 2012
 *
 * \brief
 */

#include "simulation/BamReadGenerator.h"
#include "STLExtensions.h"
#include "lookup/SAM2CRD.h"

BamReadGenerator::BamReadGenerator( CachedBAMFile const& bam,
                                        String const& regions )
{
    init();
    processBAM(bam,regions);
}

BamReadGenerator::BamReadGenerator( std::vector<CachedBAMFile> const& bams,
                                          String const& regions )
{
    init();
    typedef std::vector<CachedBAMFile>::const_iterator Itr;
    for ( Itr itr(bams.begin()), end(bams.end()); itr != end; ++itr )
        processBAM(*itr,regions);
}

void BamReadGenerator::processBAM( CachedBAMFile const& bam,
                                    String const& regions )
{
    SAM::BAMFile bamFile(bam.getPath(),regions);
    if ( !bamFile.hasReferenceDictionary() )
        FatalErr(bam.getPath() << " has no reference dictionary.");

    //std::cout << "Ref URI: " << bamFile.getUniqueRefURI() << std::endl;

    SeqGrabber<std::back_insert_iterator<vecbvec> > sg(std::back_inserter(_bvv));
    QualGrabber<std::back_insert_iterator<vecqvec> > qg(std::back_inserter(_qvv));
    RefLocGrabber<std::back_insert_iterator<RefLocusVec> > rlg(std::back_inserter(_refv),bamFile);
    PairsBuilder pb(_pm);
    copyIf(bamFile.begin(),bamFile.end(),mux(sg,mux(qg,mux(rlg,pb))),rlg);
}

void BamReadGenerator::init()
{
    size_t const RESERVE = 10000000;
    _bvv.reserve(RESERVE);
    _qvv.reserve(RESERVE);
    _refv.reserve(RESERVE);
}
