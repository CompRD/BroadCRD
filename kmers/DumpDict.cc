//////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file DumpDict.cc
 * \author tsharpe
 * \date Dec 12, 2011
 *
 * \brief
 */
#include "kmers/ReadPather.h"
#include "feudal/BinaryStream.h"

void dumpDict( KmerDict<28> const& dict )
{
    typedef KmerDict<28>::OCItr OItr;
    typedef KmerDict<28>::ICItr IItr;
    for ( OItr oitr(dict.begin()), oend(dict.end()); oitr != oend; ++oitr )
        for ( IItr iitr(oitr->begin()), iend(oitr->end()); iitr != iend; ++iitr )
            if ( iitr->getCanonicalForm() == CanonicalForm::REV )
                std::cout << static_cast<KMer<28> >(*iitr) << std::endl;
}

int main( int argc, char** argv )
{
    KmerDict<28> dict(0);
    BinaryReader::readFile("/local/scratch/tsharpe/Assembly_of_1.2_Mb_stickleback_region/data/run/frag_reads_corr.k28.ug.dict",&dict);
    dumpDict(dict);
}
